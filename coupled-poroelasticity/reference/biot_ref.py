"""
Reference implementation of 1D Biot consolidation (Phase 2).

This module does two jobs:

1. exact(): the closed-form Biot/Terzaghi solution for pressure and
   displacement, used as ground truth.
2. simulate(): a NumPy/SciPy implementation of *exactly* the discretization
   used by the OpenFOAM solver biotConsolidationFoam (collocated cell-centred
   p and w, Gauss-linear gradients, central Laplacians, implicit Euler,
   optional pressure stabilization, fixed-stress split or monolithic solve).
   The OpenFOAM solver must reproduce it to within its iteration tolerance,
   which checks the C++ implementation far more tightly than comparing with
   the exact solution alone.

Geometry and sign conventions
-----------------------------
z = 0 is the loaded, drained top; z = H is the fixed, impermeable base; z
points down and w (vertical displacement) is positive downward, so the
settlement of the top is s(t) = w(0, t). Stress is tension-positive; the
applied load is a compressive traction of magnitude sigma0.

Governing equations (uniaxial strain, quasi-static, small strain)
-----------------------------------------------------------------
    equilibrium : d/dz( Mc * dw/dz - alpha * p ) = 0
    fluid mass  : (1/M) dp/dt + alpha * d(eps)/dt = d/dz( mob * dp/dz ),  eps = dw/dz

Mc = K + 4G/3 (drained constrained modulus), M = Biot modulus,
alpha = Biot coefficient, mob = k/mu (mobility).
"""
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.special import erf, erfc


class ConvergenceError(RuntimeError):
    """A time step failed; no unconverged field may be used as a result."""


def validate_inputs(P, N, dt, write_times, beta_factor, tol, max_iter, start_time=0.0):
    values = [P.H, P.Mc, P.M, P.alpha, P.mob, P.sigma0, dt, beta_factor, tol]
    if not np.all(np.isfinite(values)) or min(values) <= 0 or P.alpha > 1:
        raise ValueError("Material values, dt and tolerances must be finite and positive; alpha <= 1")
    if int(N) != N or N < 2 or int(max_iter) != max_iter or max_iter < 1:
        raise ValueError("N must be an integer >= 2 and max_iter an integer >= 1")
    times = np.asarray(write_times, dtype=float)
    if times.ndim != 1 or len(times) == 0 or not np.all(np.isfinite(times)):
        raise ValueError("Output times must be a nonempty finite 1D list")
    if not np.isfinite(start_time) or start_time < 0 or times[0] <= start_time or np.any(np.diff(times) <= 0):
        raise ValueError("Output times must increase strictly after start_time")
    steps = (times-start_time)/dt
    if not np.allclose(steps, np.rint(steps), rtol=0, atol=1e-7):
        raise ValueError("Every output time must be on the fixed time-step grid")


class Params:
    """Material, loading and geometry parameters (SI units)."""

    def __init__(self, H=1.0, Mc=1e8, M=2e8, alpha=0.8, mob=1.14e-10, sigma0=1e6):
        self.H, self.Mc, self.M, self.alpha, self.mob, self.sigma0 = H, Mc, M, alpha, mob, sigma0

    @property
    def Mu(self):
        """Undrained constrained modulus."""
        return self.Mc + self.alpha ** 2 * self.M

    @property
    def p0(self):
        """Instantaneous (undrained) pore-pressure response to the load."""
        return self.alpha * self.M * self.sigma0 / self.Mu

    @property
    def cv(self):
        """Consolidation coefficient: mobility / constrained storage."""
        return self.mob / (1.0 / self.M + self.alpha ** 2 / self.Mc)

    @property
    def s0(self):
        """Undrained (instantaneous) settlement."""
        return self.sigma0 * self.H / self.Mu

    @property
    def s_inf(self):
        """Final (drained) settlement."""
        return self.sigma0 * self.H / self.Mc

    def with_cv(self, cv):
        """Copy with the mobility adjusted to give the requested cv."""
        q = Params(self.H, self.Mc, self.M, self.alpha, self.mob, self.sigma0)
        q.mob = cv * (1.0 / q.M + q.alpha ** 2 / q.Mc)
        return q


# ---------------------------------------------------------------------------
# Exact solution
# ---------------------------------------------------------------------------

def exact(P, z, t, n_terms=400):
    """Exact p(z, t) and w(z, t) for t > 0.

    Uses Terzaghi's Fourier series, except at very early times (Tv < 1e-3),
    where the series would need thousands of terms to resolve the thin
    boundary layer at the drained face (Gibbs ringing). There the
    semi-infinite solution is used instead; its error from ignoring the base
    is of order erfc(1/sqrt(Tv)), below 1e-13 for Tv < 1e-3.
    """
    z = np.asarray(z, dtype=float)
    Tv = P.cv * t / P.H ** 2
    if t < 0:
        raise ValueError("t must be nonnegative")
    if t == 0:
        return np.where(z == 0, 0.0, P.p0), (P.H-z)*P.sigma0/P.Mu
    if Tv < 1e-3:
        d = 2.0 * np.sqrt(P.cv * t)                    # diffusion length
        x = z / d
        p = P.p0 * erf(x)
        ierfc = np.exp(-x ** 2) / np.sqrt(np.pi) - x * erfc(x)   # integral of erfc from x to inf
        int_p = P.p0 * ((P.H - z) - d * ierfc)                   # integral of p from z to H
        w = (P.sigma0 * (P.H - z) - P.alpha * int_p) / P.Mc
        return p, w
    sp = np.zeros_like(z)
    sw = np.zeros_like(z)
    for m in range(n_terms):
        Mm = 0.5 * np.pi * (2 * m + 1)
        e = np.exp(-Mm ** 2 * Tv)
        sp += (2.0 / Mm) * np.sin(Mm * z / P.H) * e
        sw += (2.0 / Mm ** 2) * np.cos(Mm * z / P.H) * e
    p = P.p0 * sp
    w = (P.sigma0 * (P.H - z) - P.alpha * P.p0 * P.H * sw) / P.Mc
    return p, w


def exact_settlement(P, t):
    return exact(P, np.array([0.0]), t)[1][0]


# ---------------------------------------------------------------------------
# Discrete operators, N uniform cells (mirrors the OpenFOAM discretization)
# ---------------------------------------------------------------------------

def operators(P, N):
    h = P.H / N
    z = (np.arange(N) + 0.5) * h

    # Mechanics, cell balance sigma_{i+1/2} - sigma_{i-1/2} = 0 with face
    # stress sigma_f = Mc * snGrad(w) - alpha * p_f (p_f linearly interpolated):
    #     A_w w + B_p p = b_w
    A_w = np.zeros((N, N)); B_p = np.zeros((N, N)); b_w = np.zeros(N)
    for i in range(N):
        if i < N - 1:                                   # lower internal face
            A_w[i, i + 1] += P.Mc / h; A_w[i, i] -= P.Mc / h
            B_p[i, i] -= P.alpha / 2; B_p[i, i + 1] -= P.alpha / 2
        else:                                           # base: w = 0, dp/dz = 0
            A_w[i, i] -= P.Mc / (h / 2)
            B_p[i, i] -= P.alpha
        if i > 0:                                       # upper internal face
            A_w[i, i] -= P.Mc / h; A_w[i, i - 1] += P.Mc / h
            B_p[i, i] += P.alpha / 2; B_p[i, i - 1] += P.alpha / 2
        else:                                           # loaded top: sigma = -sigma0
            b_w[i] = -P.sigma0

    # Volumetric strain, Gauss-linear: eps_i = (w_f,lower - w_f,upper) / h
    E_w = np.zeros((N, N)); e_c = np.zeros(N)
    for i in range(N):
        if i < N - 1:
            E_w[i, i] += 0.5 / h; E_w[i, i + 1] += 0.5 / h
        if i > 0:
            E_w[i, i] -= 0.5 / h; E_w[i, i - 1] -= 0.5 / h
        else:  # top-face value from the fixed gradient: w_b = w_0 + (sigma0/Mc) h/2
            E_w[i, i] -= 1.0 / h
            e_c[i] -= 0.5 * P.sigma0 / P.Mc

    # Discrete Laplacian with the pressure BCs (p = 0 at top, half-cell away;
    # zero gradient at base)
    Lh = np.zeros((N, N))
    for i in range(N):
        if i < N - 1:
            Lh[i, i + 1] += 1 / h ** 2; Lh[i, i] -= 1 / h ** 2
        if i > 0:
            Lh[i, i - 1] += 1 / h ** 2; Lh[i, i] -= 1 / h ** 2
        else:
            Lh[i, i] -= 2 / h ** 2
    return dict(h=h, z=z, A_w=A_w, B_p=B_p, b_w=b_w, E_w=E_w, e_c=e_c, Lh=Lh)


def simulate(P, N, dt, write_times, scheme="fixed-stress", stab=True,
             beta_factor=1.0, tol=1e-10, max_iter=1000,
             equation_tol=1e-8, initial_state=None, start_time=0.0):
    """Run from the undrained state; return snapshots at write_times.

    Returns dict(z, t, p, w, settlement, iters), where p and w are arrays of
    shape (len(write_times), N) and iters holds fixed-stress iteration counts
    per time step (all 1 for the monolithic scheme).
    """
    validate_inputs(P, N, dt, write_times, beta_factor, tol, max_iter, start_time)
    if scheme not in ("fixed-stress", "monolithic"):
        raise ValueError("Unknown solution scheme")
    if not np.isfinite(equation_tol) or equation_tol <= 0:
        raise ValueError("Invalid equation tolerance")
    O = operators(P, N)
    n = N
    I = np.eye(n)
    h = O["h"]
    tau = P.alpha ** 2 * h ** 2 / (4 * P.Mc) if stab else 0.0
    beta = beta_factor * P.alpha ** 2 / P.Mc
    L = P.mob * O["Lh"]

    mech = lu_factor(O["A_w"])

    def solve_mech(p):
        return lu_solve(mech, O["b_w"] - O["B_p"] @ p)

    def strain(w):
        return O["E_w"] @ w + O["e_c"]

    def settlement(w):
        return w[0] + 0.5 * h * P.sigma0 / P.Mc

    if initial_state is None:
        if start_time != 0:
            raise ValueError("A restart requires pressure and displacement")
        p = np.full(n, P.p0)
        w = solve_mech(p)
    else:
        p, w = (np.array(initial_state[key], dtype=float, copy=True) for key in ("p", "w"))
        if p.shape != (n,) or w.shape != (n,) or not np.all(np.isfinite([p, w])):
            raise ValueError("Invalid restart fields")
    eps = strain(w)
    p_start, eps_start = p.copy(), eps.copy()
    content_scale = (1/P.M + P.alpha**2/P.Mc)*P.p0
    drained = 0.0
    diagnostics = []

    def residuals(p, w, p_n, eps_n):
        mass = ((p-p_n)/P.M + P.alpha*(strain(w)-eps_n)
                - dt*(L@p) - tau*(O["Lh"]@(p-p_n)))
        mech = O["A_w"]@w + O["B_p"]@p - O["b_w"]
        # A_w is integrated across one cell per unit cross-sectional area.
        return np.max(np.abs(mass))/content_scale, np.max(np.abs(mech))/(h*P.sigma0/P.H)

    if scheme == "monolithic":
        A = np.block([[I / (P.M * dt) - L - tau / dt * O["Lh"], P.alpha / dt * O["E_w"]],
                      [O["B_p"], O["A_w"]]])
        mono = lu_factor(A)
    else:
        flow = lu_factor((1 / P.M + beta) / dt * I - L - tau / dt * O["Lh"])

    write_times = sorted(write_times)
    n_steps = int(round((write_times[-1]-start_time) / dt))
    out_p, out_w, out_s, out_t, iters = [], [], [], [], []
    wi = 0
    for step in range(1, n_steps + 1):
        p_n, eps_n = p.copy(), eps.copy()
        stab_src = tau / dt * (O["Lh"] @ p_n)
        if scheme == "monolithic":
            rhs = np.concatenate([p_n / (P.M * dt) + P.alpha / dt * (eps_n - O["e_c"]) - stab_src,
                                  O["b_w"]])
            sol = lu_solve(mono, rhs)
            p, w = sol[:n], sol[n:]
            eps = strain(w)
            iters.append(1)
        else:
            k = 0
            while True:
                k += 1
                p_k = p.copy()
                rhs = (p_n / (P.M * dt) + beta / dt * p_k
                       - P.alpha / dt * (eps - eps_n) - stab_src)
                p = lu_solve(flow, rhs)
                w = solve_mech(p)
                eps = strain(w)
                change = np.max(np.abs(p - p_k)) / P.p0
                rm, rw = residuals(p, w, p_n, eps_n)
                if not np.all(np.isfinite([change, rm, rw])):
                    raise FloatingPointError("Nonfinite coupled solution")
                if change < tol and rm < equation_tol and rw < equation_tol:
                    break
                if k >= max_iter:
                    raise ConvergenceError(f"COUPLING_NOT_CONVERGED at t={start_time+step*dt:g}: "
                                           f"iterations={k}, change={change:g}, mass={rm:g}, mechanics={rw:g}")
            iters.append(k)
        t = start_time + step * dt
        rm, rw = residuals(p, w, p_n, eps_n)
        if not np.all(np.isfinite([p, w])) or not np.all(np.isfinite([rm, rw])):
            raise FloatingPointError("Nonfinite solution")
        if max(rm, rw) >= equation_tol:
            raise ConvergenceError(f"Coupled equation residual exceeds tolerance at t={t:g}")
        # Per unit cross-sectional area; outward Darcy flux at the drained top.
        drained += dt*P.mob*2*p[0]/h/(content_scale*P.H)
        physical = h*np.sum((p-p_start)/P.M + P.alpha*(eps-eps_start))/(content_scale*P.H)
        stab_content = -tau*h*np.sum(O["Lh"]@(p-p_start))/(content_scale*P.H)
        diagnostics.append([t, rm, rw, physical, stab_content, drained, physical+stab_content+drained])
        while wi < len(write_times) and abs(t - write_times[wi]) < 1e-9 * max(1.0, t):
            out_t.append(t); out_p.append(p.copy()); out_w.append(w.copy())
            out_s.append(settlement(w)); wi += 1
    if wi != len(write_times):
        raise RuntimeError("Missing output times")
    return dict(z=O["z"], t=np.array(out_t), p=np.array(out_p), w=np.array(out_w),
                settlement=np.array(out_s), iters=np.array(iters), diagnostics=np.array(diagnostics))
