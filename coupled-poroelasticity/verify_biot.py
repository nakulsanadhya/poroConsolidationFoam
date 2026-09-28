"""
Verification driver for biotConsolidationFoam (Phase 2: coupled 1D Biot).

Every study can run on two interchangeable backends:

    --solver biotConsolidationFoam   the OpenFOAM solver (default; build it first)
    --solver reference               the Python reference implementation of the
                                     identical discretization (reference/biot_ref.py)

Commands
--------
    python3 verify_biot.py check         pass/fail: exact solution, reference agreement,
                                         fixed-stress convergence
    python3 verify_biot.py compare       p, w and settlement vs the exact solution
    python3 verify_biot.py convergence   mesh and time-step refinement (p and w)
    python3 verify_biot.py oscillations  early-time pressure oscillations, with and
                                         without stabilization
    python3 verify_biot.py fixedstress   fixed-stress iteration counts vs time step
                                         and vs the fixed-stress parameter
    python3 verify_biot.py reliability   restart, monolithic agreement and failure tests
    python3 verify_biot.py all           check, reliability, and full studies

With the OpenFOAM backend, `compare` and `check` also run the reference
implementation and report the largest solver-vs-reference difference. Since
both use the same discretization, they should agree to roughly the
fixed-stress tolerance; a larger difference points to an implementation bug
even when both are close to the exact solution.

On foam-extend, run Python without the foam-extend LD_LIBRARY_PATH (see the
top-level README) so numpy imports; OF_BASHRC is then sourced for the
OpenFOAM commands.
"""

import argparse
import json
import os
import platform
import time
import shlex
import hashlib
from datetime import datetime, timezone
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TEMPLATE = HERE / "case"
RUNS = HERE / "runs"
OUT = HERE / "validation"
sys.path.insert(0, str(HERE / "reference"))
sys.path.insert(0, str(ROOT / "pressure-diffusion"))
from biot_ref import Params, exact, exact_settlement, simulate, ConvergenceError, validate_inputs  # noqa: E402
import verify as v1  # noqa: E402  (pressure-diffusion helpers: environment, time dirs, field and mesh readers)

SOLVER = "biotConsolidationFoam"


# --------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------

def _read_prop(text, key):
    m = re.search(rf"^\s*{key}\s+{key}\s+\[[^\]]*\]\s+([^;\s]+)\s*;", text, re.M)
    if not m:
        m = re.search(rf"^\s*{key}\s+([^;\s]+)\s*;", text, re.M)
    if not m:
        raise KeyError(key)
    return m.group(1)


def base_params():
    text = (TEMPLATE / "constant" / "biotProperties").read_text()
    return Params(H=1.0, Mc=float(_read_prop(text, "Mc")), M=float(_read_prop(text, "M")),
                  alpha=float(_read_prop(text, "alpha")),
                  mob=float(_read_prop(text, "mobility")),
                  sigma0=float(_read_prop(text, "sigma0")))


def _set_prop(text, key, value):
    new, n = re.subn(rf"^(\s*{key}\s+{key}\s+\[[^\]]*\]\s+)[^;]+;", rf"\g<1>{value};", text,
                     count=1, flags=re.M)
    if n == 0:
        new, n = re.subn(rf"^(\s*{key}\s+)[^;]+;", rf"\g<1>{value};", text, count=1, flags=re.M)
    if n != 1:
        raise RuntimeError(f"could not set '{key}' in biotProperties")
    return new


# --------------------------------------------------------------------------
# Backends
# --------------------------------------------------------------------------

def _run_openfoam(name, P, N, dt, write_times, stab, beta_factor, tol, max_iter):
    dst = RUNS / name
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(TEMPLATE, dst)

    bmd = dst / "constant" / "polyMesh" / "blockMeshDict"
    text, n = re.subn(r"(hex\s*\([^)]*\)\s*)\(1 1 \d+\)", rf"\g<1>(1 1 {N})", bmd.read_text(), count=1)
    if n != 1:
        raise RuntimeError("could not set cell count in blockMeshDict")
    bmd.write_text(text)
    shutil.copy(bmd, dst / "system" / "blockMeshDict")

    write_interval = min(np.diff([0.0] + list(write_times)))
    cd = dst / "system" / "controlDict"
    text = cd.read_text()
    for key, val in [("application", SOLVER), ("deltaT", f"{dt:.12g}"),
                     ("endTime", f"{max(write_times):.12g}"), ("writeInterval", f"{write_interval:.12g}")]:
        text = v1._set_entry(text, key, val)
    cd.write_text(text)

    bp = dst / "constant" / "biotProperties"
    text = bp.read_text()
    for key, val in [("Mc", f"{P.Mc:.12g}"), ("M", f"{P.M:.12g}"), ("alpha", f"{P.alpha:.12g}"),
                     ("mobility", f"{P.mob:.12g}"), ("sigma0", f"{P.sigma0:.12g}"),
                     ("stabilization", "yes" if stab else "no"),
                     ("fixedStressFactor", f"{beta_factor:.12g}"),
                     ("fixedStressTolerance", f"{tol:.12g}"),
                     ("maxFixedStressIterations", str(int(max_iter)))]:
        text = _set_prop(text, key, val)
    bp.write_text(text)

    prefix = v1._shell_prefix()
    logs = {}
    for cmd in ("blockMesh", SOLVER):
        res = subprocess.run(["bash", "-c", f"{prefix}{shlex.quote(cmd)} -case {shlex.quote(str(dst))}"], capture_output=True, text=True)
        log = res.stdout + res.stderr
        (dst / f"log.{cmd}").write_text(log)
        if "COUPLING_NOT_CONVERGED" in log:
            raise ConvergenceError(f"COUPLING_NOT_CONVERGED; see {dst / f'log.{cmd}'}")
        if "WARNING: fixed-stress split not converged" in log:
            raise RuntimeError("Old solver continued after failure; rebuild the reliability revision")
        if res.returncode != 0 or "FOAM FATAL" in log:
            raise RuntimeError(f"{cmd} failed, see {dst / f'log.{cmd}'}")
        logs[cmd] = log

    if "Reliability checks v1:" not in logs[SOLVER]:
        raise RuntimeError("Rebuild biotConsolidationFoam: reliability diagnostics missing")
    z = v1.cell_centres_z(dst)
    if not np.allclose(z, (np.arange(N)+.5)*P.H/N, rtol=0, atol=1e-9*P.H):
        raise RuntimeError("Unexpected mesh geometry or cell ordering")
    ps, ws, ss = [], [], []
    for t in write_times:
        ps.append(v1.read_field(dst, t, "p"))
        ws.append(v1.read_field(dst, t, "w"))
        ss.append(_read_patch_value(v1.time_dir(dst, t) / "w", "drained"))
    iters = np.array([int(k) for k in re.findall(r"Fixed-stress iterations = (\d+)", logs[SOLVER])])
    return dict(z=z, t=np.array(write_times, dtype=float), p=np.array(ps), w=np.array(ws),
                settlement=np.array(ss), iters=iters, diagnostics=parse_diagnostics(logs[SOLVER]))


def _read_patch_value(path, patch):
    text = path.read_text()
    block = re.search(rf"\b{patch}\s*\{{(.*?)\}}", text, re.S)
    if not block:
        raise ValueError(f"patch {patch} not found in {path}")
    m = re.search(r"value\s+uniform\s+([-+0-9.eE]+)\s*;", block.group(1))
    if not m:
        m = re.search(r"value\s+nonuniform\s+List<scalar>\s*1\s*\(\s*([-+0-9.eE]+)\s*\)", block.group(1))
    if not m:
        raise ValueError(f"Missing scalar patch value in {path}")
    return float(m.group(1))


def parse_diagnostics(log):
    rows = re.findall(r"^DIAGNOSTICS (.+)$", log, re.M)
    values = np.array([[float(x) for x in row.split()] for row in rows])
    if values.ndim != 2 or values.shape[1] != 7 or not np.all(np.isfinite(values)):
        raise RuntimeError("Missing, malformed or nonfinite per-step diagnostics")
    return values


def validate_result(r, N, dt, times, start_time=0.0):
    steps = int(round((times[-1]-start_time)/dt))
    shapes = {"z": (N,), "t": (len(times),), "p": (len(times),N),
              "w": (len(times),N), "settlement": (len(times),),
              "iters": (steps,), "diagnostics": (steps,7)}
    for key, shape in shapes.items():
        if np.shape(r[key]) != shape or not np.all(np.isfinite(r[key])):
            raise RuntimeError(f"Invalid/missing/nonfinite {key}: expected {shape}")
    if not np.allclose(r["t"], times, rtol=0, atol=1e-9):
        raise RuntimeError("Missing requested output times")
    expected = start_time+dt*np.arange(1,steps+1)
    if not np.allclose(r["diagnostics"][:,0], expected, rtol=1e-10, atol=1e-10):
        raise RuntimeError("Missing time-step diagnostics")
    if np.min(r["iters"]) < 1 or np.any(r["diagnostics"][:,1:3] >= 1e-8):
        raise RuntimeError("Coupled equation residual gate failed")
    if np.max(np.abs(r["diagnostics"][:,6])) > 1e-6:
        raise RuntimeError("Cumulative fluid balance gate failed")
    return r


def run(backend, name, P, N, dt, write_times, stab=True, beta_factor=1.0, tol=1e-10, max_iter=1000):
    validate_inputs(P, N, dt, write_times, beta_factor, tol, max_iter)
    if backend == "reference":
        r = simulate(P, N, dt, write_times, "fixed-stress", stab, beta_factor, tol, max_iter)
    else:
        if P.H != 1.0:
            raise ValueError("OpenFOAM verification template is a 1 m column")
        r = _run_openfoam(name, P, N, dt, write_times, stab, beta_factor, tol, max_iter)
    return validate_result(r, N, dt, write_times)


def scan_run(*args, **kwargs):
    # Only deliberately probed coupling nonconvergence is captured. All other
    # errors (compile, bad mesh, nonfinite fields, missing files) still abort.
    try:
        return run(*args, **kwargs)["iters"].astype(float)
    except ConvergenceError as exc:
        print(f"NC: {exc}")
        return np.array([np.nan])


def iteration_text(values):
    return "NC" if np.isnan(values).any() else f"{values.mean():.1f}"


def _l2(a):
    return float(np.sqrt(np.mean(np.asarray(a) ** 2)))


# --------------------------------------------------------------------------
# compare
# --------------------------------------------------------------------------

TABLE_TIMES = [5, 10, 20, 40, 80]


def compare(backend, quiet=False):
    P = base_params()
    N, dt = 100, 0.25
    times = [5.0 * k for k in range(1, 17)]
    r = run(backend, "compare", P, N, dt, times)
    z = r["z"]

    rows = []
    for t in TABLE_TIMES:
        j = int(np.argmin(abs(r["t"] - t)))
        pe, we = exact(P, z, t)
        rows.append(dict(Tv=P.cv * t / P.H ** 2,
                         ep=np.max(abs(r["p"][j] - pe)) / P.p0,
                         ew=np.max(abs(r["w"][j] - we)) / P.s_inf,
                         es=(r["settlement"][j] - exact_settlement(P, t)) / P.s_inf))

    agreement = None
    if backend != "reference":
        ref = simulate(P, N, dt, times, "fixed-stress", True, 1.0, 1e-10, 1000)
        agreement = (np.max(abs(r["p"] - ref["p"])) / P.p0, np.max(abs(r["w"] - ref["w"])) / P.s_inf)

    # Figure: pressure and displacement profiles, settlement history
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    zf = np.linspace(0, P.H, 300)
    for t in TABLE_TIMES:
        j = int(np.argmin(abs(r["t"] - t)))
        pe, we = exact(P, zf, t)
        line, = axes[0].plot(pe / P.sigma0, zf, lw=1.4, label=f"Tv = {P.cv * t:g}")
        axes[0].plot(r["p"][j][::5] / P.sigma0, z[::5], "o", ms=3.5, color=line.get_color())
        axes[1].plot(we / P.s_inf, zf, lw=1.4, color=line.get_color(), label=f"Tv = {P.cv * t:g}")
        axes[1].plot(r["w"][j][::5] / P.s_inf, z[::5], "o", ms=3.5, color=line.get_color())
    for ax, xl, title in [(axes[0], r"$p/\sigma_0$", "Pore pressure"),
                          (axes[1], r"$w/s_\infty$", "Vertical displacement")]:
        ax.invert_yaxis(); ax.set_xlabel(xl); ax.set_ylabel(r"$z/H$"); ax.set_title(title)
        ax.grid(alpha=0.3); ax.legend(fontsize=8)
    tf = np.linspace(0.05, 80, 400)
    se = np.array([exact_settlement(P, t) for t in tf])
    ax = axes[2]
    ax.plot(np.concatenate([[0], tf]), np.concatenate([[P.s0], se]) * 1e3, "k-", lw=1.4, label="exact")
    ax.plot(np.concatenate([[0], r["t"]]), np.concatenate([[P.s0], r["settlement"]]) * 1e3, "o",
            color="C3", ms=4, label="numerical")
    ax.axhline(P.s0 * 1e3, color="gray", ls=":", lw=1)
    ax.axhline(P.s_inf * 1e3, color="gray", ls="--", lw=1)
    ax.text(79, P.s0 * 1e3, "undrained $s_0$", ha="right", va="bottom", fontsize=8, color="gray")
    ax.text(79, P.s_inf * 1e3, r"drained $s_\infty$", ha="right", va="top", fontsize=8, color="gray")
    ax.set_xlabel("t [s]"); ax.set_ylabel("settlement [mm]"); ax.set_title("Settlement of the loaded top")
    ax.grid(alpha=0.3); ax.legend(fontsize=8, loc="center right")
    fig.suptitle(f"{_label(backend)}: markers = numerical, lines = exact "
                 f"({N} cells, $\\Delta t$ = {dt:g} s)", fontsize=11)
    plt.tight_layout()
    out = OUT / "biot_vs_exact.png"
    plt.savefig(out, dpi=140); plt.close(fig)

    if not quiet:
        print(f"\n{_label(backend)}, {N} cells, dt = {dt:g} s   "
              f"(p0 = {P.p0 / P.sigma0:.4f} sigma0, s0 = {P.s0 * 1e3:.3f} mm, s_inf = {P.s_inf * 1e3:.3f} mm)")
        print(f"{'Tv':>6} {'Linf(p)/p0':>11} {'Linf(w)/s_inf':>14} {'settlement err/s_inf':>21}")
        for row in rows:
            print(f"{row['Tv']:>6g} {row['ep']:>11.2e} {row['ew']:>14.2e} {row['es']:>+21.2e}")
        print(f"Fixed-stress iterations per step: mean {r['iters'].mean():.2f}, max {r['iters'].max()}")
        if agreement:
            print(f"Solver vs reference implementation: max |dp|/p0 = {agreement[0]:.2e}, "
                  f"max |dw|/s_inf = {agreement[1]:.2e}")
        print(f"Plot: {out}")
    return rows, r, agreement


# --------------------------------------------------------------------------
# convergence
# --------------------------------------------------------------------------

def convergence(backend):
    P = base_params()
    t_end = 20.0
    Tv = P.cv * t_end / P.H ** 2

    # Mesh study: implicit Euler's first-order time error would swamp the
    # second-order space error on fine meshes, so it is removed by Richardson
    # extrapolation in time, u = 2 u(dt/2) - u(dt), leaving an O(dt^2) remainder.
    Ns, dt_c = [10, 20, 40, 80, 160], 2e-3
    ep_s, ew_s = [], []
    for N in Ns:
        a = run(backend, f"space_N{N}_dt{dt_c:g}", P, N, dt_c, [t_end])
        b = run(backend, f"space_N{N}_dt{dt_c / 2:g}", P, N, dt_c / 2, [t_end])
        p_x = 2 * b["p"][-1] - a["p"][-1]
        w_x = 2 * b["w"][-1] - a["w"][-1]
        pe, we = exact(P, a["z"], t_end)
        ep_s.append(_l2(p_x - pe) / P.p0); ew_s.append(_l2(w_x - we) / P.s_inf)

    dts, N_fine = [2.0, 1.0, 0.5, 0.25, 0.125], 800
    ep_t, ew_t = [], []
    for dt in dts:
        r = run(backend, f"time_dt{dt:g}", P, N_fine, dt, [t_end])
        pe, we = exact(P, r["z"], t_end)
        ep_t.append(_l2(r["p"][-1] - pe) / P.p0); ew_t.append(_l2(r["w"][-1] - we) / P.s_inf)

    h = np.array([P.H / N for N in Ns])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    ax1.loglog(h, ep_s, "o-", label=r"pressure, $\|e_p\|_2/p_0$")
    ax1.loglog(h, ew_s, "s-", label=r"displacement, $\|e_w\|_2/s_\infty$")
    ax1.loglog(h, ep_s[0] * (h / h[0]) ** 2, "k--", lw=0.8, label="slope 2")
    ax1.set_xlabel(r"$\Delta z$ [m]"); ax1.set_ylabel("relative L2 error")
    ax1.set_title(f"Mesh refinement (Richardson in time, $\\Delta t$ = {dt_c:g}, {dt_c / 2:g} s; Tv = {Tv:g})",
                  fontsize=10)
    ax1.grid(alpha=0.3, which="both"); ax1.legend(fontsize=8)
    d = np.array(dts)
    ax2.loglog(d, ep_t, "o-", label=r"pressure")
    ax2.loglog(d, ew_t, "s-", label=r"displacement")
    ax2.loglog(d, ep_t[0] * d / d[0], "k--", lw=0.8, label="slope 1")
    ax2.set_xlabel(r"$\Delta t$ [s]")
    ax2.set_title(f"Time-step refinement ({N_fine} cells, Tv = {Tv:g})")
    ax2.grid(alpha=0.3, which="both"); ax2.legend(fontsize=8)
    fig.suptitle(_label(backend), fontsize=11)
    plt.tight_layout()
    out = OUT / "biot_convergence.png"
    plt.savefig(out, dpi=140); plt.close(fig)

    def orders(e):
        return [None] + [np.log2(e[i] / e[i + 1]) for i in range(len(e) - 1)]

    print(f"\nMesh refinement ({_label(backend)}, Richardson in time from dt = {dt_c:g} and {dt_c / 2:g} s, Tv = {Tv:g})")
    print(f"{'cells':>6} {'L2(p)/p0':>10} {'order':>6} {'L2(w)/s_inf':>12} {'order':>6}")
    for N, a, pa, b, pb in zip(Ns, ep_s, orders(ep_s), ew_s, orders(ew_s)):
        print(f"{N:>6} {a:>10.3e} {'' if pa is None else f'{pa:.2f}':>6} {b:>12.3e} "
              f"{'' if pb is None else f'{pb:.2f}':>6}")
    print(f"\nTime-step refinement ({_label(backend)}, {N_fine} cells, Tv = {Tv:g})")
    print(f"{'dt [s]':>6} {'L2(p)/p0':>10} {'order':>6} {'L2(w)/s_inf':>12} {'order':>6}")
    for dt, a, pa, b, pb in zip(dts, ep_t, orders(ep_t), ew_t, orders(ew_t)):
        print(f"{dt:>6g} {a:>10.3e} {'' if pa is None else f'{pa:.2f}':>6} {b:>12.3e} "
              f"{'' if pb is None else f'{pb:.2f}':>6}")
    print(f"Plot: {out}")
    for errors, low, high in [(ep_s,1.8,2.2),(ew_s,1.8,2.2),(ep_t,.85,1.15),(ew_t,.85,1.15)]:
        rates = np.log2(np.array(errors[:-1])/np.array(errors[1:]))
        if not np.all(np.isfinite(rates)) or np.any((rates < low)|(rates > high)):
            raise RuntimeError(f"Refinement-order gate failed: {rates}")
    print("REFINEMENT GATES PASSED")


# --------------------------------------------------------------------------
# oscillations
# --------------------------------------------------------------------------

def _tv_excess(p):
    return float(np.sum(np.abs(np.diff(p))) - (p.max() - p.min()))


def oscillations(backend):
    """Nearly incompressible constituents, one time step after loading."""
    P = Params(Mc=1e8, M=1e13, alpha=1.0, sigma0=1e6).with_cv(0.01)
    N = 50
    h = P.H / N
    dts = [1e-2, 1e-3, 1e-4]
    res = {}
    for dt in dts:
        for stab in (False, True):
            r = run(backend, f"osc_dt{dt:g}_{'stab' if stab else 'nostab'}", P, N, dt, [dt],
                    stab=stab, tol=1e-10, max_iter=5000)
            res[(dt, stab)] = r

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6))
    dt_show = 1e-4
    zf = np.linspace(0, 0.3 * P.H, 400)
    ax1.plot(exact(P, zf, dt_show)[0] / P.p0, zf, "k-", lw=1.3, label="exact")
    for stab, mk, col, lab in [(False, "o-", "C3", "no stabilization"), (True, "s-", "C0", "stabilized")]:
        r = res[(dt_show, stab)]
        sel = r["z"] <= 0.3 * P.H
        ax1.plot(r["p"][0][sel] / P.p0, r["z"][sel], mk, color=col, ms=4, lw=0.8, label=lab)
    ax1.axvline(1.0, color="gray", ls=":", lw=1)
    ax1.invert_yaxis(); ax1.set_xlabel(r"$p/p_0$"); ax1.set_ylabel(r"$z/H$")
    ax1.set_title(f"One step after loading, $c_v\\Delta t/\\Delta z^2$ = {P.cv * dt_show / h ** 2:.1e}")
    ax1.grid(alpha=0.3); ax1.legend(fontsize=8)

    x = [P.cv * dt / h ** 2 for dt in dts]
    for stab, mk, col, lab in [(False, "o-", "C3", "no stabilization"), (True, "s-", "C0", "stabilized")]:
        y = [max(_tv_excess(res[(dt, stab)]["p"][0]) / P.p0, 1e-16) for dt in dts]
        ax2.loglog(x, y, mk, color=col, label=lab)
    ax2.set_xlabel(r"$c_v\Delta t/\Delta z^2$"); ax2.set_ylabel(r"total-variation excess $/p_0$")
    ax2.set_title("Size of the non-physical oscillation")
    ax2.grid(alpha=0.3, which="both"); ax2.legend(fontsize=8)
    fig.suptitle(f"{_label(backend)}: nearly incompressible constituents "
                 f"($\\alpha$ = 1, M = {P.M:.0e} Pa), {N} cells", fontsize=11)
    plt.tight_layout()
    out = OUT / "biot_oscillations.png"
    plt.savefig(out, dpi=140); plt.close(fig)

    print(f"\nEarly-time oscillations ({_label(backend)}, alpha = 1, M = {P.M:.0e} Pa, {N} cells, one step)")
    print(f"{'c dt/dz^2':>10} {'stab':>5} {'max(p)/p0 - 1':>14} {'TV excess/p0':>13} {'FS iters':>9}")
    for dt in dts:
        for stab in (False, True):
            r = res[(dt, stab)]
            p = r["p"][0]
            print(f"{P.cv * dt / h ** 2:>10.1e} {('yes' if stab else 'no'):>5} "
                  f"{p.max() / P.p0 - 1:>+14.2e} {_tv_excess(p) / P.p0:>13.2e} {r['iters'].max():>9}")
    print(f"Plot: {out}")
    for dt in dts:
        p = res[(dt,True)]["p"][0]/P.p0
        if p.max()-1 > 1e-7 or p.min() < -1e-7 or _tv_excess(p) > 1e-7:
            raise RuntimeError("Stabilized oscillation gate failed")


# --------------------------------------------------------------------------
# fixed-stress convergence
# --------------------------------------------------------------------------

def fixedstress(backend):
    P = Params(Mc=1e8, M=1e13, alpha=1.0, sigma0=1e6).with_cv(0.01)
    N = 50
    h = P.H / N
    tol, max_iter = 1e-8, 3000
    dts = [1.0, 1e-1, 1e-2, 1e-3, 1e-4]
    it_dt = {}
    for dt in dts:
        for stab in (False, True):
            it_dt[(dt, stab)] = scan_run(backend, f"fs_dt{dt:g}_{'stab' if stab else 'nostab'}", P, N, dt, [3 * dt],
                    stab=stab, tol=tol, max_iter=max_iter)

    factors = [0.5, 0.6, 0.75, 1.0, 1.5, 2.0]
    dts_f = [1.0, 1e-2, 1e-4]
    it_f = {}
    for dt in dts_f:
        for f in factors:
            it_f[(dt, f)] = scan_run(backend, f"fs_beta{f:g}_dt{dt:g}", P, N, dt, [3 * dt], stab=True,
                    beta_factor=f, tol=tol, max_iter=max_iter)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6))
    x = [P.cv * dt / h ** 2 for dt in dts]
    for stab, mk, col, lab in [(False, "o-", "C3", "no stabilization"), (True, "s-", "C0", "stabilized")]:
        ax1.loglog(x, [it_dt[(dt, stab)].mean() for dt in dts], mk, color=col, label=lab)
    ax1.axhline(max_iter, color="gray", ls=":", lw=1)
    ax1.text(x[0], max_iter * 0.8, "iteration cap", ha="right", fontsize=8, color="gray")
    ax1.set_xlabel(r"$c_v\Delta t/\Delta z^2$"); ax1.set_ylabel("fixed-stress iterations per step")
    ax1.set_title(r"Iterations vs time step ($\beta = \alpha^2/M_c$)")
    ax1.grid(alpha=0.3, which="both"); ax1.legend(fontsize=8)
    for dt, col in zip(dts_f, ["C0", "C1", "C2"]):
        ax2.semilogy(factors, [it_f[(dt, f)].mean() for f in factors], "o-", color=col,
                     label=f"$c_v\\Delta t/\\Delta z^2$ = {P.cv * dt / h ** 2:.1e}")
    for dt, col in zip(dts_f, ["C0", "C1", "C2"]):
        failed = [f for f in factors if np.isnan(it_f[(dt,f)]).any()]
        if failed:
            ax2.plot(failed, [max_iter]*len(failed), "x", color=col, ms=9)
    ax2.text(.98,.98,"x = NC (iteration limit)\nNC fields are not used", transform=ax2.transAxes,
             ha="right", va="top", fontsize=8)
    ax2.axvline(0.5, color="gray", ls=":", lw=1)
    y_top = np.nanmax([it_f[k].mean() for k in it_f])
    ax2.text(0.52, 0.6 * y_top, r"$\beta = \alpha^2/(2M_c)$", fontsize=8, color="gray")
    ax2.set_xlabel(r"$\beta\,/\,(\alpha^2/M_c)$"); ax2.set_ylabel("fixed-stress iterations per step")
    ax2.set_title("Iterations vs fixed-stress parameter (stabilized)")
    ax2.grid(alpha=0.3, which="both"); ax2.legend(fontsize=8)
    fig.suptitle(f"{_label(backend)}: $\\alpha$ = 1, M = {P.M:.0e} Pa, {N} cells, tolerance {tol:g}",
                 fontsize=11)
    plt.tight_layout()
    out = OUT / "biot_fixed_stress.png"
    plt.savefig(out, dpi=140); plt.close(fig)

    print(f"\nFixed-stress iterations per step ({_label(backend)}, alpha = 1, M = {P.M:.0e} Pa, "
          f"tol {tol:g}, cap {max_iter})")
    print(f"{'c dt/dz^2':>10} {'no stabilization':>17} {'stabilized':>11}")
    for dt in dts:
        print(f"{P.cv * dt / h ** 2:>10.1e} {iteration_text(it_dt[(dt, False)]):>17} {iteration_text(it_dt[(dt, True)]):>11}")
    print(f"\nStabilized, beta = f * alpha^2/Mc")
    print(f"{'c dt/dz^2':>10} " + " ".join(f"{'f=' + format(f, 'g'):>7}" for f in factors))
    for dt in dts_f:
        print(f"{P.cv * dt / h ** 2:>10.1e} " + " ".join(f"{iteration_text(it_f[(dt, f)]):>7}" for f in factors))
    print(f"Plot: {out}")
    if any(np.isnan(v).any() for v in it_dt.values()) or any(np.isnan(it_f[(dt,1.0)]).any() for dt in dts_f):
        raise RuntimeError("Default fixed-stress factor failed in the tested scan")
    print("Parameter scan complete; NC entries are failed solves, not verified solutions.")


# --------------------------------------------------------------------------
# check
# --------------------------------------------------------------------------

def check(backend, tol_p=5e-3, tol_w=2e-3, tol_ref=1e-6):
    rows, r, agreement = compare(backend)
    ok = True
    for row in rows:
        pressure_limit = 8e-3 if row["Tv"] < .1*(1-1e-9) else tol_p
        if row["ep"] > pressure_limit or row["ew"] > tol_w or abs(row["es"]) > tol_w:
            print(f"FAIL at Tv = {row['Tv']:g}: Linf(p)/p0 = {row['ep']:.2e} (tol {pressure_limit:g}), "
                  f"Linf(w)/s_inf = {row['ew']:.2e}, settlement err = {row['es']:+.2e} (tol {tol_w:g})")
            ok = False
    if agreement and max(agreement) > tol_ref:
        print(f"FAIL: solver and reference implementation differ by {max(agreement):.2e} (tol {tol_ref:g})")
        ok = False
    d = r["diagnostics"]
    balance = np.max(np.abs(d[:,6]))
    print(f"Max coupled residuals: mass={d[:,1].max():.3e}, mechanics={d[:,2].max():.3e}")
    print(f"Final normalized content: physical={d[-1,3]:+.8e}, stabilization={d[-1,4]:+.8e}, "
          f"outward drainage={d[-1,5]:+.8e}; max balance error={balance:.3e}")
    np.savetxt(OUT/"balance.csv", d, delimiter=",",
               header="time,mass_residual,mechanics_residual,physical_content,stabilization_content,drainage,balance", comments="")
    if balance > 1e-6:
        print("FAIL: cumulative fluid balance exceeds 1e-6 of initial constrained content")
        ok = False
    print("\nCHECK PASSED" if ok else "\nCHECK FAILED")
    if not ok:
        raise RuntimeError("Verification check failed")


def _invoke_case(dst, command, logfile):
    res = subprocess.run(["bash", "-c", f"{v1._shell_prefix()}{shlex.quote(command)} -case {shlex.quote(str(dst))}"],
                         capture_output=True, text=True)
    text = res.stdout+res.stderr
    (dst/logfile).write_text(text)
    return res.returncode, text


def reliability(backend):
    """Small integration gates, including deliberate failures and a real restart."""
    P, N, dt = base_params(), 20, .25
    times = [.5,1.0,1.5,2.0]
    full = run(backend,"reliability_full",P,N,dt,times)
    mono = simulate(P,N,dt,times,scheme="monolithic")
    dp = np.max(np.abs(full["p"]-mono["p"]))/P.p0
    dw = np.max(np.abs(full["w"]-mono["w"]))/P.s_inf
    if max(dp,dw) > 1e-6:
        raise RuntimeError("Fixed-stress vs monolithic gate failed")
    print(f"PASS fixed-stress vs monolithic: dp/p0={dp:.3e}, dw/s_inf={dw:.3e}")

    first = run(backend,"reliability_restart",P,N,dt,[.5,1.0])
    if backend == "reference":
        resumed = simulate(P,N,dt,[1.5,2.0],initial_state={"p":first["p"][-1],"w":first["w"][-1]},start_time=1.0)
        validate_result(resumed,N,dt,[1.5,2.0],start_time=1.0)
    else:
        dst = RUNS/"reliability_restart"
        cd = dst/"system/controlDict"
        text = v1._set_entry(cd.read_text(),"startFrom","latestTime")
        cd.write_text(v1._set_entry(text,"endTime","2"))
        code, log = _invoke_case(dst,SOLVER,"log.restart")
        if code != 0 or "FOAM FATAL" in log or "Undrained initial pressure" in log:
            raise RuntimeError(f"Restart failed; see {dst/'log.restart'}")
        resumed = dict(z=v1.cell_centres_z(dst),t=np.array([1.5,2.0]),
                       p=np.array([v1.read_field(dst,t,"p") for t in [1.5,2.0]]),
                       w=np.array([v1.read_field(dst,t,"w") for t in [1.5,2.0]]),
                       settlement=np.array([_read_patch_value(v1.time_dir(dst,t)/"w","drained") for t in [1.5,2.0]]),
                       iters=np.array([int(k) for k in re.findall(r"Fixed-stress iterations = (\d+)",log)]),
                       diagnostics=parse_diagnostics(log))
        validate_result(resumed,N,dt,[1.5,2.0],start_time=1.0)
    dp = np.max(np.abs(resumed["p"]-full["p"][-2:]))/P.p0
    dw = np.max(np.abs(resumed["w"]-full["w"][-2:]))/P.s_inf
    ds = np.max(np.abs(resumed["settlement"]-full["settlement"][-2:]))/P.s_inf
    if max(dp,dw,ds) > 1e-8:
        raise RuntimeError(f"Restart equivalence failed: {dp}, {dw}, {ds}")
    print(f"PASS restart equivalence: dp/p0={dp:.3e}, dw/s_inf={dw:.3e}, ds/s_inf={ds:.3e}")

    try:
        run(backend,"expected_nonconvergence",P,N,dt,[dt],max_iter=1)
    except ConvergenceError:
        print("PASS iteration-cap failure is explicit (expected failure)")
    else:
        raise RuntimeError("Iteration cap was not rejected")

    # Validate guards with corrupted in-memory results, not only valid data.
    for mutation in ("nan","missing_time","missing_step","balance"):
        bad = {k:v.copy() for k,v in full.items()}
        if mutation == "nan": bad["p"][0,0] = np.nan
        if mutation == "missing_time": bad["t"] = bad["t"][:-1]
        if mutation == "missing_step": bad["diagnostics"] = bad["diagnostics"][:-1]
        if mutation == "balance": bad["diagnostics"][-1,6] = 1e-3
        try:
            validate_result(bad,N,dt,times)
        except RuntimeError:
            print(f"PASS rejects {mutation}")
        else:
            raise RuntimeError(f"Guard failed for {mutation}")

    if backend != "reference":
        # Test the C++ geometry gate, not just Python's template assumptions.
        dst = RUNS/"expected_bad_mesh"
        if dst.exists(): shutil.rmtree(dst)
        shutil.copytree(RUNS/"reliability_full",dst)
        for name in ("constant/polyMesh/blockMeshDict","system/blockMeshDict"):
            path = dst/name
            text,n = re.subn(r"simpleGrading\s*\(1 1 1\)","simpleGrading (1 1 2)",path.read_text())
            if n != 1: raise RuntimeError("Could not make graded-mesh negative test")
            path.write_text(text)
        code,log = _invoke_case(dst,"blockMesh","log.blockMesh")
        if code: raise RuntimeError("Negative-test mesh could not be generated")
        code,log = _invoke_case(dst,SOLVER,"log.rejected_mesh")
        if not code or not any(x in log for x in ("Nonuniform/non-column", "Nonorthogonal/nonuniform")):
            raise RuntimeError("C++ graded-mesh rejection failed")
        print("PASS C++ rejects graded mesh outside stabilization assumptions")
    else:
        print("C++ graded-mesh/restart I/O gates require --solver biotConsolidationFoam")
    print("RELIABILITY CHECKS PASSED")


def _label(backend):
    return "reference implementation" if backend == "reference" else backend


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["check", "reliability", "compare", "convergence", "oscillations", "fixedstress", "all"])
    ap.add_argument("--solver", default=SOLVER, choices=[SOLVER,"reference"], help=f"{SOLVER} (default) or reference")
    args = ap.parse_args()
    OUT = OUT / args.solver
    OUT.mkdir(parents=True, exist_ok=True)
    RUNS = RUNS / args.solver
    RUNS.mkdir(parents=True, exist_ok=True)
    manifest = {"timestamp_utc": datetime.now(timezone.utc).isoformat(), "argv": sys.argv,
                "backend": args.solver, "python": sys.version, "numpy": np.__version__,
                "platform": platform.platform(), "WM_PROJECT_VERSION": os.getenv("WM_PROJECT_VERSION"),
                "OF_BASHRC": os.getenv("OF_BASHRC"), "status": "started",
                "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in sorted(HERE.rglob("*")) if p.is_file() and p.suffix in (".py",".C",".H")
                    and not any(x in p.parts for x in ("runs","__pycache__"))}}
    executable = shutil.which(SOLVER) if args.solver != "reference" else None
    manifest["executable"] = executable
    if executable:
        manifest["executable_sha256"] = hashlib.sha256(Path(executable).read_bytes()).hexdigest()
    try:
        manifest["git_commit"] = subprocess.check_output(["git","-C",str(ROOT),"rev-parse","HEAD"],stderr=subprocess.DEVNULL,text=True).strip()
    except (subprocess.CalledProcessError,FileNotFoundError):
        manifest["git_commit"] = None
    manifest_path = OUT / (args.command+"-manifest.json")
    manifest_path.write_text(json.dumps(manifest,indent=2)+"\n")
    started = time.monotonic()
    try:
        if args.command == "all":
            for f in (check, reliability, convergence, oscillations, fixedstress):
                f(args.solver)
        else:
            {"check": check, "reliability": reliability, "compare": compare, "convergence": convergence,
             "oscillations": oscillations, "fixedstress": fixedstress}[args.command](args.solver)
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["error"] = str(exc)
        raise
    else:
        manifest["status"] = "completed"
    finally:
        manifest["elapsed_seconds"] = time.monotonic()-started
        manifest_path.write_text(json.dumps(manifest,indent=2)+"\n")
