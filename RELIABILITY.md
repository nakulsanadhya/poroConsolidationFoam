# Phase 2 reliability revision — 2026-09-27

## Evidence and status

The baseline C++ solver compiled and ran on Nakul's foam-extend 4.1 workstation.
Its supplied `biot/phase2-all.log` and four figures directly in `biot/validation/`
are preserved unchanged. They verify the baseline, not this revised C++ source.

This revision's Python `check`, `reliability` and full `all` study passed locally.
See `biot/validation/reference/reliability-all.log`. The revised C++ solver subsequently compiled and passed check, reliability,
and all studies on Nakul's foam-extend 4.1 workstation on 2026-09-27.
See `biot/reliability-all.log` and `biot/validation/biotConsolidationFoam/`.

| Evidence | Result |
| --- | --- |
| Baseline C++ / reference normalized differences | pressure 1.66e-10; displacement 6.03e-11 |
| Revised Python default coupled residual maxima | mass 3.99e-12; mechanics 6.47e-12 |
| Revised Python maximum cumulative normalized fluid imbalance | 1.56e-13 |
| Revised Python split / monolithic difference (small integration test) | pressure 7.60e-12; displacement 1.06e-13 |
| Revised Python restart / continuous run | identical output arrays in the tested case |
| Revised Python refinement studies | approximately order 2 in space and 1 in time |
| Failure injection | iteration limit, NaN, missing times, missing diagnostics and corrupted balance rejected |
| Revised C++ build, restart I/O and mesh rejection | Passed on foam-extend 4.1 |
| Revised C++ default coupled residual maxima | mass 5.890e-12; mechanics 1.720e-10 |
| Revised C++ maximum cumulative normalized fluid imbalance | 2.328e-10 |
| Revised C++ refinement studies | approximately order 2 in space and 1 in time |

## What changed

- A coupling step must satisfy both iterate-change and coupled-equation residual
  criteria. Nonconvergence is a fatal error in C++ and an exception in Python.
- The fixed-stress parameter scan catches only the designated nonconvergence
  exception, records `NC`, and does not use failed fields. Other errors still
  abort. Default-factor scan failures fail the study. `all` completing means the
  specified gates passed; it does not mean every exploratory parameter succeeded.
- Per-step diagnostics separately report physical content, stabilization content,
  drainage, and balance. Cumulative accounting starts afresh at each invocation;
  restarted segments are checked from their own starting state.
- Undrained initialization occurs only at time zero. On a nonzero-time restart,
  pressure and displacement are read, retained, and strain is recomputed from w.
- The C++ solver rejects unsupported geometry and boundary types: serial, >=2
  uniform cells in a constant-area orthogonal z column, empty sides, one zero-
  pressure loaded face in -z, and one impermeable fixed face in +z. The original
  default discretization schemes are enforced; arbitrary meshes/scheme overrides
  and parallel execution are intentionally unsupported by this benchmark solver.
- Verification rejects nonfinite values, incomplete time histories and malformed
  fields. The earliest analytical snapshot now has an explicit 0.8% pressure gate.
- Each backend has its own result directory. JSON manifests record command,
  environment/version information, source hashes, executable hash when available,
  and completion/failure status. Existing manifests reflect source at run time.
- `tools/run-python` replaces the personal `py` alias, preserves how to source
  OpenFOAM for child commands, and keeps its libraries out of Python's process.

## Residuals and conservation definition

Let S = 1/M + alpha^2/Mc and C = S*p0. The mass residual is the infinity norm of

```
(p - p_old)/M + alpha*(eps - eps_old)
    - dt*div(mobility*grad(p))
    - div(tau*grad(p)) + div(tau*grad(p_old))
```

divided by C. The iteration-only beta term is absent. This is a normalized
per-step fluid-content defect, not a linear solver's residual. Mechanics uses
`Laplacian(Mc,w) - alpha*grad(p)_z`, normalized by `sigma0/H`.
The default gate is 1e-8 for each; the pressure-change gate is separate.

For a segment starting at state (p_start, eps_start), integrated quantities are:

```
physical = integral[(p-p_start)/M + alpha*(eps-eps_start)] dV
stabilization = -integral div[tau*grad(p-p_start)] dV
outward drainage = sum_steps dt * integral_boundary[-mobility*grad(p).n] dA
balance = physical + stabilization + outward drainage
```

All three are divided by C times total volume. This accounts for the modified
mass equation: the stabilization contribution is a numerical contribution and
must not be hidden inside a claim of conservation of physical storage alone.
For the revised default reference case at t=80 s, these are approximately
-0.883908114, -0.00280309441 and +0.886711208, respectively.
`balance.csv` contains per-step values for the `check` case.

## Commands (from repository root)

Requires foam-extend 4.1 for C++, plus Python 3 with numpy, scipy and matplotlib
(`requirements.txt`). Source the existing foam-extend environment first. Set
`OF_BASHRC` explicitly only if it differs from `$WM_PROJECT_DIR/etc/bashrc`.

```bash
set -o pipefail
(cd biot/solver && wmake) 2>&1 | tee biot/reliability-build.log
# Stop here if compilation fails.
bash tools/run-python biot/verify_biot.py check 2>&1 | tee biot/reliability-check.log
bash tools/run-python biot/verify_biot.py reliability 2>&1 | tee biot/reliability-tests.log
# Run after those quick gates pass:
bash tools/run-python biot/verify_biot.py all 2>&1 | tee biot/reliability-all.log
```

Without OpenFOAM:

```bash
bash tools/run-python biot/verify_biot.py check --solver reference
bash tools/run-python biot/verify_biot.py reliability --solver reference
bash tools/run-python biot/verify_biot.py all --solver reference
```

Quick Python gates take seconds in the available environment. The full study
includes 150,000 small time steps in the mesh study plus the time/parameter
studies; allow several minutes or longer, especially with extra C++ diagnostics.
This is an estimate, not a measured revised-C++ runtime. The revised C++ computes
extra residual fields and therefore can be slower than the baseline.

## Acceptance gates

| Check | Required |
| --- | --- |
| Analytical pressure error / p0 | <= 8e-3 at Tv=0.05; <= 5e-3 at the remaining tabulated times |
| Analytical displacement / settlement error / s_inf | <= 2e-3 |
| C++ / Python field differences, normalized | <= 1e-6 |
| Split / monolithic differences, normalized | <= 1e-6 |
| Restart / uninterrupted p, w and settlement differences | <= 1e-8 |
| Normalized coupled residuals | < 1e-8 per step |
| Normalized cumulative balance | <= 1e-6 |
| Spatial observed orders, Richardson in time | 1.8–2.2 |
| Temporal observed orders | 0.85–1.15 |
| Stabilized early-time overshoot, undershoot, TV excess / p0 | <= 1e-7 |
| Failure injection | explicit rejection, not silent continuation |

These are explicit benchmark acceptance tolerances, not general accuracy
certificates. The solver is linear, serial, uniform-grid 1D consolidation;
no fracture propagation, contact, nonlinear elasticity, or general 2D/3D
stabilization is implemented.

## Report corrections to carry into the separate narrative

- Both the baseline and revised Phase 2 C++ solver have compiled and passed
  their specified verification checks on foam-extend 4.1.
- In the supplied early-time C++ table, stabilization required up to 26 iterations;
  the different default-factor scan had mean counts up to 20. Do not merge them.
- Iteration-limit entries are nonconvergence, not successful iteration counts.
- At fixed cv and dt, refining h increases cv*dt/h^2.
- Interior averaging annihilates an alternating pattern, but boundary rows mean
  the finite-domain highest-frequency mode is weakly coupled, not an exact
  global checkerboard null mode.
- Describe agreement with laplacianFoam as a framework consistency check. The
  separate Python implementation also shares assumptions; analytical checks and
  observed refinement orders supply additional evidence.
- Describe numerical verification under tested conditions, not a general proof
  that the solver is correct. Do not claim a new method or a universal sharp
  fixed-stress threshold. Publisher-level reference verification and a novelty
  literature review remain outstanding; this update does not claim to do them.

After these gates pass, a separate 2D Mandel benchmark is a reasonable next
milestone. It is not part of this update.
