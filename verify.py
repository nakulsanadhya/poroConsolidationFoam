"""
Verification driver for poroConsolidationFoam (Terzaghi 1D consolidation).

Runs the example case, reads the numerical excess-pore-pressure field, and
compares it with the exact series solution evaluated at the solver's own cell
centres and output times. Also runs mesh- and time-step-refinement studies
and reports observed orders of accuracy.

Usage
-----
    python verify.py compare      [--solver NAME]   # isochrones + error table
    python verify.py convergence  [--solver NAME]   # observed orders of accuracy
    python verify.py check        [--solver NAME]   # pass/fail regression test

--solver defaults to poroConsolidationFoam (build it first with wmake).
--solver laplacianFoam runs OpenFOAM's stock laplacianFoam on an equivalent
case (field u -> T, cv -> DT). It solves the same equation,
ddt(T) = laplacian(DT, T), with the same discretization, so it cross-checks
the case setup and the numerics independently of the custom solver.

Requires an OpenFOAM environment. If WM_PROJECT_DIR isn't set, the script
tries to source a known bashrc (OF_BASHRC env var, else the Debian path).
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
TEMPLATE = ROOT / "case"
RUNS = ROOT / "runs"
OUT = ROOT / "validation"
sys.path.insert(0, str(OUT))
from terzaghi_analytical import excess_pressure, degree_of_consolidation  # noqa: E402

H = 1.0      # drainage length [m], must match blockMeshDict
CV = 0.01    # coefficient of consolidation [m^2/s], must match transportProperties


# --------------------------------------------------------------------------
# Case setup and execution
# --------------------------------------------------------------------------

def _shell_prefix():
    if os.environ.get("WM_PROJECT_DIR"):
        return ""
    for cand in (os.environ.get("OF_BASHRC"), "/usr/share/openfoam/etc/bashrc"):
        if cand and Path(cand).exists():
            return f"source {cand} >/dev/null 2>&1; "
    return ""


def _set_entry(text, key, value):
    new, n = re.subn(rf"^(\s*{key}\s+)[^;]+;", rf"\g<1>{value};", text, count=1, flags=re.M)
    if n != 1:
        raise RuntimeError(f"could not set '{key}'")
    return new


def make_case(name, n_cells, dt, end_time, write_interval, solver):
    dst = RUNS / f"{name}_{solver}"
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(TEMPLATE, dst)

    # blockMeshDict lives in constant/polyMesh/, where foam-extend (and
    # OpenFOAM <= 3.x) look for it. Newer OpenFOAM versions look in system/,
    # so a copy is placed there too.
    bmd = dst / "constant" / "polyMesh" / "blockMeshDict"
    # Only touch the cell-count list right after the hex vertex list; a looser
    # pattern would also rewrite "simpleGrading (1 1 1)" into a graded mesh.
    text, n = re.subn(r"(hex\s*\([^)]*\)\s*)\(1 1 \d+\)", rf"\g<1>(1 1 {n_cells})",
                      bmd.read_text(), count=1)
    if n != 1:
        raise RuntimeError("could not set cell count in blockMeshDict")
    bmd.write_text(text)
    shutil.copy(bmd, dst / "system" / "blockMeshDict")

    cd = dst / "system" / "controlDict"
    text = cd.read_text()
    text = _set_entry(text, "application", solver)
    text = _set_entry(text, "deltaT", f"{dt:g}")
    text = _set_entry(text, "endTime", f"{end_time:g}")
    text = _set_entry(text, "writeInterval", f"{write_interval:g}")
    cd.write_text(text)

    field = "u"
    if solver == "laplacianFoam":
        field = "T"
        u = dst / "0" / "u"
        u.write_text(re.sub(r"object\s+u;", "object      T;", u.read_text()))
        u.rename(dst / "0" / "T")
        tp = dst / "constant" / "transportProperties"
        tp.write_text(re.sub(r"^\s*cv\s+cv\s", "DT              DT ", tp.read_text(), flags=re.M))
    return dst, field


def run_case(dst, solver):
    prefix = _shell_prefix()
    for cmd in ("blockMesh", solver):
        log = dst / f"log.{cmd}"
        res = subprocess.run(
            ["bash", "-c", f"{prefix}{cmd} -case {dst}"],
            capture_output=True, text=True,
        )
        log.write_text(res.stdout + res.stderr)
        if res.returncode != 0 or "FOAM FATAL" in res.stdout + res.stderr:
            raise RuntimeError(f"{cmd} failed, see {log}")


# --------------------------------------------------------------------------
# Reading results
# --------------------------------------------------------------------------

def time_dir(dst, t):
    for d in dst.iterdir():
        try:
            if d.is_dir() and abs(float(d.name) - t) < 1e-9:
                return d
        except ValueError:
            continue
    raise FileNotFoundError(f"no time directory for t={t} in {dst}")


def read_field(dst, t, field):
    text = (time_dir(dst, t) / field).read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\((.*?)\)\s*;", text, re.S)
    if not m:
        raise ValueError(f"could not parse internalField in {field} at t={t}")
    vals = np.array(m.group(2).split(), dtype=float)
    assert len(vals) == int(m.group(1))
    return vals


def _foam_list(path):
    text = path.read_text()
    text = text[text.index("}", text.index("FoamFile")) + 1:]
    m = re.search(r"(\d+)\s*\(", text)
    return int(m.group(1)), text[m.end():]


def cell_centres_z(dst):
    """Cell-centre z coordinates computed from constant/polyMesh.

    Reads the actual mesh rather than assuming uniform spacing or a cell
    ordering, so a grading or ordering mistake shows up as a wrong answer
    instead of being silently masked. For the box-shaped cells here, the
    mean of a cell's face centres equals its centroid.
    """
    mesh = dst / "constant" / "polyMesh"
    n, body = _foam_list(mesh / "points")
    pts = np.array([p.split() for p in re.findall(r"\(([^()]+)\)", body)[:n]], dtype=float)
    n, body = _foam_list(mesh / "faces")
    faces = [list(map(int, f.split())) for f in re.findall(r"\d+\(([^()]+)\)", body)[:n]]
    n, body = _foam_list(mesh / "owner")
    owner = np.array(body.split(")")[0].split(), dtype=int)[:n]
    n, body = _foam_list(mesh / "neighbour")
    neigh = np.array(body.split(")")[0].split(), dtype=int)[:n]

    fz = np.array([pts[f, 2].mean() for f in faces])
    ncell = owner.max() + 1
    acc, cnt = np.zeros(ncell), np.zeros(ncell)
    np.add.at(acc, owner, fz); np.add.at(cnt, owner, 1)
    np.add.at(acc, neigh, fz[: len(neigh)]); np.add.at(cnt, neigh, 1)
    return acc / cnt


def errors(num, exact):
    diff = num - exact
    return np.sqrt(np.mean(diff ** 2)), np.max(np.abs(diff))


# --------------------------------------------------------------------------
# Commands
# --------------------------------------------------------------------------

COMPARE_TIMES = [5, 10, 20, 40, 80]   # Tv = 0.05, 0.1, 0.2, 0.4, 0.8


def compare(solver, n_cells=100, dt=0.25):
    dst, field = make_case("compare", n_cells, dt, 80, 5, solver)
    run_case(dst, solver)
    z = cell_centres_z(dst)

    rows = []
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.8))
    for t in COMPARE_TIMES:
        Tv = CV * t / H ** 2
        num = read_field(dst, t, field)
        exact = excess_pressure(z, Tv)
        l2, linf = errors(num, exact)
        U_num = 1.0 - num.mean()
        U_exact = float(degree_of_consolidation(Tv))
        rows.append((Tv, l2, linf, U_num, U_exact))

        line, = ax1.plot(exact, z, lw=1.5, label=f"Tv = {Tv:g}")
        ax1.plot(num[::5], z[::5], "o", ms=3.5, color=line.get_color())
        ax2.plot(num - exact, z, lw=1.3, color=line.get_color(), label=f"Tv = {Tv:g}")

    ax1.invert_yaxis()
    ax1.set_xlabel(r"$u/u_0$")
    ax1.set_ylabel(r"$z/H$")
    ax1.set_title(f"{solver}: markers = numerical, lines = exact")
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.3)
    ax2.invert_yaxis()
    ax2.set_xlabel(r"numerical $-$ exact")
    ax2.set_title(f"Pointwise error ({n_cells} cells, $\\Delta t$ = {dt:g} s)")
    ax2.axvline(0, color="k", lw=0.6)
    ax2.grid(alpha=0.3)
    plt.tight_layout()
    out = OUT / "solver_vs_analytical.png"
    plt.savefig(out, dpi=140)

    print(f"\n{solver}, {n_cells} cells, dt = {dt:g} s")
    print(f"{'Tv':>6} {'L2 error':>11} {'Linf error':>11} {'U numerical':>12} {'U exact':>9}")
    for Tv, l2, linf, Un, Ue in rows:
        print(f"{Tv:>6g} {l2:>11.2e} {linf:>11.2e} {Un:>12.5f} {Ue:>9.5f}")
    print(f"Plot: {out}")
    return rows


def _orders(errs):
    return [np.log2(errs[i] / errs[i + 1]) for i in range(len(errs) - 1)]


def convergence(solver, t_end=20):
    Tv = CV * t_end / H ** 2

    # Spatial study: small dt so time error stays well below space error.
    Ns, dt_fine = [10, 20, 40, 80, 160], 2e-4
    e_space = []
    for n in Ns:
        dst, field = make_case(f"space_N{n}", n, dt_fine, t_end, t_end, solver)
        run_case(dst, solver)
        l2, _ = errors(read_field(dst, t_end, field), excess_pressure(cell_centres_z(dst), Tv))
        e_space.append(l2)

    # Temporal study: fine mesh so space error stays well below time error.
    dts, N_fine = [2.0, 1.0, 0.5, 0.25, 0.125], 800
    e_time = []
    for dt in dts:
        dst, field = make_case(f"time_dt{dt:g}", N_fine, dt, t_end, t_end, solver)
        run_case(dst, solver)
        l2, _ = errors(read_field(dst, t_end, field), excess_pressure(cell_centres_z(dst), Tv))
        e_time.append(l2)

    dz = [H / n for n in Ns]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    ax1.loglog(dz, e_space, "o-", label="L2 error")
    ax1.loglog(dz, e_space[0] * (np.array(dz) / dz[0]) ** 2, "k--", lw=0.8, label="slope 2")
    ax1.set_xlabel(r"$\Delta z$ [m]")
    ax1.set_ylabel("L2 error")
    ax1.set_title(f"Mesh refinement ($\\Delta t$ = {dt_fine:g} s, Tv = {Tv:g})")
    ax1.legend()
    ax1.grid(alpha=0.3, which="both")
    ax2.loglog(dts, e_time, "s-", color="C1", label="L2 error")
    ax2.loglog(dts, e_time[0] * (np.array(dts) / dts[0]), "k--", lw=0.8, label="slope 1")
    ax2.set_xlabel(r"$\Delta t$ [s]")
    ax2.set_title(f"Time-step refinement ({N_fine} cells, Tv = {Tv:g})")
    ax2.legend()
    ax2.grid(alpha=0.3, which="both")
    plt.tight_layout()
    out = OUT / "convergence.png"
    plt.savefig(out, dpi=140)

    print(f"\nMesh refinement ({solver}, dt = {dt_fine:g} s, Tv = {Tv:g})")
    print(f"{'cells':>6} {'L2 error':>11} {'observed order':>15}")
    oz = [None] + _orders(e_space)
    for n, e, p in zip(Ns, e_space, oz):
        print(f"{n:>6} {e:>11.3e} {('' if p is None else f'{p:.2f}'):>15}")
    print(f"\nTime-step refinement ({solver}, {N_fine} cells, Tv = {Tv:g})")
    print(f"{'dt [s]':>6} {'L2 error':>11} {'observed order':>15}")
    ot = [None] + _orders(e_time)
    for dt, e, p in zip(dts, e_time, ot):
        print(f"{dt:>6g} {e:>11.3e} {('' if p is None else f'{p:.2f}'):>15}")
    print(f"Plot: {out}")


def check(solver, tol_linf=5e-3, tol_U=2e-3):
    """Regression test: fails (exit 1) if the solution drifts from the exact one."""
    rows = compare(solver)
    ok = True
    for Tv, _, linf, Un, Ue in rows:
        if Tv < 0.1:
            continue  # earliest snapshot is dominated by the initial jump; not a pass/fail criterion
        if linf > tol_linf or abs(Un - Ue) > tol_U:
            print(f"FAIL at Tv = {Tv:g}: Linf = {linf:.2e} (tol {tol_linf:g}), "
                  f"|dU| = {abs(Un - Ue):.2e} (tol {tol_U:g})")
            ok = False
    print("\nCHECK PASSED" if ok else "\nCHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["compare", "convergence", "check"])
    ap.add_argument("--solver", default="poroConsolidationFoam")
    args = ap.parse_args()
    {"compare": compare, "convergence": convergence, "check": check}[args.command](args.solver)
