"""
Terzaghi 1D consolidation: exact analytical solution.

Governing equation (excess pore pressure u, drainage path H):
    du/dt = cv * d2u/dz2         0 <= z <= H

Boundary conditions:
    u(0, t)      = 0             (drained boundary, z=0)
    du/dz(H, t)  = 0             (impermeable boundary, z=H)

Initial condition:
    u(z, 0) = u0                 (instantaneous applied load)

Series solution (Terzaghi, 1943):
    u(z,t)/u0 = sum_{m=0}^inf  (2/M) * sin(M*z/H) * exp(-M^2 * Tv)
    M  = (pi/2) * (2m + 1)
    Tv = cv * t / H^2            (dimensionless time factor)

Average degree of consolidation:
    U(Tv) = 1 - sum_{m=0}^inf (2/M^2) * exp(-M^2 * Tv)

This module computes both and, when run directly, plots the classic
isochrones (pressure profile at several Tv) and the U-vs-Tv consolidation
curve. verify.py imports excess_pressure() and degree_of_consolidation()
as the exact reference for the numerical solution.
"""

import numpy as np
import matplotlib.pyplot as plt

H = 1.0        # drainage path length [m]
u0 = 1.0       # initial excess pore pressure (normalized to 1)
N_TERMS = 200  # series truncation; converges quickly except at very small Tv


def excess_pressure(z, Tv, n_terms=N_TERMS):
    """u(z,t)/u0 at dimensionless depth z/H and time factor Tv."""
    z = np.asarray(z, dtype=float)
    total = np.zeros_like(z)
    for m in range(n_terms):
        M = (np.pi / 2.0) * (2 * m + 1)
        total += (2.0 / M) * np.sin(M * z / H) * np.exp(-(M ** 2) * Tv)
    return total


def degree_of_consolidation(Tv, n_terms=N_TERMS):
    """Average degree of consolidation U(Tv), scalar or array input."""
    Tv = np.asarray(Tv, dtype=float)
    total = np.zeros_like(Tv)
    for m in range(n_terms):
        M = (np.pi / 2.0) * (2 * m + 1)
        total += (2.0 / M ** 2) * np.exp(-(M ** 2) * Tv)
    return 1.0 - total


if __name__ == "__main__":
    z_over_H = np.linspace(0, 1, 200)
    Tv_values = [0.05, 0.1, 0.2, 0.4, 0.8, 1.5]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Left: isochrones (pressure profile at several times)
    ax = axes[0]
    for Tv in Tv_values:
        u = excess_pressure(z_over_H * H, Tv)
        ax.plot(u, z_over_H, label=f"Tv = {Tv}")
    ax.invert_yaxis()
    ax.set_xlabel(r"$u / u_0$  (excess pore pressure)")
    ax.set_ylabel(r"$z / H$")
    ax.set_title("Terzaghi consolidation: isochrones")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # Right: average degree of consolidation vs time factor
    ax = axes[1]
    Tv_range = np.linspace(0.001, 2.0, 300)
    U = degree_of_consolidation(Tv_range)
    ax.plot(Tv_range, U, color="black")
    ax.set_xlabel(r"$T_v$  (time factor)")
    ax.set_ylabel(r"$U$  (average degree of consolidation)")
    ax.set_title("Consolidation curve")
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig("terzaghi_analytical.png", dpi=150)
    print("Plot written to terzaghi_analytical.png")

