#!/usr/bin/env python3
"""
plot_flux_phase_plane.py

Construct a phase-plane "flux" plot for the Ras–CSC–microenvironment model.

We project the 6D dynamics onto the (C, A) plane (CSC fraction vs angiogenesis)
at a fixed Ras input f_RAS. For each grid point (C, A), we:

  1) Compute quasi-steady values of T, R, L, M given C, A and the model
     parameters (assuming those four variables relax quickly).
  2) Evaluate the full 6D ODE right-hand side at [C, A, T*, R*, L*, M*].
  3) Keep only dC/dt and dA/dt to obtain a 2D vector field (flux) in (C, A).

We then overlay trajectories for two initial conditions (benign-like and
malignant-like) to show how the system flows toward the two attractors
at that Ras level.

Results:
  - results/flux_field_fras_1.0.csv      (flux field on a grid)
  - figures/phase_plane_flux_fras_1.0.png (quiver + trajectories)

This script uses the same "hysteresis-tuned" parameters as run_ras_hysteresis.py
to keep the bistable structure visible. It does not modify the global model
or the H1 test script.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from ras_csc_model import (
    RasCSCModelParams,
    ras_csc_rhs,
    hill,
    load_model_params_from_yaml,
)


def quasi_steady_TRLM(
    C: float,
    A: float,
    params: RasCSCModelParams,
) -> Tuple[float, float, float, float]:
    """
    Compute quasi-steady T, R, L, M given C and A.

    We assume:
      dT/dt = 0, dR/dt = 0, dL/dt = 0, dM/dt = 0
    and solve the linear fixed-point equations for T, R, L, M:

      T* = (alpha_T * A + gamma_T * C) / delta_T
      H_TR = hill(T*; K_R, n_R)
      R* = [eta_R * H_TR] / [eta_R * H_TR + delta_R]
      L* = (mu_L / delta_L) * A
      H_LRM = hill(L* R*; K_M, n_M)
      M* = [eta_M * H_LRM] / [eta_M * H_LRM + delta_M]

    This gives a consistent 1D manifold (T*,R*,L*,M*) for each (C,A)
    compatible with the full 6D model.
    """
    # Ensure C, A are floats and non-negative
    C_val = max(float(C), 0.0)
    A_val = max(float(A), 0.0)

    # Quasi-steady T from dT/dt = 0
    T_star = (params.alpha_T * A_val + params.gamma_T * C_val) / max(
        params.delta_T, 1e-12
    )

    # Hill(T; K_R, n_R) to drive LEPR
    H_TR = hill(T_star, params.K_R, params.n_R)

    # Quasi-steady R from linear fixed point of dR/dt = 0
    # dR/dt = eta_R * H_TR * (1 - R) - delta_R * R
    #       = eta_R * H_TR - (eta_R * H_TR + delta_R) * R
    # => R* = [eta_R * H_TR] / [eta_R * H_TR + delta_R]
    num_R = params.eta_R * H_TR
    denom_R = num_R + params.delta_R
    if denom_R <= 0.0:
        R_star = 0.0
    else:
        R_star = num_R / denom_R

    # Quasi-steady L from dL/dt = 0
    L_star = (params.mu_L * A_val) / max(params.delta_L, 1e-12)

    # Hill(L R; K_M, n_M) to drive mTOR
    LR_signal = L_star * R_star
    H_LRM = hill(LR_signal, params.K_M, params.n_M)

    # Quasi-steady M from dM/dt = 0
    # dM/dt = eta_M * H_LRM * (1 - M) - delta_M * M
    #       = eta_M * H_LRM - (eta_M * H_LRM + delta_M) * M
    # => M* = [eta_M * H_LRM] / [eta_M * H_LRM + delta_M]
    num_M = params.eta_M * H_LRM
    denom_M = num_M + params.delta_M
    if denom_M <= 0.0:
        M_star = 0.0
    else:
        M_star = num_M / denom_M

    return T_star, R_star, L_star, M_star


def compute_flux_field(
    base_params: RasCSCModelParams,
    f_ras: float,
    C_vals: np.ndarray,
    A_vals: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray,
           np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute the (dC/dt, dA/dt) flux field on a (C, A) grid at fixed Ras.

    For each grid point (C_i, A_j), we:
      - Build a parameter set with f_ras fixed.
      - Compute quasi-steady T, R, L, M.
      - Evaluate ras_csc_rhs at y = [C, A, T*, R*, L*, M*].
      - Store dC/dt and dA/dt.

    Returns:
        C_grid, A_grid, dC_dt, dA_dt, T_star, R_star, L_star, M_star
        where all arrays have shape (len(A_vals), len(C_vals)).
    """
    # Ensure 1D arrays
    C_vals = np.asarray(C_vals, dtype=float)
    A_vals = np.asarray(A_vals, dtype=float)

    # Build meshgrid for plotting
    C_grid, A_grid = np.meshgrid(C_vals, A_vals)

    # Arrays for flux components and quasi-steady variables
    dC_dt = np.zeros_like(C_grid)
    dA_dt = np.zeros_like(A_grid)
    T_star = np.zeros_like(C_grid)
    R_star = np.zeros_like(C_grid)
    L_star = np.zeros_like(C_grid)
    M_star = np.zeros_like(C_grid)

    # Use a copy of params with f_ras fixed at the desired value
    params_fixed = replace(base_params, f_ras=float(f_ras))

    # Loop over grid points
    n_rows, n_cols = C_grid.shape
    for i in range(n_rows):
        for j in range(n_cols):
            C = C_grid[i, j]
            A = A_grid[i, j]

            # Compute quasi-steady T, R, L, M at this (C, A)
            T_qs, R_qs, L_qs, M_qs = quasi_steady_TRLM(C, A, params_fixed)

            # Store quasi-steady values
            T_star[i, j] = T_qs
            R_star[i, j] = R_qs
            L_star[i, j] = L_qs
            M_star[i, j] = M_qs

            # Evaluate full RHS at this point
            y_vec = np.array([C, A, T_qs, R_qs, L_qs, M_qs], dtype=float)
            dydt = ras_csc_rhs(t=0.0, y=y_vec, params=params_fixed)

            # Keep only dC/dt and dA/dt
            dC_dt[i, j] = dydt[0]
            dA_dt[i, j] = dydt[1]

    return C_grid, A_grid, dC_dt, dA_dt, T_star, R_star, L_star, M_star


def simulate_trajectory(
    params: RasCSCModelParams,
    f_ras: float,
    y0: np.ndarray,
    t_end: float = 400.0,
    n_points: int = 800,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Simulate a single full 6D trajectory at fixed Ras and return
    (time, trajectory) for plotting.

    Args:
        params:
            Base parameter set.
        f_ras:
            Fixed Ras input for this run.
        y0:
            Initial state vector [C, A, T, R, L, M].
        t_end:
            Final time for integration.
        n_points:
            Number of evaluation points along the trajectory.
    """
    # Fix Ras input
    params_fixed = replace(params, f_ras=float(f_ras))

    # Time span and evaluation grid
    t_span = (0.0, float(t_end))
    t_eval = np.linspace(t_span[0], t_span[1], n_points)

    # RHS wrapper
    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        return ras_csc_rhs(t, y, params_fixed)

    sol = solve_ivp(
        fun=rhs,
        t_span=t_span,
        y0=np.asarray(y0, dtype=float),
        t_eval=t_eval,
        method="RK45",
        rtol=1e-6,
        atol=1e-9,
    )

    if not sol.success:
        raise RuntimeError(
            f"ODE integration failed in simulate_trajectory: {sol.message}"
        )

    t_out = sol.t
    y_out = sol.y.T
    return t_out, y_out


def main() -> None:
    """
    Build the (C, A) flux field and plot trajectories for a fixed Ras.

    Uses hysteresis-tuned parameters so that both benign and malignant
    attractors coexist at f_RAS ≈ 1.0, making the flux structure clear.
    """
    # ------------------------------------------------------------------
    # Locate project root and YAML
    # ------------------------------------------------------------------
    project_root = Path(__file__).resolve().parents[1]
    config_path = project_root / "config" / "model_parameters.yaml"

    # ------------------------------------------------------------------
    # Load base parameters
    # ------------------------------------------------------------------
    try:
        params_base = load_model_params_from_yaml(config_path)
    except Exception as exc:
        raise RuntimeError(
            f"Failed to load model parameters from {config_path}"
        ) from exc

    # ------------------------------------------------------------------
    # Hysteresis / flux tuning (same spirit as run_ras_hysteresis.py)
    #
    # Make Ras more influential and the loop slightly less self-sustaining:
    #   - rho_C, alpha_A up (Ras-driven seeding into C and A)
    #   - eta_C, eta_CM down (feedback strengths)
    # ------------------------------------------------------------------
    params_flux = replace(
        params_base,
        rho_C=params_base.rho_C * 6.0,
        alpha_A=params_base.alpha_A * 5.0,
        eta_C=params_base.eta_C * 0.6,
        eta_CM=params_base.eta_CM * 0.5,
    )

    # Choose Ras level for flux plot (middle of bistable window)
    f_ras = 1.0

    # ------------------------------------------------------------------
    # Build (C, A) grid and compute flux field
    # ------------------------------------------------------------------
    C_vals = np.linspace(0.0, 1.0, 41)
    A_vals = np.linspace(0.0, 1.0, 41)

    (
        C_grid,
        A_grid,
        dC_dt,
        dA_dt,
        T_star,
        R_star,
        L_star,
        M_star,
    ) = compute_flux_field(
        base_params=params_flux,
        f_ras=f_ras,
        C_vals=C_vals,
        A_vals=A_vals,
    )

    # ------------------------------------------------------------------
    # Simulate two trajectories at the same Ras:
    #   - Benign-like initial condition
    #   - Malignant-like initial condition
    # ------------------------------------------------------------------
    y0_benign = np.array(
        [0.01, 0.01, 0.01, 0.0, 0.01, 0.0],
        dtype=float,
    )

    y0_malignant = np.array(
        [0.7, 0.7, 10.0, 0.5, 1.0, 0.5],
        dtype=float,
    )

    _, traj_benign = simulate_trajectory(
        params=params_flux,
        f_ras=f_ras,
        y0=y0_benign,
        t_end=400.0,
        n_points=800,
    )

    _, traj_malignant = simulate_trajectory(
        params=params_flux,
        f_ras=f_ras,
        y0=y0_malignant,
        t_end=400.0,
        n_points=800,
    )

    # Project trajectories onto (C, A)
    C_benign = traj_benign[:, 0]
    A_benign = traj_benign[:, 1]

    C_malignant = traj_malignant[:, 0]
    A_malignant = traj_malignant[:, 1]

    # ------------------------------------------------------------------
    # Dump flux field to CSV for the paper
    # ------------------------------------------------------------------
    results_dir = project_root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    csv_path = results_dir / "flux_field_fras_1.0.csv"

    with csv_path.open("w", encoding="utf-8") as fh:
        fh.write(
            "C,A,dC_dt,dA_dt,T_star,R_star,L_star,M_star,f_ras\n"
        )
        n_rows, n_cols = C_grid.shape
        for i in range(n_rows):
            for j in range(n_cols):
                fh.write(
                    f"{C_grid[i, j]:.6g},{A_grid[i, j]:.6g},"
                    f"{dC_dt[i, j]:.6g},{dA_dt[i, j]:.6g},"
                    f"{T_star[i, j]:.6g},{R_star[i, j]:.6g},"
                    f"{L_star[i, j]:.6g},{M_star[i, j]:.6g},"
                    f"{f_ras:.6g}\n"
                )

    print(f"[INFO] Flux field CSV written to {csv_path}")

    # ------------------------------------------------------------------
    # Plot phase-plane flux with trajectories
    # ------------------------------------------------------------------
    figures_dir = project_root / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    fig_path = figures_dir / "phase_plane_flux_fras_1.0.png"

    fig, ax = plt.subplots(figsize=(7, 6))

    # Normalize arrows a bit for readability
    # Avoid zero division if all arrows are tiny
    mag = np.sqrt(dC_dt**2 + dA_dt**2)
    max_mag = mag.max()
    if max_mag > 0:
        U = dC_dt / max_mag
        V = dA_dt / max_mag
    else:
        U = dC_dt
        V = dA_dt

    ax.quiver(
        C_grid,
        A_grid,
        U,
        V,
        angles="xy",
        scale=None,
        scale_units="xy",
        alpha=0.6,
        width=0.003,
    )

    # Overlay trajectories
    ax.plot(
        C_benign,
        A_benign,
        "-", linewidth=2.0, label="Benign IC trajectory",
    )
    ax.plot(
        C_malignant,
        A_malignant,
        "-", linewidth=2.0, label="Malignant IC trajectory",
    )

    # Mark final steady states
    ax.plot(
        C_benign[-1],
        A_benign[-1],
        "o",
        markersize=8,
        label="Benign attractor",
    )
    ax.plot(
        C_malignant[-1],
        A_malignant[-1],
        "s",
        markersize=8,
        label="Malignant attractor",
    )

    ax.set_xlabel("CSC fraction C")
    ax.set_ylabel("Angiogenesis A")
    ax.set_title(
        f"Phase-plane flux at Ras input f_RAS = {f_ras:.1f}\n"
        "Ras–CSC–microenvironment model"
    )
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.grid(alpha=0.3)
    ax.legend(loc="best")

    plt.tight_layout()
    plt.savefig(fig_path, dpi=300)
    plt.close(fig)

    print(f"[INFO] Flux phase-plane figure saved to {fig_path}")


if __name__ == "__main__":
    main()
