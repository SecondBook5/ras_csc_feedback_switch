#!/usr/bin/env python3
"""
analyze_hysteresis_and_flux.py

Analyze Ras-driven hysteresis and the phase-plane flux structure of
the Ras–CSC–microenvironment model.

This addresses Step 8 ("Analyze the system") and contributes to
Step 9 ("Make predictions for new scenarios") of the 10-step recipe.

Part A: Ras hysteresis
    - Sweep Ras input f_ras upward from low to high, integrating to
      steady state at each step and using the previous endpoint as the
      next initial condition (up-sweep).
    - Sweep f_ras downward from high to low starting from a malignant
      steady state, again using continuation for the down-sweep.
    - Save a CSV with C* and A* for both sweeps and a figure with
      two panels: C* vs f_ras and A* vs f_ras.

Part B: Phase-plane flux and trajectories
    - For low Ras (benign) and high Ras (malignant), compute a
      vector field in the (C, A) plane, holding T, R, L, M fixed at
      their steady values for that Ras.
    - Overlay trajectories from benign and malignant initial
      conditions, and mark steady states.

"""

from __future__ import annotations

# Standard library imports
from pathlib import Path
from dataclasses import replace
from typing import Tuple

# Numerical imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Model imports
from ras_csc_model import RasCSCModelParams, ras_csc_rhs, load_model_params_from_yaml


def integrate_to_steady(
    params: RasCSCModelParams,
    y0: np.ndarray,
    t_final: float = 300.0,
    n_points: int = 600,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Integrate the model to approximate steady state.

    Args:
        params:
            Parameter set with chosen f_ras.
        y0:
            Initial condition [C, A, T, R, L, M].
        t_final:
            Final integration time.
        n_points:
            Number of evaluation points.

    Returns:
        (t, y) trajectory over time.
    """
    # Build time span and evaluation grid
    t_span = (0.0, float(t_final))
    t_eval = np.linspace(t_span[0], t_span[1], n_points)

    # Wrap RHS with fixed params
    def rhs_wrapped(t: float, y: np.ndarray) -> np.ndarray:
        return ras_csc_rhs(t, y, params)

    # Integrate ODE system
    sol = solve_ivp(
        fun=rhs_wrapped,
        t_span=t_span,
        y0=np.asarray(y0, dtype=float),
        t_eval=t_eval,
        method="RK45",
        rtol=1e-6,
        atol=1e-9,
    )

    # Check for solver success
    if not sol.success:
        raise RuntimeError(
            f"Integration failed in integrate_to_steady: {sol.message}"
        )

    # Return trajectory with shape (n_points, 6)
    return sol.t, sol.y.T


def hysteresis_sweeps(
    params_base: RasCSCModelParams,
    f_ras_min: float = 0.0,
    f_ras_max: float = 1.5,
    n_steps: int = 40,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Perform Ras up-sweep and down-sweep and collect C* and A*.

    Args:
        params_base:
            Baseline parameter set, which will be cloned with different
            f_ras values.
        f_ras_min:
            Minimum Ras input for the sweep.
        f_ras_max:
            Maximum Ras input for the sweep.
        n_steps:
            Number of points between min and max.

    Returns:
        (f_values, C_up, A_up, C_down, A_down) packed into arrays.
        For convenience we return a single array for f_values and
        two arrays of shape (n_steps,) for C and A from up and down
        sweeps respectively.
    """
    # Create Ras grid
    f_values = np.linspace(f_ras_min, f_ras_max, n_steps)

    # Prepare arrays to store steady states
    C_up = np.zeros_like(f_values)
    A_up = np.zeros_like(f_values)
    C_down = np.zeros_like(f_values)
    A_down = np.zeros_like(f_values)

    # Define benign-like starting state for up-sweep
    y_up = np.array([0.01, 0.01, 0.01, 0.0, 0.01, 0.0], dtype=float)

    # Up-sweep: walk Ras from low to high
    for i, f_val in enumerate(f_values):
        params = replace(params_base, f_ras=float(f_val))
        t, y = integrate_to_steady(params, y_up)
        y_up = y[-1, :].copy()
        C_up[i] = y_up[0]
        A_up[i] = y_up[1]

    # Define malignant-like starting state for down-sweep
    y_down = np.array([0.7, 0.7, 10.0, 0.5, 1.0, 0.5], dtype=float)

    # Down-sweep: walk Ras from high to low
    for j, f_val in enumerate(f_values[::-1]):
        params = replace(params_base, f_ras=float(f_val))
        t, y = integrate_to_steady(params, y_down)
        y_down = y[-1, :].copy()
        C_down_val = y_down[0]
        A_down_val = y_down[1]
        idx = n_steps - 1 - j
        C_down[idx] = C_down_val
        A_down[idx] = A_down_val

    # Return the Ras values and steady-state curves
    return f_values, C_up, A_up, C_down, A_down


def phase_plane_flux(
    params: RasCSCModelParams,
    y_star: np.ndarray,
    grid_size: int = 20,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute a 2D vector field in (C, A) by holding the other
    variables T, R, L, M fixed at a chosen steady state.

    Args:
        params:
            Parameter set with a fixed f_ras.
        y_star:
            Reference steady-state [C*, A*, T*, R*, L*, M*].
        grid_size:
            Number of grid points in each dimension.

    Returns:
        (C_grid, A_grid, dC, dA) arrays of shape (grid_size, grid_size).
    """
    # Extract T, R, L, M from the reference steady state
    T_star = float(y_star[2])
    R_star = float(y_star[3])
    L_star = float(y_star[4])
    M_star = float(y_star[5])

    # Define grid ranges for C and A
    C_vals = np.linspace(0.0, 1.0, grid_size)
    A_vals = np.linspace(0.0, 1.0, grid_size)

    # Create meshgrid for plotting
    C_grid, A_grid = np.meshgrid(C_vals, A_vals)

    # Allocate arrays for derivatives
    dC = np.zeros_like(C_grid)
    dA = np.zeros_like(A_grid)

    # Loop over grid points and evaluate RHS in 6D, then keep dC and dA
    for i in range(grid_size):
        for j in range(grid_size):
            C_val = float(C_grid[i, j])
            A_val = float(A_grid[i, j])
            y = np.array(
                [C_val, A_val, T_star, R_star, L_star, M_star],
                dtype=float,
            )
            rhs = ras_csc_rhs(0.0, y, params)
            dC[i, j] = rhs[0]
            dA[i, j] = rhs[1]

    # Return vector field components
    return C_grid, A_grid, dC, dA


def main() -> None:
    """
    Run hysteresis analysis and phase-plane flux plots, writing CSVs
    and figures into results/ and figures/.
    """
    # Determine project root and config path
    project_root = Path(__file__).resolve().parents[1]
    config_path = project_root / "config" / "model_parameters.yaml"

    # Load baseline parameter set
    params_base = load_model_params_from_yaml(config_path)

    # Perform Ras hysteresis sweeps
    f_vals, C_up, A_up, C_down, A_down = hysteresis_sweeps(
        params_base,
        f_ras_min=0.0,
        f_ras_max=1.5,
        n_steps=40,
    )

    # Create results directory
    results_dir = project_root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    # Path for hysteresis CSV
    hyst_csv = results_dir / "ras_hysteresis.csv"

    # Write hysteresis data to CSV
    with hyst_csv.open("w", encoding="utf-8") as fh:
        fh.write("f_ras,C_up,A_up,C_down,A_down\n")
        for f, Cu, Au, Cd, Ad in zip(f_vals, C_up, A_up, C_down, A_down):
            fh.write(
                f"{f:.6g},{Cu:.6g},{Au:.6g},{Cd:.6g},{Ad:.6g}\n"
            )

    # Make hysteresis figure
    figures_dir = project_root / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    hyst_fig = figures_dir / "ras_hysteresis.png"

    # Build figure with two stacked panels
    plt.figure(figsize=(8, 8))

    # Panel 1: C* vs f_ras
    plt.subplot(2, 1, 1)
    plt.plot(f_vals, C_up, "o-", label="Up-sweep (benign → malignant)")
    plt.plot(f_vals, C_down, "s-", label="Down-sweep (malignant → benign)")
    plt.ylabel("CSC steady state C*")
    plt.title("Ras hysteresis: CSC fraction and angiogenesis")
    plt.grid(alpha=0.3)
    plt.legend()

    # Panel 2: A* vs f_ras
    plt.subplot(2, 1, 2)
    plt.plot(f_vals, A_up, "o-", label="Up-sweep (benign → malignant)")
    plt.plot(f_vals, A_down, "s-", label="Down-sweep (malignant → benign)")
    plt.xlabel("Ras input f_RAS")
    plt.ylabel("Angiogenesis steady state A*")
    plt.grid(alpha=0.3)
    plt.legend()

    # Save hysteresis figure
    plt.tight_layout()
    plt.savefig(hyst_fig, dpi=300)

    # Now build phase-plane flux plots for low and high Ras
    # First, get benign and malignant steady states with the same helper
    benign_params = replace(params_base, f_ras=0.2)
    malignant_params = replace(params_base, f_ras=1.0)

    # Start benign from low IC
    y0_benign = np.array([0.01, 0.01, 0.01, 0.0, 0.01, 0.0], dtype=float)
    _, traj_benign = integrate_to_steady(benign_params, y0_benign)
    y_star_benign = traj_benign[-1, :].copy()

    # Start malignant from high IC
    y0_malignant = np.array([0.7, 0.7, 10.0, 0.5, 1.0, 0.5], dtype=float)
    _, traj_malignant = integrate_to_steady(malignant_params, y0_malignant)
    y_star_malignant = traj_malignant[-1, :].copy()

    # Compute vector fields
    Cg_low, Ag_low, dC_low, dA_low = phase_plane_flux(
        benign_params, y_star_benign, grid_size=25
    )
    Cg_high, Ag_high, dC_high, dA_high = phase_plane_flux(
        malignant_params, y_star_malignant, grid_size=25
    )

    # Create phase-plane figure
    flux_fig = figures_dir / "phase_plane_flux_low_high_ras.png"
    plt.figure(figsize=(10, 5))

    # Panel left: low Ras
    plt.subplot(1, 2, 1)
    plt.quiver(
        Cg_low,
        Ag_low,
        dC_low,
        dA_low,
        color="gray",
        angles="xy",
        scale_units="xy",
        scale=1.0,
        alpha=0.6,
    )
    plt.plot(traj_benign[:, 0], traj_benign[:, 1], color="tab:blue",
             label="Trajectory")
    plt.scatter(
        y_star_benign[0],
        y_star_benign[1],
        color="tab:green",
        s=60,
        label="Attractor",
    )
    plt.xlabel("CSC fraction C")
    plt.ylabel("Angiogenesis A")
    plt.title("Phase-plane flux at low Ras")
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
    plt.grid(alpha=0.3)
    plt.legend()

    # Panel right: high Ras
    plt.subplot(1, 2, 2)
    plt.quiver(
        Cg_high,
        Ag_high,
        dC_high,
        dA_high,
        color="gray",
        angles="xy",
        scale_units="xy",
        scale=1.0,
        alpha=0.6,
    )
    plt.plot(traj_malignant[:, 0], traj_malignant[:, 1],
             color="tab:orange",
             label="Trajectory")
    plt.scatter(
        y_star_malignant[0],
        y_star_malignant[1],
        color="tab:red",
        s=60,
        label="Attractor",
    )
    plt.xlabel("CSC fraction C")
    plt.ylabel("Angiogenesis A")
    plt.title("Phase-plane flux at high Ras")
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
    plt.grid(alpha=0.3)
    plt.legend()

    # Save phase-plane figure
    plt.tight_layout()
    plt.savefig(flux_fig, dpi=300)

    # Print locations to console
    print(f"[INFO] Ras hysteresis CSV written to {hyst_csv}")
    print(f"[INFO] Ras hysteresis figure saved to {hyst_fig}")
    print(f"[INFO] Phase-plane flux figure saved to {flux_fig}")


if __name__ == "__main__":
    main()
