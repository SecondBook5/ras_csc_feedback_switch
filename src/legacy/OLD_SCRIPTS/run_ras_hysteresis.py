#!/usr/bin/env python3
"""
run_ras_hysteresis.py

Construct a Ras hysteresis loop for the Ras–CSC–microenvironment model.

We compute steady-state CSC fraction (C*) and angiogenesis (A*) as a
function of Ras input f_RAS, using two branches:

  1) Up-sweep: start from a benign-like state at low Ras and increase f_RAS.
  2) Down-sweep: start from a malignant-like state at high Ras and decrease f_RAS.

Path dependence (different limits for the two branches at the same f_RAS)
indicates hysteresis. We also shade the Ras window where the two branches
differ (bistable region).

This script uses the same base YAML parameters as the main model, but applies
a small "hysteresis tuning" only here to make Ras a stronger gate and the
loop slightly less self-sustaining. The H1 test script remains unchanged.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from dataclasses import replace
from pathlib import Path

from ras_csc_model import (
    RasCSCModelParams,
    ras_csc_rhs,
    load_model_params_from_yaml,
)


def _simulate_to_steady(
    params: RasCSCModelParams,
    y0: np.ndarray,
    t_start: float = 0.0,
    t_end: float = 400.0,
    n_points: int = 400,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Integrate the ODEs from (t_start, y0) to t_end and return the full
    time grid and trajectory. The final row of the trajectory is used
    as the steady-state approximation for hysteresis scanning.
    """
    t_span = (float(t_start), float(t_end))
    t_eval = np.linspace(t_start, t_end, n_points)

    def rhs_wrapped(t: float, y: np.ndarray) -> np.ndarray:
        return ras_csc_rhs(t, y, params)

    sol = solve_ivp(
        fun=rhs_wrapped,
        t_span=t_span,
        y0=np.asarray(y0, dtype=float),
        t_eval=t_eval,
        method="RK45",
        rtol=1e-6,
        atol=1e-9,
    )

    if not sol.success:
        raise RuntimeError(
            f"ODE integration failed in _simulate_to_steady: {sol.message}"
        )

    t_out = sol.t
    y_out = sol.y.T
    return t_out, y_out


def _run_branch(
    base_params: RasCSCModelParams,
    ras_values: np.ndarray,
    y0_initial: np.ndarray,
) -> np.ndarray:
    """
    Run a sequence of simulations along a list of Ras values, always using
    the steady-state from the previous simulation as the initial condition
    for the next one.

    Args:
        base_params:
            Parameter set that will be copied and modified for each f_RAS.
        ras_values:
            1D array of Ras inputs f_RAS (scan order as given).
        y0_initial:
            Initial state [C, A, T, R, L, M] for the first Ras value.

    Returns:
        states:
            Array of shape (len(ras_values), 6) with steady states
            [C*, A*, T*, R*, L*, M*] for each Ras value in order.
    """
    n = len(ras_values)
    states = np.zeros((n, 6), dtype=float)

    y_prev = np.asarray(y0_initial, dtype=float).copy()

    for i, f_ras in enumerate(ras_values):
        # Use a copy of the base params with updated Ras input
        params_i = replace(base_params, f_ras=float(f_ras))

        _, traj = _simulate_to_steady(params=params_i, y0=y_prev)

        # Last row is the steady state for this f_RAS
        y_ss = traj[-1, :].copy()
        states[i, :] = y_ss

        # Use this steady state as the initial condition for the next Ras
        y_prev = y_ss

    return states


def main() -> None:
    """
    Build the full Ras hysteresis loop and generate:

      - A CSV file with steady states for both branches.
      - A figure showing C* and A* vs Ras with the bistable window shaded.
    """
    project_root = Path(__file__).resolve().parents[1]
    config_path = project_root / "config" / "model_parameters.yaml"

    # Load base parameters from YAML
    try:
        params_base = load_model_params_from_yaml(config_path)
    except Exception as exc:
        raise RuntimeError(
            f"Failed to load model parameters from {config_path}"
        ) from exc

    # ------------------------------------------------------------------
    # Hysteresis tuning:
    #
    # Make Ras a stronger gate:
    #   - increase rho_C (Ras -> C) and alpha_A (Ras -> A)
    # Slightly weaken self-sustaining feedback:
    #   - decrease eta_C (T  -> C) and eta_CM (M  -> C)
    #
    # This preserves the existence of a malignant high state but avoids
    # having it stable across the entire Ras axis.
    # ------------------------------------------------------------------
    params_hyst = replace(
        params_base,
        rho_C=params_base.rho_C * 6.0,
        alpha_A=params_base.alpha_A * 5.0,
        eta_C=params_base.eta_C * 0.6,
        eta_CM=params_base.eta_CM * 0.5,
    )

    # Ras range for scanning
    f_min = 0.0
    f_max = 2.0
    n_ras = 60

    ras_grid_up = np.linspace(f_min, f_max, n_ras)
    ras_grid_down_desc = np.linspace(f_max, f_min, n_ras)

    # ------------------------------------------------------------------
    # Up-sweep branch (benign -> malignant)
    # ------------------------------------------------------------------
    # Benign-like starting state:
    # [C, A, T, R, L, M]
    y0_up = np.array([0.01, 0.01, 0.01, 0.0, 0.01, 0.0], dtype=float)

    states_up = _run_branch(
        base_params=params_hyst,
        ras_values=ras_grid_up,
        y0_initial=y0_up,
    )

    # ------------------------------------------------------------------
    # Down-sweep branch (malignant -> benign)
    # ------------------------------------------------------------------
    # Malignant-like starting state:
    # [C, A, T, R, L, M]
    y0_down = np.array(
        [0.7, 0.7, 10.0, 0.5, 1.0, 0.5],
        dtype=float,
    )

    # Run branch in descending Ras
    states_down_desc = _run_branch(
        base_params=params_hyst,
        ras_values=ras_grid_down_desc,
        y0_initial=y0_down,
    )

    # Reverse to align with ascending Ras for plotting / CSV
    ras_grid_down = ras_grid_up.copy()
    states_down = states_down_desc[::-1, :]

    # ------------------------------------------------------------------
    # Save steady-state CSV
    # ------------------------------------------------------------------
    results_dir = project_root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    csv_path = results_dir / "ras_hysteresis_steady_states.csv"

    with csv_path.open("w", encoding="utf-8") as fh:
        fh.write("branch,f_ras,C_star,A_star,T_star,R_star,L_star,M_star\n")
        for f, y in zip(ras_grid_up, states_up):
            fh.write(
                f"up,{f:.6g},{y[0]:.6g},{y[1]:.6g},{y[2]:.6g},"
                f"{y[3]:.6g},{y[4]:.6g},{y[5]:.6g}\n"
            )
        for f, y in zip(ras_grid_down, states_down):
            fh.write(
                f"down,{f:.6g},{y[0]:.6g},{y[1]:.6g},{y[2]:.6g},"
                f"{y[3]:.6g},{y[4]:.6g},{y[5]:.6g}\n"
            )

    print(f"[INFO] Steady-state CSV written to {csv_path}")

    # ------------------------------------------------------------------
    # Build hysteresis figure
    # ------------------------------------------------------------------
    C_up = states_up[:, 0]
    A_up = states_up[:, 1]
    C_down = states_down[:, 0]
    A_down = states_down[:, 1]

    figures_dir = project_root / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    fig_path = figures_dir / "ras_hysteresis.png"

    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    axC, axA = axes

    # Identify bistable window based on separation between branches
    eps = 1e-3
    diff_mask = np.abs(C_up - C_down) > eps
    if diff_mask.any():
        f_bi_min = float(ras_grid_up[diff_mask].min())
        f_bi_max = float(ras_grid_up[diff_mask].max())
        for ax in axes:
            ax.axvspan(
                f_bi_min,
                f_bi_max,
                color="lightblue",
                alpha=0.25,
                zorder=0,
            )

    # Top panel: CSC fraction C*
    axC.plot(
        ras_grid_up,
        C_up,
        "o-",
        label="Up-sweep (benign → malignant)",
        linewidth=2.0,
        markersize=4,
    )
    axC.plot(
        ras_grid_down,
        C_down,
        "s-",
        label="Down-sweep (malignant → benign)",
        linewidth=2.0,
        markersize=4,
    )
    axC.set_ylabel("CSC steady state C*")
    axC.set_title("Ras hysteresis: CSC fraction and angiogenesis")
    axC.grid(alpha=0.3)
    axC.legend()

    # Bottom panel: Angiogenesis A*
    axA.plot(
        ras_grid_up,
        A_up,
        "o-",
        label="Up-sweep (benign → malignant)",
        linewidth=2.0,
        markersize=4,
    )
    axA.plot(
        ras_grid_down,
        A_down,
        "s-",
        label="Down-sweep (malignant → benign)",
        linewidth=2.0,
        markersize=4,
    )
    axA.set_xlabel("Ras input f_RAS")
    axA.set_ylabel("Angiogenesis steady state A*")
    axA.grid(alpha=0.3)
    axA.legend()

    plt.tight_layout()
    plt.savefig(fig_path, dpi=300)
    plt.close(fig)

    print(f"[INFO] Hysteresis figure saved to {fig_path}")


if __name__ == "__main__":
    main()
