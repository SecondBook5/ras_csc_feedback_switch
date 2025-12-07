#!/usr/bin/env python3
"""
sensitivity_scan.py

Perform simple one-parameter sensitivity analysis around the malignant
steady state of the Ras–CSC–microenvironment model.

This script addresses Step 7 ("Refine the model, estimate parameters")
by quantifying how the malignant steady-state CSC fraction C* responds
to perturbations in a selected subset of parameters:

    - eta_CM (M → C feedback strength)
    - eta_C  (T → C strength)
    - beta_A (C → A strength)
    - alpha_T (A → T strength)
    - eta_R  (T → R strength)

For each parameter p, we:

  1) Compute a baseline malignant steady state at f_ras = 1.0.
  2) Scale p by factors [0.5, 1.0, 1.5] and re-integrate.
  3) Record C* for each setting.
  4) Compute local sensitivities (delta C* for ±50%).

Outputs:

  - results/sensitivity_scan.csv with C* for each parameter and factor.
  - figures/sensitivity_curves.png showing C* vs scaling factor.
  - figures/sensitivity_tornado.png showing |ΔC*| for the ±50% changes.

"""

from __future__ import annotations

# Standard imports
from pathlib import Path
from dataclasses import replace
from typing import Dict, List, Tuple

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
) -> np.ndarray:
    """
    Integrate the model and return the final state as an approximate
    steady state.

    Args:
        params:
            Parameter set including f_ras.
        y0:
            Initial condition [C, A, T, R, L, M].
        t_final:
            Final integration time.
        n_points:
            Number of evaluation points.

    Returns:
        Final state y(t_final) as a NumPy array of length 6.
    """
    # Define time span and evaluation grid
    t_span = (0.0, float(t_final))
    t_eval = np.linspace(t_span[0], t_span[1], n_points)

    # Wrap RHS
    def rhs_wrapped(t: float, y: np.ndarray) -> np.ndarray:
        return ras_csc_rhs(t, y, params)

    # Integrate
    sol = solve_ivp(
        fun=rhs_wrapped,
        t_span=t_span,
        y0=np.asarray(y0, dtype=float),
        t_eval=t_eval,
        method="RK45",
        rtol=1e-6,
        atol=1e-9,
    )

    # Check success
    if not sol.success:
        raise RuntimeError(
            f"Integration failed in integrate_to_steady: {sol.message}"
        )

    # Return final state
    return sol.y[:, -1]


def main() -> None:
    """
    Perform sensitivity scans and save CSV and figures.
    """
    # Resolve paths
    project_root = Path(__file__).resolve().parents[1]
    config_path = project_root / "config" / "model_parameters.yaml"

    # Load baseline params
    params_base = load_model_params_from_yaml(config_path)

    # Make sure we are in malignant regime (high Ras)
    params_malig = replace(params_base, f_ras=1.0)

    # Start from malignant-like initial condition
    y0_malig = np.array([0.7, 0.7, 10.0, 0.5, 1.0, 0.5], dtype=float)

    # Compute baseline steady state
    y_star_base = integrate_to_steady(params_malig, y0_malig)
    C_star_base = float(y_star_base[0])

    # Parameters to scan and their attribute names on RasCSCModelParams
    targets: Dict[str, str] = {
        "eta_CM": "eta_CM",
        "eta_C": "eta_C",
        "beta_A": "beta_A",
        "alpha_T": "alpha_T",
        "eta_R": "eta_R",
    }

    # Scaling factors
    factors = np.array([0.5, 1.0, 1.5], dtype=float)

    # Storage for CSV and tornado
    records: List[Tuple[str, float, float]] = []
    tornado_vals: Dict[str, float] = {}

    # Loop over each target parameter
    for label, attr in targets.items():
        # Get baseline value
        base_val = getattr(params_malig, attr)

        # Storage for this parameter's C* values
        C_values_for_param: List[float] = []

        # Loop over scaling factors
        for factor in factors:
            # Build modified parameter set
            new_val = base_val * float(factor)
            params_mod = replace(params_malig, **{attr: new_val})

            # Integrate from the baseline steady state to speed convergence
            y_star_mod = integrate_to_steady(params_mod, y_star_base.copy())
            C_star_mod = float(y_star_mod[0])

            # Record C* for this (parameter, factor)
            C_values_for_param.append(C_star_mod)
            records.append((label, float(factor), C_star_mod))

        # For tornado: approximate sensitivity as |C(1.5) - C(0.5)|
        delta = abs(C_values_for_param[-1] - C_values_for_param[0])
        tornado_vals[label] = delta

    # Write CSV
    results_dir = project_root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    csv_path = results_dir / "sensitivity_scan.csv"

    with csv_path.open("w", encoding="utf-8") as fh:
        fh.write("parameter,scale_factor,C_star\n")
        for param, factor, C_star in records:
            fh.write(f"{param},{factor:.3g},{C_star:.6g}\n")

    # Build sensitivity curve figure
    figures_dir = project_root / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    curves_fig = figures_dir / "sensitivity_curves.png"

    plt.figure(figsize=(8, 6))
    for label in targets.keys():
        # Collect C* for this parameter
        C_param = [
            r[2] for r in records if r[0] == label
        ]
        plt.plot(
            factors,
            C_param,
            marker="o",
            label=label,
        )

    plt.axhline(
        y=C_star_base,
        color="black",
        linestyle=":",
        linewidth=1.0,
        label="Baseline C*",
    )
    plt.xlabel("Parameter scaling factor")
    plt.ylabel("CSC steady state C*")
    plt.title("One-parameter sensitivity of malignant CSC steady state")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(curves_fig, dpi=300)

    # Build tornado-style bar plot
    tornado_fig = figures_dir / "sensitivity_tornado.png"
    labels = list(tornado_vals.keys())
    deltas = [tornado_vals[k] for k in labels]

    # Sort by magnitude
    order = np.argsort(deltas)[::-1]
    labels_sorted = [labels[i] for i in order]
    deltas_sorted = [deltas[i] for i in order]

    plt.figure(figsize=(6, 4))
    y_pos = np.arange(len(labels_sorted))
    plt.barh(y_pos, deltas_sorted)
    plt.yticks(y_pos, labels_sorted)
    plt.xlabel("|ΔC*| for ±50% parameter change")
    plt.title("Local sensitivity of malignant CSC steady state")
    plt.grid(alpha=0.3, axis="x")
    plt.tight_layout()
    plt.savefig(tornado_fig, dpi=300)

    # Console report
    print(f"[INFO] Sensitivity CSV written to {csv_path}")
    print(f"[INFO] Sensitivity curves figure saved to {curves_fig}")
    print(f"[INFO] Sensitivity tornado figure saved to {tornado_fig}")


if __name__ == "__main__":
    main()
