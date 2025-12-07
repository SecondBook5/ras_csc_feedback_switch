#!/usr/bin/env python3
"""
prediction_perturbations.py

Run in silico perturbation experiments on the Ras–CSC–microenvironment
model to mirror key interventions in the Yuan et al. study:

  - Baseline (Lepr^ctrl, no drug)
  - LEPR knockout (Lepr^null)
  - mTOR inhibition (rapamycin-like)
  - Anti-angiogenic therapy (VEGFA inhibition)
  - Systemic leptin pump (high leptin)

For each condition, we:

  1) Specify a modified parameter set that approximates the biological
     intervention.
  2) Integrate from an SCC-like initial condition to steady state at
     f_ras = 1.0.
  3) Compute a simple "growth index" G* = C* × A* as a proxy for
     tumour burden (CSC fraction times angiogenesis).
  4) Save a CSV with C*, A*, M*, and G* for each condition.
  5) Save a bar plot showing growth indices relative to baseline.

This script feeds directly into Steps 9–10: making predictions for new
scenarios and comparing model fold-changes with experimental γ's.

"""

from __future__ import annotations

# Standard imports
from pathlib import Path
from dataclasses import replace
from typing import Dict, Tuple

# Numerical imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# YAML for mTOR/VEGFA/leptin fractions
import yaml

# Model imports
from ras_csc_model import RasCSCModelParams, ras_csc_rhs, load_model_params_from_yaml


def integrate_to_steady(
    params: RasCSCModelParams,
    y0: np.ndarray,
    t_final: float = 300.0,
    n_points: int = 600,
) -> np.ndarray:
    """
    Integrate the model and return the final state.

    Args:
        params:
            Parameter set including f_ras and any perturbation changes.
        y0:
            Initial condition.
        t_final:
            Final integration time.
        n_points:
            Number of evaluation points.

    Returns:
        Final state y(t_final) as a NumPy array length 6.
    """
    # Define time span and eval grid
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


def load_effect_fractions(config_path: Path) -> Dict[str, float]:
    """
    Load a few inhibitor / pump fractions from YAML so that the
    in silico perturbations share the same scale as the data.

    Returns:
        Dictionary with keys:
          - mtor_residual_fraction
          - gamma_V_growth_pump
          - lambda_lep_high
    """
    # Open YAML
    with config_path.open("r", encoding="utf-8") as fh:
        config = yaml.safe_load(fh) or {}

    # Extract relevant subsections
    pi3k_mtor = config.get("pi3k_mtor_effects", {}) or {}
    systemic_vegfa = config.get("systemic_vegfa_coupling", {}) or {}
    systemic_leptin = config.get("systemic_leptin_effects", {}) or {}

    # Read fractions with defaults
    mtor_residual = float(pi3k_mtor.get("mtor_residual_fraction", 0.32))
    gamma_V_pump = float(systemic_vegfa.get("gamma_V_growth_pump", 2.48))
    lambda_lep_high = float(systemic_leptin.get("lambda_lep_high", 2.86))

    return {
        "mtor_residual_fraction": mtor_residual,
        "gamma_V_growth_pump": gamma_V_pump,
        "lambda_lep_high": lambda_lep_high,
    }


def main() -> None:
    """
    Run perturbation simulations and save CSV and figures.
    """
    # Resolve project paths
    project_root = Path(__file__).resolve().parents[1]
    config_path = project_root / "config" / "model_parameters.yaml"

    # Load baseline params
    params_base = load_model_params_from_yaml(config_path)

    # Set Ras to high (SCC regime)
    params_base = replace(params_base, f_ras=1.0)

    # Load effect fractions from YAML
    fractions = load_effect_fractions(config_path)
    mtor_resid = fractions["mtor_residual_fraction"]

    # Define malignant-like starting condition
    y0 = np.array(
        [0.7, 0.7, 10.0, 0.5, 1.0, 0.5],
        dtype=float,
    )

    # Dictionary to store final states and growth indices
    results: Dict[str, Tuple[np.ndarray, float]] = {}

    # Condition 1: Baseline (Lepr^ctrl, no drug)
    y_baseline = integrate_to_steady(params_base, y0)
    G_baseline = float(y_baseline[0] * y_baseline[1])
    results["baseline"] = (y_baseline, G_baseline)

    # Condition 2: LEPR knockout (approximate Lepr^null)
    params_lepr_null = replace(params_base, eta_R=0.0)
    y_lepr_null = integrate_to_steady(params_lepr_null, y0)
    G_lepr_null = float(y_lepr_null[0] * y_lepr_null[1])
    results["lepr_null"] = (y_lepr_null, G_lepr_null)

    # Condition 3: mTOR inhibition (rapamycin-like)
    # Scale eta_M by the residual fraction to mimic partial inhibition
    eta_M_base = params_base.eta_M
    params_rapa = replace(params_base, eta_M=eta_M_base * mtor_resid)
    y_rapa = integrate_to_steady(params_rapa, y0)
    G_rapa = float(y_rapa[0] * y_rapa[1])
    results["rapamycin_like"] = (y_rapa, G_rapa)

    # Condition 4: Anti-angiogenic therapy
    # Decrease Ras→A and C→A couplings to mimic VEGFA blockade
    params_anti_vegfa = replace(
        params_base,
        alpha_A=params_base.alpha_A * 0.5,
        beta_A=params_base.beta_A * 0.5,
    )
    y_anti_vegfa = integrate_to_steady(params_anti_vegfa, y0)
    G_anti_vegfa = float(y_anti_vegfa[0] * y_anti_vegfa[1])
    results["anti_vegfa"] = (y_anti_vegfa, G_anti_vegfa)

    # Condition 5: High systemic leptin (pump)
    # Increase A→L coupling as a crude way to boost local leptin
    params_lep_high = replace(params_base, mu_L=params_base.mu_L * 2.0)
    y_lep_high = integrate_to_steady(params_lep_high, y0)
    G_lep_high = float(y_lep_high[0] * y_lep_high[1])
    results["leptin_high"] = (y_lep_high, G_lep_high)

    # Write CSV with all conditions
    results_dir = project_root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    csv_path = results_dir / "perturbation_predictions.csv"

    with csv_path.open("w", encoding="utf-8") as fh:
        fh.write("condition,C_star,A_star,T_star,R_star,L_star,M_star,G_star,fold_vs_baseline\n")
        for name, (y_star, G_star) in results.items():
            fold = G_star / max(G_baseline, 1e-8)
            fh.write(
                f"{name},{y_star[0]:.6g},{y_star[1]:.6g},{y_star[2]:.6g},"
                f"{y_star[3]:.6g},{y_star[4]:.6g},{y_star[5]:.6g},"
                f"{G_star:.6g},{fold:.6g}\n"
            )

    # Build bar plot of fold-changes in growth index
    labels = list(results.keys())
    folds = [
        results[name][1] / max(G_baseline, 1e-8)
        for name in labels
    ]

    figures_dir = project_root / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    fig_path = figures_dir / "perturbation_growth_folds.png"

    x = np.arange(len(labels))
    plt.figure(figsize=(8, 5))
    plt.bar(x, folds)
    plt.xticks(x, labels, rotation=45, ha="right")
    plt.ylabel("Growth index G* / G*_baseline")
    plt.title("In silico perturbation fold-changes in growth index")
    plt.grid(alpha=0.3, axis="y")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=300)

    # Report locations
    print(f"[INFO] Perturbation predictions CSV written to {csv_path}")
    print(f"[INFO] Perturbation growth figure saved to {fig_path}")


if __name__ == "__main__":
    main()
