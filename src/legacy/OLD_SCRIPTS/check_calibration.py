#!/usr/bin/env python3
"""
check_calibration.py

Check whether the Ras–CSC–microenvironment model reproduces the
experimental papilloma vs SCC fold-changes at two Ras levels.

This script plays the role of Step 6 in the "10-step recipe":
we test model performance against available data before leaning
on it for predictions.

The procedure is:

  1) Load the baseline parameter set from config/model_parameters.yaml.
  2) Define a low-Ras condition (papilloma-like) and a high-Ras
     condition (SCC-like) by setting f_ras to two values.
  3) For each Ras value, integrate the ODE system to steady state
     from an appropriate initial condition.
  4) Extract C*, A*, T*, and L* at each Ras level.
  5) Compute model SCC/papilloma ratios and compare them to the
     experimental ratios encoded in the YAML state_variables
     (C_scc / C_pap, etc.).
  6) Save a CSV and generate a simple bar plot for the paper.

"""

from __future__ import annotations

# Import standard library
from pathlib import Path
from dataclasses import replace
from typing import Tuple, Dict, Any

# Import numerical libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Import YAML to read experimental fold-changes
import yaml

# Import model components
from ras_csc_model import RasCSCModelParams, ras_csc_rhs, load_model_params_from_yaml


def integrate_to_steady_state(
    params: RasCSCModelParams,
    y0: np.ndarray,
    t_final: float = 300.0,
    n_points: int = 600,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Integrate the model to (approximate) steady state for a given
    parameter set and initial condition.

    Args:
        params:
            Parameter set including f_ras.
        y0:
            Initial state [C, A, T, R, L, M].
        t_final:
            Final integration time.
        n_points:
            Number of evaluation points for the trajectory.

    Returns:
        (t, y) where t is the time grid and y is the trajectory
        with shape (n_points, 6).
    """
    # Build time span and evaluation grid
    t_span = (0.0, float(t_final))
    t_eval = np.linspace(t_span[0], t_span[1], n_points)

    # Wrap the RHS so solve_ivp only sees (t, y)
    def rhs_wrapped(t: float, y: np.ndarray) -> np.ndarray:
        # Delegate to the shared ras_csc_rhs function
        return ras_csc_rhs(t, y, params)

    # Run the solver
    sol = solve_ivp(
        fun=rhs_wrapped,
        t_span=t_span,
        y0=np.asarray(y0, dtype=float),
        t_eval=t_eval,
        method="RK45",
        rtol=1e-6,
        atol=1e-9,
    )

    # Check success and raise if integration failed
    if not sol.success:
        raise RuntimeError(
            f"Integration failed in integrate_to_steady_state: {sol.message}"
        )

    # Return time and trajectory (transpose y to shape (n_points, 6))
    return sol.t, sol.y.T


def load_state_variables(config_path: Path) -> Dict[str, Any]:
    """
    Load state_variables from the YAML file so we can compute
    experimental fold-changes.

    Args:
        config_path:
            Path to model_parameters.yaml.

    Returns:
        The state_variables mapping as a dictionary.

    Raises:
        RuntimeError if the expected mapping is missing.
    """
    # Open YAML file
    with config_path.open("r", encoding="utf-8") as fh:
        config = yaml.safe_load(fh) or {}

    # Get the state_variables section
    state_vars = config.get("state_variables", None)

    # Ensure it is a dictionary
    if not isinstance(state_vars, dict):
        raise RuntimeError(
            "Expected a 'state_variables' mapping in the YAML for "
            "experimental fold-change comparison."
        )

    # Return the mapping
    return state_vars


def main() -> None:
    """
    Entry point for the calibration check.

    This function:

      1) Loads model parameters.
      2) Runs low-Ras and high-Ras simulations to steady state.
      3) Computes SCC/papilloma ratios for C, A, T, and L.
      4) Compares them with experimental ratios from the YAML.
      5) Saves a CSV and bar plot in results/ and figures/.
    """
    # Resolve project root as parent of src/
    project_root = Path(__file__).resolve().parents[1]

    # Build the path to the YAML configuration
    config_path = project_root / "config" / "model_parameters.yaml"

    # Load parameter set using existing helper
    params_base = load_model_params_from_yaml(config_path)

    # Load experimental state variables from YAML
    state_vars = load_state_variables(config_path)

    # Extract experimental papilloma and SCC values (all > 0)
    C_pap_exp = float(state_vars.get("C_pap", 1.0))
    C_scc_exp = float(state_vars.get("C_scc", 1.0))
    A_pap_exp = float(state_vars.get("A_pap", 1.0))
    A_scc_exp = float(state_vars.get("A_scc", 1.0))
    T_pap_exp = float(state_vars.get("T_pap", 1.0))
    T_scc_exp = float(state_vars.get("T_scc", 1.0))
    L_pap_exp = float(state_vars.get("L_pap", 1.0))
    L_scc_exp = float(state_vars.get("L_scc", 1.0))

    # Compute experimental SCC/papilloma ratios
    gamma_C_exp = C_scc_exp / C_pap_exp
    gamma_A_exp = A_scc_exp / A_pap_exp
    gamma_T_exp = T_scc_exp / T_pap_exp
    gamma_L_exp = L_scc_exp / L_pap_exp

    # Define Ras levels for pap and SCC
    f_ras_pap = 0.1
    f_ras_scc = 1.0

    # Set up benign-like initial condition
    y0_benign = np.array(
        [0.01, 0.01, 0.01, 0.0, 0.01, 0.0],
        dtype=float,
    )

    # Integrate low-Ras condition
    params_low = replace(params_base, f_ras=f_ras_pap)
    t_low, y_low = integrate_to_steady_state(params_low, y0_benign)

    # Extract steady state at end of trajectory
    C_pap_mod = float(y_low[-1, 0])
    A_pap_mod = float(y_low[-1, 1])
    T_pap_mod = float(y_low[-1, 2])
    L_pap_mod = float(y_low[-1, 4])

    # Use endpoint of low-Ras as initial condition for high-Ras run
    y0_high = y_low[-1, :].copy()

    # Integrate high-Ras condition
    params_high = replace(params_base, f_ras=f_ras_scc)
    t_high, y_high = integrate_to_steady_state(params_high, y0_high)

    # Extract high-Ras steady state
    C_scc_mod = float(y_high[-1, 0])
    A_scc_mod = float(y_high[-1, 1])
    T_scc_mod = float(y_high[-1, 2])
    L_scc_mod = float(y_high[-1, 4])

    # Compute model SCC/papilloma ratios
    gamma_C_mod = C_scc_mod / max(C_pap_mod, 1e-8)
    gamma_A_mod = A_scc_mod / max(A_pap_mod, 1e-8)
    gamma_T_mod = T_scc_mod / max(T_pap_mod, 1e-8)
    gamma_L_mod = L_scc_mod / max(L_pap_mod, 1e-8)

    # Create results directory if needed
    results_dir = project_root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    # Path for calibration CSV
    csv_path = results_dir / "calibration_check.csv"

    # Write summary CSV
    with csv_path.open("w", encoding="utf-8") as fh:
        fh.write("quantity,gamma_exp,gamma_model\n")
        fh.write(f"C,{gamma_C_exp:.6g},{gamma_C_mod:.6g}\n")
        fh.write(f"A,{gamma_A_exp:.6g},{gamma_A_mod:.6g}\n")
        fh.write(f"T,{gamma_T_exp:.6g},{gamma_T_mod:.6g}\n")
        fh.write(f"L,{gamma_L_exp:.6g},{gamma_L_mod:.6g}\n")

    # Prepare bar plot comparing experimental vs model ratios
    labels = ["C (CSC)", "A (Angio)", "T (TGFb)", "L (Leptin)"]
    exp_vals = [gamma_C_exp, gamma_A_exp, gamma_T_exp, gamma_L_exp]
    mod_vals = [gamma_C_mod, gamma_A_mod, gamma_T_mod, gamma_L_mod]

    # Choose x positions for bars
    x = np.arange(len(labels))
    width = 0.35

    # Create figures directory
    figures_dir = project_root / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Build figure
    plt.figure(figsize=(8, 5))
    plt.bar(x - width / 2, exp_vals, width=width, label="Experiment")
    plt.bar(x + width / 2, mod_vals, width=width, label="Model")

    # Label axes and ticks
    plt.xticks(x, labels, rotation=0)
    plt.ylabel("SCC / papilloma fold-change")
    plt.title("Calibration check: model steady states vs experimental ratios")
    plt.legend()
    plt.grid(alpha=0.3, axis="y")

    # Save figure
    fig_path = figures_dir / "calibration_check.png"
    plt.tight_layout()
    plt.savefig(fig_path, dpi=300)

    # Print summary to console
    print(f"[INFO] Calibration CSV written to {csv_path}")
    print(f"[INFO] Calibration figure saved to {fig_path}")


if __name__ == "__main__":
    main()
