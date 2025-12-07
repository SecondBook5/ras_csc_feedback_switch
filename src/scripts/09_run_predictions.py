#!/usr/bin/env python3
"""
run_predictions.py

Prediction and falsification experiments for the Ras–CSC feedback model.

This script corresponds to the "prediction" and "falsification" steps
in the 10-step systems biology workflow. It uses the calibrated Ras–CSC
ODE model to simulate specific perturbations that are biologically
interpretable and testable, then summarizes both numerical and graphical
outputs.

Why this script exists:
    - The model calibration step only ensures that the ODE system can
      match the training data at a few conditions.
    - To be scientifically useful, the model must also make concrete,
      testable predictions about new perturbations or regimes that were
      not used for calibration.
    - This script implements those perturbations in a single place,
      saves reproducible outputs, and prints a concise textual summary
      of the main findings.

What this script does:
    1. Loads calibrated parameters from results/calibration.

    2. Computes steady-state predictions under:
         - High Ras (SCC-like).
         - Lepr loss-of-function (LeprKO-like) at high Ras.
         - TGFβ neutralization at high Ras.
         - mTOR inhibition at high Ras.
         - Anti-angiogenic therapy at high Ras.

       It saves a CSV of steady-state values and prints clear
       statements about how each perturbation affects C, A, T, R, L, M.

    3. Runs two time-course "therapy" experiments:
         - Ras inhibition: start from high-Ras steady state, then drop
           Ras to a low value and integrate with time to test Ras
           independence of CSCs.
         - Anti-angiogenic therapy: start from high-Ras steady state
           and increase the angiogenesis degradation rate.

       It saves each time course to CSV and generates publication-ready
       figures for C(t), A(t), T(t), R(t), L(t), M(t).

    4. Optional bimodality check:
         - Loads CSC signature scores from the scRNA pipeline.
         - Fits one- and two-component Gaussian mixtures.
         - Compares BICs and prints whether there is evidence for a
           bimodal CSC distribution.
         - Generates a figure with the empirical density and the
           fitted two-component mixture, if scikit-learn is available.

Outputs
-------
Numerical:
    results/model/prediction_steady_states_CORRECTED.csv
    results/model/prediction_timecourse_ras_inhibition_CORRECTED.csv
    results/model/prediction_timecourse_antiangiogenic_CORRECTED.csv
    results/model/prediction_bimodality_summary_CORRECTED.json

Figures:
    figures/main/prediction_steady_state_bars_CORRECTED.png
    figures/main/prediction_ras_inhibition_timecourse_CORRECTED.png
    figures/main/prediction_antiangiogenic_timecourse_CORRECTED.png
    figures/main/prediction_csc_bimodality_CORRECTED.png   (optional)

Console:
    A short textual summary of the main qualitative predictions.
"""

from __future__ import annotations

# Import typing utilities for clearer type hints
from typing import Dict, Any, List, Tuple

# Import standard library utilities
import sys
import json
from pathlib import Path

# Import numerical libraries
import numpy as np
import pandas as pd

# Import plotting library
import matplotlib.pyplot as plt

# Import ODE solver
from scipy.integrate import odeint  # type: ignore

# Try to import GaussianMixture for the optional bimodality check
try:
    # Import GaussianMixture for 1D mixture fitting
    from sklearn.mixture import GaussianMixture  # type: ignore

    # Flag that mixture modeling is available
    HAS_SKLEARN: bool = True
except Exception:
    # If import fails, disable bimodality modeling
    HAS_SKLEARN = False

# ----------------------------------------------------------------------
# Ensure src/ is on the Python path and import core model utilities
# ----------------------------------------------------------------------

#   Get the directory of this script file
CURRENT_DIR: Path = Path(__file__).resolve().parent
#   The src root is the parent of the scripts directory
SRC_ROOT: Path = CURRENT_DIR.parent
#   The project root is the parent of src/
ROOT: Path = SRC_ROOT.parent

#   Add src root to sys.path if not already present
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

#   Import the ODE right-hand side and steady-state integrator
try:
    #   ras_csc_ode is the ODE RHS; simulate_steady_state integrates it
    from ras_csc_core import ras_csc_ode, simulate_steady_state, get_ras_input  # type: ignore
except ImportError as exc:  # pragma: no cover - defensive
    #   Raise a descriptive error if the core module cannot be imported
    raise ImportError(
        "Could not import 'ras_csc_ode', 'simulate_steady_state', or 'get_ras_input' "
        "from 'ras_csc_core'. Ensure ras_csc_core.py is in src/ and "
        "this script lives in src/scripts/."
    ) from exc

# ======================================================================
# PATH CONSTANTS
# ======================================================================

# Comment: determine script, src, and project root directories
CURRENT_DIR: Path = Path(__file__).resolve().parent
SRC_ROOT: Path = CURRENT_DIR.parent
ROOT: Path = SRC_ROOT.parent

# Comment: define locations of calibration results and processed data
CALIBRATION_DIR: Path = ROOT / "results" / "calibration"
PARAM_JSON: Path = CALIBRATION_DIR / "optimized_parameters_CORRECTED.json"

PROCESSED_DIR: Path = ROOT / "data" / "processed" / "omics_summaries"
RESULTS_MODEL_DIR: Path = ROOT / "results" / "model"
RESULTS_PRED_DIR: Path = ROOT / "results" / "predictions"
FIG_DIR: Path = ROOT / "figures" / "main"

# Comment: ensure that output directories exist
RESULTS_MODEL_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_PRED_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


# ======================================================================
# PARAMETER LOADING
# ======================================================================

def load_parameters() -> Dict[str, float]:
    """
    Load calibrated parameters from JSON and return as a float-valued dict.

    This centralizes parameter loading so that prediction experiments
    share exactly the same parameter set as the calibration and
    temporal scripts. The function checks that:
      - The JSON file exists.
      - The JSON is a non-empty dictionary.
      - Each value can be safely converted to a float.

    Returns
    -------
    Dict[str, float]
        Mapping from parameter names to numeric values.

    Raises
    ------
    FileNotFoundError
        If the calibration file does not exist.
    ValueError
        If the JSON is empty, malformed, or contains non-numeric values.
    """
    # Comment: check that the calibration JSON file exists
    if not PARAM_JSON.exists():
        raise FileNotFoundError(
            f"Calibrated parameters not found at {PARAM_JSON}. "
            "Run run_model_calibration.py before running predictions."
        )

    # Comment: open the JSON file and parse its contents
    with PARAM_JSON.open("r", encoding="utf-8") as fh:
        raw: Any = json.load(fh)

    # Comment: validate that the parsed structure is a non-empty dict
    if not isinstance(raw, dict) or not raw:
        raise ValueError(
            f"Parameter file {PARAM_JSON} is empty or not a JSON object."
        )

    # Comment: convert each parameter value to float defensively
    params: Dict[str, float] = {}
    for key, value in raw.items():
        try:
            params[str(key)] = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Parameter '{key}' with value '{value}' could not be converted to float."
            ) from exc

    # Comment: report how many parameters were successfully loaded
    print(
        f"[INFO] Loaded {len(params)} calibrated parameters from {PARAM_JSON}")
    return params


# ======================================================================
# TIME-COURSE SIMULATION
# ======================================================================

def simulate_time_course(
    f_ras: float,
    params: Dict[str, float],
    y0: List[float],
    t_max: float = 200.0,
    n_steps: int = 400,
    atol: float = 1e-6,
    rtol: float = 1e-6,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Integrate the Ras–CSC ODE system over time for a fixed Ras input.

    This function wraps scipy.integrate.odeint and the central ODE
    function from ras_csc_core. It returns the full trajectory so
    that temporal therapy responses can be analyzed.

    Parameters
    ----------
    f_ras : float
        Effective Ras input level for this simulation, typically
        between 0 and 1.
    params : Dict[str, float]
        Calibrated parameter dictionary.
    y0 : List[float]
        Initial state vector [C, A, T, R, L, M] at time 0.
    t_max : float, optional
        Final time point for integration.
    n_steps : int, optional
        Number of time points from 0 to t_max.
    atol : float, optional
        Absolute tolerance for the ODE solver.
    rtol : float, optional
        Relative tolerance for the ODE solver.

    Returns
    -------
    t : np.ndarray
        Time grid of shape (n_steps,).
    traj : np.ndarray
        State trajectory of shape (n_steps, 6) for [C, A, T, R, L, M].

    Notes
    -----
    This function does not attempt to detect steady state. It always
    returns the full time course for visualization and qualitative
    interpretation.
    """
    # Comment: build the time grid from 0 to t_max
    t: np.ndarray = np.linspace(0.0, float(t_max), int(n_steps))

    # Comment: cast the initial condition to a 1D float array
    y0_arr: np.ndarray = np.asarray(y0, dtype=float)
    if y0_arr.shape != (6,):
        raise ValueError(
            f"simulate_time_course expected y0 of length 6, got shape {y0_arr.shape}"
        )

    # Comment: integrate the ODEs using odeint with defensive error handling
    try:
        sol: np.ndarray = odeint(
            func=ras_csc_ode,
            y0=y0_arr,
            t=t,
            args=(float(f_ras), dict(params)),
            atol=float(atol),
            rtol=float(rtol),
        )
    except Exception as exc:  # pragma: no cover - defensive
        print(
            f"[ERROR] ODE integration failed in simulate_time_course "
            f"for f_ras={f_ras:.3f}: {exc}"
        )
        # Comment: if integration fails, return a flat trajectory
        fallback = np.tile(y0_arr, (t.size, 1))
        return t, fallback

    # Comment: handle any NaNs or infinities and clip states to [0, 1]
    if not np.all(np.isfinite(sol)):
        print(
            f"[WARN] Non-finite values in trajectory for f_ras={f_ras:.3f}; "
            "replacing NaNs and clipping."
        )
    sol = np.nan_to_num(sol, nan=0.0, posinf=1.0, neginf=0.0)
    sol = np.clip(sol, 0.0, 1.0)

    # Comment: return time grid and cleaned trajectory
    return t, sol


# ======================================================================
# PERTURBATION PARAMETER HELPERS
# ======================================================================

def make_lepr_ko_params(params: Dict[str, float]) -> Dict[str, float]:
    """
    Return a parameter set approximating Lepr loss-of-function.

    The Lepr knockout primarily affects receptor induction and
    stability. To capture this effect at the ODE level, this function
    scales down the TGFβ driven LEPR induction term and slightly
    increases LEPR decay.

    Parameters
    ----------
    params : Dict[str, float]
        Baseline calibrated parameters.

    Returns
    -------
    Dict[str, float]
        Modified parameter dictionary representing LeprKO.
    """
    # Comment: copy the original parameters to avoid mutating in place
    perturbed = dict(params)

    # Comment: reduce LEPR induction strength to a small fraction
    eta_R = float(perturbed.get("eta_R", 0.6))
    perturbed["eta_R"] = 0.05 * eta_R

    # Comment: increase LEPR degradation modestly
    delta_R = float(perturbed.get("delta_R", 0.6))
    perturbed["delta_R"] = 1.5 * delta_R

    # Comment: return the modified parameter set
    return perturbed


def make_tgfb_neutralization_params(params: Dict[str, float]) -> Dict[str, float]:
    """
    Return a parameter set approximating TGFβ neutralization.

    TGFβ neutralization is modeled here by reducing the production
    terms that generate TGFβ from angiogenesis and CSCs.

    Parameters
    ----------
    params : Dict[str, float]
        Baseline calibrated parameters.

    Returns
    -------
    Dict[str, float]
        Modified parameter dictionary representing TGFβ blockade.
    """
    # Comment: copy the parameter dictionary to avoid side effects
    perturbed = dict(params)

    # Comment: reduce angiogenesis and CSC contributions to TGFβ
    alpha_T = float(perturbed.get("alpha_T", 0.8))
    gamma_T = float(perturbed.get("gamma_T", 0.15))
    perturbed["alpha_T"] = 0.2 * alpha_T
    perturbed["gamma_T"] = 0.2 * gamma_T

    # Comment: return the modified parameter set
    return perturbed


def make_mtor_inhibition_params(params: Dict[str, float]) -> Dict[str, float]:
    """
    Return a parameter set approximating mTOR inhibition.

    mTOR inhibition is represented by weakening the leptin-LEPR driven
    activation of mTOR and slightly increasing its decay. This is a
    simple proxy for pharmacologic mTOR blockade.

    Parameters
    ----------
    params : Dict[str, float]
        Baseline calibrated parameters.

    Returns
    -------
    Dict[str, float]
        Modified parameter dictionary representing mTOR inhibition.
    """
    # Comment: copy parameters so we do not mutate the original dict
    perturbed = dict(params)

    # Comment: reduce mTOR activation sensitivity
    eta_M = float(perturbed.get("eta_M", 1.0))
    perturbed["eta_M"] = 0.3 * eta_M

    # Comment: increase mTOR degradation moderately
    delta_M = float(perturbed.get("delta_M", 0.6))
    perturbed["delta_M"] = 1.5 * delta_M

    # Comment: return the modified parameters
    return perturbed


def make_antiangiogenic_params(params: Dict[str, float]) -> Dict[str, float]:
    """
    Return a parameter set approximating anti-angiogenic therapy.

    Anti-angiogenic therapy is encoded as an increased degradation
    rate of the angiogenesis variable. This is a coarse model of
    vessel pruning and reduced vascular support.

    Parameters
    ----------
    params : Dict[str, float]
        Baseline calibrated parameters.

    Returns
    -------
    Dict[str, float]
        Modified parameter dictionary representing anti-angiogenic therapy.
    """
    # Comment: copy input parameters
    perturbed = dict(params)

    # Comment: increase angiogenesis degradation substantially
    delta_A = float(perturbed.get("delta_A", 0.6))
    perturbed["delta_A"] = 3.0 * delta_A

    # Comment: return the modified parameter set
    return perturbed


# ======================================================================
# STEADY-STATE PREDICTIONS
# ======================================================================

def compute_steady_state_predictions(params: Dict[str, float]) -> pd.DataFrame:
    """
    Compute steady-state predictions for several perturbation scenarios.

    Scenarios implemented:
        - baseline_high_ras: calibrated parameters at high Ras.
        - lepr_ko: Lepr loss-of-function at high Ras.
        - tgfb_neut: TGFβ neutralization at high Ras.
        - mtor_inhib: mTOR inhibition at high Ras.
        - anti_angio: anti-angiogenic therapy at high Ras.

    For each scenario, the function uses simulate_steady_state from
    ras_csc_core to obtain the approximate steady-state values of
    [C, A, T, R, L, M] and returns a tidy DataFrame.

    Parameters
    ----------
    params : Dict[str, float]
        Baseline calibrated parameters.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns:
            scenario, f_ras, C, A, T, R, L, M
    """
    # Comment: use a Ras level near the SCC mapping from get_ras_input
    f_ras_high: float = float(get_ras_input("SCC"))

    # Comment: benign-like initial condition for all steady-state runs
    y0_benign: List[float] = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

    # Comment: define all perturbation parameter sets
    scenarios: List[Tuple[str, Dict[str, float]]] = [
        ("baseline_high_ras", params),
        ("lepr_ko", make_lepr_ko_params(params)),
        ("tgfb_neut", make_tgfb_neutralization_params(params)),
        ("mtor_inhib", make_mtor_inhibition_params(params)),
        ("anti_angio", make_antiangiogenic_params(params)),
    ]

    # Comment: accumulate steady-state rows here
    rows: List[Dict[str, Any]] = []

    # Comment: compute steady state for each scenario
    for name, pset in scenarios:
        # Comment: integrate to approximate steady state
        y_ss = simulate_steady_state(
            f_ras=f_ras_high,
            params=pset,
            y0=y0_benign,
            t_max=200.0,
            n_steps=2000,
        )

        # Comment: append results as a tidy row
        rows.append(
            {
                "scenario": name,
                "f_ras": f_ras_high,
                "C": float(y_ss[0]),
                "A": float(y_ss[1]),
                "T": float(y_ss[2]),
                "R": float(y_ss[3]),
                "L": float(y_ss[4]),
                "M": float(y_ss[5]),
            }
        )

    # Comment: construct DataFrame from accumulated rows
    df = pd.DataFrame(rows)

    # Comment: save steady-state predictions to disk
    out_path = RESULTS_MODEL_DIR / "prediction_steady_states_CORRECTED.csv"
    df.to_csv(out_path, index=False)
    print(f"[SAVED] Steady-state predictions to {out_path}")

    # Comment: return DataFrame to caller for further use
    return df


def plot_steady_state_bars(df: pd.DataFrame) -> None:
    """
    Generate a bar plot comparing steady states across scenarios.

    The plot shows each variable (C, A, T, R, L, M) for each scenario
    as grouped bars so that qualitative differences between perturbations
    are visually obvious.

    Parameters
    ----------
    df : pd.DataFrame
        Steady-state prediction DataFrame from compute_steady_state_predictions.
    """
    # Comment: define variable order for plotting
    variables: List[str] = ["C", "A", "T", "R", "L", "M"]

    # Comment: extract scenario names in a stable order
    scenarios = list(df["scenario"].values)

    # Comment: set up bar positions
    n_vars = len(variables)
    n_scen = len(scenarios)
    x = np.arange(n_vars, dtype=float)
    width = 0.8 / max(n_scen, 1)

    # Comment: create the matplotlib figure and axes
    fig, ax = plt.subplots(figsize=(7.0, 4.5))

    # Comment: iterate over scenarios and draw bars
    for i, scen in enumerate(scenarios):
        # Comment: get the row for this scenario
        row = df[df["scenario"] == scen].iloc[0]
        values = [float(row[v]) for v in variables]
        # Comment: horizontal offset for this scenario
        offset = (i - (n_scen - 1) / 2.0) * width
        # Comment: plot the bars for this scenario
        ax.bar(x + offset, values, width=width, label=scen)

    # Comment: configure x-axis labels and limits
    ax.set_xticks(x)
    ax.set_xticklabels(variables)
    ax.set_ylabel("Steady-state value")
    ax.set_title("Steady-state predictions under model perturbations")

    # Comment: add legend and tight layout
    ax.legend(title="Scenario", fontsize=8)
    fig.tight_layout()

    # Comment: save the figure to disk
    fig_path = FIG_DIR / "prediction_steady_state_bars_CORRECTED.png"
    fig.savefig(fig_path, dpi=400)
    plt.close(fig)
    print(f"[SAVED] Steady-state prediction figure to {fig_path}")


# ======================================================================
# TIME-COURSE THERAPY EXPERIMENTS
# ======================================================================

def run_ras_inhibition_therapy(
    params: Dict[str, float],
) -> pd.DataFrame:
    """
    Simulate Ras inhibition therapy and save the resulting time course.

    The experiment:
        - Compute the high-Ras steady state.
        - Use that state as initial condition.
        - Lower Ras from high (SCC-like) to low (benign-like).
        - Integrate the ODEs and record the full trajectory.

    The main question:
        Does C(t) collapse to a low value or remain high despite
        Ras inhibition, indicating Ras independence?

    Parameters
    ----------
    params : Dict[str, float]
        Baseline calibrated parameters.

    Returns
    -------
    pd.DataFrame
        Time-course DataFrame with columns:
            time, C, A, T, R, L, M, f_ras
    """
    # Comment: get Ras inputs for SCC and Normal conditions
    f_ras_high: float = float(get_ras_input("SCC"))
    f_ras_low: float = float(get_ras_input("Normal"))

    # Comment: compute high-Ras steady state to represent pre-therapy tumor
    y_ss_high = simulate_steady_state(
        f_ras=f_ras_high,
        params=params,
        y0=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
        t_max=200.0,
        n_steps=2000,
    )

    # Comment: simulate time course after Ras inhibition using low Ras
    t, traj = simulate_time_course(
        f_ras=f_ras_low,
        params=params,
        y0=list(y_ss_high.tolist()),
        t_max=200.0,
        n_steps=400,
    )

    # Comment: build the time-course DataFrame
    df = pd.DataFrame(
        {
            "time": t.astype(float),
            "C": traj[:, 0].astype(float),
            "A": traj[:, 1].astype(float),
            "T": traj[:, 2].astype(float),
            "R": traj[:, 3].astype(float),
            "L": traj[:, 4].astype(float),
            "M": traj[:, 5].astype(float),
            "f_ras": float(f_ras_low),
        }
    )

    # Comment: save the time course to disk
    out_path = RESULTS_MODEL_DIR / "prediction_timecourse_ras_inhibition_CORRECTED.csv"
    df.to_csv(out_path, index=False)
    print(f"[SAVED] Ras inhibition time course to {out_path}")

    # Comment: return the time-course DataFrame
    return df


def plot_ras_inhibition_therapy(df: pd.DataFrame) -> None:
    """
    Plot the Ras inhibition therapy time course for key variables.

    The figure contains a panel showing the trajectories of C, A, T, R,
    L, and M over time after Ras is lowered. The qualitative shape
    indicates whether the CSC population collapses or persists.

    Parameters
    ----------
    df : pd.DataFrame
        Time-course DataFrame from run_ras_inhibition_therapy.
    """
    # Comment: create a figure and axes
    fig, ax = plt.subplots(figsize=(7.0, 4.5))

    # Comment: plot each variable as a time series
    ax.plot(df["time"], df["C"], label="C (CSC fraction)")
    ax.plot(df["time"], df["A"], label="A (Angio)")
    ax.plot(df["time"], df["T"], label="T (TGFβ)")
    ax.plot(df["time"], df["R"], label="R (Lepr)")
    ax.plot(df["time"], df["L"], label="L (Leptin)")
    ax.plot(df["time"], df["M"], label="M (mTOR)")

    # Comment: label axes and add title
    ax.set_xlabel("Time")
    ax.set_ylabel("State value")
    ax.set_title("Ras inhibition therapy: response of Ras–CSC system")

    # Comment: add legend and layout
    ax.legend(fontsize=8)
    fig.tight_layout()

    # Comment: save figure to disk
    fig_path = FIG_DIR / "prediction_ras_inhibition_timecourse_CORRECTED.png"
    fig.savefig(fig_path, dpi=400)
    plt.close(fig)
    print(f"[SAVED] Ras inhibition time-course figure to {fig_path}")


def run_antiangiogenic_therapy(
    params: Dict[str, float],
) -> pd.DataFrame:
    """
    Simulate anti-angiogenic therapy and save the resulting time course.

    The experiment:
        - Compute the high-Ras steady state under baseline parameters.
        - Use that state as initial condition.
        - Switch to a parameter set with increased angiogenesis
          degradation (anti-angiogenic).
        - Integrate the ODEs forward in time with high Ras kept on.

    The main question:
        Is vessel pruning sufficient to collapse the CSC compartment
        despite persistent Ras activation?

    Parameters
    ----------
    params : Dict[str, float]
        Baseline calibrated parameters.

    Returns
    -------
    pd.DataFrame
        Time-course DataFrame with columns:
            time, C, A, T, R, L, M, f_ras
    """
    # Comment: high Ras input for SCC-like regime
    f_ras_high: float = float(get_ras_input("SCC"))

    # Comment: pre-therapy high-Ras steady state
    y_ss_high = simulate_steady_state(
        f_ras=f_ras_high,
        params=params,
        y0=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
        t_max=200.0,
        n_steps=2000,
    )

    # Comment: perturbation parameters representing anti-angiogenic therapy
    params_anti = make_antiangiogenic_params(params)

    # Comment: simulate time course with anti-angiogenic parameters
    t, traj = simulate_time_course(
        f_ras=f_ras_high,
        params=params_anti,
        y0=list(y_ss_high.tolist()),
        t_max=200.0,
        n_steps=400,
    )

    # Comment: assemble the time-course DataFrame
    df = pd.DataFrame(
        {
            "time": t.astype(float),
            "C": traj[:, 0].astype(float),
            "A": traj[:, 1].astype(float),
            "T": traj[:, 2].astype(float),
            "R": traj[:, 3].astype(float),
            "L": traj[:, 4].astype(float),
            "M": traj[:, 5].astype(float),
            "f_ras": float(f_ras_high),
        }
    )

    # Comment: save the time course to CSV
    out_path = RESULTS_MODEL_DIR / "prediction_timecourse_antiangiogenic_CORRECTED.csv"
    df.to_csv(out_path, index=False)
    print(f"[SAVED] Anti-angiogenic time course to {out_path}")

    # Comment: return DataFrame
    return df


def plot_antiangiogenic_therapy(df: pd.DataFrame) -> None:
    """
    Plot the anti-angiogenic therapy time course for key variables.

    The figure shows how C, A, T, R, L, and M evolve after an increase
    in angiogenesis degradation while Ras remains high.

    Parameters
    ----------
    df : pd.DataFrame
        Time-course DataFrame from run_antiangiogenic_therapy.
    """
    # Comment: create figure and axes
    fig, ax = plt.subplots(figsize=(7.0, 4.5))

    # Comment: plot each variable as a time series
    ax.plot(df["time"], df["C"], label="C (CSC fraction)")
    ax.plot(df["time"], df["A"], label="A (Angio)")
    ax.plot(df["time"], df["T"], label="T (TGFβ)")
    ax.plot(df["time"], df["R"], label="R (Lepr)")
    ax.plot(df["time"], df["L"], label="L (Leptin)")
    ax.plot(df["time"], df["M"], label="M (mTOR)")

    # Comment: label axes and add title
    ax.set_xlabel("Time")
    ax.set_ylabel("State value")
    ax.set_title("Anti-angiogenic therapy: response of Ras–CSC system")

    # Comment: add legend and tight layout
    ax.legend(fontsize=8)
    fig.tight_layout()

    # Comment: save figure to disk
    fig_path = FIG_DIR / "prediction_antiangiogenic_timecourse_CORRECTED.png"
    fig.savefig(fig_path, dpi=400)
    plt.close(fig)
    print(f"[SAVED] Anti-angiogenic time-course figure to {fig_path}")


# ======================================================================
# OPTIONAL CSC BIMODALITY CHECK
# ======================================================================

def run_csc_bimodality_check() -> Dict[str, Any]:
    """
    Optional step: test for bimodality in CSC signature scores.

    This function:
        - Loads per-cell module scores and metadata from the scRNA
          pipeline output.
        - Extracts the Yuan CSC signature score.
        - Fits one- and two-component Gaussian mixtures when
          scikit-learn is available.
        - Compares BIC between the one- and two-component fits.
        - Saves a small JSON summary and an optional figure with the
          empirical density and the fitted two-component mixture.

    Returns
    -------
    Dict[str, Any]
        Dictionary summarizing:
            - has_sklearn
            - n_cells
            - bic_1
            - bic_2
            - delta_bic
            - evidence_for_bimodality
    """
    # Comment: define paths to scRNA module scores and metadata
    scores_path = PROCESSED_DIR / "scc_scRNA_module_scores_per_cell.csv"
    meta_path = PROCESSED_DIR / "scc_scRNA_K14pos_metadata_with_CSC_labels.csv"

    # Comment: check that required files exist
    if not scores_path.exists() or not meta_path.exists():
        print(
            "[WARN] scRNA outputs not found for bimodality check. "
            "Skipping CSC bimodality analysis."
        )
        return {
            "has_sklearn": HAS_SKLEARN,
            "n_cells": 0,
            "bic_1": None,
            "bic_2": None,
            "delta_bic": None,
            "evidence_for_bimodality": False,
        }

    # Comment: load both tables
    scores_df = pd.read_csv(scores_path)
    meta_df = pd.read_csv(meta_path)

    # Comment: merge by cell_id to align cells
    if "cell_id" not in scores_df.columns or "cell_id" not in meta_df.columns:
        print(
            "[WARN] Missing 'cell_id' column in scRNA outputs. "
            "Skipping CSC bimodality analysis."
        )
        return {
            "has_sklearn": HAS_SKLEARN,
            "n_cells": 0,
            "bic_1": None,
            "bic_2": None,
            "delta_bic": None,
            "evidence_for_bimodality": False,
        }

    merged = pd.merge(scores_df, meta_df, on="cell_id", how="inner")

    # Comment: verify CSC signature column exists
    if "CSC_signature_score" not in merged.columns:
        print(
            "[WARN] CSC_signature_score not found in scRNA outputs. "
            "Skipping CSC bimodality analysis."
        )
        return {
            "has_sklearn": HAS_SKLEARN,
            "n_cells": 0,
            "bic_1": None,
            "bic_2": None,
            "delta_bic": None,
            "evidence_for_bimodality": False,
        }

    # Comment: extract CSC signature values and drop non-finite entries
    x = merged["CSC_signature_score"].to_numpy(dtype=float)
    x = x[np.isfinite(x)]
    n_cells = int(x.size)

    # Comment: handle trivial cases with very few cells
    if n_cells < 50:
        print(
            f"[WARN] Only {n_cells} cells for CSC bimodality check. "
            "Too few for a stable mixture fit."
        )
        return {
            "has_sklearn": HAS_SKLEARN,
            "n_cells": n_cells,
            "bic_1": None,
            "bic_2": None,
            "delta_bic": None,
            "evidence_for_bimodality": False,
        }

    # Comment: if sklearn is not available, we cannot fit mixtures
    if not HAS_SKLEARN:
        print(
            "[WARN] scikit-learn not available. "
            "Cannot perform Gaussian mixture-based bimodality check."
        )
        return {
            "has_sklearn": False,
            "n_cells": n_cells,
            "bic_1": None,
            "bic_2": None,
            "delta_bic": None,
            "evidence_for_bimodality": False,
        }

    # Comment: reshape to column vector for GaussianMixture
    X = x.reshape(-1, 1)

    # Comment: fit a one-component Gaussian mixture
    gm1 = GaussianMixture(
        n_components=1, covariance_type="full", random_state=1)
    gm1.fit(X)
    bic_1 = float(gm1.bic(X))

    # Comment: fit a two-component Gaussian mixture
    gm2 = GaussianMixture(
        n_components=2, covariance_type="full", random_state=1)
    gm2.fit(X)
    bic_2 = float(gm2.bic(X))

    # Comment: compute BIC difference (positive means 2-component preferred)
    delta_bic = bic_1 - bic_2

    # Comment: simple heuristic for "evidence for bimodality"
    evidence = bool(delta_bic > 10.0)

    # Comment: summarize results in a dictionary
    summary = {
        "has_sklearn": True,
        "n_cells": n_cells,
        "bic_1": bic_1,
        "bic_2": bic_2,
        "delta_bic": delta_bic,
        "evidence_for_bimodality": evidence,
    }

    # Comment: save JSON summary
    json_path = RESULTS_PRED_DIR / "prediction_bimodality_summary_CORRECTED.json"
    with json_path.open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[SAVED] CSC bimodality summary to {json_path}")

    # Comment: generate and save a density plus mixture figure
    try:
        # Comment: create a dense grid for plotting
        xs = np.linspace(x.min(), x.max(), 400).reshape(-1, 1)
        pdf2 = np.exp(gm2.score_samples(xs))

        fig, ax = plt.subplots(figsize=(5.0, 4.0))

        # Comment: plot empirical histogram
        ax.hist(
            x,
            bins=40,
            density=True,
            alpha=0.4,
            label="Empirical CSC signature",
        )

        # Comment: overlay fitted two-component mixture density
        ax.plot(xs[:, 0], pdf2, linewidth=2.0, label="Two-component GMM")

        # Comment: label axes and add title
        ax.set_xlabel("CSC signature score (Yuan)")
        ax.set_ylabel("Density")
        ax.set_title("CSC signature distribution and two-component GMM")

        # Comment: add legend and layout
        ax.legend(fontsize=8)
        fig.tight_layout()

        # Comment: save figure to disk
        fig_path = FIG_DIR / "prediction_csc_bimodality_CORRECTED.png"
        fig.savefig(fig_path, dpi=400)
        plt.close(fig)
        print(f"[SAVED] CSC bimodality figure to {fig_path}")
    except Exception as exc:  # pragma: no cover - defensive
        print(f"[WARN] Failed to generate CSC bimodality figure: {exc}")

    # Comment: return summary dictionary to caller
    return summary


# ======================================================================
# TEXTUAL SUMMARY OF PREDICTIONS
# ======================================================================

def summarize_predictions_console(
    steady_df: pd.DataFrame,
    ras_df: pd.DataFrame,
    anti_df: pd.DataFrame,
    bimodal_summary: Dict[str, Any],
) -> None:
    """
    Print a concise textual summary of the main model predictions.

    This function translates numerical differences into human-readable
    statements that can be copied into the manuscript as the narrative
    around the prediction and falsification steps.

    Parameters
    ----------
    steady_df : pd.DataFrame
        Steady-state prediction DataFrame.
    ras_df : pd.DataFrame
        Ras inhibition therapy time-course DataFrame.
    anti_df : pd.DataFrame
        Anti-angiogenic therapy time-course DataFrame.
    bimodal_summary : Dict[str, Any]
        Summary dictionary from run_csc_bimodality_check.
    """
    # Comment: helper to compute qualitative change between scenarios
    def describe_change(base: float, pert: float, name: str) -> str:
        # Comment: compute absolute difference
        diff = pert - base
        # Comment: threshold for calling something "changed"
        if abs(diff) < 0.02:
            return f"{name}: essentially unchanged (Δ ≈ {diff:.2f})"
        direction = "increased" if diff > 0.0 else "decreased"
        return f"{name}: {direction} by {diff:.2f} (from {base:.2f} to {pert:.2f})"

    print("\n" + "=" * 70)
    print("SUMMARY OF MODEL PREDICTIONS")
    print("=" * 70)

    # Comment: extract baseline steady state
    base = steady_df[steady_df["scenario"] == "baseline_high_ras"].iloc[0]

    # Comment: report each perturbation relative to baseline for C and M
    for scen in ["lepr_ko", "tgfb_neut", "mtor_inhib", "anti_angio"]:
        row = steady_df[steady_df["scenario"] == scen].iloc[0]
        print(f"\nScenario: {scen}")
        print(
            "  " + describe_change(float(base["C"]), float(row["C"]), "CSC fraction C"))
        print(
            "  " + describe_change(float(base["M"]), float(row["M"]), "mTOR activity M"))
        print(
            "  " + describe_change(float(base["A"]), float(row["A"]), "Angiogenesis A"))
        print(
            "  " + describe_change(float(base["T"]), float(row["T"]), "TGFβ T"))

    # Comment: interpret Ras inhibition therapy trajectory for C
    C_start = float(ras_df["C"].iloc[0])
    C_end = float(ras_df["C"].iloc[-1])
    print("\nRas inhibition therapy:")
    print(f"  C starts at {C_start:.2f} and ends at {C_end:.2f}.")
    if C_end > 0.5 * C_start:
        print("  Interpretation: CSCs remain partially high after Ras inhibition "
              "(model predicts partial Ras independence).")
    else:
        print("  Interpretation: CSCs collapse strongly after Ras inhibition "
              "(model predicts Ras dependence).")

    # Comment: interpret anti-angiogenic therapy trajectory for C
    C_start_anti = float(anti_df["C"].iloc[0])
    C_end_anti = float(anti_df["C"].iloc[-1])
    print("\nAnti-angiogenic therapy:")
    print(f"  C starts at {C_start_anti:.2f} and ends at {C_end_anti:.2f}.")
    if C_end_anti < 0.5 * C_start_anti:
        print("  Interpretation: CSCs collapse when angiogenesis is pruned "
              "even with Ras kept high.")
    else:
        print("  Interpretation: CSCs do not fully collapse under "
              "anti-angiogenic therapy alone.")

    # Comment: summarize bimodality evidence if available
    print("\nCSC bimodality check (optional):")
    if bimodal_summary.get("n_cells", 0) == 0:
        print("  Bimodality could not be evaluated due to missing or insufficient data.")
    elif not bimodal_summary.get("has_sklearn", False):
        print("  Bimodality analysis requires scikit-learn, which is not available.")
    else:
        delta_bic = bimodal_summary.get("delta_bic", 0.0)
        evidence = bimodal_summary.get("evidence_for_bimodality", False)
        print(f"  Cells used: {bimodal_summary['n_cells']}")
        print(f"  BIC(1-component) = {bimodal_summary['bic_1']:.1f}")
        print(f"  BIC(2-component) = {bimodal_summary['bic_2']:.1f}")
        print(f"  ΔBIC = BIC1 - BIC2 = {delta_bic:.1f}")
        if evidence:
            print("  Interpretation: strong evidence for a bimodal CSC signature "
                  "distribution in scRNA data.")
        else:
            print("  Interpretation: no strong evidence for a clean bimodal CSC "
                  "signature distribution by this simple GMM criterion.")


# ======================================================================
# MAIN ORCHESTRATION
# ======================================================================

def main() -> None:
    """
    Orchestrate prediction experiments and optional bimodality check.

    Steps:
        1. Load calibrated parameters.
        2. Compute steady-state predictions under several perturbations.
        3. Run Ras inhibition therapy time course and plot it.
        4. Run anti-angiogenic therapy time course and plot it.
        5. Generate a grouped bar plot for steady-state scenarios.
        6. Perform optional CSC bimodality analysis from scRNA data.
        7. Print a concise textual summary of all findings.
    """
    # Comment: section header for the console
    print("\n" + "=" * 70)
    print("RUNNING RAS–CSC MODEL PREDICTION EXPERIMENTS")
    print("=" * 70)

    # Comment: load calibrated parameters
    params = load_parameters()

    # Comment: compute and save steady-state predictions
    print("\n[STEP] Computing steady-state predictions...")
    steady_df = compute_steady_state_predictions(params)
    plot_steady_state_bars(steady_df)

    # Comment: simulate Ras inhibition therapy
    print("\n[STEP] Simulating Ras inhibition therapy...")
    ras_df = run_ras_inhibition_therapy(params)
    plot_ras_inhibition_therapy(ras_df)

    # Comment: simulate anti-angiogenic therapy
    print("\n[STEP] Simulating anti-angiogenic therapy...")
    anti_df = run_antiangiogenic_therapy(params)
    plot_antiangiogenic_therapy(anti_df)

    # Comment: run optional CSC bimodality analysis
    print("\n[STEP] Optional CSC bimodality analysis from scRNA data...")
    bimodal_summary = run_csc_bimodality_check()

    # Comment: print a concise textual summary of the predictions
    summarize_predictions_console(
        steady_df=steady_df,
        ras_df=ras_df,
        anti_df=anti_df,
        bimodal_summary=bimodal_summary,
    )

    # Comment: footer message
    print("\n" + "=" * 70)
    print("PREDICTION AND FALSIFICATION STEP COMPLETE")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    # Comment: invoke main when the script is executed directly
    main()
