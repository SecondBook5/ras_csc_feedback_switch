#!/usr/bin/env python3
"""
fit_ras_csc_model.py

Calibrate Ras–CSC–microenvironment feedback model parameters to
RNA-Seq-derived module scores.

This script:
  1) Loads module scores per sample from:
       data/processed/rnaseq/module_scores_by_sample.csv
  2) Aggregates to dataset × condition × module means.
  3) Interprets the modules as:
       - TGFb_bulk  ~ T (TGFβ)
       - Angio_bulk ~ A (angiogenesis)
       - CSC_bulk   ~ C (CSC activity)
       - mTOR_bulk  ~ M (mTOR activity)
  4) Uses a coarse random search over a limited parameter space to
     find a single global parameter set that roughly matches the
     relative pattern of module means across conditions.
  5) Exports:
       data/processed/model_fits/model_calibration_results.json
     with keys:
       - model_A_params: dict of RasCSCParams (feedback intact)
       - model_B_params: same but k_C_M = 0 (feedback removed)

Important:
  - This is NOT a full Bayesian or MCMC calibration.
  - It is a reproducible, explicit quantitative fit that you can
    tighten or replace later.
"""

from __future__ import annotations

# Import standard libraries for filesystem and types
import json
from dataclasses import asdict
from pathlib import Path
from typing import Dict, Tuple

# Import numerical libraries
import numpy as np
import pandas as pd

# Import model definitions from the calibration model file
from ras_csc_model_calib import RasCSCParams, ras_csc_rhs, simulate_trajectory


#--------------------------------------------------------------------
# Paths and constants
#--------------------------------------------------------------------

# Create a Path object for the project root relative paths
PROJECT_ROOT: Path = Path(".").resolve()

# Path to the module scores CSV produced by 03_rnaseq_compute_module_scores.R
MODULE_SCORES_PATH: Path = PROJECT_ROOT / "data" / "processed" / "rnaseq" / "module_scores_by_sample.csv"

# Path to the output JSON for model calibration
CALIB_OUTPUT_DIR: Path = PROJECT_ROOT / "data" / "processed" / "model_fits"
CALIB_OUTPUT_PATH: Path = CALIB_OUTPUT_DIR / "model_calibration_results.json"

# Fixed time span and step for steady-state simulation
T_SPAN: Tuple[float, float] = (0.0, 100.0)
DT: float = 0.1

# Number of random parameter samples to try in this coarse search
N_SAMPLES: int = 500  # You can increase to 2000+ for a tighter fit


#--------------------------------------------------------------------
# Utility / helper functions
#--------------------------------------------------------------------

def load_module_scores(path: Path) -> pd.DataFrame:
    """
    Load module scores and apply basic sanity checks.

    The CSV is expected to have columns:
      - sample_id
      - dataset
      - condition
      - TGFb_bulk
      - mTOR_bulk
      - Angio_bulk
      - CSC_bulk
    """
    # Ensure the file exists before trying to read
    if not path.exists():
        # Raise a clear error if the file is missing
        raise FileNotFoundError(f"[ERROR] Module scores file not found: {path}")

    # Read the CSV into a DataFrame
    df = pd.read_csv(path)

    # Define the required columns we expect
    required_cols = {
        "sample_id",
        "dataset",
        "condition",
        "TGFb_bulk",
        "mTOR_bulk",
        "Angio_bulk",
        "CSC_bulk",
    }

    # Check that all required columns are present
    missing = required_cols - set(df.columns)
    if missing:
        # Raise an error if some columns are missing
        raise ValueError(
            f"[ERROR] module_scores_by_sample.csv missing columns: {sorted(missing)}"
        )

    # Return the loaded DataFrame
    return df


def summarize_module_scores(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate module scores to dataset × condition × module means.

    The output DataFrame has columns:
      - dataset
      - condition
      - module   (one of TGFb_bulk, mTOR_bulk, Angio_bulk, CSC_bulk)
      - n_samples
      - mean_score
    """
    # Melt the wide module columns into long format
    df_long = df.melt(
        id_vars=["sample_id", "dataset", "condition"],
        value_vars=["TGFb_bulk", "mTOR_bulk", "Angio_bulk", "CSC_bulk"],
        var_name="module",
        value_name="score",
    )

    # Group by dataset, condition, and module and compute mean and count
    summary = (
        df_long.groupby(["dataset", "condition", "module"], as_index=False)
        .agg(
            n_samples=("score", "size"),
            mean_score=("score", "mean"),
        )
    )

    # Return the summary table
    return summary


def condition_config(dataset: str, condition: str) -> Dict[str, float]:
    """
    Provide condition-specific configuration such as f_ras and
    a flag for Lepr knockout based on dataset and condition.

    This is a hard-coded mapping for the specific GSE190411 design.

    Returns a dictionary with keys:
      - f_ras: float
      - lepr_knockout: float (0.0 or 1.0; used as multiplier)
    """
    # Normalize the strings to avoid case issues
    d = dataset
    c = condition

    # Default values if no explicit mapping found
    f_ras = 1.0
    lepr_knockout = 0.0  # 0 = WT, 1 = KO

    # Map Bl6 Normal vs SCC
    if d == "Bl6" and c == "Normal":
        # Mild Ras activation in benign state
        f_ras = 0.3
        lepr_knockout = 0.0
    elif d == "Bl6" and c == "SCC":
        # Strong Ras activation in SCC
        f_ras = 1.0
        lepr_knockout = 0.0

    # Map PAP_SCC Papilloma vs SCC
    elif d == "PAP_SCC" and c == "Papilloma":
        # Intermediate Ras level in papilloma
        f_ras = 0.7
        lepr_knockout = 0.0
    elif d == "PAP_SCC" and c == "SCC":
        # High Ras in invasive SCC
        f_ras = 1.0
        lepr_knockout = 0.0

    # Map PDV WT vs LeprKO
    elif d == "PDV" and c == "PDV_WT":
        # Ras on, Lepr intact
        f_ras = 1.0
        lepr_knockout = 0.0
    elif d == "PDV" and c == "PDV_LeprKO":
        # Ras on, but Lepr is knocked out
        f_ras = 1.0
        lepr_knockout = 1.0

    # Return the configuration as a dict
    return {"f_ras": f_ras, "lepr_knockout": lepr_knockout}


def simulate_steady_state(
    base_params: RasCSCParams,
    f_ras: float,
    lepr_knockout: float,
    t_span: Tuple[float, float],
    dt: float,
) -> Dict[str, float]:
    """
    Simulate the ODE to steady state for a specific condition.

    This function:
      1) Clones the base RasCSCParams.
      2) Applies condition-specific f_ras.
      3) Applies a crude Lepr knockout by scaling k_R_prod and k_M_max.
      4) Integrates from a fixed y0 to t_end.
      5) Returns a dictionary with approximate steady-state
         values for C, A, T, R, and M.

    The M value is recomputed from the final state using
    the same algebra as in ras_csc_rhs.
    """
    # Create a copy of the base parameters (dataclasses are mutable)
    params = RasCSCParams(**asdict(base_params))

    # Set the condition-specific Ras input
    params.f_ras = float(f_ras)

    # Apply Lepr knockout scaling (1 = full KO, 0 = WT)
    # Here we simply scale k_R_prod and k_M_max by (1 - lepr_knockout).
    # This is crude but captures the idea that Lepr→mTOR axis is disabled.
    scale = 1.0 - float(lepr_knockout)
    params.k_R_prod *= scale
    params.k_M_max *= scale

    # Define a fixed initial condition near the benign state
    y0 = np.array([0.01, 0.01, 0.01, 0.01], dtype=float)

    # Run the RK4 integrator
    times, traj = simulate_trajectory(
        rhs=ras_csc_rhs,
        y0=y0,
        params=params,
        t_span=t_span,
        dt=dt,
    )

    # Extract the final state
    C, A, T, R = traj[-1, :]

    # Recompute M at steady state to align with the algebraic definition
    L = params.L_sys + params.k_L_A * A
    S = L * R
    S_q = S ** params.q_M
    K_Mq = params.K_M ** params.q_M
    denom_M = S_q + K_Mq
    if denom_M > 1e-12:
        M = params.k_M_max * (S_q / denom_M)
    else:
        M = 0.0

    # Return a dictionary of the steady-state values
    return {"C": C, "A": A, "T": T, "R": R, "M": M}


def build_observation_matrix(summary: pd.DataFrame) -> pd.DataFrame:
    """
    Construct a "wide" observation matrix indexed by dataset+condition,
    with columns corresponding to the four modules interpreted as
    C, A, T, M.

    Output columns:
      - dataset
      - condition
      - obs_C   (CSC_bulk)
      - obs_A   (Angio_bulk)
      - obs_T   (TGFb_bulk)
      - obs_M   (mTOR_bulk)
    """
    # Pivot the long summary table into wide format: one row per dataset+condition
    obs_wide = summary.pivot_table(
        index=["dataset", "condition"],
        columns="module",
        values="mean_score",
    ).reset_index()

    # Rename module columns to map onto state variables
    rename_map = {
        "CSC_bulk": "obs_C",
        "Angio_bulk": "obs_A",
        "TGFb_bulk": "obs_T",
        "mTOR_bulk": "obs_M",
    }

    # Apply the renaming, ignoring any missing keys
    obs_wide = obs_wide.rename(columns=rename_map)

    # Check that the required mapped columns exist
    required_obs = {"obs_C", "obs_A", "obs_T", "obs_M"}
    missing = required_obs - set(obs_wide.columns)
    if missing:
        raise ValueError(
            f"[ERROR] Observation matrix missing expected module columns: {sorted(missing)}"
        )

    # Return the observation matrix
    return obs_wide


def evaluate_param_set(
    base_params: RasCSCParams,
    obs_wide: pd.DataFrame,
    t_span: Tuple[float, float],
    dt: float,
) -> float:
    """
    Compute a scalar error for a candidate parameter set.

    Steps:
      1) For each dataset+condition in obs_wide:
           - Get condition_config (f_ras, lepr_knockout).
           - Simulate steady state for that condition.
      2) Build two matrices:
           - observed (rows = conditions, cols = C,A,T,M)
           - simulated (same shape).
      3) Convert each matrix to z-scores across conditions
         for each variable to compare *patterns* rather than
         absolute scales.
      4) Compute sum of squared differences between
         z_sim and z_obs across all variables and conditions.

    Returns:
      - A non-negative scalar error; smaller is better.
    """
    # Lists to accumulate rows
    obs_rows = []
    sim_rows = []

    # Loop over each row (dataset+condition) in the observation matrix
    for _, row in obs_wide.iterrows():
        # Extract dataset and condition
        dataset = str(row["dataset"])
        condition = str(row["condition"])

        # Construct condition-specific configuration
        cfg = condition_config(dataset, condition)

        # Simulate steady state for this condition
        sim = simulate_steady_state(
            base_params=base_params,
            f_ras=cfg["f_ras"],
            lepr_knockout=cfg["lepr_knockout"],
            t_span=t_span,
            dt=dt,
        )

        # Append observed values [C,A,T,M] from module means
        obs_rows.append(
            [
                row["obs_C"],
                row["obs_A"],
                row["obs_T"],
                row["obs_M"],
            ]
        )

        # Append simulated values [C,A,T,M] from the ODE
        sim_rows.append(
            [
                sim["C"],
                sim["A"],
                sim["T"],
                sim["M"],
            ]
        )

    # Convert lists to numpy arrays
    obs_arr = np.asarray(obs_rows, dtype=float)
    sim_arr = np.asarray(sim_rows, dtype=float)

    # Z-score each column (variable) across conditions for both
    # observed and simulated matrices, guarding against zero variance.
    def zscore_mat(mat: np.ndarray) -> np.ndarray:
        # Compute mean across rows
        mu = mat.mean(axis=0)
        # Compute standard deviation across rows
        sigma = mat.std(axis=0)
        # Avoid division by zero: replace 0 with 1
        sigma_safe = np.where(sigma < 1e-8, 1.0, sigma)
        # Return the z-scored matrix
        return (mat - mu) / sigma_safe

    # Z-score observed and simulated
    z_obs = zscore_mat(obs_arr)
    z_sim = zscore_mat(sim_arr)

    # Compute sum of squared differences across all entries
    diff = z_sim - z_obs
    error = float(np.sum(diff ** 2))

    # Return the scalar error
    return error


def sample_base_params(rng: np.random.Generator) -> RasCSCParams:
    """
    Sample a candidate RasCSCParams from simple uniform ranges.

    These ranges are intentionally narrow and biologically plausible
    in order of magnitude. They can be tuned later.

    For now, we:
      - Fix time-scale-like parameters (decays) near 1.
      - Randomize key gains and thresholds in [0.1, 2.0] ranges.
    """
    # Sample helper: uniform in [low, high]
    def u(low: float, high: float) -> float:
        return float(rng.uniform(low, high))

    # Build parameter object with sampled values
    params = RasCSCParams(
        # CSC equation
        f_ras=1.0,                   # Overridden per condition
        k_C_ras=u(0.01, 0.3),
        k_C_TGFb=u(0.1, 1.0),
        K_C=u(0.1, 1.0),
        n_C=u(2.0, 5.0),
        k_C_M=u(0.1, 1.5),
        k_C_deg=u(0.3, 1.5),

        # Angiogenesis equation
        k_A_ras=u(0.01, 0.3),
        k_A_C=u(0.2, 1.5),
        k_A_deg=u(0.3, 1.5),

        # TGFβ equation
        k_T_A=u(0.2, 2.0),
        k_T_C=u(0.0, 0.5),
        k_T_deg=u(0.3, 1.5),

        # Lepr equation
        k_R_prod=u(0.2, 2.0),
        K_R=u(0.1, 1.0),
        p_R=u(2.0, 5.0),
        k_R_deg=u(0.3, 1.5),

        # Leptin algebraic
        L_sys=u(0.0, 0.3),
        k_L_A=u(0.2, 1.5),

        # mTOR algebraic
        k_M_max=u(0.2, 2.0),
        K_M=u(0.1, 1.0),
        q_M=u(2.0, 5.0),

        # Numerical option
        clip_state=True,
    )

    # Return sampled parameters
    return params


def run_random_search(
    obs_wide: pd.DataFrame,
    n_samples: int,
    t_span: Tuple[float, float],
    dt: float,
    seed: int = 123,
) -> RasCSCParams:
    """
    Run a simple random search over RasCSCParams.

    This procedure:
      1) Draws n_samples random parameter sets from sample_base_params.
      2) Evaluates each with evaluate_param_set.
      3) Keeps the parameter set with minimal error.

    Returns:
      - RasCSCParams instance with the best (lowest-error) parameters.
    """
    # Initialize random generator with a fixed seed for reproducibility
    rng = np.random.default_rng(seed)

    # Initialize best variables
    best_params: RasCSCParams | None = None
    best_error: float = np.inf

    # Loop over number of samples to try
    for i in range(n_samples):
        # Sample a candidate parameter set
        candidate = sample_base_params(rng)

        try:
            # Compute the error for this candidate
            err = evaluate_param_set(
                base_params=candidate,
                obs_wide=obs_wide,
                t_span=t_span,
                dt=dt,
            )
        except Exception as exc:
            # If simulation fails, skip this candidate but log the issue
            # (in a real project we might use logging instead of print)
            print(f"[WARN] Candidate {i} failed during evaluation: {exc}")
            continue

        # If this is the best so far, store it
        if err < best_error:
            best_error = err
            best_params = candidate
            print(f"[INFO] New best candidate at i={i}, error={best_error:.4f}")

    # If no candidate succeeded, raise an error
    if best_params is None:
        raise RuntimeError("[ERROR] Random search failed; no valid parameter sets found.")

    # Print final best error for transparency
    print(f"[INFO] Best error after {n_samples} samples: {best_error:.4f}")

    # Return the best parameter set
    return best_params


#--------------------------------------------------------------------
# Main
#--------------------------------------------------------------------

def main() -> None:
    """
    Main entry point:
      - Load module scores.
      - Summarize to dataset×condition×module means.
      - Build observation matrix.
      - Run random search to calibrate model_A_params.
      - Construct model_B_params by zeroing k_C_M (remove feedback).
      - Write results JSON for downstream use.
    """
    # Ensure output directory exists
    CALIB_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load per-sample module scores
    print(f"[STEP] Loading module scores from: {MODULE_SCORES_PATH}")
    df_scores = load_module_scores(MODULE_SCORES_PATH)

    # Summarize to dataset×condition×module
    print("[STEP] Summarizing module scores")
    df_summary = summarize_module_scores(df_scores)

    # Build observation matrix mapping modules to C,A,T,M
    print("[STEP] Building observation matrix")
    obs_wide = build_observation_matrix(df_summary)

    # Run random search to find best parameter set for model A
    print("[STEP] Running random search for model_A_params")
    best_params_A = run_random_search(
        obs_wide=obs_wide,
        n_samples=N_SAMPLES,
        t_span=T_SPAN,
        dt=DT,
        seed=123,
    )

    # Construct model B as "feedback off": copy A and set k_C_M = 0
    print("[STEP] Constructing model_B_params (feedback off: k_C_M=0)")
    params_A_dict = asdict(best_params_A)
    params_B_dict = dict(params_A_dict)
    params_B_dict["k_C_M"] = 0.0

    # Assemble calibration results for JSON
    calib_results = {
        "model_A_params": params_A_dict,
        "model_B_params": params_B_dict,
        "meta": {
            "description": "Coarse random-search calibration from RNA-Seq module scores",
            "n_samples": N_SAMPLES,
            "t_span": list(T_SPAN),
            "dt": DT,
        },
    }

    # Write JSON to disk
    with open(CALIB_OUTPUT_PATH, "w", encoding="utf-8") as f:
        json.dump(calib_results, f, indent=2)

    print(f"[DONE] Wrote calibration results to: {CALIB_OUTPUT_PATH}")


if __name__ == "__main__":
    main()
