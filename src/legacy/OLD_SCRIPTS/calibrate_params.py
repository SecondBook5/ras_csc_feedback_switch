#!/usr/bin/env python3
"""
Analytic / heuristic calibration of Ras–CSC feedback model parameters
from bulk RNA-Seq calibration targets.

This script is designed to be consistent with:
  - ras_csc_model.RasCSCParams
  - fit_ras_csc_model.py / evaluate_ras_csc_fit.py

It performs a coarse, transparent calibration:

  1) Load bulk targets:
       data/processed/rnaseq/ras_csc_calibration_targets.csv

     Required columns:
       - dataset
       - condition
       - C_target, A_target, T_target, M_target

  2) Rescale each target variable (C, A, T, M) to [0, 1] across all
     dataset/condition combinations, so we can reason in a
     dimensionless space.

  3) Use the rescaled CSC targets for SCC-like conditions
     (SCC, PDV_WT) to set the effective activation needed in the
     CSC equation:

         dC/dt = activation * (1 - C) - k_C_deg * C

     At steady state:

         C* = activation / (activation + k_C_deg)

     Given C* and k_C_deg we can solve for "activation".

     We then split activation into three pieces:

         activation = k_C_ras * f_ras + k_C_TGFb * Hill_T + k_C_M * M

     and choose positive k_C_ras, k_C_TGFb, k_C_M that sum to the
     required activation, assuming Hill_T ~ 1 and M ~ 1 in the SCC-like
     state after rescaling.

  4) Set the remaining parameters to biologically plausible values
     on order 0.1–1.0. These are documented explicitly below.

  5) Construct:
       - model_A_params: feedback ON (k_C_M > 0)
       - model_B_params: feedback OFF (k_C_M = 0)

  6) Write results to:
       data/processed/model_fits/model_calibration_results.json

Schema:

  {
    "model_A_params": { RasCSCParams fields ... },
    "model_B_params": { RasCSCParams fields ... },
    "meta": {
        "description": "...",
        "source_csv": "...",
        "normalization": "minmax_per_variable",
        "scc_like_conditions": [...],
        "t_span": [0.0, 100.0],
        "dt": 0.1
    }
  }

This script does not attempt full Bayesian or MCMC inference. It gives a
transparent, analytic starting point consistent with the ODE structure.
For rigorous fitting, use fit_ras_csc_model.py on top of this.
"""

from __future__ import annotations

# Import standard libraries
import json
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Dict, List, Tuple

# Import third-party libraries
import numpy as np
import pandas as pd

# Import your ODE parameter dataclass
from ras_csc_model import RasCSCParams


# ---------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------

def get_project_root() -> Path:
    """
    Infer the project root as the parent of the 'src' directory that
    contains this script.

    Returns:
        Path: Absolute path to the project root directory.

    Raises:
        RuntimeError: If the inferred root does not exist.
    """
    # Get the absolute path to this script
    script_path: Path = Path(__file__).resolve()
    # Move one level up from src/ to the project root
    root_dir: Path = script_path.parent.parent

    # Check that the root directory exists
    if not root_dir.is_dir():
        # Raise a clear error if the layout is unexpected
        raise RuntimeError(
            f"[ERROR] Inferred project root does not exist: {root_dir}"
        )

    # Return the verified project root
    return root_dir


# ---------------------------------------------------------------------
# Data loading and normalization
# ---------------------------------------------------------------------

def load_calibration_targets(root_dir: Path) -> pd.DataFrame:
    """
    Load the bulk RNA-Seq calibration targets that define the target
    values for (C, A, T, M) in z-score space.

    The expected file is:
      data/processed/rnaseq/ras_csc_calibration_targets.csv

    Required columns:
      - dataset
      - condition
      - C_target
      - A_target
      - T_target
      - M_target

    Args:
        root_dir: Project root directory.

    Returns:
        DataFrame: Cleaned targets table.

    Raises:
        FileNotFoundError: If the CSV is missing.
        ValueError: If required columns are absent or data is empty.
    """
    # Build the path to the calibration targets CSV
    csv_path: Path = root_dir / "data" / "processed" / \
        "rnaseq" / "ras_csc_calibration_targets.csv"

    # Check that the file exists before reading
    if not csv_path.exists():
        # Prompt the user to run the RNA pipeline first
        raise FileNotFoundError(
            f"[ERROR] Calibration targets CSV not found at {csv_path}. "
            f"Run 05_export_ras_csc_targets.R first."
        )

    # Read the CSV into a DataFrame
    df: pd.DataFrame = pd.read_csv(csv_path)

    # Define required columns
    required_cols: List[str] = [
        "dataset", "condition",
        "C_target", "A_target", "T_target", "M_target",
    ]

    # Identify any missing columns
    missing: List[str] = [c for c in required_cols if c not in df.columns]
    if missing:
        # Raise an informative error if the CSV is not as expected
        raise ValueError(
            f"[ERROR] Calibration targets CSV is missing required columns: {missing}"
        )

    # Drop rows with any missing target values
    df_clean: pd.DataFrame = df.dropna(
        subset=["C_target", "A_target", "T_target", "M_target"]
    )

    # Ensure that we still have data
    if df_clean.empty:
        raise ValueError(
            "[ERROR] No usable rows remain in calibration targets after "
            "dropping rows with missing values."
        )

    # Return the cleaned DataFrame
    return df_clean


def minmax_normalize_column(col: pd.Series) -> Tuple[pd.Series, float, float]:
    """
    Perform simple min–max normalization of a numeric column.

    For each value x in col:

        x_norm = (x - min_x) / (max_x - min_x)

    If max_x == min_x, the column has zero range and all normalized
    values are set to 0.5 for lack of a better choice, with a warning.

    Args:
        col: Pandas Series of numeric values.

    Returns:
        Tuple of:
          - normalized Series
          - min value (float)
          - max value (float)
    """
    # Compute min and max
    min_x: float = float(col.min())
    max_x: float = float(col.max())
    range_x: float = max_x - min_x

    # Guard against zero range
    if range_x <= 1e-12:
        # All values identical; return 0.5 everywhere
        norm = pd.Series(0.5, index=col.index)
        return norm, min_x, max_x

    # Compute normalized values
    norm = (col - min_x) / range_x

    return norm, min_x, max_x


def normalize_targets(df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, Dict[str, float]]]:
    """
    Normalize C_target, A_target, T_target, and M_target to [0, 1]
    across all dataset/condition combinations.

    Args:
        df: DataFrame with raw targets.

    Returns:
        Tuple of:
          - DataFrame with added columns:
              C_norm, A_norm, T_norm, M_norm
          - Dict with min/max per variable for bookkeeping:
              {
                "C": {"min": ..., "max": ...},
                ...
              }
    """
    # Initialize dict to store scaling info
    scales: Dict[str, Dict[str, float]] = {}

    # Copy the input DataFrame so we do not modify in place
    out: pd.DataFrame = df.copy()

    # Normalize each variable and record min/max
    for var, col_name in [("C", "C_target"),
                          ("A", "A_target"),
                          ("T", "T_target"),
                          ("M", "M_target")]:
        # Perform min–max normalization
        norm_col, min_x, max_x = minmax_normalize_column(out[col_name])

        # Store normalized column
        out[f"{var}_norm"] = norm_col

        # Record scaling parameters
        scales[var] = {"min": min_x, "max": max_x}

    return out, scales


# ---------------------------------------------------------------------
# Analytic calibration of RasCSCParams
# ---------------------------------------------------------------------

def pick_scc_like_rows(df_norm: pd.DataFrame) -> pd.DataFrame:
    """
    Select rows that represent "high malignancy / high feedback" states.

    We treat these as SCC-like conditions that should be near the
    upper fixed point in the dimensionless model.

    Currently uses:
      - condition in {"SCC", "PDV_WT"}

    Args:
        df_norm: DataFrame with normalized targets.

    Returns:
        DataFrame subset with SCC-like rows. If none, falls back to
        all rows.
    """
    # Define which condition labels we treat as SCC-like
    scc_like_labels = {"SCC", "PDV_WT"}

    # Filter based on the condition column
    mask = df_norm["condition"].astype(str).isin(scc_like_labels)
    scc_like = df_norm.loc[mask].copy()

    # If no rows match, fall back to using the full DataFrame
    if scc_like.empty:
        scc_like = df_norm.copy()

    return scc_like


def calibrate_ras_csc_params(df_norm: pd.DataFrame) -> RasCSCParams:
    """
    Construct a RasCSCParams instance using the normalized calibration
    targets.

    Strategy:

      1) Use SCC-like rows to estimate a target CSC steady state C* in
         the normalized space (C_norm between 0 and 1).

      2) Set k_C_deg = 1.0. Solve for the effective activation:

             C* = activation / (activation + k_C_deg)

         so:

             activation = k_C_deg * C* / (1 - C* + eps)

      3) Split activation into contributions from Ras, TGFb, and mTOR:

             activation = k_C_ras * f_ras + k_C_TGFb * Hill_T + k_C_M * M

         For calibration we assume Hill_T ≈ 1 and M ≈ 1 in SCC-like
         conditions after normalization. We set:

             k_C_ras  = 0.3 * activation
             k_C_TGFb = 0.3 * activation
             k_C_M    = 0.4 * activation

         This guarantees they sum to activation and are all positive.

      4) Set angiogenesis and TGFb parameters so that A and T track C:

             dA/dt = k_A_ras * f_ras + k_A_C * C - k_A_deg * A
             dT/dt = k_T_A * A + k_T_C * C - k_T_deg * T

         We choose:
             k_A_deg = k_T_deg = 1.0
             k_A_ras small (0.05)
             k_A_C   moderate (1.0)
             k_T_A   moderate (1.0)
             k_T_C   small (0.1)

      5) Lepr dynamics:

             dR/dt = k_R_prod * Hill_T(T) - k_R_deg * R

         We reuse the same Hill form as for CSC with separate K_R, p_R.
         We choose K_R = 0.5, p_R = 3.0, k_R_deg = 1.0, k_R_prod = 1.0.

      6) Leptin and mTOR QSS:

             L = L_sys + k_L_A * A
             M = k_M_max * (S^q_M / (K_M^q_M + S^q_M)),
                 S = L * R

         We pick L_sys small (0.1) and k_L_A = 1.0 so that in high-A
         states, L is of order 1–few. We choose K_M = 0.5, q_M = 2.0,
         k_M_max = 1.0.

    This is intentionally simple and transparent. It gives a consistent
    RasCSCParams object that reflects the data ranges but does not claim
    to be statistically optimal.
    """
    # Identify SCC-like rows
    scc_like = pick_scc_like_rows(df_norm)

    # Compute mean normalized CSC level in SCC-like conditions
    C_star: float = float(scc_like["C_norm"].mean())

    # Clip C_star away from 0 and 1 to avoid infinite activation
    C_star = max(min(C_star, 0.95), 0.05)

    # Set CSC decay rate
    k_C_deg: float = 1.0

    # Small epsilon to avoid division by zero
    eps: float = 1e-8

    # Solve for activation so that C* is a fixed point
    activation: float = k_C_deg * C_star / (1.0 - C_star + eps)

    # Set Ras input to "on" (1.0)
    f_ras: float = 1.0

    # Split activation into Ras, TGFb, and mTOR components
    k_C_ras: float = 0.3 * activation
    k_C_TGFb: float = 0.3 * activation
    k_C_M: float = 0.4 * activation

    # Choose Hill parameters for TGFb -> CSC
    K_C: float = 0.5
    n_C: float = 3.0

    # Angiogenesis parameters
    k_A_deg: float = 1.0
    k_A_ras: float = 0.05
    k_A_C: float = 1.0

    # TGFb parameters
    k_T_deg: float = 1.0
    k_T_A: float = 1.0
    k_T_C: float = 0.1

    # Lepr parameters
    k_R_deg: float = 1.0
    k_R_prod: float = 1.0
    K_R: float = 0.5
    p_R: float = 3.0

    # Leptin QSS parameters
    L_sys: float = 0.1
    k_L_A: float = 1.0

    # mTOR QSS parameters
    k_M_max: float = 1.0
    K_M: float = 0.5
    q_M: float = 2.0

    # Construct RasCSCParams
    params = RasCSCParams(
        f_ras=f_ras,
        k_C_ras=k_C_ras,
        k_C_TGFb=k_C_TGFb,
        K_C=K_C,
        n_C=n_C,
        k_C_M=k_C_M,
        k_C_deg=k_C_deg,
        k_A_ras=k_A_ras,
        k_A_C=k_A_C,
        k_A_deg=k_A_deg,
        k_T_A=k_T_A,
        k_T_C=k_T_C,
        k_T_deg=k_T_deg,
        k_R_prod=k_R_prod,
        K_R=K_R,
        p_R=p_R,
        k_R_deg=k_R_deg,
        L_sys=L_sys,
        k_L_A=k_L_A,
        k_M_max=k_M_max,
        K_M=K_M,
        q_M=q_M,
        clip_state=True,
    )

    return params


def make_feedback_off_variant(params: RasCSCParams) -> RasCSCParams:
    """
    Create a second RasCSCParams instance representing Model B
    (feedback off) by copying Model A and forcing k_C_M = 0.0.

    Args:
        params: RasCSCParams for Model A (feedback on).

    Returns:
        RasCSCParams for Model B (feedback off).
    """
    # Convert the dataclass to a dict
    p_dict: Dict[str, object] = asdict(params)

    # Force feedback coefficient to zero
    p_dict["k_C_M"] = 0.0

    # Construct a new RasCSCParams from the modified dict
    params_b = RasCSCParams(**p_dict)

    return params_b


# ---------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------

def main() -> None:
    """
    Run the analytic calibration workflow and write a JSON file with:

      - model_A_params: RasCSCParams for feedback model
      - model_B_params: RasCSCParams with k_C_M = 0
      - meta: information about normalization and source data

    The output path is:
      data/processed/model_fits/model_calibration_results.json
    """
    # Infer project root
    root_dir: Path = get_project_root()

    # Load raw calibration targets
    targets_df: pd.DataFrame = load_calibration_targets(root_dir)

    # Normalize to [0, 1]
    targets_norm, scales = normalize_targets(targets_df)

    # Calibrate Model A parameters from normalized data
    params_A: RasCSCParams = calibrate_ras_csc_params(targets_norm)

    # Build Model B parameters (feedback off)
    params_B: RasCSCParams = make_feedback_off_variant(params_A)

    # Prepare output directory
    out_dir: Path = root_dir / "data" / "processed" / "model_fits"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Output path
    out_path: Path = out_dir / "model_calibration_results.json"

    # Build meta information
    meta = {
        "description": "Analytic / heuristic calibration from bulk RNA-Seq targets "
                       "with min–max normalization and SCC-like anchoring.",
        "source_csv": str(root_dir / "data" / "processed" / "rnaseq" / "ras_csc_calibration_targets.csv"),
        "normalization": "minmax_per_variable",
        "scales": scales,
        "scc_like_conditions": ["SCC", "PDV_WT"],
        "t_span": [0.0, 100.0],
        "dt": 0.1,
    }

    # Build JSON payload in the schema expected by evaluate_ras_csc_fit.py
    payload = {
        "model_A_params": asdict(params_A),
        "model_B_params": asdict(params_B),
        "meta": meta,
    }

    # Write JSON file
    with out_path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    # Print a concise summary to stdout
    print(f"[OK] Wrote analytic calibration parameters to: {out_path}")
    print("\n[INFO] Model A (feedback ON) key parameters:")
    print(f"  k_C_ras  = {params_A.k_C_ras:.4f}")
    print(f"  k_C_TGFb = {params_A.k_C_TGFb:.4f}")
    print(f"  k_C_M    = {params_A.k_C_M:.4f}")
    print(f"  k_C_deg  = {params_A.k_C_deg:.4f}")
    print("\n[INFO] Model B (feedback OFF) sets k_C_M = 0.0")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        # Print any error clearly and exit with non-zero status
        print(f"[ERROR] calibrate_params.py failed: {exc}", file=sys.stderr)
        sys.exit(1)
