#!/usr/bin/env python3
"""
Evaluate the calibrated Ras–CSC ODE model against RNA-Seq calibration targets.

This script is meant to sit *after* you have:
  1) Built ras_csc_calibration_targets.csv from RNA module scores:
       data/processed/rnaseq/ras_csc_calibration_targets.csv
  2) Calibrated the Ras–CSC model parameters via fit_ras_csc_model.py:
       data/processed/model_fits/model_calibration_results.json

The goal of this script is to:
  - Load the calibrated parameters (model_A_params) from the JSON file.
  - For each (dataset, condition) row in ras_csc_calibration_targets.csv:
      * Simulate the Ras–CSC ODE from a common benign initial condition.
      * Run long enough to approximate the steady state.
      * Compute model-predicted values for (C, A, T, M).
  - Compare these to the corresponding target values
      (C_target, A_target, T_target, M_target).
  - Write a CSV summarizing model vs data and per-condition residuals:
      data/processed/model_fits/model_vs_data.csv
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple, List

import json
import sys

import numpy as np
import pandas as pd

from ras_csc_model import RasCSCParams, ras_csc_rhs, simulate_trajectory


def compute_mtor_from_state(
    C: float,
    A: float,
    T: float,
    R: float,
    params: RasCSCParams,
) -> float:
    """
    Compute mTOR activity M from the quasi-steady-state definitions
    given a state vector (C, A, T, R) and RasCSCParams.

    Mirrors the algebra inside ras_csc_rhs so we can read out M at the
    final state for comparison to M_target.
    """
    # Leptin: L = L_sys + k_L_A * A
    L = params.L_sys + params.k_L_A * A

    # Signal S = L * R
    S = L * R

    # Hill transform
    S_q = S ** params.q_M
    K_Mq = params.K_M ** params.q_M
    denom = S_q + K_Mq

    if denom <= 1e-12:
        return 0.0

    M = params.k_M_max * (S_q / denom)
    return float(M)


def load_calibrated_params(json_path: Path) -> RasCSCParams:
    """
    Load the calibrated Ras–CSC model parameters from the JSON file
    produced by fit_ras_csc_model.py.
    """
    if not json_path.exists():
        raise FileNotFoundError(
            f"Calibration JSON not found at {json_path}. "
            "Run fit_ras_csc_model.py first."
        )

    with json_path.open("r", encoding="utf-8") as f:
        try:
            data = json.load(f)
        except json.JSONDecodeError as e:
            raise ValueError(
                f"Failed to parse JSON at {json_path}: {e}"
            ) from e

    if "model_A_params" not in data:
        raise KeyError(
            "JSON is missing 'model_A_params'. "
            "Check fit_ras_csc_model.py output."
        )

    param_dict = data["model_A_params"]

    try:
        params = RasCSCParams(**param_dict)
    except TypeError as e:
        raise ValueError(
            "Failed to construct RasCSCParams from model_A_params. "
            "Check that all required fields are present and names match "
            "the RasCSCParams dataclass."
        ) from e

    return params


def load_calibration_targets(csv_path: Path) -> pd.DataFrame:
    """
    Load the bulk RNA calibration targets that specify, for each
    (dataset, condition), the target values for (C, A, T, M).

    Requires columns:
      dataset, condition, C_target, A_target, T_target, M_target
    """
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Calibration targets CSV not found at {csv_path}. "
            "Ensure 05_export_ras_csc_targets.R has been run."
        )

    df = pd.read_csv(csv_path)

    required_cols = ["dataset", "condition",
                     "C_target", "A_target", "T_target", "M_target"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"Calibration targets CSV is missing required columns: {missing}"
        )

    df_clean = df.dropna(
        subset=["C_target", "A_target", "T_target", "M_target"])
    if df_clean.empty:
        raise ValueError(
            "Calibration targets CSV has no usable rows after dropping NAs. "
            "Inspect data/processed/rnaseq/ras_csc_calibration_targets.csv."
        )

    return df_clean


def evaluate_model_fit(
    root_dir: Path,
    t_span: Tuple[float, float] = (0.0, 100.0),
    dt: float = 0.1,
    y0: Tuple[float, float, float, float] = (0.01, 0.01, 0.01, 0.01),
) -> None:
    """
    Run the calibrated Ras–CSC model for each bulk RNA calibration target
    and quantify how well the model reproduces the target (C, A, T, M) values.

    Writes:
      data/processed/model_fits/model_vs_data.csv
    """
    calib_json_path = root_dir / "data" / "processed" / \
        "model_fits" / "model_calibration_results.json"
    targets_csv_path = root_dir / "data" / "processed" / \
        "rnaseq" / "ras_csc_calibration_targets.csv"
    out_csv_path = root_dir / "data" / "processed" / \
        "model_fits" / "model_vs_data.csv"

    params = load_calibrated_params(calib_json_path)
    targets_df = load_calibration_targets(targets_csv_path)

    y0_arr = np.array(y0, dtype=float)
    records: List[Dict[str, float]] = []

    for _, row in targets_df.iterrows():
        dataset = str(row["dataset"])
        condition = str(row["condition"])

        C_target = float(row["C_target"])
        A_target = float(row["A_target"])
        T_target = float(row["T_target"])
        M_target = float(row["M_target"])

        try:
            times, traj = simulate_trajectory(
                rhs=ras_csc_rhs,
                y0=y0_arr,
                params=params,
                t_span=t_span,
                dt=dt,
            )
        except Exception as e:
            raise RuntimeError(
                f"Simulation failed for dataset={dataset}, condition={condition}: {e}"
            ) from e

        C_model = float(traj[-1, 0])
        A_model = float(traj[-1, 1])
        T_model = float(traj[-1, 2])
        R_model = float(traj[-1, 3])

        M_model = compute_mtor_from_state(
            C=C_model,
            A=A_model,
            T=T_model,
            R=R_model,
            params=params,
        )

        C_resid = C_model - C_target
        A_resid = A_model - A_target
        T_resid = T_model - T_target
        M_resid = M_model - M_target

        records.append(
            {
                "dataset": dataset,
                "condition": condition,
                "C_target": C_target,
                "A_target": A_target,
                "T_target": T_target,
                "M_target": M_target,
                "C_model": C_model,
                "A_model": A_model,
                "T_model": T_model,
                "M_model": M_model,
                "C_resid": C_resid,
                "A_resid": A_resid,
                "T_resid": T_resid,
                "M_resid": M_resid,
            }
        )

    results_df = pd.DataFrame.from_records(records)

    # Squared errors
    results_df["C_sqerr"] = results_df["C_resid"] ** 2
    results_df["A_sqerr"] = results_df["A_resid"] ** 2
    results_df["T_sqerr"] = results_df["T_resid"] ** 2
    results_df["M_sqerr"] = results_df["M_resid"] ** 2

    sse_C = float(results_df["C_sqerr"].sum())
    sse_A = float(results_df["A_sqerr"].sum())
    sse_T = float(results_df["T_sqerr"].sum())
    sse_M = float(results_df["M_sqerr"].sum())

    n = len(results_df)
    rmse_C = float(np.sqrt(sse_C / n))
    rmse_A = float(np.sqrt(sse_A / n))
    rmse_T = float(np.sqrt(sse_T / n))
    rmse_M = float(np.sqrt(sse_M / n))

    print("\n[SUMMARY] Model vs data fit metrics:")
    print(f"  N conditions: {n}")
    print(
        f"  SSE:  C={sse_C:.4f}, A={sse_A:.4f}, T={sse_T:.4f}, M={sse_M:.4f}")
    print(
        f"  RMSE: C={rmse_C:.4f}, A={rmse_A:.4f}, T={rmse_T:.4f}, M={rmse_M:.4f}")

    results_to_save = results_df.drop(
        columns=["C_sqerr", "A_sqerr", "T_sqerr", "M_sqerr"])

    out_csv_path.parent.mkdir(parents=True, exist_ok=True)
    results_to_save.to_csv(out_csv_path, index=False)
    print(f"[OK] Wrote model vs data table to: {out_csv_path}")


if __name__ == "__main__":
    script_path = Path(__file__).resolve()
    root_dir = script_path.parent.parent

    try:
        evaluate_model_fit(root_dir=root_dir)
    except Exception as exc:
        print(f"[ERROR] evaluate_ras_csc_fit failed: {exc}", file=sys.stderr)
        sys.exit(1)
