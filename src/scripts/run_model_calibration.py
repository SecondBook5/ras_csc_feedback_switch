#!/usr/bin/env python3
"""
run_model_calibration.py

Calibrate the Ras–CSC–microenvironment model to RNA-seq-derived targets.

This script:
    1. Loads calibration targets from:
           data/processed/rnaseq/ras_csc_calibration_targets.csv
       where targets are already in z-score space for modules:
           C, A, T, M
    2. Uses the shared ras_csc_core model to simulate steady states
       across conditions.
    3. Applies a z-score transform to model outputs across conditions
       to ensure consistent comparison to z-scored data.
    4. Runs non-linear least squares (scipy.optimize.least_squares)
       to find a parameter set that minimizes residuals in z-score units.
    5. Saves:
           - results/calibration/optimized_parameters_CORRECTED.json
           - results/calibration/model_vs_data_CORRECTED.csv
    6. Generates publication-quality calibration figures:
           - figures/main/calibration_model_vs_data_CORRECTED.png
           - figures/main/calibration_residuals_CORRECTED.png
    7. Prints human-readable summaries and sanity checks.

All ODE logic is delegated to src/ras_csc_core.py.
"""

from __future__ import annotations

# ----------------------------------------------------------------------
# Standard library imports
# ----------------------------------------------------------------------
from typing import Dict, Any, Tuple, List
import json
from pathlib import Path
import sys

# ----------------------------------------------------------------------
# Third-party imports
# ----------------------------------------------------------------------
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import seaborn as sns

# ----------------------------------------------------------------------
# Global plotting style for professional-looking figures
# ----------------------------------------------------------------------
sns.set_style("whitegrid")

plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["font.size"] = 11
plt.rcParams["axes.labelsize"] = 12
plt.rcParams["axes.titlesize"] = 14
plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 10
plt.rcParams["legend.fontsize"] = 10
plt.rcParams["figure.titlesize"] = 16
plt.rcParams["figure.dpi"] = 300

# ----------------------------------------------------------------------
# Ensure src/ is on the Python path so that ras_csc_core can be imported
# ----------------------------------------------------------------------
CURRENT_DIR: Path = Path(__file__).resolve().parent
SRC_ROOT: Path = CURRENT_DIR.parent

if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from ras_csc_core import simulate_steady_state, get_ras_input  # type: ignore


# ======================================================================
# PATH CONSTANTS
# ======================================================================

ROOT: Path = SRC_ROOT.parent

TARGETS_CSV: Path = ROOT / "data" / "processed" / "rnaseq" / "ras_csc_calibration_targets.csv"

CALIBRATION_DIR: Path = ROOT / "results" / "calibration"
FIG_MAIN_DIR: Path = ROOT / "figures" / "main"

CALIBRATION_DIR.mkdir(parents=True, exist_ok=True)
FIG_MAIN_DIR.mkdir(parents=True, exist_ok=True)


# ======================================================================
# LOADING TARGETS
# ======================================================================

def load_targets() -> Dict[str, Dict[str, Any]]:
    """
    Load z-scored RNA-seq calibration targets from CSV.

    Expected CSV columns:
        - dataset
        - condition
        - C_target
        - A_target
        - T_target
        - M_target

    Returns:
        Dict[str, Dict[str, Any]]:
            Dictionary keyed by "dataset_condition" with values:
                {
                    "C": float,
                    "A": float,
                    "T": float,
                    "M": float,
                    "condition": str,
                    "dataset": str,
                }

    Raises:
        FileNotFoundError:
            If the targets CSV is missing.
        ValueError:
            If required columns are missing or the file is empty.
    """
    if not TARGETS_CSV.exists():
        raise FileNotFoundError(
            f"Calibration targets file not found at: {TARGETS_CSV}"
        )

    df: pd.DataFrame = pd.read_csv(TARGETS_CSV)

    if df.empty:
        raise ValueError(
            f"Calibration targets file at {TARGETS_CSV} is empty."
        )

    required_cols = [
        "dataset",
        "condition",
        "C_target",
        "A_target",
        "T_target",
        "M_target",
    ]
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(
            f"Calibration targets missing required columns: {missing}"
        )

    targets: Dict[str, Dict[str, Any]] = {}
    for _, row in df.iterrows():
        key: str = f"{row['dataset']}_{row['condition']}"
        targets[key] = {
            "C": float(row["C_target"]),
            "A": float(row["A_target"]),
            "T": float(row["T_target"]),
            "M": float(row["M_target"]),
            "condition": str(row["condition"]),
            "dataset": str(row["dataset"]),
        }

    print(f"[INFO] Loaded {len(targets)} calibration targets from {TARGETS_CSV}")
    return targets


# ======================================================================
# PARAMETER HANDLING
# ======================================================================

def build_initial_parameters() -> Tuple[Dict[str, float], np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    """
    Construct initial parameter guesses, fixed parameters, and bounds.

    Returns:
        Tuple containing:
            - fixed_params: dict of parameters held constant
            - x0: np.ndarray of initial guesses for optimized parameters
            - bounds: (lower_bounds, upper_bounds) as np.ndarray pairs
    """
    fixed_params: Dict[str, float] = {
        "eta_C": 1.2,    # TGFβ → CSC contribution
        "alpha_A": 0.4,  # Ras → angiogenesis
        "gamma_T": 0.15,  # CSC → TGFβ
        "eta_R": 0.6,    # TGFβ → LEPR
        "K_R": 0.35,
        "p_R": 2.5,
        "delta_R": 0.6,
        "mu_L": 0.9,     # angiogenesis → leptin
        "delta_L": 0.6,
        "K_M": 0.25,
        "q_M": 3.5,
    }

    x0: np.ndarray = np.array(
        [
            0.15,  # rho_C
            2.5,   # eta_C_M
            0.25,  # K_C
            3.5,   # n_C
            0.8,   # delta_C
            1.0,   # beta_A
            0.6,   # delta_A
            0.8,   # alpha_T
            1.2,   # delta_T
            1.0,   # eta_M
            0.6,   # delta_M
        ],
        dtype=float,
    )

    lower_bounds: np.ndarray = np.array(
        [
            0.05,  # rho_C
            0.5,   # eta_C_M
            0.1,   # K_C
            2.0,   # n_C
            0.3,   # delta_C
            0.3,   # beta_A
            0.3,   # delta_A
            0.3,   # alpha_T
            0.5,   # delta_T
            0.3,   # eta_M
            0.3,   # delta_M
        ],
        dtype=float,
    )

    upper_bounds: np.ndarray = np.array(
        [
            0.5,   # rho_C
            5.0,   # eta_C_M
            0.5,   # K_C
            5.0,   # n_C
            2.0,   # delta_C
            2.0,   # beta_A
            1.5,   # delta_A
            2.0,   # alpha_T
            3.0,   # delta_T
            2.0,   # eta_M
            1.5,   # delta_M
        ],
        dtype=float,
    )

    bounds: Tuple[np.ndarray, np.ndarray] = (lower_bounds, upper_bounds)
    return fixed_params, x0, bounds


# ======================================================================
# OBJECTIVE FUNCTION
# ======================================================================

def objective_function(
    param_vec: np.ndarray,
    targets: Dict[str, Dict[str, Any]],
    fixed_params: Dict[str, float],
) -> np.ndarray:
    """
    Compute residuals between z-scored model outputs and z-scored targets.
    """
    params: Dict[str, float] = dict(fixed_params)
    params.update(
        {
            "rho_C": float(param_vec[0]),
            "eta_C_M": float(param_vec[1]),
            "K_C": float(param_vec[2]),
            "n_C": float(param_vec[3]),
            "delta_C": float(param_vec[4]),
            "beta_A": float(param_vec[5]),
            "delta_A": float(param_vec[6]),
            "alpha_T": float(param_vec[7]),
            "delta_T": float(param_vec[8]),
            "eta_M": float(param_vec[9]),
            "delta_M": float(param_vec[10]),
        }
    )

    model_raw: Dict[str, Dict[str, float]] = {}

    for key, target in targets.items():
        condition: str = target["condition"]
        f_ras: float = get_ras_input(condition)
        local_params: Dict[str, float] = dict(params)

        if "LeprKO" in condition:
            local_params["mu_L"] = local_params["mu_L"] * 0.2

        ss = simulate_steady_state(
            f_ras=f_ras,
            params=local_params,
            y0=None,
            t_max=200.0,
            n_steps=2000,
        )

        model_raw[key] = {
            "C_raw": float(ss[0]),
            "A_raw": float(ss[1]),
            "T_raw": float(ss[2]),
            "M_raw": float(ss[5]),
        }

    model_z: Dict[str, Dict[str, float]] = {k: {} for k in model_raw.keys()}
    for module in ["C", "A", "T", "M"]:
        raw_key = f"{module}_raw"
        values = np.array(
            [model_raw[k][raw_key] for k in model_raw.keys()],
            dtype=float,
        )
        mean_val: float = float(values.mean())
        std_val: float = float(values.std())
        if std_val < 1e-6:
            std_val = 1.0

        for idx, key in enumerate(model_raw.keys()):
            z_value: float = (values[idx] - mean_val) / std_val
            model_z[key][f"{module}_z"] = z_value

    residuals: List[float] = []
    for key, target in targets.items():
        residuals.append(model_z[key]["C_z"] - float(target["C"]))
        residuals.append(model_z[key]["A_z"] - float(target["A"]))
        residuals.append(model_z[key]["T_z"] - float(target["T"]))
        residuals.append(model_z[key]["M_z"] - float(target["M"]))

    return np.asarray(residuals, dtype=float)


# ======================================================================
# FITTING WRAPPER
# ======================================================================

def fit_model(
    targets: Dict[str, Dict[str, Any]],
) -> Tuple[Dict[str, float], float, float]:
    """
    Run non-linear least squares to fit the model to z-scored targets.
    """
    print("\n" + "=" * 70)
    print("RAS–CSC MODEL CALIBRATION (Z-SCORE SPACE)")
    print("=" * 70)

    fixed_params, x0, bounds = build_initial_parameters()

    print(f"[INFO] Number of targets: {len(targets)}")
    print(f"[INFO] Initial parameter vector: {x0}")
    print(f"[INFO] Lower bounds: {bounds[0]}")
    print(f"[INFO] Upper bounds: {bounds[1]}")

    try:
        result = least_squares(
            fun=objective_function,
            x0=x0,
            args=(targets, fixed_params),
            bounds=bounds,
            method="trf",
            verbose=2,
            max_nfev=500,
            ftol=1e-6,
            xtol=1e-6,
        )
    except Exception as exc:
        raise RuntimeError(
            f"Least-squares optimization failed with error: {exc}"
        ) from exc

    opt_params: Dict[str, float] = dict(fixed_params)
    opt_params.update(
        {
            "rho_C": float(result.x[0]),
            "eta_C_M": float(result.x[1]),
            "K_C": float(result.x[2]),
            "n_C": float(result.x[3]),
            "delta_C": float(result.x[4]),
            "beta_A": float(result.x[5]),
            "delta_A": float(result.x[6]),
            "alpha_T": float(result.x[7]),
            "delta_T": float(result.x[8]),
            "eta_M": float(result.x[9]),
            "delta_M": float(result.x[10]),
        }
    )

    residuals_opt: np.ndarray = objective_function(
        param_vec=result.x,
        targets=targets,
        fixed_params=fixed_params,
    )

    rss: float = float(np.sum(residuals_opt ** 2))
    rmse: float = float(np.sqrt(rss / residuals_opt.size))

    print(f"\n[SUCCESS] Optimization finished with status {result.status}")
    print(f"  Message: {result.message}")
    print(f"  RSS  = {rss:.3f}")
    print(f"  RMSE = {rmse:.3f} z-score units")

    print("\n[SANITY CHECK] Parameter ratios:")
    print(
        f"  rho_C / delta_C   = {opt_params['rho_C'] / opt_params['delta_C']:.3f}"
    )
    print(
        f"  beta_A / delta_A  = {opt_params['beta_A'] / opt_params['delta_A']:.3f}"
    )
    print(
        f"  eta_M / delta_M   = {opt_params['eta_M'] / opt_params['delta_M']:.3f}"
    )

    return opt_params, rss, rmse


# ======================================================================
# VALIDATION: MODEL VS DATA TABLE
# ======================================================================

def build_model_vs_data_table(
    params: Dict[str, float],
    targets: Dict[str, Dict[str, Any]],
) -> pd.DataFrame:
    """
    Build a tidy DataFrame comparing model outputs to data per condition.

    Columns:
        - dataset
        - condition
        - module
        - data_zscore
        - model_zscore
        - model_raw
        - residual
    """
    raw: Dict[str, Dict[str, float]] = {}

    for key, target in targets.items():
        condition: str = target["condition"]
        f_ras: float = get_ras_input(condition)
        local_params: Dict[str, float] = dict(params)

        if "LeprKO" in condition:
            local_params["mu_L"] = local_params["mu_L"] * 0.2

        ss = simulate_steady_state(
            f_ras=f_ras,
            params=local_params,
            y0=None,
            t_max=200.0,
            n_steps=2000,
        )

        raw[key] = {
            "C": float(ss[0]),
            "A": float(ss[1]),
            "T": float(ss[2]),
            "M": float(ss[5]),
        }

    model_z: Dict[str, Dict[str, float]] = {k: {} for k in raw.keys()}

    for module in ["C", "A", "T", "M"]:
        values = np.array(
            [raw[k][module] for k in raw.keys()],
            dtype=float,
        )
        mean_val: float = float(values.mean())
        std_val: float = float(values.std())
        if std_val < 1e-6:
            std_val = 1.0

        for idx, key in enumerate(raw.keys()):
            z_value: float = (values[idx] - mean_val) / std_val
            model_z[key][module] = z_value

    rows: List[Dict[str, Any]] = []
    for key, target in targets.items():
        for module in ["C", "A", "T", "M"]:
            rows.append(
                {
                    "dataset": target["dataset"],
                    "condition": target["condition"],
                    "module": module,
                    "data_zscore": float(target[module]),
                    "model_zscore": float(model_z[key][module]),
                    "model_raw": float(raw[key][module]),
                    "residual": float(model_z[key][module] - target[module]),
                }
            )

    df: pd.DataFrame = pd.DataFrame(rows)
    return df


# ======================================================================
# FIGURE GENERATION
# ======================================================================

def plot_calibration_results(df: pd.DataFrame) -> None:
    """
    Generate publication-quality calibration figures.

    Figures:
        1. 2×2 panel of C, A, T, M:
             - x-axis: condition
             - y-axis: z-score
             - lines/points for data vs model
        2. Residual distribution by module:
             - boxplot + jitter of residuals

    Saves:
        - figures/main/calibration_model_vs_data_CORRECTED.png
        - figures/main/calibration_residuals_CORRECTED.png
    """
    df_plot = df.copy()

    preferred_order = ["Normal", "Papilloma", "SCC", "PDV_LeprKO", "PDV_WT"]
    existing_conditions = df_plot["condition"].unique().tolist()
    ordered_conditions: List[str] = [
        c for c in preferred_order if c in existing_conditions
    ]
    for c in sorted(existing_conditions):
        if c not in ordered_conditions:
            ordered_conditions.append(c)

    df_plot["condition"] = pd.Categorical(
        df_plot["condition"], categories=ordered_conditions, ordered=True
    )

    df_long = pd.melt(
        df_plot,
        id_vars=["dataset", "condition", "module"],
        value_vars=["data_zscore", "model_zscore"],
        var_name="source",
        value_name="zscore",
    )

    source_map = {
        "data_zscore": "RNA-seq (z-score)",
        "model_zscore": "Model (z-score)",
    }
    df_long["source"] = df_long["source"].map(source_map)

    modules = ["C", "A", "T", "M"]
    pretty_names = {
        "C": "CSC module",
        "A": "Angiogenesis module",
        "T": "TGFβ module",
        "M": "mTOR module",
    }

    fig, axes = plt.subplots(2, 2, figsize=(10, 6), sharey=False)
    axes = axes.flatten()

    palette = sns.color_palette("colorblind", n_colors=2)

    for ax, module in zip(axes, modules):
        sub = df_long[df_long["module"] == module].copy()
        if sub.empty:
            ax.set_visible(False)
            continue

        ax.axhline(0.0, color="0.8", linewidth=1.0, zorder=0)

        sns.pointplot(
            data=sub,
            x="condition",
            y="zscore",
            hue="source",
            dodge=0.4,
            markers=["o", "s"],
            linestyles=["-", "--"],
            palette=palette,
            errorbar="se",
            ax=ax,
        )

        ax.set_title(pretty_names.get(module, module))
        ax.set_xlabel("")
        ax.set_ylabel("z-score")
        ax.tick_params(axis="x", rotation=30)
        ax.grid(alpha=0.3)

        ax.legend_.remove()

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 1.02),
    )

    plt.tight_layout(rect=[0, 0, 1, 0.92])
    out_model_vs_data = FIG_MAIN_DIR / "calibration_model_vs_data_CORRECTED.png"
    fig.savefig(out_model_vs_data, bbox_inches="tight")
    plt.close(fig)
    print(f"[SAVED] Calibration model-vs-data figure to {out_model_vs_data}")

    fig_res, ax_res = plt.subplots(figsize=(6, 4))

    ax_res.axhline(0.0, color="0.8", linewidth=1.0, zorder=0)

    sns.boxplot(
        data=df_plot,
        x="module",
        y="residual",
        whis=1.5,
        width=0.6,
        fliersize=0,
        ax=ax_res,
    )

    sns.stripplot(
        data=df_plot,
        x="module",
        y="residual",
        color="black",
        size=4,
        alpha=0.6,
        jitter=0.2,
        ax=ax_res,
    )

    ax_res.set_xlabel("Module")
    ax_res.set_ylabel("Residual (model – data, z-score)")
    ax_res.set_title("Residuals by module")
    ax_res.grid(alpha=0.3)

    plt.tight_layout()
    out_res = FIG_MAIN_DIR / "calibration_residuals_CORRECTED.png"
    fig_res.savefig(out_res, bbox_inches="tight")
    plt.close(fig_res)
    print(f"[SAVED] Calibration residuals figure to {out_res}")


# ======================================================================
# MAIN
# ======================================================================

def main() -> Tuple[Dict[str, float], pd.DataFrame]:
    """
    Entry point for Ras–CSC model calibration.
    """
    print("\n" + "=" * 70)
    print("RUNNING RAS–CSC MODEL CALIBRATION PIPELINE")
    print("=" * 70)

    targets = load_targets()
    opt_params, rss, rmse = fit_model(targets)

    param_path: Path = CALIBRATION_DIR / "optimized_parameters_CORRECTED.json"
    with param_path.open("w", encoding="utf-8") as f:
        json.dump(opt_params, f, indent=2)
    print(f"\n[SAVED] Optimized parameters to {param_path}")

    df_results: pd.DataFrame = build_model_vs_data_table(
        params=opt_params,
        targets=targets,
    )

    results_csv: Path = CALIBRATION_DIR / "model_vs_data_CORRECTED.csv"
    df_results.to_csv(results_csv, index=False)
    print(f"[SAVED] Model vs data comparison table to {results_csv}")

    print("\n[SANITY CHECK] Model raw outputs by condition (C and M):")
    for condition in df_results["condition"].unique():
        subset = df_results[
            (df_results["condition"] == condition)
            & (df_results["module"].isin(["C", "M"]))
        ]
        for module in ["C", "M"]:
            row = subset[subset["module"] == module].head(1)
            if not row.empty:
                value = float(row["model_raw"].iloc[0])
                print(f"  {condition:15s} {module}: {value:.3f}")

    try:
        plot_calibration_results(df_results)
    except Exception as exc:
        print(f"[WARN] Failed to generate calibration figures: {exc}")

    print("\n" + "=" * 70)
    print("CALIBRATION COMPLETE")
    print("=" * 70)

    return opt_params, df_results


if __name__ == "__main__":
    main()
