#!/usr/bin/env python3
"""
run_sensitivity_analyses.py

Run sensitivity and robustness analyses for the Ras–CSC feedback model.

This script assumes that:
    - The model has already been calibrated using run_model_calibration.py.
    - Optimized parameters are stored in:
          results/calibration/optimized_parameters_CORRECTED.json
    - Core ODE / steady-state logic is defined in:
          src/ras_csc_core.py

Analyses performed (numerical + figures):

    0. Calibration RSS summary:
         - Load per-point residuals from calibration.
         - Compute global RSS and per-module RSS.
         - Save:
               results/sensitivity/calibration_rss_CORRECTED.csv
         - Plot:
               figures/sensitivity/calibration_rss_by_module_CORRECTED.png
         - Print a human-readable summary of which modules are
           best / worst fit.

    1. Hysteresis sensitivity:
         - For a set of key parameters, scale them by factors
           {0.8, 1.0, 1.2} around the calibrated value.
         - For each scaled set, compute the maximum CSC gap between
           Ras ↑ and Ras ↓ branches.
         - Save:
               results/sensitivity/hysteresis_sensitivity_CORRECTED.csv
         - Plot:
               figures/sensitivity/hysteresis_sensitivity_CORRECTED.png

    2. Perturbation-effect sensitivity:
         - For the same scaled parameter sets, re-run a high-Ras
           LEPR / mTOR perturbation experiment (mu_L × 0.3, eta_M × 0.3).
         - Save:
               results/sensitivity/perturbation_sensitivity_CORRECTED.csv
         - Plot:
               figures/sensitivity/perturbation_sensitivity_CORRECTED.png

    3. LEPR–mTOR perturbation grid:
         - For a grid of (mu_L_scale, eta_M_scale) values, simulate
           high-Ras steady states and record C, M, and ΔC vs baseline.
         - Save:
               results/perturbations/perturbation_grid_CORRECTED.csv
         - Plot:
               figures/sensitivity/perturbation_grid_deltaC_CORRECTED.png

    4. Hysteresis robustness to initial conditions:
         - Repeat Ras ↑ / Ras ↓ sweeps for several initial-condition
           pairs (benign- and malignant-like seeds).
         - Save:
               results/model/ras_hysteresis_initial_seeds_CORRECTED.csv
         - Plot:
               figures/model/ras_hysteresis_initial_seeds_CSC_CORRECTED.png

    5. Ras-derivative profiles:
         - Perform a Ras ↑ sweep on the calibrated parameters
           (following the benign trajectory).
         - Approximate dX/dRas for X ∈ {C, A, T, R, L, M} using centered
           finite differences along the sweep.
         - Save:
               results/model/ras_derivatives_CORRECTED.csv
         - Plot:
               figures/model/ras_derivatives_CORRECTED.png
"""

from __future__ import annotations

# Import typing utilities for clearer type hints
from typing import Dict, Any, List, Tuple

# Import standard library modules
import json
import sys
from pathlib import Path

# Import numerical libraries
import numpy as np
import pandas as pd

# Import plotting libraries for figure generation
import matplotlib.pyplot as plt
import seaborn as sns


# ----------------------------------------------------------------------
# Ensure src/ is on the Python path and import core model utilities
# ----------------------------------------------------------------------

CURRENT_DIR: Path = Path(__file__).resolve().parent
SRC_ROOT: Path = CURRENT_DIR.parent

if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

try:
    from ras_csc_core import simulate_steady_state  # type: ignore
except ImportError as exc:  # pragma: no cover - defensive
    raise ImportError(
        "Could not import 'simulate_steady_state' from 'ras_csc_core'. "
        "Ensure that ras_csc_core.py is located in src/ and that this "
        "script is in src/scripts/."
    ) from exc


# ======================================================================
# GLOBAL STYLE (FIGURES)
# ======================================================================

sns.set_style("whitegrid")

plt.rcParams.update(
    {
        "figure.dpi": 300,
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    }
)

COLOR_BASELINE = "#0072B2"
COLOR_LEPR = "#D55E00"
COLOR_MTOR = "#009E73"


# ======================================================================
# PATH CONSTANTS
# ======================================================================

ROOT: Path = SRC_ROOT.parent

RESULTS_DIR: Path = ROOT / "results"
CALIBRATION_DIR: Path = RESULTS_DIR / "calibration"
SENSITIVITY_DIR: Path = RESULTS_DIR / "sensitivity"
PERTURBATION_DIR: Path = RESULTS_DIR / "perturbations"
MODEL_RESULTS_DIR: Path = RESULTS_DIR / "model"

FIG_ROOT: Path = ROOT / "figures"
FIG_SENSITIVITY_DIR: Path = FIG_ROOT / "sensitivity"
FIG_MODEL_DIR: Path = FIG_ROOT / "model"

SENSITIVITY_DIR.mkdir(parents=True, exist_ok=True)
PERTURBATION_DIR.mkdir(parents=True, exist_ok=True)
MODEL_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_SENSITIVITY_DIR.mkdir(parents=True, exist_ok=True)
FIG_MODEL_DIR.mkdir(parents=True, exist_ok=True)

PARAM_JSON: Path = CALIBRATION_DIR / "optimized_parameters_CORRECTED.json"


# ======================================================================
# PARAMETER LOADING
# ======================================================================

def load_parameters() -> Dict[str, float]:
    """
    Load optimized parameters from JSON and return a float-valued dict.

    Returns
    -------
    Dict[str, float]
        Dictionary mapping parameter names to float values.

    Raises
    ------
    FileNotFoundError
        If the parameter JSON is missing.
    ValueError
        If the JSON is empty, malformed, or contains non-numeric values.
    """
    if not PARAM_JSON.exists():
        raise FileNotFoundError(
            f"Optimized parameters JSON not found at: {PARAM_JSON}. "
            "Run run_model_calibration.py first."
        )

    with PARAM_JSON.open("r", encoding="utf-8") as f:
        raw_params: Any = json.load(f)

    if not isinstance(raw_params, dict) or not raw_params:
        raise ValueError(
            f"Parameter file {PARAM_JSON} is empty or malformed."
        )

    clean_params: Dict[str, float] = {}
    for key, value in raw_params.items():
        try:
            clean_params[str(key)] = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Parameter '{key}' with value '{value}' could not be "
                "converted to float."
            ) from exc

    print(f"[INFO] Loaded {len(clean_params)} parameters from {PARAM_JSON}")
    return clean_params


# ======================================================================
# ANALYSIS 0: CALIBRATION RSS SUMMARY (RESIDUALS PER MODULE)
# ======================================================================

def _find_calibration_residuals_file() -> Path:
    """
    Try to locate a calibration residuals CSV in the calibration directory.

    This is intentionally defensive: it tries a few reasonable filenames
    and raises a clear error if none are found.
    """
    candidate_names = [
        # Preferred: explicit residuals CSV written by run_model_calibration.py
        "model_calibration_residuals_CORRECTED.csv",
        "model_calibration_residuals.csv",
        # Fallbacks: full model-vs-data tables that contain a 'residual' column
        "model_vs_data_CORRECTED.csv",
        "model_vs_data.csv",
        # Older names if you ever use them
        "model_calibration_fits_CORRECTED.csv",
        "model_calibration_fits.csv",
    ]

    for name in candidate_names:
        candidate = CALIBRATION_DIR / name
        if candidate.exists():
            return candidate

    raise FileNotFoundError(
        "Could not find a calibration residuals CSV in "
        f"{CALIBRATION_DIR}. Expected one of: "
        + ", ".join(candidate_names)
        + ".\nPlease ensure that run_model_calibration.py writes a "
        "per-point residual table and update the filename here if needed."
    )


def load_calibration_residuals() -> pd.DataFrame:
    """
    Load per-point calibration residuals from CSV.

    The file is expected to contain at least:
        - a column with residuals (named 'residual' or 'resid').
        - optionally a 'module' or 'variable' column to group by.

    Returns
    -------
    pd.DataFrame
        Residual table with standardized column names:
            'module' (if available) and 'residual'.
    """
    resid_path = _find_calibration_residuals_file()
    print(f"[INFO] Using calibration residuals file: {resid_path}")

    df = pd.read_csv(resid_path)

    # Normalise residual column name
    if "residual" in df.columns:
        resid_col = "residual"
    elif "resid" in df.columns:
        resid_col = "resid"
    else:
        raise ValueError(
            f"Calibration residuals file {resid_path} does not contain a "
            "'residual' or 'resid' column."
        )
    df = df.rename(columns={resid_col: "residual"})

    # Normalise module/variable name if present
    if "module" in df.columns:
        df = df.rename(columns={"module": "module"})
    elif "variable" in df.columns:
        df = df.rename(columns={"variable": "module"})
    else:
        df["module"] = "ALL"

    # Drop any non-finite residuals
    df = df[np.isfinite(df["residual"].to_numpy())].copy()
    if df.empty:
        raise ValueError("Calibration residuals table is empty after filtering.")

    return df


def run_calibration_rss_analysis() -> pd.DataFrame:
    """
    Compute RSS per module and globally from calibration residuals.

    Outputs
    -------
    - results/sensitivity/calibration_rss_CORRECTED.csv
    - figures/sensitivity/calibration_rss_by_module_CORRECTED.png
    - figures/sensitivity/calibration_residuals_distribution_CORRECTED.png
    - Printed textual summary of fit quality per module.

    Returns
    -------
    pd.DataFrame
        Table with columns:
            module, rss, n_points, mse
    """
    # Load per-point residuals (module + residual columns)
    df_resid = load_calibration_residuals()

    # Group residuals by module
    group = df_resid.groupby("module", as_index=False)

    # Compute RSS for each module
    rss_df = group["residual"].apply(
        lambda s: float(np.sum(np.square(s.to_numpy(dtype=float))))
    ).rename(columns={"residual": "rss"})

    # Number of points per module
    size_df = group.size()
    rss_df["n_points"] = size_df["size"].astype(int)

    # MSE per module
    rss_df["mse"] = rss_df["rss"] / rss_df["n_points"].replace(0, np.nan)

    # Global metrics
    global_rss = float(
        np.sum(np.square(df_resid["residual"].to_numpy(dtype=float))))
    global_n = int(df_resid.shape[0])
    global_mse = global_rss / global_n if global_n > 0 else float("nan")

    global_row = pd.DataFrame(
        [
            {
                "module": "GLOBAL",
                "rss": global_rss,
                "n_points": global_n,
                "mse": global_mse,
            }
        ]
    )

    out_df = pd.concat([rss_df, global_row], ignore_index=True)

    # Save RSS summary
    out_csv = SENSITIVITY_DIR / "calibration_rss_CORRECTED.csv"
    out_df.to_csv(out_csv, index=False)
    print(f"[SAVED] Calibration RSS summary to {out_csv}")

    # ------------------------------------------------------------------
    # Figure 1: RSS per module (fancier barplot with annotations)
    # ------------------------------------------------------------------
    plot_df = rss_df.sort_values("rss", ascending=False).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(4.8, 3.2), constrained_layout=False)

    # Use a gradient palette so higher RSS visually stands out
    palette = sns.color_palette("Blues", n_colors=plot_df.shape[0])

    sns.barplot(
        data=plot_df,
        x="module",
        y="rss",
        palette=palette,
        ax=ax,
        edgecolor="black",
    )

    ax.set_ylabel("RSS (sum of squared residuals)")
    ax.set_xlabel("Module")
    ax.set_title("Calibration fit per module (RSS)")
    ax.tick_params(axis="x", rotation=45)
    ax.grid(axis="y", alpha=0.3)

    # Annotate bars with RSS values (compact, scientific if needed)
    for patch, (_, row) in zip(ax.patches, plot_df.iterrows()):
        height = patch.get_height()
        if height <= 0 or not np.isfinite(height):
            continue
        label = f"{height:.2g}" if height > 100 else f"{height:.2f}"
        ax.annotate(
            label,
            xy=(patch.get_x() + patch.get_width() / 2.0, height),
            xytext=(0, 2),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=7,
            rotation=90,
        )

    fig.subplots_adjust(left=0.18, right=0.98, top=0.88, bottom=0.30)
    out_png_rss = FIG_SENSITIVITY_DIR / "calibration_rss_by_module_CORRECTED.png"
    fig.savefig(out_png_rss)
    plt.close(fig)
    print(f"[SAVED] Calibration RSS figure to {out_png_rss}")

    # ------------------------------------------------------------------
    # Figure 2: Fancier residual plots
    #   Left: global distribution (hist + KDE + mean ± SD lines)
    #   Right: violin + box overlay per module + jittered points
    # ------------------------------------------------------------------
    fig2, axes = plt.subplots(
        1,
        2,
        figsize=(8.0, 3.4),
        constrained_layout=False,
        sharey=False,
    )

    # Global statistics
    resid_vals = df_resid["residual"].to_numpy(dtype=float)
    resid_mean = float(np.mean(resid_vals))
    resid_std = float(np.std(resid_vals))

    # Left: histogram + KDE + mean / ±1 SD markers
    sns.histplot(
        df_resid["residual"],
        bins=15,
        kde=True,
        ax=axes[0],
        stat="count",
        edgecolor="black",
        linewidth=0.5,
    )
    axes[0].axvline(0.0, color="grey", linestyle="--",
                    linewidth=0.8, label="0")
    axes[0].axvline(
        resid_mean,
        color="black",
        linestyle="-",
        linewidth=1.0,
        label="mean",
    )
    axes[0].axvline(
        resid_mean + resid_std,
        color="black",
        linestyle=":",
        linewidth=0.9,
        label="mean ± SD",
    )
    axes[0].axvline(
        resid_mean - resid_std,
        color="black",
        linestyle=":",
        linewidth=0.9,
    )

    axes[0].set_xlabel("Residual (model − data, z-score)")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Global residual distribution")
    axes[0].grid(alpha=0.3)

    # Compact text box with summary stats
    text_str = f"mean = {resid_mean:.3f}\nSD   = {resid_std:.3f}\nN    = {len(resid_vals)}"
    axes[0].text(
        0.97,
        0.97,
        text_str,
        transform=axes[0].transAxes,
        ha="right",
        va="top",
        fontsize=7,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8),
    )

    axes[0].legend(frameon=False, fontsize=7)

    # Right: violin + boxplot overlay + jitter
    # Sort modules by median residual absolute value so ordering is informative
    median_abs = (
        df_resid.groupby("module")["residual"]
        .apply(lambda s: float(np.median(np.abs(s.to_numpy(dtype=float)))))
        .sort_values(ascending=False)
    )
    ordered_modules = median_abs.index.tolist()
    df_resid["module"] = pd.Categorical(
        df_resid["module"], categories=ordered_modules, ordered=True
    )

    # Violin for shape
    sns.violinplot(
        data=df_resid,
        x="module",
        y="residual",
        inner=None,
        cut=0,
        scale="width",
        linewidth=0.7,
        ax=axes[1],
        palette="pastel",
    )

    # Narrow boxplot overlaid for median / IQR
    sns.boxplot(
        data=df_resid,
        x="module",
        y="residual",
        whis=1.5,
        width=0.3,
        fliersize=0,
        boxprops={"facecolor": "white", "zorder": 3},
        medianprops={"color": "black", "linewidth": 1.2},
        whiskerprops={"linewidth": 0.8},
        capprops={"linewidth": 0.8},
        ax=axes[1],
    )

    # Jittered points (low alpha so it does not clutter)
    sns.stripplot(
        data=df_resid,
        x="module",
        y="residual",
        color="black",
        size=2.5,
        alpha=0.4,
        jitter=0.2,
        ax=axes[1],
        zorder=4,
    )

    axes[1].axhline(0.0, color="grey", linestyle="--", linewidth=0.8)
    axes[1].set_xlabel("Module")
    axes[1].set_ylabel("Residual (model − data, z-score)")
    axes[1].set_title("Residuals by module")
    axes[1].tick_params(axis="x", rotation=45)
    axes[1].grid(alpha=0.3, axis="y")

    fig2.subplots_adjust(
        left=0.08,
        right=0.99,
        top=0.88,
        bottom=0.24,
        wspace=0.35,
    )
    out_png_resid = (
        FIG_SENSITIVITY_DIR / "calibration_residuals_distribution_CORRECTED.png"
    )
    fig2.savefig(out_png_resid)
    plt.close(fig2)
    print(f"[SAVED] Calibration residuals figure to {out_png_resid}")

    # ------------------------------------------------------------------
    # Textual summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("CALIBRATION RSS SUMMARY")
    print("=" * 70)
    print(
        f"Global RSS: {global_rss:.3f} over {global_n} data points "
        f"(MSE = {global_mse:.4f})"
    )
    print("\nPer-module RSS (sorted, highest to lowest):")
    for _, row in plot_df.iterrows():
        module = str(row["module"])
        rss_val = float(row["rss"])
        n_val = int(row["n_points"])
        mse_val = float(row["mse"])
        print(
            f"  {module:15s}  RSS = {rss_val:8.3f}  "
            f"n = {n_val:5d}  MSE = {mse_val:7.4f}"
        )
    print("=" * 70 + "\n")

    return out_df


# ======================================================================
# HELPER: RAS HYSTERESIS GAP COMPUTATION
# ======================================================================

def compute_ras_hysteresis_gap(
    params: Dict[str, float],
    n_points: int = 50,
    f_ras_min: float = 0.0,
    f_ras_max: float = 1.0,
) -> float:
    """
    Compute the maximum CSC gap between Ras ↑ and Ras ↓ branches.

    Returns
    -------
    float
        Maximum absolute |C_down - C_up| in central Ras range.
    """
    ras_up: np.ndarray = np.linspace(f_ras_min, f_ras_max, int(n_points))
    ras_down: np.ndarray = np.linspace(f_ras_max, f_ras_min, int(n_points))

    C_up_vals: List[float] = []
    C_down_vals: List[float] = []
    f_up_vals: List[float] = []
    f_down_vals: List[float] = []

    y_current_up: List[float] = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    for f_ras in ras_up:
        ss_up: np.ndarray = simulate_steady_state(
            f_ras=float(f_ras),
            params=params,
            y0=y_current_up,
            t_max=150.0,
            n_steps=1500,
        )
        y_current_up = ss_up.tolist()
        f_up_vals.append(float(f_ras))
        C_up_vals.append(float(ss_up[0]))

    y_current_down: List[float] = [0.8, 0.8, 0.5, 0.5, 0.5, 0.8]
    for f_ras in ras_down:
        ss_down: np.ndarray = simulate_steady_state(
            f_ras=float(f_ras),
            params=params,
            y0=y_current_down,
            t_max=150.0,
            n_steps=1500,
        )
        y_current_down = ss_down.tolist()
        f_down_vals.append(float(f_ras))
        C_down_vals.append(float(ss_down[0]))

    f_up_arr: np.ndarray = np.asarray(f_up_vals, dtype=float)
    f_down_arr: np.ndarray = np.asarray(f_down_vals, dtype=float)
    C_up_arr: np.ndarray = np.asarray(C_up_vals, dtype=float)
    C_down_arr: np.ndarray = np.asarray(C_down_vals, dtype=float)

    central_ras: np.ndarray = np.linspace(0.2, 0.8, 25)
    max_gap: float = 0.0

    for fr in central_ras:
        up_mask: np.ndarray = np.abs(f_up_arr - fr) < 0.02
        down_mask: np.ndarray = np.abs(f_down_arr - fr) < 0.02

        C_up_mean: float = float(
            C_up_arr[up_mask].mean()) if np.any(up_mask) else np.nan
        C_down_mean: float = float(
            C_down_arr[down_mask].mean()) if np.any(down_mask) else np.nan

        if np.isnan(C_up_mean) or np.isnan(C_down_mean):
            continue

        gap: float = float(abs(C_down_mean - C_up_mean))
        if gap > max_gap:
            max_gap = gap

    return max_gap


# ======================================================================
# HELPER: HIGH-RAS PERTURBATION EXPERIMENT
# ======================================================================

def compute_high_ras_perturbation_effects(
    params: Dict[str, float],
    f_ras: float = 0.9,
    mu_L_scale: float = 0.3,
    eta_M_scale: float = 0.3,
) -> Dict[str, float]:
    """
    Run a high-Ras LEPR / mTOR perturbation experiment and quantify ΔC/ΔM.

    Returns
    -------
    Dict[str, float]
        Baseline and perturbed C/M values and their differences.
    """
    y0_high: List[float] = [0.8, 0.8, 0.5, 0.5, 0.5, 0.8]

    ss_base: np.ndarray = simulate_steady_state(
        f_ras=float(f_ras),
        params=params,
        y0=y0_high,
        t_max=200.0,
        n_steps=2000,
    )

    C_base: float = float(ss_base[0])
    M_base: float = float(ss_base[5])

    params_lepr: Dict[str, float] = dict(params)
    params_lepr["mu_L"] = params_lepr.get("mu_L", 0.9) * float(mu_L_scale)
    ss_lepr: np.ndarray = simulate_steady_state(
        f_ras=float(f_ras),
        params=params_lepr,
        y0=y0_high,
        t_max=200.0,
        n_steps=2000,
    )
    C_lepr: float = float(ss_lepr[0])
    M_lepr: float = float(ss_lepr[5])

    params_mtor: Dict[str, float] = dict(params)
    params_mtor["eta_M"] = params_mtor.get("eta_M", 1.0) * float(eta_M_scale)
    ss_mtor: np.ndarray = simulate_steady_state(
        f_ras=float(f_ras),
        params=params_mtor,
        y0=y0_high,
        t_max=200.0,
        n_steps=2000,
    )
    C_mtor: float = float(ss_mtor[0])
    M_mtor: float = float(ss_mtor[5])

    delta_C_lepr: float = C_base - C_lepr
    delta_M_lepr: float = M_base - M_lepr
    delta_C_mtor: float = C_base - C_mtor
    delta_M_mtor: float = M_base - M_mtor

    return {
        "C_base": C_base,
        "M_base": M_base,
        "C_lepr": C_lepr,
        "M_lepr": M_lepr,
        "delta_C_lepr": delta_C_lepr,
        "delta_M_lepr": delta_M_lepr,
        "C_mtor": C_mtor,
        "M_mtor": M_mtor,
        "delta_C_mtor": delta_C_mtor,
        "delta_M_mtor": delta_M_mtor,
    }


# ======================================================================
# ANALYSIS 1: HYSTERESIS SENSITIVITY
# ======================================================================

def run_hysteresis_sensitivity(
    params: Dict[str, float],
    scales: Tuple[float, float, float] = (0.8, 1.0, 1.2),
    n_points: int = 50,
) -> pd.DataFrame:
    """
    Evaluate how Ras hysteresis strength changes with parameter scaling.

    Returns
    -------
    pd.DataFrame
        Columns: param, baseline_value, scale, scaled_value,
                 max_gap, delta_gap_vs_baseline
    """
    scale_list: List[float] = [float(s) for s in scales]
    rows: List[Dict[str, Any]] = []

    baseline_gaps: Dict[str, float] = {}

    for param_name, base_value in params.items():
        if abs(base_value) < 1e-8:
            continue

        baseline_gap: float = np.nan

        for scale in scale_list:
            local_params: Dict[str, float] = dict(params)
            scaled_value: float = base_value * float(scale)
            local_params[param_name] = scaled_value

            try:
                max_gap: float = compute_ras_hysteresis_gap(
                    params=local_params,
                    n_points=int(n_points),
                    f_ras_min=0.0,
                    f_ras_max=1.0,
                )
            except Exception as exc:  # pragma: no cover - defensive
                print(
                    f"[WARN] Hysteresis gap computation failed for "
                    f"param={param_name}, scale={scale:.3f}: {exc}"
                )
                max_gap = float("nan")

            if abs(scale - 1.0) < 1e-9:
                baseline_gap = max_gap
                baseline_gaps[param_name] = max_gap

            if not np.isnan(baseline_gap):
                delta_gap: float = max_gap - baseline_gap
            else:
                delta_gap = float("nan")

            rows.append(
                {
                    "param": param_name,
                    "baseline_value": float(base_value),
                    "scale": float(scale),
                    "scaled_value": float(scaled_value),
                    "max_gap": float(max_gap),
                    "delta_gap_vs_baseline": float(delta_gap),
                }
            )

    df: pd.DataFrame = pd.DataFrame(rows)
    out_csv: Path = SENSITIVITY_DIR / "hysteresis_sensitivity_CORRECTED.csv"
    df.to_csv(out_csv, index=False)
    print(f"[SAVED] Hysteresis sensitivity table to {out_csv}")
    return df


def plot_hysteresis_sensitivity(df: pd.DataFrame) -> None:
    """
    Plot hysteresis sensitivity for the most influential parameters.

    - Ranks parameters by max |delta_gap_vs_baseline| across non-1 scales.
    - Plots top 6 in a 2×3 grid, max_gap vs scale.
    """
    if df.empty:
        print("[WARN] Hysteresis sensitivity DataFrame is empty; skipping plot.")
        return

    non_baseline = df[df["scale"] != 1.0].copy()
    non_baseline["abs_delta"] = non_baseline["delta_gap_vs_baseline"].abs()
    rank_df = (
        non_baseline.groupby("param", as_index=False)["abs_delta"]
        .max()
        .sort_values("abs_delta", ascending=False)
    )
    top_params = rank_df["param"].head(6).tolist()

    plot_df = df[df["param"].isin(top_params)].copy()
    plot_df["scale_str"] = plot_df["scale"].map(lambda x: f"{x:.1f}")

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(9.0, 4.0),
        constrained_layout=False,
        sharex=True,
    )
    axes = axes.ravel()

    for ax, param in zip(axes, top_params):
        sub = plot_df[plot_df["param"] == param].sort_values("scale")
        ax.plot(
            sub["scale"].values,
            sub["max_gap"].values,
            marker="o",
            linewidth=1.6,
        )
        ax.set_xticks(sub["scale"].values)
        ax.set_xlabel("Scale")
        ax.set_ylabel("Max CSC gap")
        ax.set_title(param)

        ax.axvline(1.0, color="grey", linestyle="--", linewidth=0.8)
        ax.grid(alpha=0.3)

    for k in range(len(top_params), len(axes)):
        fig.delaxes(axes[k])

    fig.suptitle("Hysteresis sensitivity for top parameters", y=0.98)
    fig.subplots_adjust(
        left=0.07,
        right=0.98,
        top=0.9,
        bottom=0.12,
        wspace=0.35,
        hspace=0.45,
    )

    out_png: Path = FIG_SENSITIVITY_DIR / "hysteresis_sensitivity_CORRECTED.png"
    fig.savefig(out_png)
    plt.close(fig)
    print(f"[SAVED] Hysteresis sensitivity figure to {out_png}")


# ======================================================================
# ANALYSIS 2: PERTURBATION SENSITIVITY
# ======================================================================

def run_perturbation_sensitivity(
    params: Dict[str, float],
    scales: Tuple[float, float, float] = (0.8, 1.0, 1.2),
    f_ras: float = 0.9,
    mu_L_scale: float = 0.3,
    eta_M_scale: float = 0.3,
) -> pd.DataFrame:
    """
    Evaluate how LEPR/mTOR perturbation effects change with parameter scaling.

    Returns
    -------
    pd.DataFrame
        Columns: param, baseline_value, scale, scaled_value,
                 C_base, M_base, C_lepr, M_lepr, delta_C_lepr, delta_M_lepr,
                 C_mtor, M_mtor, delta_C_mtor, delta_M_mtor
    """
    scale_list: List[float] = [float(s) for s in scales]
    rows: List[Dict[str, Any]] = []

    for param_name, base_value in params.items():
        if abs(base_value) < 1e-8:
            continue

        for scale in scale_list:
            local_params: Dict[str, float] = dict(params)
            scaled_value: float = base_value * float(scale)
            local_params[param_name] = scaled_value

            try:
                result: Dict[str, float] = compute_high_ras_perturbation_effects(
                    params=local_params,
                    f_ras=float(f_ras),
                    mu_L_scale=float(mu_L_scale),
                    eta_M_scale=float(eta_M_scale),
                )
            except Exception as exc:  # pragma: no cover - defensive
                print(
                    f"[WARN] Perturbation sensitivity failed for "
                    f"param={param_name}, scale={scale:.3f}: {exc}"
                )
                result = {
                    "C_base": float("nan"),
                    "M_base": float("nan"),
                    "C_lepr": float("nan"),
                    "M_lepr": float("nan"),
                    "delta_C_lepr": float("nan"),
                    "delta_M_lepr": float("nan"),
                    "C_mtor": float("nan"),
                    "M_mtor": float("nan"),
                    "delta_C_mtor": float("nan"),
                    "delta_M_mtor": float("nan"),
                }

            row: Dict[str, Any] = {
                "param": param_name,
                "baseline_value": float(base_value),
                "scale": float(scale),
                "scaled_value": float(scaled_value),
            }
            row.update(result)
            rows.append(row)

    df: pd.DataFrame = pd.DataFrame(rows)
    out_csv: Path = SENSITIVITY_DIR / "perturbation_sensitivity_CORRECTED.csv"
    df.to_csv(out_csv, index=False)
    print(f"[SAVED] Perturbation sensitivity table to {out_csv}")
    return df


def plot_perturbation_sensitivity(df: pd.DataFrame) -> None:
    """
    Plot perturbation sensitivity focusing on ΔC for LEPR↓ and mTOR↓.

    - Ranks parameters by max |ΔC_lepr| and |ΔC_mtor| across scales.
    - Plots top 5 parameters in a 1×2 panel (LEPR, mTOR).
    """
    if df.empty:
        print("[WARN] Perturbation sensitivity DataFrame is empty; skipping plot.")
        return

    lepr_rank = (
        df.groupby("param", as_index=False)["delta_C_lepr"]
        .apply(lambda s: s.abs().max())
        .rename(columns={"delta_C_lepr": "max_abs_delta_lepr"})
        .sort_values("max_abs_delta_lepr", ascending=False)
    )
    mtor_rank = (
        df.groupby("param", as_index=False)["delta_C_mtor"]
        .apply(lambda s: s.abs().max())
        .rename(columns={"delta_C_mtor": "max_abs_delta_mtor"})
        .sort_values("max_abs_delta_mtor", ascending=False)
    )

    top_lepr = lepr_rank["param"].head(5).tolist()
    top_mtor = mtor_rank["param"].head(5).tolist()

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(9.0, 3.2),
        constrained_layout=False,
        sharey=False,
    )

    lepr_df = df[df["param"].isin(top_lepr)].copy()
    mtor_df = df[df["param"].isin(top_mtor)].copy()

    for ax, sub, title, col, color in [
        (axes[0], lepr_df, "ΔC under LEPR perturbation", "delta_C_lepr", COLOR_LEPR),
        (axes[1], mtor_df, "ΔC under mTOR perturbation", "delta_C_mtor", COLOR_MTOR),
    ]:
        for param in sub["param"].unique():
            p_df = sub[sub["param"] == param].sort_values("scale")
            ax.plot(
                p_df["scale"].values,
                p_df[col].values,
                marker="o",
                linewidth=1.6,
                label=param,
            )
        ax.axhline(0.0, color="grey", linestyle="--", linewidth=0.8)
        ax.set_xlabel("Scale")
        ax.set_ylabel("ΔCSC (C_base − C_perturbed)")
        ax.set_title(title)
        ax.grid(alpha=0.3)

    axes[0].legend(
        title="Top LEPR-sensitive params",
        fontsize=7,
        loc="upper right",
        frameon=False,
    )
    axes[1].legend(
        title="Top mTOR-sensitive params",
        fontsize=7,
        loc="upper right",
        frameon=False,
    )

    fig.subplots_adjust(
        left=0.07,
        right=0.99,
        top=0.88,
        bottom=0.15,
        wspace=0.35,
    )
    fig.suptitle("Sensitivity of LEPR/mTOR perturbations to parameter scaling", y=0.98)

    out_png: Path = FIG_SENSITIVITY_DIR / "perturbation_sensitivity_CORRECTED.png"
    fig.savefig(out_png)
    plt.close(fig)
    print(f"[SAVED] Perturbation sensitivity figure to {out_png}")


# ======================================================================
# ANALYSIS 3: LEPR–MTOR PERTURBATION GRID
# ======================================================================

def run_perturbation_grid(
    params: Dict[str, float],
    f_ras: float = 0.9,
    mu_scales: Tuple[float, ...] = (0.3, 0.5, 0.7, 1.0),
    eta_scales: Tuple[float, ...] = (0.3, 0.5, 0.7, 1.0),
) -> pd.DataFrame:
    """
    Generate a 2D grid of LEPR/mTOR perturbations and record steady states.

    Returns
    -------
    pd.DataFrame
        One row per (mu_L_scale, eta_M_scale) with module values and ΔC.
    """
    y0_high: List[float] = [0.8, 0.8, 0.5, 0.5, 0.5, 0.8]

    ss_base: np.ndarray = simulate_steady_state(
        f_ras=float(f_ras),
        params=params,
        y0=y0_high,
        t_max=200.0,
        n_steps=2000,
    )

    C_base: float = float(ss_base[0])
    rows: List[Dict[str, Any]] = []

    for mu_scale in mu_scales:
        for eta_scale in eta_scales:
            local_params: Dict[str, float] = dict(params)
            local_params["mu_L"] = local_params.get("mu_L", 0.9) * float(mu_scale)
            local_params["eta_M"] = local_params.get("eta_M", 1.0) * float(eta_scale)

            try:
                ss: np.ndarray = simulate_steady_state(
                    f_ras=float(f_ras),
                    params=local_params,
                    y0=y0_high,
                    t_max=200.0,
                    n_steps=2000,
                )
            except Exception as exc:  # pragma: no cover - defensive
                print(
                    f"[WARN] Perturbation grid simulation failed for "
                    f"mu_L_scale={mu_scale:.3f}, eta_M_scale={eta_scale:.3f}: {exc}"
                )
                ss = np.array([np.nan] * 6, dtype=float)

            C_val: float = float(ss[0])
            A_val: float = float(ss[1])
            T_val: float = float(ss[2])
            R_val: float = float(ss[3])
            L_val: float = float(ss[4])
            M_val: float = float(ss[5])

            delta_C: float = C_base - C_val

            rows.append(
                {
                    "f_ras": float(f_ras),
                    "mu_L_scale": float(mu_scale),
                    "eta_M_scale": float(eta_scale),
                    "C": C_val,
                    "A": A_val,
                    "T": T_val,
                    "R": R_val,
                    "L": L_val,
                    "M": M_val,
                    "C_base": C_base,
                    "delta_C_vs_base": delta_C,
                }
            )

    df: pd.DataFrame = pd.DataFrame(rows)
    out_csv: Path = PERTURBATION_DIR / "perturbation_grid_CORRECTED.csv"
    df.to_csv(out_csv, index=False)
    print(f"[SAVED] LEPR–mTOR perturbation grid to {out_csv}")
    return df


def plot_perturbation_grid(df: pd.DataFrame) -> None:
    """
    Plot ΔC vs baseline across the LEPR–mTOR perturbation grid as a heatmap.
    """
    if df.empty:
        print("[WARN] Perturbation grid DataFrame is empty; skipping plot.")
        return

    pivot = df.pivot_table(
        index="mu_L_scale",
        columns="eta_M_scale",
        values="delta_C_vs_base",
    )

    fig, ax = plt.subplots(figsize=(4.2, 3.5), constrained_layout=False)
    im = ax.imshow(
        pivot.values,
        origin="lower",
        aspect="auto",
        cmap="viridis",
    )

    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([f"{v:.2f}" for v in pivot.columns])
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels([f"{v:.2f}" for v in pivot.index])

    ax.set_xlabel("eta_M scale")
    ax.set_ylabel("mu_L scale")
    ax.set_title("ΔCSC (C_base − C_perturbed) across LEPR–mTOR grid")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("ΔCSC vs baseline")

    fig.subplots_adjust(
        left=0.12,
        right=0.98,
        top=0.9,
        bottom=0.15,
    )
    out_png: Path = FIG_SENSITIVITY_DIR / "perturbation_grid_deltaC_CORRECTED.png"
    fig.savefig(out_png)
    plt.close(fig)
    print(f"[SAVED] Perturbation grid heatmap to {out_png}")


# ======================================================================
# ANALYSIS 4: HYSTERESIS UNDER MULTIPLE INITIAL CONDITIONS
# ======================================================================

def run_hysteresis_initial_condition_sweeps(
    params: Dict[str, float],
    n_points: int = 50,
) -> pd.DataFrame:
    """
    Repeat Ras hysteresis sweeps under multiple initial-condition pairs.

    Returns
    -------
    pd.DataFrame
        Columns: seed_label, direction, f_ras, C, A, T, R, L, M
    """
    seed_configs: List[Dict[str, Any]] = [
        {
            "label": "default",
            "y0_up": [0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
            "y0_down": [0.8, 0.8, 0.5, 0.5, 0.5, 0.8],
        },
        {
            "label": "slightly_more_benign",
            "y0_up": [0.05, 0.05, 0.05, 0.05, 0.05, 0.05],
            "y0_down": [0.8, 0.8, 0.6, 0.6, 0.6, 0.8],
        },
        {
            "label": "slightly_more_malignant",
            "y0_up": [0.2, 0.2, 0.1, 0.1, 0.1, 0.2],
            "y0_down": [0.9, 0.9, 0.6, 0.6, 0.6, 0.9],
        },
    ]

    ras_up: np.ndarray = np.linspace(0.0, 1.0, int(n_points))
    ras_down: np.ndarray = np.linspace(1.0, 0.0, int(n_points))

    rows: List[Dict[str, Any]] = []

    for cfg in seed_configs:
        label: str = str(cfg["label"])
        y_up: List[float] = list(cfg["y0_up"])
        y_down: List[float] = list(cfg["y0_down"])

        for f_ras in ras_up:
            ss_up: np.ndarray = simulate_steady_state(
                f_ras=float(f_ras),
                params=params,
                y0=y_up,
                t_max=150.0,
                n_steps=1500,
            )
            y_up = ss_up.tolist()
            rows.append(
                {
                    "seed_label": label,
                    "direction": "up",
                    "f_ras": float(f_ras),
                    "C": float(ss_up[0]),
                    "A": float(ss_up[1]),
                    "T": float(ss_up[2]),
                    "R": float(ss_up[3]),
                    "L": float(ss_up[4]),
                    "M": float(ss_up[5]),
                }
            )

        for f_ras in ras_down:
            ss_down: np.ndarray = simulate_steady_state(
                f_ras=float(f_ras),
                params=params,
                y0=y_down,
                t_max=150.0,
                n_steps=1500,
            )
            y_down = ss_down.tolist()
            rows.append(
                {
                    "seed_label": label,
                    "direction": "down",
                    "f_ras": float(f_ras),
                    "C": float(ss_down[0]),
                    "A": float(ss_down[1]),
                    "T": float(ss_down[2]),
                    "R": float(ss_down[3]),
                    "L": float(ss_down[4]),
                    "M": float(ss_down[5]),
                }
            )

    df: pd.DataFrame = pd.DataFrame(rows)
    out_csv: Path = MODEL_RESULTS_DIR / "ras_hysteresis_initial_seeds_CORRECTED.csv"
    df.to_csv(out_csv, index=False)
    print(f"[SAVED] Hysteresis initial-condition sweeps to {out_csv}")
    return df


def plot_hysteresis_initial_conditions(df: pd.DataFrame) -> None:
    """
    Plot CSC hysteresis curves for each initial-condition seed.

    - 1×N panels, one per seed_label.
    - Ras ↑ and Ras ↓ trajectories for C only.
    """
    if df.empty:
        print("[WARN] Hysteresis initial-condition DataFrame empty; skipping plot.")
        return

    seeds = list(df["seed_label"].unique())
    n_seeds = len(seeds)

    fig, axes = plt.subplots(
        1,
        n_seeds,
        figsize=(3.2 * n_seeds, 2.8),
        constrained_layout=False,
        sharey=True,
    )
    if n_seeds == 1:
        axes = [axes]

    for ax, seed in zip(axes, seeds):
        sub = df[df["seed_label"] == seed].copy()
        up = sub[sub["direction"] == "up"].sort_values("f_ras")
        down = sub[sub["direction"] == "down"].sort_values("f_ras")

        ax.plot(
            up["f_ras"].values,
            up["C"].values,
            label="Ras ↑",
            linewidth=1.6,
            marker="o",
            markersize=3,
        )
        ax.plot(
            down["f_ras"].values,
            down["C"].values,
            label="Ras ↓",
            linewidth=1.6,
            marker="o",
            markersize=3,
        )

        ax.set_xlabel("Ras input")
        ax.set_ylabel("CSC steady state (C)")
        ax.set_title(seed.replace("_", " "))
        ax.grid(alpha=0.3)

    axes[0].legend(frameon=False, loc="upper left")

    fig.subplots_adjust(
        left=0.08,
        right=0.99,
        top=0.88,
        bottom=0.16,
        wspace=0.35,
    )
    fig.suptitle("Robustness of hysteresis to initial conditions", y=0.98)

    out_png: Path = FIG_MODEL_DIR / "ras_hysteresis_initial_seeds_CSC_CORRECTED.png"
    fig.savefig(out_png)
    plt.close(fig)
    print(f"[SAVED] Hysteresis initial-condition figure to {out_png}")


# ======================================================================
# ANALYSIS 5: RAS-DERIVATIVE PROFILES (UP BRANCH)
# ======================================================================

def run_ras_derivative_profiles(
    params: Dict[str, float],
    n_points: int = 200,
) -> pd.DataFrame:
    """
    Compute dX/dRas along the benign → malignant (up) branch.

    Returns
    -------
    pd.DataFrame
        Columns: module, f_ras, value, dX_dRas
    """
    ras_up: np.ndarray = np.linspace(0.0, 1.0, int(n_points))
    y_current: List[float] = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

    states: List[np.ndarray] = []
    for f_ras in ras_up:
        ss: np.ndarray = simulate_steady_state(
            f_ras=float(f_ras),
            params=params,
            y0=y_current,
            t_max=150.0,
            n_steps=1500,
        )
        y_current = ss.tolist()
        states.append(ss)

    states_arr: np.ndarray = np.vstack(states).astype(float)
    rows: List[Dict[str, Any]] = []

    module_names: List[str] = ["C", "A", "T", "R", "L", "M"]
    for idx, module in enumerate(module_names):
        values: np.ndarray = states_arr[:, idx].copy()
        derivs: np.ndarray = np.zeros_like(values, dtype=float)

        for i in range(1, len(ras_up) - 1):
            delta_f: float = float(ras_up[i + 1] - ras_up[i - 1])
            if abs(delta_f) < 1e-12:
                derivs[i] = 0.0
            else:
                derivs[i] = float(
                    (values[i + 1] - values[i - 1]) / delta_f
                )

        if len(ras_up) >= 2:
            delta_f_left: float = float(ras_up[1] - ras_up[0])
            derivs[0] = (
                0.0
                if abs(delta_f_left) < 1e-12
                else float((values[1] - values[0]) / delta_f_left)
            )

            delta_f_right: float = float(ras_up[-1] - ras_up[-2])
            derivs[-1] = (
                0.0
                if abs(delta_f_right) < 1e-12
                else float((values[-1] - values[-2]) / delta_f_right)
            )

        for f_ras, val, dv in zip(ras_up, values, derivs):
            rows.append(
                {
                    "module": module,
                    "f_ras": float(f_ras),
                    "value": float(val),
                    "dX_dRas": float(dv),
                }
            )

    df: pd.DataFrame = pd.DataFrame(rows)
    out_csv: Path = MODEL_RESULTS_DIR / "ras_derivatives_CORRECTED.csv"
    df.to_csv(out_csv, index=False)
    print(f"[SAVED] Ras-derivative profiles to {out_csv}")
    return df


def plot_ras_derivative_profiles(df: pd.DataFrame) -> None:
    """
    Plot dX/dRas profiles for all modules in a 2×3 panel figure.
    """
    if df.empty:
        print("[WARN] Ras-derivative DataFrame empty; skipping plot.")
        return

    modules = ["C", "A", "T", "R", "L", "M"]
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(9.0, 4.0),
        constrained_layout=False,
        sharex=True,
    )
    axes = axes.ravel()

    for ax, mod in zip(axes, modules):
        sub = df[df["module"] == mod].sort_values("f_ras")
        ax.plot(
            sub["f_ras"].values,
            sub["dX_dRas"].values,
            linewidth=1.6,
        )
        ax.axhline(0.0, color="grey", linestyle="--", linewidth=0.8)
        ax.set_xlabel("Ras input")
        ax.set_ylabel("d" + mod + "/dRas")
        ax.set_title(mod)
        ax.grid(alpha=0.3)

    fig.subplots_adjust(
        left=0.07,
        right=0.98,
        top=0.9,
        bottom=0.12,
        wspace=0.35,
        hspace=0.45,
    )
    fig.suptitle("Sensitivity of each module to Ras input (dX/dRas)", y=0.98)

    out_png: Path = FIG_MODEL_DIR / "ras_derivatives_CORRECTED.png"
    fig.savefig(out_png)
    plt.close(fig)
    print(f"[SAVED] Ras-derivative figure to {out_png}")


# ======================================================================
# MAIN ORCHESTRATION
# ======================================================================

def main() -> None:
    """
    Entry point for all Ras–CSC sensitivity and robustness analyses.

    Orchestrates:
        - Parameter loading.
        - Calibration RSS summary + figure.
        - Hysteresis sensitivity + figure.
        - Perturbation-effect sensitivity + figure.
        - LEPR–mTOR perturbation grid + heatmap.
        - Hysteresis sweeps under multiple initial conditions + figure.
        - Ras-derivative profile computation + figure.
    """
    print("\n" + "=" * 70)
    print("RUNNING RAS–CSC SENSITIVITY AND ROBUSTNESS ANALYSES")
    print("=" * 70)

    params: Dict[str, float] = load_parameters()

    print("\n[STEP] Calibration RSS summary...")
    try:
        run_calibration_rss_analysis()
    except Exception as exc:
        print(
            "[WARN] Calibration RSS summary failed; "
            "continuing with other analyses. "
            f"Details: {exc}"
        )

    print("\n[STEP] Hysteresis sensitivity...")
    df_hyst = run_hysteresis_sensitivity(params=params)
    plot_hysteresis_sensitivity(df_hyst)

    print("\n[STEP] Perturbation-effect sensitivity...")
    df_pert = run_perturbation_sensitivity(params=params)
    plot_perturbation_sensitivity(df_pert)

    print("\n[STEP] LEPR–mTOR perturbation grid...")
    df_grid = run_perturbation_grid(params=params)
    plot_perturbation_grid(df_grid)

    print("\n[STEP] Hysteresis initial-condition sweeps...")
    df_ic = run_hysteresis_initial_condition_sweeps(params=params)
    plot_hysteresis_initial_conditions(df_ic)

    print("\n[STEP] Ras-derivative profiles...")
    df_deriv = run_ras_derivative_profiles(params=params)
    plot_ras_derivative_profiles(df_deriv)

    print("\n" + "=" * 70)
    print("SENSITIVITY AND ROBUSTNESS ANALYSES COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
