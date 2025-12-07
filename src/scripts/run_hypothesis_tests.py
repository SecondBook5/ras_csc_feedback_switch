#!/usr/bin/env python3
"""
run_hypothesis_tests.py

Hypothesis testing for the Ras–CSC feedback model.

This script uses the *already calibrated* parameters to test:
    1. Bistability via Ras hysteresis:
         - Sweep Ras up from low to high (benign start).
         - Sweep Ras down from high to low (malignant start).
         - Compare CSC steady states at the same Ras input.

    2. LEPR / mTOR axis perturbations:
         - High-Ras malignant state with no perturbation.
         - Reduce leptin signaling (mu_L scaled down).
         - Reduce mTOR gain (eta_M scaled down).
         - Evaluate whether CSC fraction collapses.

All ODE logic and steady-state integration are delegated to:
    - src/ras_csc_core.py
"""

from __future__ import annotations

from typing import Dict, Any, Tuple

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ----------------------------------------------------------------------
# Ensure src/ is on the Python path
# ----------------------------------------------------------------------

CURRENT_DIR: Path = Path(__file__).resolve().parent
SRC_ROOT: Path = CURRENT_DIR.parent

if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from ras_csc_core import simulate_steady_state  # type: ignore  # noqa: E402


# ======================================================================
# GLOBAL STYLE SETTINGS
# ======================================================================

# Basic seaborn style
sns.set_style("whitegrid")

# Matplotlib rcParams tuned for journal-style figures
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

COLOR_UP = "#0072B2"     # blue – benign start
COLOR_DOWN = "#E69F00"   # orange – malignant start
COLOR_BASELINE = "#0072B2"
COLOR_LEPR = "#D55E00"
COLOR_MTOR = "#009E73"


# ======================================================================
# PATH CONSTANTS
# ======================================================================

ROOT: Path = SRC_ROOT.parent

CALIBRATION_DIR: Path = ROOT / "results" / "calibration"
HYPOTHESIS_DIR: Path = ROOT / "results" / "hypothesis_tests"
FIG_MAIN_DIR: Path = ROOT / "figures" / "main"

HYPOTHESIS_DIR.mkdir(parents=True, exist_ok=True)
FIG_MAIN_DIR.mkdir(parents=True, exist_ok=True)

PARAM_JSON: Path = CALIBRATION_DIR / "optimized_parameters_CORRECTED.json"


# ======================================================================
# PARAMETER LOADING
# ======================================================================

def load_parameters() -> Dict[str, float]:
    """
    Load optimized parameters for hypothesis testing.

    Returns
    -------
    Dict[str, float]
        Dictionary of parameter values.

    Raises
    ------
    FileNotFoundError
        If the parameter JSON file does not exist.
    ValueError
        If the JSON file is empty or malformed.
    """
    if not PARAM_JSON.exists():
        raise FileNotFoundError(
            f"Optimized parameters JSON not found at: {PARAM_JSON}. "
            "Run run_model_calibration.py first."
        )

    with PARAM_JSON.open("r", encoding="utf-8") as f:
        params_raw = json.load(f)

    if not isinstance(params_raw, dict) or not params_raw:
        raise ValueError(
            f"Parameter file {PARAM_JSON} is empty or malformed."
        )

    clean_params: Dict[str, float] = {}
    for key, value in params_raw.items():
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
# TEST 1: RAS HYSTERESIS (BISTABILITY)
# ======================================================================

def run_ras_hysteresis_test(
    params: Dict[str, float],
    n_points: int = 50,
    f_ras_min: float = 0.0,
    f_ras_max: float = 1.0,
) -> Tuple[pd.DataFrame, bool, float]:
    """
    Sweep Ras up and down to detect hysteresis in CSC steady state.

    Returns
    -------
    df : pd.DataFrame
        Table with columns direction, f_ras, C, A, T, R, L, M.
    hysteresis_detected : bool
        True if CSC up/down branches differ by > 0.1 in the central Ras range.
    max_gap : float
        Maximum absolute difference between up and down CSC branches.
    """
    print("\n" + "=" * 70)
    print("TEST 1: RAS HYSTERESIS (BISTABILITY)")
    print("=" * 70)

    ras_up = np.linspace(f_ras_min, f_ras_max, n_points)
    ras_down = np.linspace(f_ras_max, f_ras_min, n_points)

    records: list[Dict[str, Any]] = []

    y_current_up = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    print("[INFO] Sweeping Ras UP (benign → malignant)...")
    for f_ras in ras_up:
        ss_up = simulate_steady_state(
            f_ras=float(f_ras),
            params=params,
            y0=y_current_up,
            t_max=150.0,
            n_steps=1500,
        )
        y_current_up = ss_up.tolist()
        records.append(
            {
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

    y_current_down = [0.8, 0.8, 0.5, 0.5, 0.5, 0.8]
    print("[INFO] Sweeping Ras DOWN (malignant → benign)...")
    for f_ras in ras_down:
        ss_down = simulate_steady_state(
            f_ras=float(f_ras),
            params=params,
            y0=y_current_down,
            t_max=150.0,
            n_steps=1500,
        )
        y_current_down = ss_down.tolist()
        records.append(
            {
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

    df = pd.DataFrame(records)

    hysteresis_detected: bool = False
    max_gap: float = 0.0

    central_ras = np.linspace(0.2, 0.8, 25)
    for fr in central_ras:
        C_up = df[
            (df["direction"] == "up")
            & (np.abs(df["f_ras"] - fr) < 0.02)
        ]["C"].mean()

        C_down = df[
            (df["direction"] == "down")
            & (np.abs(df["f_ras"] - fr) < 0.02)
        ]["C"].mean()

        if np.isnan(C_up) or np.isnan(C_down):
            continue

        gap = float(abs(C_down - C_up))
        if gap > max_gap:
            max_gap = gap
        if gap > 0.1:
            hysteresis_detected = True

    print(f"\n[RESULT] Hysteresis detected: {hysteresis_detected}")
    print(f"  Maximum CSC gap between branches: {max_gap:.3f}")

    hysteresis_csv: Path = ROOT / "results" / "model" / "ras_hysteresis_CORRECTED.csv"
    hysteresis_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(hysteresis_csv, index=False)
    print(f"[SAVED] Hysteresis data to {hysteresis_csv}")

    return df, hysteresis_detected, max_gap


def plot_ras_hysteresis(
    df: pd.DataFrame,
    output_path: Path,
) -> None:
    """
    Plot Ras hysteresis curves for key modules and save as a PNG.

    The layout is 1×4 panels (C, A, T, M) with panel labels A–D and a
    shared legend below the axes.
    """
    modules = [("C", "CSC"), ("A", "Angiogenesis"),
               ("T", "TGFβ"), ("M", "mTOR")]
    panel_labels = ["A", "B", "C", "D"]

    fig, axes = plt.subplots(
        1,
        4,
        figsize=(10.0, 2.8),
        constrained_layout=False,
        sharex=True,
    )

    df_up = df[df["direction"] == "up"].sort_values("f_ras")
    df_down = df[df["direction"] == "down"].sort_values("f_ras")

    for ax, (mod, name), label in zip(axes, modules, panel_labels):
        ax.plot(
            df_up["f_ras"].values,
            df_up[mod].values,
            color=COLOR_UP,
            linewidth=1.6,
            marker="o",
            markersize=3,
            label="Ras ↑ (benign start)",
        )
        ax.plot(
            df_down["f_ras"].values,
            df_down[mod].values,
            color=COLOR_DOWN,
            linewidth=1.6,
            marker="o",
            markersize=3,
            label="Ras ↓ (malignant start)",
        )

        ax.set_xlabel("Ras input")
        ax.set_ylabel(f"{name} steady state")
        ax.set_title(name)

        ax.grid(alpha=0.3)

        # Panel label
        ax.text(
            -0.18,
            1.07,
            label,
            transform=ax.transAxes,
            fontsize=11,
            fontweight="bold",
            ha="left",
            va="bottom",
        )

    # Shared legend below all panels
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, -0.02),
    )

    fig.subplots_adjust(
        left=0.07,
        right=0.99,
        top=0.9,
        bottom=0.22,
        wspace=0.35,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)
    print(f"[SAVED] Hysteresis figure to {output_path}")


# ======================================================================
# TEST 2: LEPR / MTOR PERTURBATIONS
# ======================================================================

def run_perturbation_tests(
    params: Dict[str, float],
    f_ras: float = 0.9,
    mu_L_scale: float = 0.3,
    eta_M_scale: float = 0.3,
    threshold: float = 0.2,
) -> Dict[str, Any]:
    """
    Test whether perturbations to the Leptin–mTOR axis can collapse CSCs.

    Returns
    -------
    Dict[str, Any]
        Baseline and perturbed steady states and summary metrics.
    """
    print("\n" + "=" * 70)
    print("TEST 2: LEPR / MTOR PERTURBATIONS")
    print("=" * 70)

    y0_high = [0.8, 0.8, 0.5, 0.5, 0.5, 0.8]

    print("\n[INFO] Baseline (no perturbation) at high Ras...")
    ss_baseline = simulate_steady_state(
        f_ras=float(f_ras),
        params=params,
        y0=y0_high,
        t_max=200.0,
        n_steps=2000,
    )

    C_base = float(ss_baseline[0])
    M_base = float(ss_baseline[5])
    print(f"  Baseline: C={C_base:.3f}, M={M_base:.3f}")

    print(
        f"\n[INFO] Perturbation 1: Reduce leptin production (mu_L × {mu_L_scale:.2f})..."
    )
    params_lepr = dict(params)
    params_lepr["mu_L"] = params_lepr.get("mu_L", 0.9) * mu_L_scale
    ss_lepr = simulate_steady_state(
        f_ras=float(f_ras),
        params=params_lepr,
        y0=y0_high,
        t_max=200.0,
        n_steps=2000,
    )
    C_lepr = float(ss_lepr[0])
    M_lepr = float(ss_lepr[5])
    print(f"  After LEPR-axis perturbation: C={C_lepr:.3f}, M={M_lepr:.3f}")

    print(
        f"\n[INFO] Perturbation 2: Reduce mTOR gain (eta_M × {eta_M_scale:.2f})..."
    )
    params_mtor = dict(params)
    params_mtor["eta_M"] = params_mtor.get("eta_M", 1.0) * eta_M_scale
    ss_mtor = simulate_steady_state(
        f_ras=float(f_ras),
        params=params_mtor,
        y0=y0_high,
        t_max=200.0,
        n_steps=2000,
    )
    C_mtor = float(ss_mtor[0])
    M_mtor = float(ss_mtor[5])
    print(f"  After mTOR perturbation: C={C_mtor:.3f}, M={M_mtor:.3f}")

    delta_C_lepr = C_base - C_lepr
    delta_M_lepr = M_base - M_lepr
    delta_C_mtor = C_base - C_mtor
    delta_M_mtor = M_base - M_mtor

    lepr_effective = bool(delta_C_lepr > threshold)
    mtor_effective = bool(delta_C_mtor > threshold)

    print("\n[RESULT] Perturbation effects on CSC and mTOR:")
    print(
        f"  LEPR axis: ΔC={delta_C_lepr:.3f}, ΔM={delta_M_lepr:.3f} "
        f"{'→ EFFECTIVE' if lepr_effective else '→ WEAK'}"
    )
    print(
        f"  mTOR axis: ΔC={delta_C_mtor:.3f}, ΔM={delta_M_mtor:.3f} "
        f"{'→ EFFECTIVE' if mtor_effective else '→ WEAK'}"
    )

    results: Dict[str, Any] = {
        "baseline": {"C": C_base, "M": M_base},
        "lepr_perturbation": {
            "C": C_lepr,
            "M": M_lepr,
            "delta_C": delta_C_lepr,
            "delta_M": delta_M_lepr,
            "effective": lepr_effective,
        },
        "mtor_perturbation": {
            "C": C_mtor,
            "M": M_mtor,
            "delta_C": delta_C_mtor,
            "delta_M": delta_M_mtor,
            "effective": mtor_effective,
        },
    }
    return results


def plot_perturbation_bars(
    perturbation_results: Dict[str, Any],
    output_path: Path,
) -> None:
    """
    Plot CSC and mTOR steady-state responses to LEPR / mTOR perturbations.

    Creates a 1×2 panel figure:
      A – CSC steady state (C)
      B – mTOR steady state (M)
    """
    labels = ["Baseline", "LEPR↓", "mTOR↓"]

    C_vals = [
        perturbation_results["baseline"]["C"],
        perturbation_results["lepr_perturbation"]["C"],
        perturbation_results["mtor_perturbation"]["C"],
    ]
    M_vals = [
        perturbation_results["baseline"]["M"],
        perturbation_results["lepr_perturbation"]["M"],
        perturbation_results["mtor_perturbation"]["M"],
    ]

    colors = [COLOR_BASELINE, COLOR_LEPR, COLOR_MTOR]

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(6.0, 2.8),
        constrained_layout=False,
        sharey=False,
    )

    # Panel A – CSC
    ax_c = axes[0]
    ax_c.bar(labels, C_vals, color=colors, width=0.6)
    ax_c.set_ylabel("CSC steady state (C)")
    ax_c.set_title("CSC response to LEPR / mTOR perturbations")
    ax_c.set_ylim(0.0, max(C_vals) * 1.1)
    ax_c.grid(axis="y", alpha=0.3)
    ax_c.text(
        -0.25,
        1.07,
        "A",
        transform=ax_c.transAxes,
        fontsize=11,
        fontweight="bold",
        ha="left",
        va="bottom",
    )

    # Panel B – mTOR
    ax_m = axes[1]
    ax_m.bar(labels, M_vals, color=colors, width=0.6)
    ax_m.set_ylabel("mTOR steady state (M)")
    ax_m.set_title("mTOR response to LEPR / mTOR perturbations")
    ax_m.set_ylim(0.0, max(M_vals) * 1.1)
    ax_m.grid(axis="y", alpha=0.3)
    ax_m.text(
        -0.25,
        1.07,
        "B",
        transform=ax_m.transAxes,
        fontsize=11,
        fontweight="bold",
        ha="left",
        va="bottom",
    )

    for ax in axes:
        ax.set_axisbelow(True)

    fig.subplots_adjust(
        left=0.08,
        right=0.99,
        top=0.9,
        bottom=0.18,
        wspace=0.4,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)
    print(f"[SAVED] Perturbation figure to {output_path}")


# ======================================================================
# SUMMARY TABLE
# ======================================================================

def build_hypothesis_summary(
    hysteresis_detected: bool,
    max_gap: float,
    perturbation_results: Dict[str, Any],
) -> pd.DataFrame:
    """
    Build a summary table of hypothesis test outcomes.
    """
    lepr_info = perturbation_results["lepr_perturbation"]
    mtor_info = perturbation_results["mtor_perturbation"]

    summary = {
        "Test": [
            "Bistability (Ras hysteresis)",
            "LEPR-axis perturbation",
            "mTOR-axis perturbation",
        ],
        "H1_prediction": [
            "Hysteresis loop exists",
            "Perturbation collapses high-C state",
            "Perturbation collapses high-C state",
        ],
        "Result": [
            "PASS" if hysteresis_detected else "GRADED",
            "PASS" if lepr_info["effective"] else "PARTIAL",
            "PASS" if mtor_info["effective"] else "PARTIAL",
        ],
        "Metric": [
            f"Max CSC gap = {max_gap:.3f}",
            f"ΔC = {lepr_info['delta_C']:.3f}",
            f"ΔC = {mtor_info['delta_C']:.3f}",
        ],
    }

    df_summary = pd.DataFrame(summary)

    print("\n" + "=" * 70)
    print("HYPOTHESIS TEST SUMMARY")
    print("=" * 70)
    print(df_summary.to_string(index=False))

    n_pass = sum(
        [
            hysteresis_detected,
            bool(lepr_info["effective"]),
            bool(mtor_info["effective"]),
        ]
    )

    print(f"\n[VERDICT] {n_pass}/3 tests fully support H1.")
    if n_pass >= 2:
        print("  Interpretation: Strong support for the feedback-based hypothesis.")
    elif n_pass == 1:
        print("  Interpretation: Moderate support; feedback is present but limited.")
    else:
        print("  Interpretation: Weak support; model behaves more like graded response.")

    summary_csv: Path = HYPOTHESIS_DIR / "hypothesis_test_summary_CORRECTED.csv"
    df_summary.to_csv(summary_csv, index=False)
    print(f"\n[SAVED] Hypothesis summary to {summary_csv}")

    return df_summary


# ======================================================================
# MAIN
# ======================================================================

def main() -> pd.DataFrame:
    """
    Entry point for Ras–CSC hypothesis testing.
    """
    print("\n" + "=" * 70)
    print("RUNNING RAS–CSC HYPOTHESIS TESTS")
    print("=" * 70)

    params = load_parameters()

    df_hyst, hysteresis_detected, max_gap = run_ras_hysteresis_test(
        params=params
    )
    hysteresis_fig: Path = FIG_MAIN_DIR / "ras_hysteresis_CORRECTED.png"
    plot_ras_hysteresis(df=df_hyst, output_path=hysteresis_fig)

    perturbation_results = run_perturbation_tests(params=params)
    perturb_fig: Path = FIG_MAIN_DIR / "lepr_mtor_perturbations_CORRECTED.png"
    plot_perturbation_bars(
        perturbation_results=perturbation_results,
        output_path=perturb_fig,
    )

    df_summary = build_hypothesis_summary(
        hysteresis_detected=hysteresis_detected,
        max_gap=max_gap,
        perturbation_results=perturbation_results,
    )

    print("\n" + "=" * 70)
    print("HYPOTHESIS TESTING COMPLETE")
    print("=" * 70)

    return df_summary


if __name__ == "__main__":
    main()
