#!/usr/bin/env python3
"""
run_temporal_analysis.py

Temporal behavior analysis for the Ras–CSC feedback model.

This script does three things:

1) MODEL TIME COURSE
   - Load calibrated parameters.
   - Simulate a Ras "switch-on" experiment:
       * Start from a benign-like state at low Ras.
       * Jump Ras to a malignant-like value (e.g., 0.9).
       * Integrate the full ODE system for a fixed duration.
   - Save the full time course C(t), A(t), T(t), R(t), L(t), M(t).
   - Compute half-rise times t_50 for each module as a compact
     temporal ordering summary.
   - Generate a multi-panel figure of the time course with half-rise
     markers.

2) DATA PSEUDOTIME PROFILES
   - Load per-cell module scores and metadata from the scRNA pipeline.
   - Use normalized pseudotime (pseudotime_norm) as a data-derived
     progression coordinate.
   - Bin pseudotime into equal-width bins.
   - For each bin and for each module, compute mean, standard
     deviation, and cell count.
   - Save these pseudotime profiles for downstream inspection.
   - Generate a multi-panel figure of module score vs. pseudotime.

Outputs
-------
results/model/model_time_course_CORRECTED.csv
results/model/model_half_times_CORRECTED.csv
results/temporal/scrna_pseudotime_profiles_CORRECTED.csv

figures/temporal/temporal_model_time_course_CORRECTED.png
figures/temporal/temporal_scrna_pseudotime_profiles_CORRECTED.png
"""

from __future__ import annotations

# Import typing utilities
from typing import Dict, Any, List, Tuple

# Import standard library modules
import json
from pathlib import Path

# Import numerical libraries
import numpy as np
import pandas as pd
import sys

# Import plotting library
import matplotlib.pyplot as plt

CURRENT_DIR: Path = Path(__file__).resolve().parent
SRC_ROOT: Path = CURRENT_DIR.parent

if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

try:
    from ras_csc_core import ras_csc_ode # type: ignore
except ImportError as exc:  # pragma: no cover - defensive
    raise ImportError(
        "Could not import 'ras_csc_ode' from 'ras_csc_core'. "
        "Ensure that ras_csc_core.py is located in src/ and that this "
        "script is in src/scripts/."
    ) from exc

# Import ODE integrator
from scipy.integrate import odeint  # type: ignore


# ----------------------------------------------------------------------
# PATH CONSTANTS
# ----------------------------------------------------------------------

#   determine script, src, and project root directories
CURRENT_DIR: Path = Path(__file__).resolve().parent
SRC_ROOT: Path = CURRENT_DIR.parent
ROOT: Path = SRC_ROOT.parent

#   define where calibration and processed data live
CALIBRATION_DIR: Path = ROOT / "results" / "calibration"
PARAM_JSON: Path = CALIBRATION_DIR / "optimized_parameters_CORRECTED.json"

PROCESSED_DIR: Path = ROOT / "data" / "processed" / "omics_summaries"
RESULTS_MODEL_DIR: Path = ROOT / "results" / "model"
RESULTS_TEMPORAL_DIR: Path = ROOT / "results" / "temporal"

#   define temporal figure directory
FIG_TEMPORAL_DIR: Path = ROOT / "figures" / "temporal"

#   ensure output directories exist
RESULTS_MODEL_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_TEMPORAL_DIR.mkdir(parents=True, exist_ok=True)
FIG_TEMPORAL_DIR.mkdir(parents=True, exist_ok=True)


# ----------------------------------------------------------------------
# PARAMETER LOADING
# ----------------------------------------------------------------------

def load_parameters() -> Dict[str, float]:
    """
    Load optimized parameters from JSON and return as a float-valued dict.

    This centralizes parameter loading so that both temporal and
    sensitivity scripts use exactly the same parameter set. The
    function enforces:
      - JSON file existence
      - Non-empty structure
      - Conversion of every value to float

    Returns
    -------
    Dict[str, float]
        Mapping from parameter name to numeric value.

    Raises
    ------
    FileNotFoundError
        If the calibration file does not exist.
    ValueError
        If the file is empty, malformed, or contains non-numeric values.
    """
    #   check the calibration JSON exists
    if not PARAM_JSON.exists():
        raise FileNotFoundError(
            f"Calibrated parameters not found at {PARAM_JSON}. "
            "Run run_model_calibration.py first."
        )

    #   read JSON safely
    with PARAM_JSON.open("r", encoding="utf-8") as fh:
        raw: Any = json.load(fh)

    #   validate structure
    if not isinstance(raw, dict) or not raw:
        raise ValueError(
            f"Parameter file {PARAM_JSON} is empty or not a JSON object."
        )

    #   convert each value to float defensively
    params: Dict[str, float] = {}
    for key, value in raw.items():
        try:
            params[str(key)] = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Parameter '{key}' with value '{value}' "
                "could not be converted to float."
            ) from exc

    #   report how many parameters were loaded
    print(
        f"[INFO] Loaded {len(params)} calibrated parameters from {PARAM_JSON}")
    return params


# ----------------------------------------------------------------------
# MODEL TIME COURSE SIMULATION
# ----------------------------------------------------------------------

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

    This wraps scipy.integrate.odeint, using the central ras_csc_ode
    function from ras_csc_core. It returns the full trajectory, not
    just the steady state, so that temporal ordering can be studied.

    Parameters
    ----------
    f_ras : float
        Effective Ras input (0–1) for this time-course simulation.
    params : Dict[str, float]
        Calibrated model parameters.
    y0 : List[float]
        Initial state [C, A, T, R, L, M].
    t_max : float, optional
        Final integration time.
    n_steps : int, optional
        Number of time points between 0 and t_max.
    atol : float, optional
        Absolute tolerance for ODE solver.
    rtol : float, optional
        Relative tolerance for ODE solver.

    Returns
    -------
    t : np.ndarray
        Time grid of shape (n_steps,).
    traj : np.ndarray
        State trajectory of shape (n_steps, 6) for [C, A, T, R, L, M].
    """
    #   build time grid
    t: np.ndarray = np.linspace(0.0, float(t_max), int(n_steps))

    #   cast initial condition to numpy array
    y0_arr: np.ndarray = np.asarray(y0, dtype=float)
    if y0_arr.shape != (6,):
        raise ValueError(
            f"simulate_time_course expected y0 of length 6, got shape {y0_arr.shape}"
        )

    #   integrate ODE with defensive try/except
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
        #   return flat trajectory as a safe fallback
        fallback = np.tile(y0_arr, (t.size, 1))
        return t, fallback

    #   clip to [0, 1] and ensure finite values
    if not np.all(np.isfinite(sol)):
        print(
            f"[WARN] Non-finite values in trajectory for f_ras={f_ras:.3f}; "
            "clipping and replacing NaNs."
        )
    sol = np.nan_to_num(sol, nan=0.0, posinf=1.0, neginf=0.0)
    sol = np.clip(sol, 0.0, 1.0)

    return t, sol


def compute_half_rise_time(t: np.ndarray, x: np.ndarray) -> float:
    """
    Compute a half-rise (or half-decay) time for a single trajectory.

    Given a time series x(t), this function defines:
        x0 = x(t=0)
        xT = x(t=t_max)
        x_half = x0 + 0.5 * (xT - x0)

    If xT > x0, it finds the earliest time t where x >= x_half.
    If xT < x0, it finds the earliest time t where x <= x_half.

    If the dynamic range is tiny or the condition is never met,
    returns NaN.

    Parameters
    ----------
    t : np.ndarray
        Time grid.
    x : np.ndarray
        Time series values.

    Returns
    -------
    float
        Half-rise time t_50 (or NaN if undefined).
    """
    #   basic shape checks
    if t.shape != x.shape:
        raise ValueError(
            f"compute_half_rise_time expects t and x with same shape, got {t.shape} vs {x.shape}"
        )

    #   compute initial and final values
    x0: float = float(x[0])
    xT: float = float(x[-1])

    #   ignore if the overall change is negligible
    if abs(xT - x0) < 1e-4:
        return float("nan")

    #   compute target half-height
    x_half: float = x0 + 0.5 * (xT - x0)

    #   choose inequality based on direction of change
    if xT > x0:
        mask = x >= x_half
    else:
        mask = x <= x_half

    #   find earliest time index satisfying condition
    idx = np.where(mask)[0]
    if idx.size == 0:
        return float("nan")

    return float(t[idx[0]])


def run_model_time_course(params: Dict[str, float]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run a Ras "switch-on" experiment and summarize temporal ordering.

    This function:
      - Starts from a benign-like initial condition.
      - Sets a high Ras input (e.g., 0.9).
      - Integrates the ODE system over time.
      - Builds a DataFrame with the full time course.
      - Computes half-rise times for each module.

    Returns
    -------
    time_course_df : pd.DataFrame
        Columns: time, C, A, T, R, L, M, f_ras.
    half_times_df : pd.DataFrame
        Columns: module, t_half.
    """
    #   define Ras input for malignant-like regime
    f_ras_high: float = 0.9

    #   use benign-like initial condition consistent with hysteresis script
    y0_benign: List[float] = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

    #   integrate time course
    t, traj = simulate_time_course(
        f_ras=f_ras_high,
        params=params,
        y0=y0_benign,
        t_max=200.0,
        n_steps=400,
    )

    #   unpack trajectory into named columns
    C = traj[:, 0]
    A = traj[:, 1]
    T = traj[:, 2]
    R = traj[:, 3]
    L = traj[:, 4]
    M = traj[:, 5]

    #   build time-course DataFrame
    time_course_df = pd.DataFrame(
        {
            "time": t.astype(float),
            "C": C.astype(float),
            "A": A.astype(float),
            "T": T.astype(float),
            "R": R.astype(float),
            "L": L.astype(float),
            "M": M.astype(float),
            "f_ras": float(f_ras_high),
        }
    )

    #   compute half-rise times for each module
    module_names: List[str] = ["C", "A", "T", "R", "L", "M"]
    half_times: List[Dict[str, Any]] = []

    for name, series in zip(module_names, [C, A, T, R, L, M]):
        #   compute half-rise time for this series
        t_half = compute_half_rise_time(t, series)
        #   store result as a dictionary
        half_times.append({"module": name, "t_half": float(t_half)})

    #   build half-rise summary DataFrame
    half_times_df = pd.DataFrame(half_times)

    #   write outputs to disk
    tc_path: Path = RESULTS_MODEL_DIR / "model_time_course_CORRECTED.csv"
    ht_path: Path = RESULTS_MODEL_DIR / "model_half_times_CORRECTED.csv"

    #   save time course and half times
    time_course_df.to_csv(tc_path, index=False)
    half_times_df.to_csv(ht_path, index=False)

    #   log save locations
    print(f"[SAVED] Model time course to {tc_path}")
    print(f"[SAVED] Model half-rise times to {ht_path}")

    return time_course_df, half_times_df


# ----------------------------------------------------------------------
# DATA PSEUDOTIME PROFILES
# ----------------------------------------------------------------------

def load_scrna_temporal_frame() -> pd.DataFrame:
    """
    Load per-cell module scores and metadata including pseudotime.

    This merges:
      - scc_scRNA_module_scores_per_cell.csv
      - scc_scRNA_K14pos_metadata_with_CSC_labels.csv

    and returns a single DataFrame keyed by cell_id.

    Returns
    -------
    pd.DataFrame
        Per-cell table containing pseudotime_norm and module scores.

    Raises
    ------
    FileNotFoundError
        If required CSVs are missing.
    """
    #   define paths to scRNA outputs
    scores_path: Path = PROCESSED_DIR / "scc_scRNA_module_scores_per_cell.csv"
    meta_path: Path = PROCESSED_DIR / "scc_scRNA_K14pos_metadata_with_CSC_labels.csv"

    #   check both files exist
    if not scores_path.exists():
        raise FileNotFoundError(
            f"Module scores file not found at {scores_path}. "
            "Run 02_process_scrna.R first."
        )
    if not meta_path.exists():
        raise FileNotFoundError(
            f"Metadata file not found at {meta_path}. "
            "Run 02_process_scrna.R first."
        )

    #   read tables
    scores_df = pd.read_csv(scores_path)
    meta_df = pd.read_csv(meta_path)

    #   basic sanity checks
    if "cell_id" not in scores_df.columns:
        raise ValueError("scores_df is missing 'cell_id' column.")
    if "cell_id" not in meta_df.columns:
        raise ValueError("meta_df is missing 'cell_id' column.")
    if "pseudotime_norm" not in meta_df.columns:
        raise ValueError(
            "meta_df is missing 'pseudotime_norm'. "
            "Update 02_process_scrna.R to export it."
        )

    #   merge on cell_id
    merged = pd.merge(scores_df, meta_df, on="cell_id", how="inner")

    #   drop cells with missing pseudotime
    merged = merged[np.isfinite(merged["pseudotime_norm"].to_numpy())].copy()

    #   log number of usable cells
    print(
        f"[INFO] Loaded {len(merged)} cells with module scores and pseudotime_norm."
    )
    return merged


def compute_pseudotime_profiles(
    df: pd.DataFrame,
    n_bins: int = 20,
) -> pd.DataFrame:
    """
    Bin pseudotime and compute module profiles along the trajectory.

    Parameters
    ----------
    df : pd.DataFrame
        Per-cell table with columns:
          - pseudotime_norm
          - Angio_module_score
          - TGFb_module_score
          - mTOR_module_score
          - CSC_module_score
          - CSC_signature_score
    n_bins : int, optional
        Number of equal-width pseudotime bins between 0 and 1.

    Returns
    -------
    pd.DataFrame
        Long-format table with columns:
          - module
          - pt_bin_start, pt_bin_end, pt_bin_center
          - mean, std, n
    """
    #   ensure required columns are present
    required_cols = [
        "pseudotime_norm",
        "Angio_module_score",
        "TGFb_module_score",
        "mTOR_module_score",
        "CSC_module_score",
        "CSC_signature_score",
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"compute_pseudotime_profiles missing columns: {missing}")

    #   set up bin edges in [0, 1]
    edges = np.linspace(0.0, 1.0, int(n_bins) + 1)
    pt = df["pseudotime_norm"].to_numpy(dtype=float)

    #   assign each cell to a bin index (0..n_bins-1)
    bin_idx = np.digitize(pt, edges, right=False) - 1
    bin_idx = np.clip(bin_idx, 0, n_bins - 1)

    #   attach bin index to DataFrame
    df = df.copy()
    df["pt_bin_index"] = bin_idx

    #   map from bin index to numeric boundaries and centers
    bin_starts = edges[:-1]
    bin_ends = edges[1:]
    bin_centers = 0.5 * (bin_starts + bin_ends)

    #   modules to summarize
    module_cols = {
        "Angio": "Angio_module_score",
        "TGFb": "TGFb_module_score",
        "mTOR": "mTOR_module_score",
        "CSC_module": "CSC_module_score",
        "CSC_signature": "CSC_signature_score",
    }

    #   assemble long-format summary
    rows: List[Dict[str, Any]] = []
    for b in range(n_bins):
        #   subset cells in this pseudotime bin
        subset = df.loc[df["pt_bin_index"] == b]

        #   skip empty bins
        if subset.empty:
            continue

        #   get bin boundaries
        start = float(bin_starts[b])
        end = float(bin_ends[b])
        center = float(bin_centers[b])

        #   summarize each module in this bin
        for module_name, col in module_cols.items():
            #   extract numeric module values
            values = subset[col].to_numpy(dtype=float)
            values = values[np.isfinite(values)]

            #   compute mean, std, n with edge-case handling
            if values.size == 0:
                mean_val = float("nan")
                std_val = float("nan")
                n_val = 0
            else:
                mean_val = float(values.mean())
                std_val = float(values.std(ddof=1)) if values.size > 1 else 0.0
                n_val = int(values.size)

            #   append summary row
            rows.append(
                {
                    "module": module_name,
                    "pt_bin_index": int(b),
                    "pt_bin_start": start,
                    "pt_bin_end": end,
                    "pt_bin_center": center,
                    "mean": mean_val,
                    "std": std_val,
                    "n": n_val,
                }
            )

    #   build summary DataFrame
    profiles_df = pd.DataFrame(rows)

    #   write to disk
    out_path: Path = RESULTS_TEMPORAL_DIR / \
        "scrna_pseudotime_profiles_CORRECTED.csv"
    profiles_df.to_csv(out_path, index=False)
    print(f"[SAVED] scRNA pseudotime profiles to {out_path}")

    return profiles_df


# ----------------------------------------------------------------------
# FIGURE GENERATION
# ----------------------------------------------------------------------

def make_model_time_course_figure(
    time_course_df: pd.DataFrame,
    half_times_df: pd.DataFrame,
) -> None:
    """
    Generate a multi-panel figure for C, A, T, R, L, M vs. time.

    The figure:
      - Uses 3x2 subplots, one per variable.
      - Shows trajectories in [0, 1].
      - Overlays a vertical dashed line at each variable's half-rise
        time when available.

    The result is saved to the temporal figures directory.
    """
    #   ensure expected columns are present
    required_cols = ["time", "C", "A", "T", "R", "L", "M"]
    for col in required_cols:
        if col not in time_course_df.columns:
            raise ValueError(
                f"make_model_time_course_figure missing column '{col}' in time_course_df."
            )

    #   create a mapping from module name to half-rise time
    half_map: Dict[str, float] = {}
    for _, row in half_times_df.iterrows():
        #   safely extract module and t_half
        module_name = str(row.get("module", ""))
        t_half_val = row.get("t_half", np.nan)
        #   store only finite half-times
        if module_name and np.isfinite(t_half_val):
            half_map[module_name] = float(t_half_val)

    #   define subplot layout and module ordering
    module_order: List[str] = ["C", "A", "T", "R", "L", "M"]
    titles: Dict[str, str] = {
        "C": "CSC fraction C(t)",
        "A": "Angiogenesis A(t)",
        "T": "TGFβ T(t)",
        "R": "LEPR R(t)",
        "L": "Leptin L(t)",
        "M": "mTOR M(t)",
    }

    #   extract shared time vector
    t = time_course_df["time"].to_numpy(dtype=float)

    #   create 3x2 panel figure
    fig, axes = plt.subplots(
        nrows=3,
        ncols=2,
        figsize=(8.0, 9.0),
        sharex=True,
        sharey=True,
    )

    #   flatten axes for easy iteration
    axes_flat = axes.ravel()

    #   iterate over modules and axes simultaneously
    for ax, module_name in zip(axes_flat, module_order):
        #   extract variable trajectory
        y = time_course_df[module_name].to_numpy(dtype=float)

        #   plot time course
        ax.plot(t, y, linewidth=2.0)

        #   add half-rise vertical line if available
        if module_name in half_map:
            t_half = half_map[module_name]
            ax.axvline(
                t_half,
                linestyle="--",
                linewidth=1.2,
            )
            #   annotate t_half numerically
            ax.text(
                t_half,
                0.05,
                f"t₅₀={t_half:.1f}",
                rotation=90,
                va="bottom",
                ha="right",
                fontsize=8,
            )

        #   set axis title and limits
        ax.set_title(titles[module_name], fontsize=11)
        ax.set_ylim(-0.02, 1.02)

        #   enable grid for visual guidance
        ax.grid(alpha=0.3, linestyle=":")

    #   label shared axes
    for ax in axes[-1, :]:
        ax.set_xlabel("Model time (a.u.)", fontsize=11)
    for ax in axes[:, 0]:
        ax.set_ylabel("Activity / fraction", fontsize=11)

    #   add overall title and adjust layout
    fig.suptitle("Ras switch-on time course: model modules", fontsize=14)
    fig.tight_layout(rect=[0, 0.02, 1, 0.97])

    #   define output path and save
    out_path: Path = FIG_TEMPORAL_DIR / "temporal_model_time_course_CORRECTED.png"
    fig.savefig(out_path, dpi=400)
    plt.close(fig)

    #   log save location
    print(f"[SAVED] Model temporal figure to {out_path}")


def make_pseudotime_profile_figure(
    profiles_df: pd.DataFrame,
) -> None:
    """
    Generate a multi-panel figure for module scores vs. pseudotime.

    The figure:
      - Uses one row of panels (up to 5 modules).
      - For each module, plots mean score vs. pseudotime bin center.
      - Shades ±1 standard deviation as a band.
      - Annotates approximate trends along pseudotime.

    The result is saved to the temporal figures directory.
    """
    #   ensure required columns are present
    required_cols = [
        "module",
        "pt_bin_center",
        "mean",
        "std",
        "n",
    ]
    for col in required_cols:
        if col not in profiles_df.columns:
            raise ValueError(
                f"make_pseudotime_profile_figure missing column '{col}' in profiles_df."
            )

    #   define the module plotting order
    module_order: List[str] = [
        "Angio",
        "TGFb",
        "mTOR",
        "CSC_module",
        "CSC_signature",
    ]

    #   filter to only modules we actually have
    available_modules = [
        m for m in module_order if m in profiles_df["module"].unique()
    ]
    if not available_modules:
        raise ValueError(
            "No expected modules found in profiles_df['module']; "
            "cannot generate pseudotime figure."
        )

    #   determine number of panels to draw
    n_mod = len(available_modules)

    #   set up figure with one row of panels
    fig, axes = plt.subplots(
        nrows=1,
        ncols=n_mod,
        figsize=(3.2 * n_mod, 3.5),
        sharex=True,
    )

    #   normalize axes variable to a list for n_mod == 1
    if n_mod == 1:
        axes_list = [axes]
    else:
        axes_list = list(axes)

    #   pretty display labels for modules
    label_map: Dict[str, str] = {
        "Angio": "Angiogenesis module",
        "TGFb": "TGFβ module",
        "mTOR": "mTOR module",
        "CSC_module": "CSC module",
        "CSC_signature": "Yuan CSC signature",
    }

    #   iterate through each module and build its panel
    for ax, module_name in zip(axes_list, available_modules):
        #   subset rows for this module
        sub = profiles_df.loc[profiles_df["module"] == module_name].copy()

        #   sort by pseudotime center
        sub = sub.sort_values("pt_bin_center")

        #   extract pseudotime and summary statistics
        x = sub["pt_bin_center"].to_numpy(dtype=float)
        y = sub["mean"].to_numpy(dtype=float)
        y_std = sub["std"].to_numpy(dtype=float)

        #   plot mean line
        ax.plot(
            x,
            y,
            linewidth=2.0,
        )

        #   plot ±1 standard deviation band
        ax.fill_between(
            x,
            y - y_std,
            y + y_std,
            alpha=0.25,
        )

        #   decorate axes
        ax.set_title(label_map.get(module_name, module_name), fontsize=11)
        ax.set_xlim(-0.02, 1.02)
        ax.grid(alpha=0.3, linestyle=":")

        #   set y-label only on first panel
        if ax is axes_list[0]:
            ax.set_ylabel("Module score (a.u.)", fontsize=11)

        #   set x-label on all panels
        ax.set_xlabel("Pseudotime (normalized)", fontsize=10)

    #   adjust layout and add figure title
    fig.suptitle("scRNA modules along pseudotime trajectory", fontsize=14)
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])

    #   define output path and save
    out_path: Path = FIG_TEMPORAL_DIR / \
        "temporal_scrna_pseudotime_profiles_CORRECTED.png"
    fig.savefig(out_path, dpi=400)
    plt.close(fig)

    #   log save location
    print(f"[SAVED] scRNA pseudotime figure to {out_path}")


# ----------------------------------------------------------------------
# MAIN ORCHESTRATION
# ----------------------------------------------------------------------

def main() -> None:
    """
    Orchestrate temporal analysis for model and scRNA data.

    Steps
    -----
    1. Load calibrated parameters.
    2. Run Ras switch-on time course and compute half-rise times.
    3. Generate model time-course figure.
    4. Load scRNA module scores + pseudotime from R pipeline.
    5. Compute binned pseudotime module profiles.
    6. Generate scRNA pseudotime profile figure.
    """
    #   section header
    print("\n" + "=" * 70)
    print("RUNNING TEMPORAL ANALYSIS: MODEL TIME COURSE + scRNA PSEUDOTIME")
    print("=" * 70)

    #   load parameters
    params = load_parameters()

    #   run model time-course analysis
    print("\n[STEP] Model Ras switch-on time course...")
    time_course_df, half_times_df = run_model_time_course(params=params)

    #   generate model time course figure
    print("\n[STEP] Generating model temporal figure...")
    make_model_time_course_figure(time_course_df, half_times_df)

    #   build scRNA pseudotime profiles
    print("\n[STEP] scRNA pseudotime profiles...")
    scrna_df = load_scrna_temporal_frame()
    profiles_df = compute_pseudotime_profiles(scrna_df, n_bins=20)

    #   generate pseudotime profile figure
    print("\n[STEP] Generating scRNA pseudotime figure...")
    make_pseudotime_profile_figure(profiles_df)

    #   final message
    print("\n" + "=" * 70)
    print("TEMPORAL ANALYSIS + FIGURE GENERATION COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
