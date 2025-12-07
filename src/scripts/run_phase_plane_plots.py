#!/usr/bin/env python3
"""
run_phase_plane_plots.py

Advanced dynamical visualizations for the Ras–CSC feedback model.

This script assumes:
    - The model has been calibrated with run_model_calibration.py.
    - Optimized parameters are stored in:
          results/calibration/optimized_parameters_CORRECTED.json
    - The core ODE and steady-state integrator live in:
          src/ras_csc_core.py

What this script does:
    1. Loads the calibrated parameter set.
    2. Picks a Ras input in the "danger zone" (user configurable).
    3. Computes low-branch and high-branch steady states at that Ras.
    4. For each branch, builds a 2D C–M phase plane:
         - Holds A, T, R, L fixed at the branch steady state.
         - Sweeps over a grid of (C, M) pairs.
         - Evaluates the ODE right-hand side to get (dC/dt, dM/dt).
         - Plots a vector field and marks the steady state.

Why this is separate:
    - Calibration and hypothesis tests establish that the model makes
      sense numerically.
    - This script gives a geometric view of the feedback loop that is
      easier to explain in figures and text.

Outputs:
    - figures/phase_plane/phase_plane_C_M_CORRECTED.png
"""

from __future__ import annotations

# Import typing utilities for clearer type hints
from typing import Dict, Any, Tuple

# Import standard library modules
# Import json to read the calibrated parameter file
import json
# Import sys to modify the Python path so src/ can be imported cleanly
import sys
# Import Path for filesystem-safe path handling
from pathlib import Path

# Import numerical and plotting libraries
# Import numpy for array operations
import numpy as np
# Import matplotlib for plotting
import matplotlib.pyplot as plt
# Import seaborn to set a publication-friendly style
import seaborn as sns


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
    from ras_csc_core import ras_csc_ode, simulate_steady_state  # type: ignore
except ImportError as exc:  # pragma: no cover - defensive
    #   Raise a descriptive error if the core module cannot be imported
    raise ImportError(
        "Could not import 'ras_csc_ode' and 'simulate_steady_state' "
        "from 'ras_csc_core'. Ensure ras_csc_core.py is in src/ and "
        "this script lives in src/scripts/."
    ) from exc


# ======================================================================
# PATH CONSTANTS AND PLOTTING STYLE
# ======================================================================

#   Path to the calibrated parameter JSON from the calibration step
PARAM_JSON: Path = ROOT / "results" / "calibration" / "optimized_parameters_CORRECTED.json"

#   Directory for advanced phase-plane figures
PHASE_PLANE_FIG_DIR: Path = ROOT / "figures" / "phase_plane"

#   Make sure the figure directory exists
PHASE_PLANE_FIG_DIR.mkdir(parents=True, exist_ok=True)

#   Set a clean whitegrid style for publication-quality figures
sns.set_style("whitegrid")
#   Increase DPI for sharper saved images
plt.rcParams["figure.dpi"] = 300
#   Set base font size for clarity without changing font family
plt.rcParams["font.size"] = 11


# ======================================================================
# PARAMETER LOADING
# ======================================================================

def load_parameters() -> Dict[str, float]:
    """
    Load calibrated parameters from JSON and return as a float-valued dict.

    This helper keeps parameter loading consistent with other scripts.
    It enforces:
        - File existence.
        - Non-empty JSON content.
        - Numeric conversion of values.

    Returns:
        Dict[str, float]:
            Mapping from parameter names to float values.

    Raises:
        FileNotFoundError:
            If the JSON parameter file does not exist.
        ValueError:
            If the JSON content is empty, malformed, or non-numeric.
    """
    #   Check for presence of the calibrated parameter file
    if not PARAM_JSON.exists():
        #   Raise an explicit error if the file is missing
        raise FileNotFoundError(
            f"Calibrated parameter file not found at: {PARAM_JSON}. "
            "Run run_model_calibration.py before generating phase-plane plots."
        )

    #   Open the JSON file and load its contents
    with PARAM_JSON.open("r", encoding="utf-8") as f:
        raw: Any = json.load(f)

    #   Validate that the JSON content is a non-empty dictionary
    if not isinstance(raw, dict) or not raw:
        #   Raise an error if the content is unusable
        raise ValueError(
            f"Parameter file {PARAM_JSON} is empty or malformed."
        )

    #   Create a clean dictionary with float-cast values
    params: Dict[str, float] = {}
    #   Iterate over all key-value pairs from the JSON
    for key, value in raw.items():
        try:
            #   Convert each value to float defensively
            params[str(key)] = float(value)
        except (TypeError, ValueError) as exc:
            #   If conversion fails, raise a detailed error
            raise ValueError(
                f"Parameter '{key}' with value '{value}' could not be "
                "converted to float."
            ) from exc

    #   Print a brief summary for sanity
    print(f"[INFO] Loaded {len(params)} calibrated parameters from {PARAM_JSON}")

    #   Return the cleaned parameter dictionary
    return params


# ======================================================================
# VECTOR FIELD HELPER FOR C–M PHASE PLANE
# ======================================================================

def compute_vector_field_C_M(
    params: Dict[str, float],
    f_ras: float,
    A_fixed: float,
    T_fixed: float,
    R_fixed: float,
    L_fixed: float,
    C_range: Tuple[float, float] = (0.0, 1.0),
    M_range: Tuple[float, float] = (0.0, 1.0),
    n_grid: int = 25,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute the projected vector field (dC/dt, dM/dt) on a C–M grid.

    Strategy:
        1. Set up a regular grid in the C–M plane.
        2. For each grid point, construct a full state:
             y = [C, A_fixed, T_fixed, R_fixed, L_fixed, M].
        3. Evaluate the ODE RHS:
             dy/dt = ras_csc_ode(y, t, f_ras, params)
        4. Extract dC/dt and dM/dt at each grid point.

    This gives a 2D slice through the full 6D system, which helps
    visualize attraction toward the low- and high-CSC states.

    Args:
        params:
            Calibrated parameter dictionary.
        f_ras:
            Ras input level at which the slice is taken.
        A_fixed:
            Fixed angiogenesis level used in the slice.
        T_fixed:
            Fixed TGFβ level used in the slice.
        R_fixed:
            Fixed LEPR level used in the slice.
        L_fixed:
            Fixed leptin level used in the slice.
        C_range:
            Tuple (C_min, C_max) defining the C axis range.
        M_range:
            Tuple (M_min, M_max) defining the M axis range.
        n_grid:
            Number of grid points along each axis.

    Returns:
        Tuple of:
            - CC: 2D array of C coordinates.
            - MM: 2D array of M coordinates.
            - dC: 2D array of dC/dt values.
            - dM: 2D array of dM/dt values.
    """
    #   Create a linearly spaced grid for C
    C_vals: np.ndarray = np.linspace(C_range[0], C_range[1], int(n_grid))
    #   Create a linearly spaced grid for M
    M_vals: np.ndarray = np.linspace(M_range[0], M_range[1], int(n_grid))

    #   Build a meshgrid so each (C,M) combination is represented
    CC, MM = np.meshgrid(C_vals, M_vals)

    #   Allocate arrays for the time derivatives
    dC: np.ndarray = np.zeros_like(CC, dtype=float)
    dM: np.ndarray = np.zeros_like(MM, dtype=float)

    #   Loop over the grid indices for C and M
    for i in range(CC.shape[0]):
        #   Inner loop over columns for each row
        for j in range(CC.shape[1]):
            #   Build the full 6D state vector at this grid point
            y: np.ndarray = np.array(
                [
                    float(CC[i, j]),   # C
                    float(A_fixed),    # A
                    float(T_fixed),    # T
                    float(R_fixed),    # R
                    float(L_fixed),    # L
                    float(MM[i, j]),   # M
                ],
                dtype=float,
            )

            #   Evaluate the ODE right-hand side at t=0
            dy_dt = ras_csc_ode(
                y=y,
                t=0.0,
                f_ras=float(f_ras),
                params=params,
            )

            #   Extract dC/dt and dM/dt from the derivative vector
            dC[i, j] = float(dy_dt[0])
            dM[i, j] = float(dy_dt[5])

    #   Return the grid and derivative arrays
    return CC, MM, dC, dM


# ======================================================================
# PHASE-PLANE PLOT FOR C–M
# ======================================================================

def plot_C_M_phase_plane(
    params: Dict[str, float],
    f_ras: float = 0.6,
    C_range: Tuple[float, float] = (0.0, 1.0),
    M_range: Tuple[float, float] = (0.0, 1.0),
    n_grid: int = 40,
) -> None:
    """
    Plot a two-panel C–M phase plane at a given Ras input with:
        - colored streamlines (direction field),
        - speed background,
        - dC/dt = 0 and dM/dt = 0 nullclines,
        - low- and high-branch steady states marked.

    Panel A: slice near the low-CSC branch.
    Panel B: slice near the high-CSC branch.
    """

    # Compute low-branch steady state from a benign-like initial condition
    y_low: np.ndarray = simulate_steady_state(
        f_ras=float(f_ras),
        params=params,
        y0=[0.05, 0.05, 0.05, 0.05, 0.05, 0.05],
        t_max=200.0,
        n_steps=2000,
    )

    # Compute high-branch steady state from a malignant-like initial condition
    y_high: np.ndarray = simulate_steady_state(
        f_ras=float(f_ras),
        params=params,
        y0=[0.9, 0.9, 0.6, 0.6, 0.6, 0.9],
        t_max=200.0,
        n_steps=2000,
    )

    # Unpack steady states for readability
    C_low, A_low, T_low, R_low, L_low, M_low = [float(v) for v in y_low]
    C_high, A_high, T_high, R_high, L_high, M_high = [float(v) for v in y_high]

    # Build a figure with two horizontally aligned panels
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(7.0, 3.0),
        sharex=True,
        sharey=True,
    )

    # Helper to style each panel (so we don't duplicate code)
    def _style_axis(ax: plt.Axes, panel_label: str, title_text: str) -> None:
        # Add bold panel label in the top-left corner
        ax.text(
            0.02,
            0.98,
            panel_label,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=12,
            fontweight="bold",
        )
        # Set axis labels only once; y-label on left panel only
        ax.set_xlabel("CSC fraction (C)")
        ax.set_xlim(C_range[0], C_range[1])
        ax.set_ylim(M_range[0], M_range[1])
        ax.set_title(title_text, fontsize=11)
        ax.grid(alpha=0.3)
        ax.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
        ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])

    # ---------- Panel A: low-CSC branch slice ----------
    ax_left = axes[0]

    # Compute vector field for low-branch slice
    CC_low, MM_low, dC_low, dM_low = compute_vector_field_C_M(
        params=params,
        f_ras=float(f_ras),
        A_fixed=A_low,
        T_fixed=T_low,
        R_fixed=R_low,
        L_fixed=L_low,
        C_range=C_range,
        M_range=M_range,
        n_grid=int(n_grid),
    )

    # Compute speed (magnitude of vector) for background and coloring
    speed_low: np.ndarray = np.sqrt(dC_low ** 2 + dM_low ** 2)

    # Normalize vectors for direction field
    mag_low = speed_low.copy()
    mag_low[mag_low == 0.0] = 1.0
    dC_low_norm = dC_low / mag_low
    dM_low_norm = dM_low / mag_low

    # Plot speed background as a soft colormap
    im_low = ax_left.pcolormesh(
        CC_low,
        MM_low,
        speed_low,
        shading="gouraud",
        cmap="viridis",
        alpha=0.6,
    )

    # Overlay streamlines colored by speed
    ax_left.streamplot(
        CC_low,
        MM_low,
        dC_low_norm,
        dM_low_norm,
        color=speed_low,
        cmap="viridis",
        density=1.5,
        linewidth=0.8,
        arrowsize=0.9,
    )

    # Add nullclines: dC/dt = 0 (solid), dM/dt = 0 (dashed)
    ax_left.contour(
        CC_low,
        MM_low,
        dC_low,
        levels=[0.0],
        colors="black",
        linewidths=1.1,
    )
    ax_left.contour(
        CC_low,
        MM_low,
        dM_low,
        levels=[0.0],
        colors="black",
        linewidths=1.1,
        linestyles="--",
    )

    # Mark low-branch steady state
    ax_left.scatter(
        [C_low],
        [M_low],
        s=35,
        color="black",
        edgecolor="white",
        linewidth=0.7,
        zorder=5,
    )

    # Style left panel
    ax_left.set_ylabel("mTOR activity (M)")
    _style_axis(ax_left, panel_label="A",
                title_text="Slice near low-CSC branch")

    # ---------- Panel B: high-CSC branch slice ----------
    ax_right = axes[1]

    # Compute vector field for high-branch slice
    CC_high, MM_high, dC_high, dM_high = compute_vector_field_C_M(
        params=params,
        f_ras=float(f_ras),
        A_fixed=A_high,
        T_fixed=T_high,
        R_fixed=R_high,
        L_fixed=L_high,
        C_range=C_range,
        M_range=M_range,
        n_grid=int(n_grid),
    )

    # Compute speed for high-branch slice
    speed_high: np.ndarray = np.sqrt(dC_high ** 2 + dM_high ** 2)

    # Normalize vectors
    mag_high = speed_high.copy()
    mag_high[mag_high == 0.0] = 1.0
    dC_high_norm = dC_high / mag_high
    dM_high_norm = dM_high / mag_high

    # Plot background speed field
    ax_right.pcolormesh(
        CC_high,
        MM_high,
        speed_high,
        shading="gouraud",
        cmap="viridis",
        alpha=0.6,
    )

    # Overlay streamlines
    ax_right.streamplot(
        CC_high,
        MM_high,
        dC_high_norm,
        dM_high_norm,
        color=speed_high,
        cmap="viridis",
        density=1.5,
        linewidth=0.8,
        arrowsize=0.9,
    )

    # Add nullclines: dC/dt = 0 (solid), dM/dt = 0 (dashed)
    ax_right.contour(
        CC_high,
        MM_high,
        dC_high,
        levels=[0.0],
        colors="black",
        linewidths=1.1,
    )
    ax_right.contour(
        CC_high,
        MM_high,
        dM_high,
        levels=[0.0],
        colors="black",
        linewidths=1.1,
        linestyles="--",
    )

    # Mark high-branch steady state
    ax_right.scatter(
        [C_high],
        [M_high],
        s=35,
        color="black",
        edgecolor="white",
        linewidth=0.7,
        zorder=5,
    )

    # Style right panel
    _style_axis(ax_right, panel_label="B",
                title_text="Slice near high-CSC branch")

    # ---------- Figure-level tweaks ----------
    # Add a modest suptitle with Ras value
    fig.suptitle(f"C–M phase plane at Ras = {f_ras:.2f}", fontsize=12)

    # Adjust layout to reduce whitespace and keep panels tight
    fig.subplots_adjust(
        left=0.09,
        right=0.98,
        bottom=0.18,
        top=0.80,
        wspace=0.25,
    )

    # Add a single colorbar for speed on the right side
    cbar = fig.colorbar(
        im_low,
        ax=axes.ravel().tolist(),
        fraction=0.035,
        pad=0.02,
    )
    cbar.set_label("Speed √[(dC/dt)² + (dM/dt)²]", fontsize=9)

    # Save the figure
    out_png: Path = PHASE_PLANE_FIG_DIR / "phase_plane_C_M_CORRECTED.png"
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)

    print(f"[SAVED] C–M phase-plane figure to {out_png}")



# ======================================================================
# MAIN ORCHESTRATION
# ======================================================================

def main() -> None:
    """
    Entry point for advanced Ras–CSC phase-plane visualizations.

    This function:
        - Loads calibrated parameters.
        - Generates a two-panel C–M phase plane at a user-specified
          Ras input (default 0.6).
        - Writes the resulting figure to figures/phase_plane/.

    In the current course timeline, this sits after:
        - Calibration (run_model_calibration.py)
        - Hypothesis tests (run_hypothesis_tests.py)
        - Sensitivity analyses (run_sensitivity_analyses.py)
    """
    #   Print a header so logs are easy to read
    print("\n" + "=" * 70)
    print("RUNNING RAS–CSC PHASE-PLANE PLOTS")
    print("=" * 70)

    #   Load the calibrated parameters from JSON
    params: Dict[str, float] = load_parameters()

    #   Choose a Ras value in the bistable region for visualization
    #   You can change f_ras here if you want a different slice
    f_ras_visual: float = 0.6

    #   Generate and save the C–M phase-plane figure
    plot_C_M_phase_plane(
        params=params,
        f_ras=f_ras_visual,
        C_range=(0.0, 1.0),
        M_range=(0.0, 1.0),
        n_grid=25,
    )

    #   Print completion footer
    print("\n" + "=" * 70)
    print("PHASE-PLANE PLOTS COMPLETE")
    print("=" * 70)


#   Execute main when the script is run directly
if __name__ == "__main__":
    main()
