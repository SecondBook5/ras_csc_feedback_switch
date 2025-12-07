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
    4. Builds a 2D C–M phase plane for each branch:
         - Holds A, T, R, L fixed at the branch steady state.
         - Sweeps over a grid of (C, M) pairs.
         - Evaluates the ODE right-hand side to get (dC/dt, dM/dt).
         - Plots a vector field, nullclines, and marks the steady state.
         - Overlays sample trajectories projected into the C–M plane.
    5. Computes Jacobians at the low- and high-branch steady states
       via finite differences and visualizes their eigenvalue spectra.

Outputs:
    - figures/phase_plane/phase_plane_C_M_CORRECTED.png
    - figures/phase_plane/phase_plane_C_M_jacobian_spectra_CORRECTED.png
"""

from __future__ import annotations

from typing import Dict, Any, Tuple, List

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ----------------------------------------------------------------------
# Ensure src/ is on the Python path and import core model utilities
# ----------------------------------------------------------------------

CURRENT_DIR: Path = Path(__file__).resolve().parent
SRC_ROOT: Path = CURRENT_DIR.parent
ROOT: Path = SRC_ROOT.parent

if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

try:
    from ras_csc_core import ras_csc_ode, simulate_steady_state  # type: ignore
except ImportError as exc:  # pragma: no cover - defensive
    raise ImportError(
        "Could not import 'ras_csc_ode' and 'simulate_steady_state' "
        "from 'ras_csc_core'. Ensure ras_csc_core.py is in src/ and "
        "this script lives in src/scripts/."
    ) from exc


# ======================================================================
# PATH CONSTANTS AND PLOTTING STYLE
# ======================================================================

PARAM_JSON: Path = ROOT / "results" / "calibration" / "optimized_parameters_CORRECTED.json"
PHASE_PLANE_FIG_DIR: Path = ROOT / "figures" / "phase_plane"
PHASE_PLANE_FIG_DIR.mkdir(parents=True, exist_ok=True)

sns.set_style("whitegrid")
plt.rcParams["figure.dpi"] = 300
plt.rcParams["font.size"] = 11

COLOR_STABLE = "#009E73"   # for locally stable steady states
COLOR_UNSTABLE = "#D55E00" # for unstable or marginal states
COLOR_TRAJ_1 = "#0072B2"
COLOR_TRAJ_2 = "#CC79A7"
COLOR_TRAJ_3 = "#F0E442"


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

    Returns
    -------
    Dict[str, float]
        Mapping from parameter names to float values.

    Raises
    ------
    FileNotFoundError
        If the JSON parameter file does not exist.
    ValueError
        If the JSON content is empty, malformed, or non-numeric.
    """
    if not PARAM_JSON.exists():
        raise FileNotFoundError(
            f"Calibrated parameter file not found at: {PARAM_JSON}. "
            "Run run_model_calibration.py before generating phase-plane plots."
        )

    with PARAM_JSON.open("r", encoding="utf-8") as f:
        raw: Any = json.load(f)

    if not isinstance(raw, dict) or not raw:
        raise ValueError(
            f"Parameter file {PARAM_JSON} is empty or malformed."
        )

    params: Dict[str, float] = {}
    for key, value in raw.items():
        try:
            params[str(key)] = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Parameter '{key}' with value '{value}' could not be "
                "converted to float."
            ) from exc

    print(f"[INFO] Loaded {len(params)} calibrated parameters from {PARAM_JSON}")
    return params


# ======================================================================
# NUMERIC JACOBIAN (FINITE DIFFERENCES)
# ======================================================================

def compute_numeric_jacobian(
    y_ss: np.ndarray,
    params: Dict[str, float],
    f_ras: float,
    eps: float = 1e-4,
) -> np.ndarray:
    """
    Compute the 6×6 Jacobian matrix at a steady state using finite differences.

    Strategy
    --------
    For each state variable y_k:
        1. Perturb y_k by ±eps (central difference).
        2. Evaluate ras_csc_ode at the perturbed states.
        3. Approximate column k of the Jacobian as:
               (f(y + e_k*eps) - f(y - e_k*eps)) / (2*eps)

    This avoids needing an analytic Jacobian in ras_csc_core.

    Parameters
    ----------
    y_ss : np.ndarray
        Steady-state vector of length 6.
    params : Dict[str, float]
        Calibrated parameter dictionary.
    f_ras : float
        Ras input level at which the Jacobian is evaluated.
    eps : float, optional
        Perturbation step size used for finite differences.

    Returns
    -------
    np.ndarray
        A 6×6 Jacobian matrix evaluated at (y_ss, f_ras, params).
    """
    y_base: np.ndarray = np.asarray(y_ss, dtype=float).copy()
    if y_base.shape != (6,):
        raise ValueError(
            f"compute_numeric_jacobian expects a 6D state vector; got shape {y_base.shape}."
        )

    f_base: np.ndarray = ras_csc_ode(
        y=y_base,
        t=0.0,
        f_ras=float(f_ras),
        params=params,
    )
    f_base = np.asarray(f_base, dtype=float)

    J: np.ndarray = np.zeros((6, 6), dtype=float)

    for k in range(6):
        y_plus: np.ndarray = y_base.copy()
        y_minus: np.ndarray = y_base.copy()

        y_plus[k] += eps
        if y_minus[k] - eps >= 0.0:
            y_minus[k] -= eps
        else:
            # If subtracting eps would go negative, use forward difference
            y_minus[k] = y_base[k]
            y_plus[k] = y_base[k] + eps

        f_plus: np.ndarray = ras_csc_ode(
            y=y_plus,
            t=0.0,
            f_ras=float(f_ras),
            params=params,
        )
        f_minus: np.ndarray = ras_csc_ode(
            y=y_minus,
            t=0.0,
            f_ras=float(f_ras),
            params=params,
        )

        f_plus = np.asarray(f_plus, dtype=float)
        f_minus = np.asarray(f_minus, dtype=float)

        denom: float = float(2.0 * eps) if y_minus[k] != y_plus[k] else float(eps)
        if denom == 0.0:
            raise ZeroDivisionError("Jacobian finite-difference denominator is zero.")

        J[:, k] = (f_plus - f_minus) / denom

    return J


def summarize_jacobian_stability(
    J: np.ndarray,
    label: str,
) -> Tuple[np.ndarray, float]:
    """
    Compute eigenvalues of the Jacobian and report the largest real part.

    Parameters
    ----------
    J : np.ndarray
        6×6 Jacobian matrix.
    label : str
        Label used for printed summaries (e.g., 'low branch').

    Returns
    -------
    Tuple[np.ndarray, float]
        - eigenvalues (complex array of shape (6,))
        - max_real: largest real part among the eigenvalues
    """
    eigvals: np.ndarray = np.linalg.eigvals(J)
    max_real: float = float(np.max(np.real(eigvals)))

    print(f"\n[INFO] Jacobian stability summary for {label}:")
    print(f"  Max Re(λ) = {max_real:.4f}")
    print("  Eigenvalues:")
    for val in eigvals:
        print(f"    λ = {val.real: .4f} + {val.imag: .4f}i")

    return eigvals, max_real


def plot_jacobian_spectra(
    eig_low: np.ndarray,
    eig_high: np.ndarray,
    f_ras: float,
) -> None:
    """
    Plot eigenvalue spectra for low- and high-branch Jacobians in the complex plane.

    Parameters
    ----------
    eig_low : np.ndarray
        Eigenvalues at the low-branch steady state.
    eig_high : np.ndarray
        Eigenvalues at the high-branch steady state.
    f_ras : float
        Ras input level for annotation.
    """
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(7.0, 3.0),
        sharex=True,
        sharey=True,
    )

    ax_low, ax_high = axes

    # Low branch spectrum
    ax_low.axvline(0.0, color="grey", linestyle="--", linewidth=0.8)
    ax_low.scatter(
        np.real(eig_low),
        np.imag(eig_low),
        s=35,
        color=COLOR_STABLE,
        edgecolor="black",
        linewidth=0.6,
    )
    ax_low.set_title("Low-branch Jacobian spectrum")
    ax_low.set_xlabel("Re(λ)")
    ax_low.set_ylabel("Im(λ)")
    ax_low.grid(alpha=0.3)

    # High branch spectrum
    ax_high.axvline(0.0, color="grey", linestyle="--", linewidth=0.8)
    ax_high.scatter(
        np.real(eig_high),
        np.imag(eig_high),
        s=35,
        color=COLOR_UNSTABLE,
        edgecolor="black",
        linewidth=0.6,
    )
    ax_high.set_title("High-branch Jacobian spectrum")
    ax_high.set_xlabel("Re(λ)")
    ax_high.grid(alpha=0.3)

    fig.suptitle(f"Jacobian eigenvalue spectra at Ras = {f_ras:.2f}", fontsize=12)
    fig.subplots_adjust(
        left=0.09,
        right=0.98,
        bottom=0.18,
        top=0.82,
        wspace=0.25,
    )

    out_png: Path = PHASE_PLANE_FIG_DIR / "phase_plane_C_M_jacobian_spectra_CORRECTED.png"
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)

    print(f"[SAVED] Jacobian spectra figure to {out_png}")


# ======================================================================
# TRAJECTORY INTEGRATION (FOR OVERLAYS)
# ======================================================================

def integrate_trajectory_euler(
    y0: List[float],
    params: Dict[str, float],
    f_ras: float,
    t_max: float = 200.0,
    dt: float = 0.5,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Integrate a trajectory in the full 6D system using explicit Euler,
    then return the projected C(t) and M(t) coordinates.

    This is intended purely for visualization: trajectories are used to
    illustrate how trajectories move in the C–M plane for a fixed Ras.

    Parameters
    ----------
    y0 : List[float]
        Initial 6D state vector.
    params : Dict[str, float]
        Calibrated parameter dictionary.
    f_ras : float
        Ras input level.
    t_max : float, optional
        Final time for integration.
    dt : float, optional
        Time step for explicit Euler.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        - C_traj: C(t) samples
        - M_traj: M(t) samples
    """
    y: np.ndarray = np.asarray(y0, dtype=float).copy()
    if y.shape != (6,):
        raise ValueError(
            f"integrate_trajectory_euler expects a 6D state; got shape {y.shape}."
        )

    n_steps: int = int(np.ceil(t_max / dt))
    C_traj: List[float] = []
    M_traj: List[float] = []

    for _ in range(n_steps):
        C_traj.append(float(y[0]))
        M_traj.append(float(y[5]))

        dy_dt: np.ndarray = ras_csc_ode(
            y=y,
            t=0.0,
            f_ras=float(f_ras),
            params=params,
        )
        dy_dt = np.asarray(dy_dt, dtype=float)

        y = y + dt * dy_dt
        y = np.maximum(y, 0.0)

    return np.asarray(C_traj, dtype=float), np.asarray(M_traj, dtype=float)


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

    Strategy
    --------
    1. Set up a regular grid in the C–M plane.
    2. For each grid point, construct a full state:
         y = [C, A_fixed, T_fixed, R_fixed, L_fixed, M].
    3. Evaluate the ODE RHS:
         dy/dt = ras_csc_ode(y, t, f_ras, params)
    4. Extract dC/dt and dM/dt at each grid point.

    This gives a 2D slice through the full 6D system, which helps
    visualize attraction toward the low- and high-CSC states.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        - CC: 2D array of C coordinates.
        - MM: 2D array of M coordinates.
        - dC: 2D array of dC/dt values.
        - dM: 2D array of dM/dt values.
    """
    C_vals: np.ndarray = np.linspace(C_range[0], C_range[1], int(n_grid))
    M_vals: np.ndarray = np.linspace(M_range[0], M_range[1], int(n_grid))

    CC, MM = np.meshgrid(C_vals, M_vals)

    dC: np.ndarray = np.zeros_like(CC, dtype=float)
    dM: np.ndarray = np.zeros_like(MM, dtype=float)

    for i in range(CC.shape[0]):
        for j in range(CC.shape[1]):
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

            dy_dt: np.ndarray = ras_csc_ode(
                y=y,
                t=0.0,
                f_ras=float(f_ras),
                params=params,
            )

            dC[i, j] = float(dy_dt[0])
            dM[i, j] = float(dy_dt[5])

    return CC, MM, dC, dM


# ======================================================================
# PHASE-PLANE PLOT FOR C–M (WITH TRAJECTORIES AND JACOBIANS)
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
        - low- and high-branch steady states marked and colored by
          local stability (from Jacobian eigenvalues),
        - sample trajectories overlaid in C–M space.

    Panel A: slice near the low-CSC branch.
    Panel B: slice near the high-CSC branch.

    Also computes Jacobians at the two steady states and writes a
    separate Jacobian eigenvalue spectrum figure.
    """
    # --- Steady states at this Ras value ---
    y_low: np.ndarray = simulate_steady_state(
        f_ras=float(f_ras),
        params=params,
        y0=[0.05, 0.05, 0.05, 0.05, 0.05, 0.05],
        t_max=200.0,
        n_steps=2000,
    )

    y_high: np.ndarray = simulate_steady_state(
        f_ras=float(f_ras),
        params=params,
        y0=[0.9, 0.9, 0.6, 0.6, 0.6, 0.9],
        t_max=200.0,
        n_steps=2000,
    )

    C_low, A_low, T_low, R_low, L_low, M_low = [float(v) for v in y_low]
    C_high, A_high, T_high, R_high, L_high, M_high = [float(v) for v in y_high]

    print(f"\n[INFO] Low-branch steady state at Ras={f_ras:.2f}: "
          f"C={C_low:.3f}, M={M_low:.3f}")
    print(f"[INFO] High-branch steady state at Ras={f_ras:.2f}: "
          f"C={C_high:.3f}, M={M_high:.3f}")

    # --- Jacobians and stability classification ---
    J_low: np.ndarray = compute_numeric_jacobian(
        y_ss=y_low,
        params=params,
        f_ras=float(f_ras),
    )
    eig_low, max_real_low = summarize_jacobian_stability(
        J=J_low,
        label="low branch",
    )

    J_high: np.ndarray = compute_numeric_jacobian(
        y_ss=y_high,
        params=params,
        f_ras=float(f_ras),
    )
    eig_high, max_real_high = summarize_jacobian_stability(
        J=J_high,
        label="high branch",
    )

    # Marker colors based on local stability
    color_low = COLOR_STABLE if max_real_low < 0.0 else COLOR_UNSTABLE
    color_high = COLOR_STABLE if max_real_high < 0.0 else COLOR_UNSTABLE

    # --- Build trajectories for overlay (same for both panels) ---
    traj_inits: List[List[float]] = [
        [0.05, 0.05, 0.05, 0.05, 0.05, 0.05],
        [0.4, 0.4, 0.3, 0.3, 0.3, 0.4],
        [0.9, 0.9, 0.6, 0.6, 0.6, 0.9],
    ]
    traj_colors: List[str] = [COLOR_TRAJ_1, COLOR_TRAJ_2, COLOR_TRAJ_3]

    traj_CM: List[Tuple[np.ndarray, np.ndarray]] = []
    for y0 in traj_inits:
        C_traj, M_traj = integrate_trajectory_euler(
            y0=y0,
            params=params,
            f_ras=float(f_ras),
            t_max=200.0,
            dt=0.5,
        )
        traj_CM.append((C_traj, M_traj))

    # --- Figure with two C–M slices ---
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(7.0, 3.0),
        sharex=True,
        sharey=True,
    )

    def _style_axis(ax: plt.Axes, panel_label: str, title_text: str) -> None:
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
        ax.set_xlabel("CSC fraction (C)")
        ax.set_xlim(C_range[0], C_range[1])
        ax.set_ylim(M_range[0], M_range[1])
        ax.set_title(title_text, fontsize=11)
        ax.grid(alpha=0.3)
        ax.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
        ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])

    # ---------- Panel A: low-CSC branch slice ----------
    ax_left = axes[0]

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

    speed_low: np.ndarray = np.sqrt(dC_low ** 2 + dM_low ** 2)
    mag_low: np.ndarray = speed_low.copy()
    mag_low[mag_low == 0.0] = 1.0
    dC_low_norm: np.ndarray = dC_low / mag_low
    dM_low_norm: np.ndarray = dM_low / mag_low

    im_low = ax_left.pcolormesh(
        CC_low,
        MM_low,
        speed_low,
        shading="gouraud",
        cmap="viridis",
        alpha=0.6,
    )

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

    # Overlay trajectories (projected into this slice)
    for (C_traj, M_traj), col in zip(traj_CM, traj_colors):
        ax_left.plot(
            C_traj,
            M_traj,
            color=col,
            linewidth=1.1,
            alpha=0.9,
        )

    ax_left.scatter(
        [C_low],
        [M_low],
        s=40,
        color=color_low,
        edgecolor="white",
        linewidth=0.8,
        zorder=5,
        label="Low-branch SS",
    )

    ax_left.set_ylabel("mTOR activity (M)")
    _style_axis(
        ax_left,
        panel_label="A",
        title_text="Slice near low-CSC branch",
    )

    # ---------- Panel B: high-CSC branch slice ----------
    ax_right = axes[1]

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

    speed_high: np.ndarray = np.sqrt(dC_high ** 2 + dM_high ** 2)
    mag_high: np.ndarray = speed_high.copy()
    mag_high[mag_high == 0.0] = 1.0
    dC_high_norm: np.ndarray = dC_high / mag_high
    dM_high_norm: np.ndarray = dM_high / mag_high

    ax_right.pcolormesh(
        CC_high,
        MM_high,
        speed_high,
        shading="gouraud",
        cmap="viridis",
        alpha=0.6,
    )

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

    for (C_traj, M_traj), col in zip(traj_CM, traj_colors):
        ax_right.plot(
            C_traj,
            M_traj,
            color=col,
            linewidth=1.1,
            alpha=0.9,
        )

    ax_right.scatter(
        [C_high],
        [M_high],
        s=40,
        color=color_high,
        edgecolor="white",
        linewidth=0.8,
        zorder=5,
        label="High-branch SS",
    )

    _style_axis(
        ax_right,
        panel_label="B",
        title_text="Slice near high-CSC branch",
    )

    # --- Figure-level tweaks ---
    fig.suptitle(f"C–M phase plane with trajectories at Ras = {f_ras:.2f}", fontsize=12)
    fig.subplots_adjust(
        left=0.09,
        right=0.98,
        bottom=0.18,
        top=0.80,
        wspace=0.25,
    )

    cbar = fig.colorbar(
        im_low,
        ax=axes.ravel().tolist(),
        fraction=0.035,
        pad=0.02,
    )
    cbar.set_label("Speed √[(dC/dt)² + (dM/dt)²]", fontsize=9)

    out_png: Path = PHASE_PLANE_FIG_DIR / "phase_plane_C_M_CORRECTED.png"
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)

    print(f"[SAVED] C–M phase-plane figure to {out_png}")

    # --- Separate Jacobian spectra figure ---
    plot_jacobian_spectra(
        eig_low=eig_low,
        eig_high=eig_high,
        f_ras=float(f_ras),
    )


# ======================================================================
# MAIN ORCHESTRATION
# ======================================================================

def main() -> None:
    """
    Entry point for advanced Ras–CSC phase-plane visualizations.

    This function:
        - Loads calibrated parameters.
        - Generates a two-panel C–M phase plane at a user-specified
          Ras input (default 0.6), including trajectories and local
          stability markers.
        - Generates a separate Jacobian eigenvalue spectrum figure.
    """
    print("\n" + "=" * 70)
    print("RUNNING RAS–CSC PHASE-PLANE PLOTS")
    print("=" * 70)

    params: Dict[str, float] = load_parameters()

    f_ras_visual: float = 0.6

    plot_C_M_phase_plane(
        params=params,
        f_ras=f_ras_visual,
        C_range=(0.0, 1.0),
        M_range=(0.0, 1.0),
        n_grid=25,
    )

    print("\n" + "=" * 70)
    print("PHASE-PLANE PLOTS COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
