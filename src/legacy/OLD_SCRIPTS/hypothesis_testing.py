#!/usr/bin/env python3
"""
hypothesis_testing.py

Test of the Alternative Hypothesis (H1):

  "Disruption of the leptin–LEPR–mTOR feedback arm (M → C) will drive
   the system from a high-CSC, high-angiogenesis steady state toward
   a lower-CSC, lower-angiogenesis state, even when oncogenic Ras is
   held constant."

The procedure uses the six-variable Ras–CSC–microenvironment model
implemented in `ras_csc_model.py` and parameters loaded from
`config/model_parameters.yaml`:

  Phase 1 (Tumour establishment under full feedback):
      - Ras input f_ras is fixed at 1.0.
      - Full feedback loop is active (η_CM > 0).
      - System starts from an SCC-like "malignant" initial condition
        with high C, high A, elevated T, R, L, and M, and is simulated
        forward until it approaches a steady state.

  Phase 2 (Feedback disruption):
      - Ras remains ON (f_ras = 1.0).
      - The M → C feedback is cut by setting η_CM = 0 in a copied
        parameter set.
      - The Phase 1 steady state is used as the initial condition, and
        the system is simulated forward again until it settles.

Outputs:
  - Console summary of pre/post steady states and % drops in C and A.
  - figures/hypothesis_test_result.png           (timecourse plot)
  - results/hypothesis_test_summary.csv          (pre/post steady states)
  - results/hypothesis_test_timeseries.csv       (full trajectories)
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from pathlib import Path
from dataclasses import replace

from ras_csc_model import (
    RasCSCModelParams,
    ras_csc_rhs,
    load_model_params_from_yaml,
)


def simulate_phase(
    params: RasCSCModelParams,
    y0: np.ndarray,
    t_start: float,
    t_end: float,
    n_points: int = 400,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Integrate the Ras–CSC model over a given time window.

    Wraps `solve_ivp` to run a single forward simulation with specified
    parameters and initial conditions, and returns a dense time grid and
    the corresponding state trajectory.
    """
    t_span = (float(t_start), float(t_end))
    t_eval = np.linspace(t_start, t_end, n_points)

    def rhs_wrapped(t, y):
        return ras_csc_rhs(t, y, params)

    sol = solve_ivp(
        fun=rhs_wrapped,
        t_span=t_span,
        y0=np.asarray(y0, dtype=float),
        t_eval=t_eval,
        method="RK45",
        rtol=1e-6,
        atol=1e-9,
    )

    if not sol.success:
        raise RuntimeError(
            f"ODE integration failed in simulate_phase: {sol.message}"
        )

    t_out = sol.t
    y_out = sol.y.T
    return t_out, y_out


def main() -> None:
    """
    Entry point for the hypothesis test.

    Steps:
      1) Load parameters from YAML and run Phase 1 with full feedback
         from an SCC-like high state.
      2) Copy the parameters, cut M → C feedback (η_CM = 0), and run
         Phase 2 starting from the Phase 1 endpoint.
      3) Compute and print pre/post steady-state values for C and A,
         along with percentage changes.
      4) Dump pre/post steady states and full trajectories to CSV.
      5) Generate and save a plot of C, A, and M over time with the
         intervention marked.
    """
    project_root = Path(__file__).resolve().parents[1]
    config_path = project_root / "config" / "model_parameters.yaml"

    try:
        params_full = load_model_params_from_yaml(config_path)
    except Exception as exc:
        raise RuntimeError(
            f"Failed to load model parameters from {config_path}"
        ) from exc

    # SCC-like malignant initial condition:
    # [C, A, T, R, L, M]
    y0_phase1 = np.array(
        [
            0.7,   # C: CSC fraction
            0.7,   # A: angiogenesis
            10.0,  # T: TGFβ activity
            0.5,   # R: LEPR competence
            1.0,   # L: local leptin
            0.5,   # M: mTOR activity
        ],
        dtype=float,
    )

    # Phase 1 time window
    t1_start = 0.0
    t1_end = 200.0

    # Phase 1: Ras ON, full loop with η_CM as loaded
    t1, y1 = simulate_phase(
        params=params_full,
        y0=y0_phase1,
        t_start=t1_start,
        t_end=t1_end,
        n_points=400,
    )

    # Pre-treatment (Phase 1 steady state)
    C_pre = float(y1[-1, 0])
    A_pre = float(y1[-1, 1])
    T_pre = float(y1[-1, 2])
    R_pre = float(y1[-1, 3])
    L_pre = float(y1[-1, 4])
    M_pre = float(y1[-1, 5])

    # Create parameter set with M → C feedback cut
    params_cut = replace(params_full, eta_CM=0.0)

    # Phase 2 time window
    t2_start = t1_end
    t2_end = t2_start + 200.0

    # Use Phase 1 endpoint as initial condition for Phase 2
    y0_phase2 = y1[-1, :].copy()

    # Phase 2: Ras ON, but η_CM = 0
    t2, y2 = simulate_phase(
        params=params_cut,
        y0=y0_phase2,
        t_start=t2_start,
        t_end=t2_end,
        n_points=400,
    )

    # Post-treatment (Phase 2 steady state)
    C_post = float(y2[-1, 0])
    A_post = float(y2[-1, 1])
    T_post = float(y2[-1, 2])
    R_post = float(y2[-1, 3])
    L_post = float(y2[-1, 4])
    M_post = float(y2[-1, 5])

    def pct_drop(before: float, after: float) -> float:
        if before == 0.0:
            return 0.0
        return (before - after) / abs(before) * 100.0

    C_drop = pct_drop(C_pre, C_post)
    A_drop = pct_drop(A_pre, A_post)

    print("\n--- Hypothesis Test: mTOR feedback disruption ---")
    print("Ras input f_RAS is held constant throughout.")
    print("Phase 1: full loop (η_CM > 0).")
    print("Phase 2: M → C feedback cut (η_CM = 0).")
    print(f"\n[Pre-treatment steady state at t = {t1_end:.1f}]")
    print(f"  CSC fraction C*: {C_pre:.4f}")
    print(f"  Angiogenesis A*: {A_pre:.4f}")
    print(f"\n[Post-treatment steady state at t = {t2_end:.1f}]")
    print(f"  CSC fraction C*: {C_post:.4f}")
    print(f"  Angiogenesis A*: {A_post:.4f}")
    print("\n[Effect size]")
    print(f"  CSC reduction:        {C_drop:.1f}%")
    print(f"  Angiogenesis reduction: {A_drop:.1f}%")

    # Simple sanity check: did Phase 1 actually end in a high state?
    if C_pre < 0.3 or A_pre < 0.3:
        print(
            "\n[NOTE] Pre-treatment state is not strongly malignant-like "
            "(C and/or A < 0.3). The current parameter set may not place "
            "the system in the high-C/high-A regime; interpret H1 test "
            "with caution."
        )

    if C_drop > 50.0 and A_drop > 50.0:
        print(
            "\n>>> CONCLUSION: For this parameter set, behaviour is "
            "consistent with the Alternative Hypothesis (H1)."
        )
        print(
            "    Disrupting the M → C feedback substantially collapses "
            "CSC fraction and angiogenesis despite Ras being fixed."
        )
    else:
        print(
            "\n>>> CONCLUSION: For this parameter set, behaviour does not "
            "strongly support H1."
        )
        print(
            "    The malignant-like state (if present) does not collapse "
            "by >50% in both C and A when M → C is cut."
        )

    # --------------------------------------------------------------
    # CSV outputs: pre/post steady states and full trajectories
    # --------------------------------------------------------------
    results_dir = project_root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    # Pre/post summary
    summary_path = results_dir / "hypothesis_test_summary.csv"
    with summary_path.open("w", encoding="utf-8") as fh:
        fh.write("phase,time,C,A,T,R,L,M\n")
        fh.write(
            f"pre,{t1_end:.6g},{C_pre:.6g},{A_pre:.6g},"
            f"{T_pre:.6g},{R_pre:.6g},{L_pre:.6g},{M_pre:.6g}\n"
        )
        fh.write(
            f"post,{t2_end:.6g},{C_post:.6g},{A_post:.6g},"
            f"{T_post:.6g},{R_post:.6g},{L_post:.6g},{M_post:.6g}\n"
        )

    # Full time series
    timeseries_path = results_dir / "hypothesis_test_timeseries.csv"
    with timeseries_path.open("w", encoding="utf-8") as fh:
        fh.write("time,phase,C,A,T,R,L,M\n")
        # Phase 1
        for t, row in zip(t1, y1):
            C, A, T, R, L, M = row
            fh.write(
                f"{t:.6g},pre,{C:.6g},{A:.6g},{T:.6g},"
                f"{R:.6g},{L:.6g},{M:.6g}\n"
            )
        # Phase 2
        for t, row in zip(t2, y2):
            C, A, T, R, L, M = row
            fh.write(
                f"{t:.6g},post,{C:.6g},{A:.6g},{T:.6g},"
                f"{R:.6g},{L:.6g},{M:.6g}\n"
            )

    # --------------------------------------------------------------
    # Timecourse figure
    # --------------------------------------------------------------
    t_full = np.concatenate([t1, t2])
    y_full = np.vstack([y1, y2])

    C_series = y_full[:, 0]
    A_series = y_full[:, 1]
    M_series = y_full[:, 5]

    figures_dir = project_root / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    out_path = figures_dir / "hypothesis_test_result.png"

    plt.figure(figsize=(10, 6))

    plt.plot(t_full, C_series, label="CSC fraction C", linewidth=2.0)
    plt.plot(t_full, A_series, label="Angiogenesis A", linewidth=2.0)
    plt.plot(
        t_full,
        M_series,
        label="mTOR activity M",
        linewidth=1.5,
        linestyle="--",
    )

    plt.axvline(
        x=t1_end,
        color="black",
        linestyle=":",
        linewidth=1.5,
        label="Cut M → C (η_CM = 0)",
    )

    plt.title("Ras–CSC–microenvironment model: mTOR feedback disruption test")
    plt.xlabel("Time (arbitrary units)")
    plt.ylabel("State variables (dimensionless)")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)

    print(f"\n[INFO] Plot saved to {out_path}")
    print(f"[INFO] Summary CSV written to {summary_path}")
    print(f"[INFO] Time-series CSV written to {timeseries_path}")


if __name__ == "__main__":
    main()
