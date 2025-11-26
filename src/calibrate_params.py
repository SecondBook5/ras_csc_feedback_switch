#!/usr/bin/env python3
"""
Calibrate Ras–CSC feedback model parameters from GSE190411 Pap/SCC module scores.

This script implements Step 4.1 and 4.2 of the plan:

  - Step 4.1: Non-dimensionalize Pap vs SCC module scores so that
              Pap -> 0 and SCC -> 1 for each module X in {C, V, R, M}.

  - Step 4.2: Choose simple decay rates and solve the steady-state
              equations for a minimal 5-variable feedback model so that
              the dimensionless Pap and SCC points

                  Pap*:  (C*, V*, L*, R*, M*) = (0, 0, 0, 0, 0)
                  SCC*:  (1, 1, 1, 1, 1)

              are fixed points of Model A (with feedback).

IMPORTANT:
    - The Pap/SCC module scores from GSE190411 are *real data* and are
      used to define the scaling between model state variables and
      log2(TPM) space.

    - The kinetic parameters (k_* and K_*) are not uniquely identifiable
      from 8 bulk RNA-seq samples. Here we choose them so that Pap* and
      SCC* are equilibria in a simple ODE system, with time scales set
      by hand (explicitly documented).

    - Model B (no feedback) is *deliberately* not calibrated to match
      both Pap and SCC for C, because the structural argument in the
      plan is that this is impossible with fixed Ras.

Outputs:
    data/processed/model_calibration_results.json
        {
          "targets": { ... Pap/SCC raw scores and scales ... },
          "model_A": { ... calibrated parameters ... },
          "model_B": { ... no-feedback parameters ... }
        }
"""

from __future__ import annotations

# Import standard libraries
import json
import os
from dataclasses import asdict, dataclass
from typing import Dict

# Import third-party libraries
import pandas as pd


# ---------------------------------------------------------------------
# Utilities for paths
# ---------------------------------------------------------------------
def get_repo_root() -> str:
    """
    Find the repository root by going up from this file's directory.

    The script is expected to live under 'src/' in the repo. We go one
    level up to reach the project root and assert that it exists.
    """
    # Get the absolute path of this script
    here: str = os.path.abspath(__file__)
    # Get directory containing this script (should be 'src/')
    src_dir: str = os.path.dirname(here)
    # Move one level up to reach repo root
    repo_root: str = os.path.dirname(src_dir)

    # Verify that repo_root exists
    if not os.path.isdir(repo_root):
        # Raise an error if the layout is not as expected
        raise RuntimeError(
            f"[ERROR] Computed repo root does not exist: {repo_root}"
        )

    # Return the verified repo root path
    return repo_root


# ---------------------------------------------------------------------
# Dataclasses for targets and parameters
# ---------------------------------------------------------------------
@dataclass
class PapSccTargets:
    """
    Store Pap vs SCC module scores and linear scaling factors.

    The idea is to make a simple linear map between model state variables
    (dimensionless) and module scores (log2(TPM+1) or similar).

    For each module X in {C, V, R, M}, we define:

        X_model = 0  ->  X_raw = X_pap
        X_model = 1  ->  X_raw = X_scc

    So:

        X_raw = X_offset + X_scale * X_model
              = X_pap + (X_scc - X_pap) * X_model
    """

    # Raw Pap (mpos) scores
    C_pap: float
    V_pap: float
    R_pap: float
    M_pap: float

    # Raw SCC (mpos) scores
    C_scc: float
    V_scc: float
    R_scc: float
    M_scc: float

    # Linear scaling (offset + scale * X_model)
    C_offset: float
    C_scale: float
    V_offset: float
    V_scale: float
    R_offset: float
    R_scale: float
    M_offset: float
    M_scale: float


@dataclass
class ModelParams:
    """
    Store a set of kinetic parameters for the Ras–CSC feedback model.

    This is kept generic as a dict-like container; ras_csc_model.py
    will later unpack these into its RasCSCParams dataclass.

    Fields:
        params: Mapping from parameter name to float.
    """
    params: Dict[str, float]


# ---------------------------------------------------------------------
# Load Pap/SCC module scores and build targets
# ---------------------------------------------------------------------
def load_pap_scc_targets() -> PapSccTargets:
    """
    Load Pap vs SCC module scores for mpos samples and compute offsets/scales.

    This function:
        1. Reads 'module_scores_pap_scc.csv' in data/processed/.
        2. Extracts the Pap mpos and SCC mpos rows.
        3. Constructs a PapSccTargets object with:
             - raw Pap/SCC scores,
             - X_offset = X_pap,
             - X_scale  = X_scc - X_pap.

    Raises:
        FileNotFoundError: If the processed module score file is missing.
        RuntimeError: If Pap/SCC mpos rows are not uniquely defined.

    Returns:
        PapSccTargets: Dataclass with raw scores and scaling factors.
    """
    # Resolve repo root and scores path
    repo_root: str = get_repo_root()
    scores_path: str = os.path.join(
        repo_root, "data", "processed", "module_scores_pap_scc.csv"
    )

    # Check that the module score file exists
    if not os.path.exists(scores_path):
        raise FileNotFoundError(
            f"[ERROR] Cannot find module_scores_pap_scc.csv at: {scores_path}\n"
            f"Run module_scores_pap_scc.py before calibrating parameters."
        )

    # Load CSV into a DataFrame
    df = pd.read_csv(scores_path)

    # Ensure required columns are present
    required_cols = {
        "condition", "stem_subset", "C_score", "V_score", "R_score", "M_score"
    }
    missing = required_cols - set(df.columns)
    if missing:
        raise RuntimeError(
            f"[ERROR] module_scores_pap_scc.csv is missing required columns: {missing}"
        )

    # Filter for Pap mpos and SCC mpos
    df_pap = df[
        (df["condition"] == "Pap") & (df["stem_subset"] == "mpos")
    ]
    df_scc = df[
        (df["condition"] == "SCC") & (df["stem_subset"] == "mpos")
    ]

    # We expect exactly one row for each
    if df_pap.shape[0] != 1 or df_scc.shape[0] != 1:
        raise RuntimeError(
            "[ERROR] Expected exactly one Pap mpos and one SCC mpos row in "
            "module_scores_pap_scc.csv, but found "
            f"{df_pap.shape[0]} Pap mpos and {df_scc.shape[0]} SCC mpos rows."
        )

    # Extract as Python dicts
    pap_row = df_pap.iloc[0].to_dict()
    scc_row = df_scc.iloc[0].to_dict()

    # Raw scores
    C_pap = float(pap_row["C_score"])
    V_pap = float(pap_row["V_score"])
    R_pap = float(pap_row["R_score"])
    M_pap = float(pap_row["M_score"])

    C_scc = float(scc_row["C_score"])
    V_scc = float(scc_row["V_score"])
    R_scc = float(scc_row["R_score"])
    M_scc = float(scc_row["M_score"])

    # Compute scales (SCC - Pap)
    C_scale = C_scc - C_pap
    V_scale = V_scc - V_pap
    R_scale = R_scc - R_pap
    M_scale = M_scc - M_pap

    # Guard against degenerate scaling
    if abs(C_scale) < 1e-8:
        raise RuntimeError("[ERROR] ΔC (SCC - Pap) is ~0; cannot scale C.")
    if abs(R_scale) < 1e-8:
        raise RuntimeError("[ERROR] ΔR (SCC - Pap) is ~0; cannot scale R.")
    if abs(M_scale) < 1e-8:
        raise RuntimeError("[ERROR] ΔM (SCC - Pap) is ~0; cannot scale M.")
    if abs(V_scale) < 1e-8:
        raise RuntimeError("[ERROR] ΔV (SCC - Pap) is ~0; cannot scale V.")

    # Build and return targets
    return PapSccTargets(
        C_pap=C_pap,
        V_pap=V_pap,
        R_pap=R_pap,
        M_pap=M_pap,
        C_scc=C_scc,
        V_scc=V_scc,
        R_scc=R_scc,
        M_scc=M_scc,
        C_offset=C_pap,
        C_scale=C_scale,
        V_offset=V_pap,
        V_scale=V_scale,
        R_offset=R_pap,
        R_scale=R_scale,
        M_offset=M_pap,
        M_scale=M_scale,
    )


# ---------------------------------------------------------------------
# Analytic calibration of Model A parameters
# ---------------------------------------------------------------------
def calibrate_model_a(targets: PapSccTargets) -> ModelParams:
    """
    Calibrate a simple 5-variable feedback model so that Pap* and SCC*
    are steady states in dimensionless space.

    We work entirely in *dimensionless* variables here. The link between
    those and the real log2(TPM+1) scores is encoded in 'targets' via
    offsets and scales.

    Model A structure (dimensionless, with Pap* = 0, SCC* = 1):

        dC/dt = k_C_M * M - k_C_decay * C
        dV/dt = k_V_C * C - k_V_deg * V
        dL/dt = k_L_V * V - k_L_deg * L
        dR/dt = k_R_TGFb * u_TGFb - k_R_deg * R
        dM/dt = k_M_act * Hill(L * R; K_M, m_M) - k_M_deg * M

    with:

        Hill(x; K, m) = x^m / (K^m + x^m)

    Assumptions:
        - Pap* is the origin: (C, V, L, R, M) = (0, 0, 0, 0, 0)
        - SCC* is (1, 1, 1, 1, 1) under "high TGFβ" (u_TGFb = 1).
        - Time scales (decay rates) are chosen to be 1.0 for simplicity.
        - Ras baseline term is set to 0 in this reduced model, so that
          the difference between Pap and SCC is carried entirely by the
          positive feedback loop and TGFβ-driven LEPR.

    From the steady-state conditions, we choose:

        k_C_decay = 1.0, k_C_M = 1.0
        k_V_deg   = 1.0, k_V_C = 1.0
        k_L_deg   = 1.0, k_L_V = 1.0
        k_R_deg   = 1.0, k_R_TGFb = 1.0, u_TGFb_SCC = 1.0
        k_M_deg   = 1.0, K_M = 0.5, m_M = 2.0
        k_M_act   = 1.0 / Hill(1; K_M, m_M)

    Returns:
        ModelParams: Wrapper around a dict of parameter values for Model A.
    """
    # Set decay rates (time scales) to 1.0
    k_C_decay: float = 1.0
    k_V_deg: float = 1.0
    k_L_deg: float = 1.0
    k_R_deg: float = 1.0
    k_M_deg: float = 1.0

    # Feedback from M to C chosen so that SCC*=(1,1,1,1,1) is a fixed point
    k_C_M: float = k_C_decay  # so that M=1, C=1 solves 0 = k_C_M*1 - k_C_decay*1

    # V dynamics chosen similarly: V tracks C at steady state
    k_V_C: float = k_V_deg

    # L dynamics: L tracks V at steady state
    k_L_V: float = k_L_deg

    # R dynamics: R tracks u_TGFb at steady state (for u_TGFb = 1 -> R=1)
    k_R_TGFb: float = k_R_deg
    u_TGFb: float = 1.0

    # M dynamics: choose Hill parameters, then solve for k_M_act such that
    # L=R=1 => M=1 is steady state.
    K_M: float = 0.5
    m_M: float = 2.0

    # Compute Hill(1; K_M, m_M)
    hill_num: float = 1.0 ** m_M
    hill_den: float = (K_M ** m_M) + (1.0 ** m_M)
    hill_1: float = hill_num / hill_den

    # Avoid division by zero if something pathological happens
    if hill_1 <= 0.0:
        raise RuntimeError(
            f"[ERROR] Hill(1; K={K_M}, m={m_M}) <= 0; cannot solve for k_M_act."
        )

    # Solve for k_M_act from steady-state condition at SCC*:
    #   0 = k_M_act * Hill(1) - k_M_deg * 1  => k_M_act = k_M_deg / Hill(1)
    k_M_act: float = k_M_deg / hill_1

    # We are not using a Hill nonlinearity on C itself in this reduced
    # calibration, but ras_csc_model.py may include such a term. For now,
    # set K_C and n_C to reasonable sigmoid-shaping values.
    K_C: float = 0.3
    n_C: float = 3.0

    # Baseline CSC production from Ras in this reduced model is set to 0.0,
    # consistent with Pap* = 0 under M=0. The structural argument about
    # Ras being fixed is handled separately in the report (Model B).
    k_C_base: float = 0.0

    # Baseline leptin production term set to 0 so Pap* has L=0.
    k_L0: float = 0.0

    # In this reduced calibration, we treat "clip_state" as True so that
    # ras_csc_model can optionally bound state variables in [0, 1].
    clip_state: float = 1.0  # stored as 1.0; ras_csc_model can cast to bool

    # Pack all parameters relevant to ras_csc_model.RasCSCParams
    params: Dict[str, float] = {
        "k_C_base": k_C_base,
        "k_C_M": k_C_M,
        "K_C": K_C,
        "n_C": n_C,
        "k_C_decay": k_C_decay,
        "k_V_C": k_V_C,
        "k_V_deg": k_V_deg,
        "k_L0": k_L0,
        "k_L_V": k_L_V,
        "k_L_deg": k_L_deg,
        "k_R_TGFb": k_R_TGFb,
        "u_TGFb": u_TGFb,
        "k_R_deg": k_R_deg,
        "k_M_act": k_M_act,
        "K_M": K_M,
        "m_M": m_M,
        "k_M_deg": k_M_deg,
        # We store clip_state as 1.0; ras_csc_model can convert to bool
        "clip_state": clip_state,
    }

    # Wrap and return as ModelParams
    return ModelParams(params=params)


def define_model_b_from_model_a(model_a: ModelParams) -> ModelParams:
    """
    Define Model B (no-feedback) parameters based on Model A.

    Model B removes the M->C feedback by setting k_C_M = 0. The point of
    Model B is *not* to match both Pap and SCC; rather, it is to serve as
    the structurally restricted comparator:

        dC/dt = Ras * k_C_base - k_C_decay * C

    Under fixed Ras in Model B:

        C* = Ras * k_C_base / k_C_decay

    so there is a single CSC steady state for a given Ras, which cannot
    simultaneously represent both Pap and SCC CSC states.

    Here, we simply copy Model A's parameters and set k_C_M = 0.0. You
    can optionally adjust k_C_base in the report if you want C* to sit
    near the Pap value, but the structural impossibility result does not
    depend on its exact value.
    """
    # Copy Model A parameters
    params_b: Dict[str, float] = dict(model_a.params)

    # Remove M->C feedback
    params_b["k_C_M"] = 0.0

    # Set a small baseline CSC production so that C* > 0 in Model B
    # when Ras is nonzero (Ras will be chosen in ras_csc_model.py).
    params_b["k_C_base"] = 0.1

    # Return wrapped object
    return ModelParams(params=params_b)


# ---------------------------------------------------------------------
# Main calibration routine
# ---------------------------------------------------------------------
def main() -> None:
    """
    Run the full calibration workflow:

        1. Load Pap/SCC module scores from GSE190411 (Pap/SCC mpos).
        2. Build PapSccTargets with raw scores and linear scaling.
        3. Calibrate Model A parameters so Pap* and SCC* are steady states.
        4. Define Model B parameters by removing M->C feedback.
        5. Write everything to 'data/processed/model_calibration_results.json'.

    The JSON file is later read by ras_csc_model.py to instantiate
    RasCSCParams for both models in a reproducible way.
    """
    # Load Pap/SCC module target information
    targets = load_pap_scc_targets()

    # Calibrate Model A (feedback)
    model_a = calibrate_model_a(targets)

    # Build Model B (no-feedback) from Model A
    model_b = define_model_b_from_model_a(model_a)

    # Prepare JSON payload
    payload = {
        "targets": asdict(targets),
        "model_A": model_a.params,
        "model_B": model_b.params,
    }

    # Resolve output path
    repo_root: str = get_repo_root()
    out_dir: str = os.path.join(repo_root, "data", "processed")
    os.makedirs(out_dir, exist_ok=True)
    out_path: str = os.path.join(out_dir, "model_calibration_results.json")

    # Write to JSON
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    # Print a concise summary
    print(f"[OK] Wrote calibrated parameters and targets to: {out_path}")
    print("\nPap (mpos) raw module scores:")
    print(f"  C = {targets.C_pap:.6f}, V = {targets.V_pap:.6f}, "
          f"R = {targets.R_pap:.6f}, M = {targets.M_pap:.6f}")
    print("SCC (mpos) raw module scores:")
    print(f"  C = {targets.C_scc:.6f}, V = {targets.V_scc:.6f}, "
          f"R = {targets.R_scc:.6f}, M = {targets.M_scc:.6f}")
    print("\nModel A parameters (feedback):")
    for k, v in model_a.params.items():
        print(f"  {k} = {v:.6f}")
    print("\nModel B parameters (no feedback, k_C_M forced to 0):")
    for k, v in model_b.params.items():
        print(f"  {k} = {v:.6f}")


if __name__ == "__main__":
    main()
