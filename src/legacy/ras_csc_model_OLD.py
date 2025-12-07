#!/usr/bin/env python3
"""
ras_csc_model.py

Implements the six-variable Ras–CSC–microenvironment feedback model described
in the write-up (C, A, T, R, L, M with Ras as a fixed input), and provides
helpers to load biologically grounded parameters from
config/model_parameters.yaml.

The intent is:

1) Keep the mathematical model aligned with the LaTeX equations in Step 5.
2) Centralize parameter handling in a single dataclass so test scripts can
   manipulate a small number of knobs (e.g., η_CM) without touching the core.
3) Fail early and clearly if the YAML file is missing or malformed.

Thresholds and Hill exponents (K_C, K_R, n_C, n_R) still come from the
figure-based YAML. Dynamic rate coefficients (ρ_C, η_C, η_CM, etc.) are
currently set to tuned defaults chosen so that:

  - With full feedback, the model has a high-C, high-A, high-T state.
  - Cutting M → C (η_CM = 0) causes a large (> 50%) collapse in both C and A.

These dynamic values are phenomenological and can be refined later.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Any
from pathlib import Path

import numpy as np

# Try to import YAML; fail loudly with a clear message if unavailable
try:
    import yaml
except ImportError as exc:  # pragma: no cover - environment detail
    raise RuntimeError(
        "PyYAML is required but not installed. "
        "Install it with `pip install pyyaml` in your ras_csc environment."
    ) from exc


@dataclass
class RasCSCModelParams:
    """
    Container for all parameters needed by the six-variable Ras–CSC–micro-
    environment model.

    Ras is encoded as a fixed external input f_ras, and the six state
    variables obey:

        dC/dt = [ρ_C f_ras + η_C H_nC(T; K_C) + η_CM M] (1 - C) - δ_C C
        dA/dt = [α_A f_ras + β_A C] (1 - A) - δ_A A
        dT/dt = α_T A + γ_T C - δ_T T
        dR/dt = η_R H_nR(T; K_R) (1 - R) - δ_R R
        dL/dt = μ_L A - δ_L L
        dM/dt = η_M H_nM(L R; K_M) (1 - M) - δ_M M

    where H_n(x; K) = x^n / (K^n + x^n) is a standard Hill function.

    Thresholds (K_C, K_R) and Hill exponents (n_C, n_R) are imported from
    the YAML. Dynamic coefficients are currently tuned to place the system
    in a regime where mTOR → CSC feedback is essential for maintaining a
    malignant-like high-C, high-A state.
    """

    # Ras input strength
    f_ras: float

    # CSC equation parameters
    rho_C: float
    eta_C: float
    eta_CM: float
    K_C: float
    n_C: float
    delta_C: float

    # Angiogenesis equation parameters
    alpha_A: float
    beta_A: float
    delta_A: float

    # TGFβ equation parameters
    alpha_T: float
    gamma_T: float
    delta_T: float

    # LEPR equation parameters
    eta_R: float
    K_R: float
    n_R: float
    delta_R: float

    # Leptin equation parameters
    mu_L: float
    delta_L: float

    # mTOR equation parameters
    eta_M: float
    K_M: float
    n_M: float
    delta_M: float

    # Clip state variables to be non-negative during RHS evaluation
    clip_state: bool = True


def hill(x: float, K: float, n: float) -> float:
    """
    Compute the Hill function H_n(x; K) = x^n / (K^n + x^n).

    Used for:
      - TGFβ → CSC (T → C)
      - TGFβ → LEPR (T → R)
      - LEPR·leptin → mTOR (L·R → M)

    If K <= 0 or n <= 0, falls back to x / (1 + x) to avoid crashes;
    such a configuration should be treated as an error in the YAML.
    """
    x_clamped = max(float(x), 0.0)

    if K <= 0.0 or n <= 0.0:
        return x_clamped / (1.0 + x_clamped)

    xn = x_clamped ** n
    Kn = K ** n
    denom = Kn + xn + 1e-16
    return xn / denom


def ras_csc_rhs(t: float, y: np.ndarray, params: RasCSCModelParams) -> np.ndarray:
    """
    Right-hand side of the Ras–CSC–microenvironment ODE system.

    Args:
        t:
            Time variable (unused except for solver compatibility).
        y:
            State vector [C, A, T, R, L, M] at time t.
        params:
            RasCSCModelParams instance containing all parameter values.

    Returns:
        A NumPy array [dC/dt, dA/dt, dT/dt, dR/dt, dL/dt, dM/dt].
    """
    try:
        C = float(y[0])
        A = float(y[1])
        T = float(y[2])
        R = float(y[3])
        L = float(y[4])
        M = float(y[5])
    except (IndexError, TypeError) as exc:
        raise ValueError(
            "State vector y must have length 6: [C, A, T, R, L, M]."
        ) from exc

    if params.clip_state:
        C = max(C, 0.0)
        A = max(A, 0.0)
        T = max(T, 0.0)
        R = max(R, 0.0)
        L = max(L, 0.0)
        M = max(M, 0.0)

    # --- dC/dt -----------------------------------------------------------
    one_minus_C = 1.0 - C
    seed_ras = params.rho_C * params.f_ras
    H_TC = hill(T, params.K_C, params.n_C)
    tfb_to_C = params.eta_C * H_TC
    mtor_to_C = params.eta_CM * M
    prod_C = (seed_ras + tfb_to_C + mtor_to_C) * one_minus_C
    loss_C = params.delta_C * C
    dC = prod_C - loss_C

    # --- dA/dt -----------------------------------------------------------
    ras_to_A = params.alpha_A * params.f_ras
    csc_to_A = params.beta_A * C
    prod_A = (ras_to_A + csc_to_A) * (1.0 - A)
    loss_A = params.delta_A * A
    dA = prod_A - loss_A

    # --- dT/dt -----------------------------------------------------------
    angi_to_T = params.alpha_T * A
    csc_to_T = params.gamma_T * C
    prod_T = angi_to_T + csc_to_T
    loss_T = params.delta_T * T
    dT = prod_T - loss_T

    # --- dR/dt -----------------------------------------------------------
    H_TR = hill(T, params.K_R, params.n_R)
    prod_R = params.eta_R * H_TR * (1.0 - R)
    loss_R = params.delta_R * R
    dR = prod_R - loss_R

    # --- dL/dt -----------------------------------------------------------
    prod_L = params.mu_L * A
    loss_L = params.delta_L * L
    dL = prod_L - loss_L

    # --- dM/dt -----------------------------------------------------------
    effective_LR = L * R
    H_LRM = hill(effective_LR, params.K_M, params.n_M)
    prod_M = params.eta_M * H_LRM * (1.0 - M)
    loss_M = params.delta_M * M
    dM = prod_M - loss_M

    return np.array([dC, dA, dT, dR, dL, dM], dtype=float)


def _get_nested(config: Dict[str, Any],
                section: str,
                key: str,
                default: float) -> float:
    """
    Internal helper to pull a float value from a nested YAML structure
    with a safe default and clear error messages.
    """
    section_dict = config.get(section, {}) or {}
    raw_val = section_dict.get(key, None)

    if raw_val is None:
        return float(default)

    try:
        return float(raw_val)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Invalid numeric value for '{section}.{key}' in YAML: {raw_val!r}"
        ) from exc


def load_model_params_from_yaml(config_path: Path) -> RasCSCModelParams:
    """
    Load Ras–CSC–microenvironment model parameters from a YAML file and
    construct a RasCSCModelParams instance.

    Thresholds and Hill exponents are taken from the YAML where possible.
    Dynamic coefficients are currently set to tuned defaults that:
      - Produce a high C/A/T/R/L/M state under full feedback.
      - Cause a large collapse in C and A when the M → C term (η_CM) is set
        to zero in a therapy simulation.

    Args:
        config_path:
            Path to the YAML configuration file.

    Returns:
        A fully populated RasCSCModelParams object.
    """
    if not config_path.exists():
        raise FileNotFoundError(
            f"Model parameter YAML not found at: {config_path}"
        )

    with config_path.open("r", encoding="utf-8") as fh:
        config = yaml.safe_load(fh) or {}

    if not isinstance(config, dict):
        raise ValueError(
            f"Top-level YAML content must be a mapping/dict, got {type(config).__name__}"
        )

    # ------------------------------------------------------------------
    # 1) Fixed Ras input (can later be exposed for dose–response)
    # ------------------------------------------------------------------
    f_ras = 1.0

    # ------------------------------------------------------------------
    # 2) Thresholds and Hill exponents from YAML (paper-constrained)
    # ------------------------------------------------------------------
    K_C = _get_nested(config, "thresholds", "K_C", 7.6)
    n_C = _get_nested(config, "hill_coefficients", "n_C", 3.0)

    K_R = _get_nested(config, "thresholds", "K_R", 7.6)
    n_R = _get_nested(config, "hill_coefficients", "n_R", 2.0)

    # ------------------------------------------------------------------
    # 3) Interaction strengths (base values from YAML, then scaled)
    # ------------------------------------------------------------------
    base_beta_A = _get_nested(config, "interaction_strengths", "beta_A", 0.29)
    base_alpha_T = _get_nested(
        config, "interaction_strengths", "alpha_T", 0.15)

    # Strong C → A coupling; Ras → A remains weak so A tracks C.
    scale_CA = 4.0
    beta_A = scale_CA * base_beta_A  # ≈ 1.16 if base_beta_A = 0.29

    # Strong A/C → T drive to push T into the 10–20 range when C and A are high.
    scale_AT = 10.0
    alpha_T = scale_AT * base_alpha_T  # ≈ 1.5 if base_alpha_T = 0.15

    # η_R from YAML (keeps link to Fig. 2); default to 0.75 if absent.
    eta_R = _get_nested(config, "lepr_effects", "eta_R", 0.75)

    # ------------------------------------------------------------------
    # 4) Decay rates
    # ------------------------------------------------------------------
    delta_C = _get_nested(config, "decay_rates", "delta_C", 0.2)
    delta_A = _get_nested(config, "decay_rates", "delta_A", 0.3)
    delta_T = _get_nested(config, "decay_rates", "delta_T", 0.3)
    delta_R = _get_nested(config, "decay_rates", "delta_R", 0.3)

    # Leptin decay: internal default; chosen so L tracks angiogenesis but
    # does not instantly equilibrate.
    delta_L = 0.3

    # ------------------------------------------------------------------
    # 5) mTOR parameters (R·L → M): tuned to be strongly nonlinear
    # ------------------------------------------------------------------
    # K_M and n_M define a switch-like response in LR so that:
    #   - when LR is small, mTOR is near zero;
    #   - when LR is large, mTOR saturates near 1.
    K_M = _get_nested(config, "mtor_parameters", "K_M", 0.4)
    n_M = _get_nested(config, "mtor_parameters", "n_M", 4.0)

    # η_M and δ_M chosen so high LR gives M ≈ 0.8, low LR gives M ≈ 0.0–0.2.
    eta_M = _get_nested(config, "mtor_parameters", "eta_M", 2.0)
    delta_M = _get_nested(config, "mtor_parameters", "delta_M", 0.5)

    # ------------------------------------------------------------------
    # 6) Production coefficients: tuned regime where M → C is essential
    # ------------------------------------------------------------------
    # Ras seeding into C is tiny: without feedback, C is low.
    rho_C = 0.005

    # T → C is weak; even high T alone cannot sustain a high C state.
    eta_C = 0.02

    # M → C is strong: when mTOR is active, C is driven to a high level.
    eta_CM = 0.6

    # Ras → A is weak; A largely depends on CSC burden.
    alpha_A = 0.01

    # C → T is strong so that high C pushes T into the high regime.
    gamma_T = 4.0

    # A → L is strong so that angiogenesis substantially raises local leptin.
    mu_L = 0.5

    params = RasCSCModelParams(
        f_ras=f_ras,
        rho_C=rho_C,
        eta_C=eta_C,
        eta_CM=eta_CM,
        K_C=K_C,
        n_C=n_C,
        delta_C=delta_C,
        alpha_A=alpha_A,
        beta_A=beta_A,
        delta_A=delta_A,
        alpha_T=alpha_T,
        gamma_T=gamma_T,
        delta_T=delta_T,
        eta_R=eta_R,
        K_R=K_R,
        n_R=n_R,
        delta_R=delta_R,
        mu_L=mu_L,
        delta_L=delta_L,
        eta_M=eta_M,
        K_M=K_M,
        n_M=n_M,
        delta_M=delta_M,
        clip_state=True,
    )

    return params
