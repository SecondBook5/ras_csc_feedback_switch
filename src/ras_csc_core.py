#!/usr/bin/env python3
"""
ras_csc_core.py

Shared core for the Ras–CSC–microenvironment ODE system.

This module is meant to be the *single source of truth* for:
- The state variables and their ordering.
- The Hill functions and parameter usage.
- The steady-state simulation routine.

All higher-level scripts (calibration, hysteresis tests, perturbations)
should import from here instead of re-implementing the ODEs independently.
"""

from __future__ import annotations

# Import typing utilities for better clarity and static checking
from typing import Dict, List, Optional, Sequence

# Import numpy for numerical operations
import numpy as np

# Import odeint for ODE integration
from scipy.integrate import odeint


# ---------------------------------------------------------------------------
# Hill function
# ---------------------------------------------------------------------------

def hill(x: float, K: float, n: float) -> float:
    """
    Compute a Hill-type activation function.

    This function encodes a saturating, sigmoidal response of the form:

        H(x) = x^n / (K^n + x^n)

    where:
      - x is the input (e.g., ligand or upstream signal),
      - K is the half-maximal constant,
      - n is the Hill coefficient controlling steepness.

    The inputs are clipped into biologically reasonable ranges to avoid
    numerical overflow or division-by-zero issues during integration.

    Args:
        x: Input signal level.
        K: Half-saturation constant.
        n: Hill coefficient controlling steepness.

    Returns:
        The Hill function value between 0 and 1.
    """
    # Clip x to a non-negative, finite range to avoid overflow
    x_clipped: float = float(np.clip(x, 0.0, 10.0))
    # Clip K away from zero so K^n is defined and finite
    K_clipped: float = float(np.clip(K, 1e-3, 10.0))
    # Clip n to a moderate range to avoid extreme exponents
    n_clipped: float = float(np.clip(n, 1.0, 6.0))

    # Compute denominator safely
    denom: float = K_clipped ** n_clipped + x_clipped ** n_clipped

    # Guard against division by zero just in case
    if denom == 0.0:
        return 0.0

    # Return the Hill response
    return (x_clipped ** n_clipped) / denom


# ---------------------------------------------------------------------------
# ODE right-hand side
# ---------------------------------------------------------------------------

def ras_csc_rhs(
    y: Sequence[float],
    t: float,
    f_ras: float,
    params: Dict[str, float],
) -> List[float]:
    """
    Right-hand side of the Ras–CSC–microenvironment ODE system.

    The state vector y has the following components:
        C: Cancer stem cell fraction
        A: Angiogenesis level
        T: TGFβ level
        R: LEPR signaling level
        L: Leptin level
        M: mTOR pathway activity

    The equations implement the positive feedback loop:

        Ras → Angio → TGFβ → LEPR → Leptin → mTOR → CSC

    plus direct Ras → CSC seeding and moderate TGFβ → CSC input.

    Args:
        y: Current state vector [C, A, T, R, L, M].
        t: Time (unused but required by the ODE solver).
        f_ras: Ras input level (0–1).
        params: Dictionary of kinetic parameters. Must contain keys such as:
            - rho_C, eta_C_M, K_C, n_C, delta_C
            - alpha_A, beta_A, delta_A
            - alpha_T, gamma_T, delta_T
            - eta_R, K_R, p_R, delta_R
            - mu_L, delta_L
            - eta_M, K_M, q_M, delta_M
            - eta_C (TGFβ → CSC)

    Returns:
        Time derivatives [dC_dt, dA_dt, dT_dt, dR_dt, dL_dt, dM_dt].
    """
    # Unpack state variables
    C, A, T, R, L, M = y

    # Clip all states into biologically reasonable ranges
    C = float(np.clip(C, 0.0, 1.0))
    A = float(np.clip(A, 0.0, 1.0))
    T = float(np.clip(T, 0.0, 1.0))
    R = float(np.clip(R, 0.0, 1.0))
    L = float(np.clip(L, 0.0, 1.0))
    M = float(np.clip(M, 0.0, 1.0))

    # Safely retrieve parameters with explicit defaults for optional ones
    rho_C: float = float(params.get("rho_C", 0.1))
    eta_C_M: float = float(params.get("eta_C_M", 1.0))
    K_C: float = float(params.get("K_C", 0.3))
    n_C: float = float(params.get("n_C", 3.0))
    delta_C: float = float(params.get("delta_C", 1.0))

    alpha_A: float = float(params.get("alpha_A", 0.4))
    beta_A: float = float(params.get("beta_A", 1.0))
    delta_A: float = float(params.get("delta_A", 1.0))

    alpha_T: float = float(params.get("alpha_T", 0.8))
    gamma_T: float = float(params.get("gamma_T", 0.15))
    delta_T: float = float(params.get("delta_T", 1.0))

    eta_R: float = float(params.get("eta_R", 0.6))
    K_R: float = float(params.get("K_R", 0.35))
    p_R: float = float(params.get("p_R", 2.5))
    delta_R: float = float(params.get("delta_R", 0.6))

    mu_L: float = float(params.get("mu_L", 0.9))
    delta_L: float = float(params.get("delta_L", 0.6))

    eta_M: float = float(params.get("eta_M", 1.0))
    K_M: float = float(params.get("K_M", 0.25))
    q_M: float = float(params.get("q_M", 3.5))
    delta_M: float = float(params.get("delta_M", 0.6))

    eta_C: float = float(params.get("eta_C", 1.2))

    # Compute dC/dt: CSC dynamics with Ras seeding + mTOR + TGFβ inputs
    dC_dt: float = (
        rho_C * f_ras
        + eta_C_M * M
        + eta_C * hill(T, K_C, n_C)
    ) * (1.0 - C) - delta_C * C

    # Compute dA/dt: Angiogenesis controlled by Ras and CSCs
    dA_dt: float = (
        alpha_A * f_ras
        + beta_A * C
    ) * (1.0 - A) - delta_A * A

    # Compute dT/dt: TGFβ from Angio and CSCs
    dT_dt: float = alpha_T * A + gamma_T * C - delta_T * T

    # Compute dR/dt: LEPR induced by TGFβ
    dR_dt: float = eta_R * hill(T, K_R, p_R) * (1.0 - R) - delta_R * R

    # Compute dL/dt: Leptin proportional to Angio
    dL_dt: float = mu_L * A - delta_L * L

    # Compute dM/dt: mTOR activated by L * R
    dM_dt: float = eta_M * hill(L * R, K_M, q_M) * (1.0 - M) - delta_M * M

    # Return time derivatives as a list
    return [dC_dt, dA_dt, dT_dt, dR_dt, dL_dt, dM_dt]


# ---------------------------------------------------------------------------
# Steady-state simulation
# ---------------------------------------------------------------------------

def simulate_steady_state(
    f_ras: float,
    params: Dict[str, float],
    y0: Optional[Sequence[float]] = None,
    t_max: float = 200.0,
    n_steps: int = 2000,
) -> np.ndarray:
    """
    Simulate the ODE system forward in time until approximate steady state.

    This function integrates the ODEs from an initial condition y0 up to
    time t_max, returning the final state vector. It applies clipping and
    basic NaN/Inf checks to avoid propagating numerical pathologies into
    downstream calibration or hypothesis tests.

    Args:
        f_ras: Ras input level used for the simulation.
        params: Dictionary of kinetic parameters.
        y0: Optional initial condition; if None, a default moderate state is used.
        t_max: Final simulation time.
        n_steps: Number of time steps for integration.

    Returns:
        A numpy array of the final state [C, A, T, R, L, M].
        If integration fails, a default moderate state is returned.
    """
    # Use a moderate default initial state if none is provided
    if y0 is None:
        y0 = [0.3, 0.3, 0.2, 0.2, 0.2, 0.3]

    # Build time grid for integration
    t: np.ndarray = np.linspace(0.0, float(t_max), int(n_steps))

    try:
        # Integrate the ODEs using odeint
        sol: np.ndarray = odeint(
            ras_csc_rhs,
            y0,
            t,
            args=(float(f_ras), params),
            atol=1e-6,
            rtol=1e-6,
        )

        # Extract the final state
        final_state: np.ndarray = sol[-1, :]

        # Check for NaNs or infinities
        if np.any(np.isnan(final_state)) or np.any(np.isinf(final_state)):
            print(f"[WARN] NaN/Inf in steady state at f_ras={f_ras}. Using default fallback.")
            return np.array([0.3, 0.3, 0.2, 0.2, 0.2, 0.3], dtype=float)

        # Clip final state into biologically reasonable range
        return np.clip(final_state, 0.0, 1.0)

    except Exception as exc:
        # Catch any integration failure and report it clearly
        print(f"[ERROR] ODE integration failed at f_ras={f_ras}: {exc}")
        # Return a safe fallback state
        return np.array([0.3, 0.3, 0.2, 0.2, 0.2, 0.3], dtype=float)
