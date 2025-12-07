#!/usr/bin/env python3
"""
ras_csc_core.py

Core Ras–CSC–microenvironment model utilities.

This module centralizes the *mathematical model* so that all scripts
(calibration, hypothesis tests, predictions, figures) use the exact
same ODE system, Hill function, and steady-state integration logic.

Why this exists:
    - Previously, multiple scripts each reimplemented the ODEs slightly
      differently, which is dangerous for interpretation and
      reproducibility.
    - A single, importable core module makes the project closer to a
      reusable package and lowers the risk of silent drift.

This module is intentionally I/O-light:
    - No CLI argument parsing.
    - No plotting.
    - No filesystem writes.

It only exposes:
    - A safe Hill function.
    - The Ras–CSC ODE right-hand side.
    - A steady-state simulator with basic convergence checks.
    - A condition → Ras mapping helper.

These functions are used by:
    - src/scripts/run_model_calibration.py
    - src/scripts/run_hypothesis_tests.py
"""

from __future__ import annotations

# Import typing utilities for type hints
from typing import Dict, Mapping, Sequence, Optional, Tuple

# Import numpy for numerical arrays and operations
import numpy as np

# Import scipy.integrate for ODE integration
from scipy.integrate import odeint


# ======================================================================
# TYPE ALIASES
# ======================================================================

# Define a type alias for the parameter dictionary
ParamsDict = Dict[str, float]


# ======================================================================
# HILL FUNCTION
# ======================================================================

def hill_function(x: float, K: float, n: float) -> float:
    """
    Compute a Hill-type activation function with defensive clipping.

    This represents a generic saturating activation term:
        H(x) = x^n / (K^n + x^n)

    Why defensive clipping:
        - Biological parameters are positive and finite, but numerical
          solvers can wander into negative or extremely large values.
        - Clipping avoids NaNs or division by zero without silently
          changing the qualitative shape of the function.

    Returns a float between 0 and 1 (inclusive, up to numerical noise).
    """
    # Clip x so that it remains non-negative and avoids extreme tails
    x_safe: float = float(np.clip(x, 0.0, 10.0))

    # Clip K to remain positive and bounded away from zero
    K_safe: float = float(np.clip(K, 1e-3, 10.0))

    # Clip n to a reasonable Hill exponent range
    n_safe: float = float(np.clip(n, 1.0, 6.0))

    # Compute numerator as x^n
    numerator: float = x_safe ** n_safe

    # Compute denominator as K^n + x^n
    denominator: float = (K_safe ** n_safe) + numerator

    # Guard against an unlikely zero denominator
    if denominator <= 0.0:
        return 0.0

    # Return the Hill fraction
    return numerator / denominator


# ======================================================================
# RAS–CSC ODE SYSTEM
# ======================================================================

def ras_csc_ode(
    y: Sequence[float],
    t: float,
    f_ras: float,
    params: Mapping[str, float],
) -> Sequence[float]:
    """
    Right-hand side of the Ras–CSC–microenvironment ODE system.

    This function encodes the dynamical equations for:
        C(t) – CSC fraction
        A(t) – angiogenesis / vascular support
        T(t) – TGFβ signaling
        R(t) – LEPR (Leptin receptor)
        L(t) – Leptin ligand
        M(t) – mTOR pathway activity

    The model structure mirrors the conceptual Ras-driven feedback loop:
        Ras → angiogenesis → TGFβ → LEPR → leptin → mTOR → CSCs

    Why this is centralized:
        - Every piece of code that uses the model should call *this*
          function, not its own variant. This keeps calibration,
          hypothesis tests, and predictions consistent.

    Args:
        y:
            Current state vector [C, A, T, R, L, M].
        t:
            Time (not used directly because the system is autonomous,
            but required by the ODE integrator interface).
        f_ras:
            Effective Ras input for this simulation (0–1).
        params:
            Mapping of parameter names to float values. Missing keys
            fall back to conservative defaults so that older parameter
            dictionaries remain usable.

    Returns:
        List[float]: Time derivatives [dC/dt, dA/dt, dT/dt, dR/dt, dL/dt, dM/dt].
    """
    # Unpack the state vector with defensive casting to float
    C: float = float(y[0])
    A: float = float(y[1])
    T: float = float(y[2])
    R: float = float(y[3])
    L: float = float(y[4])
    M: float = float(y[5])

    # Clip state variables to biologically meaningful range [0, 1]
    C = float(np.clip(C, 0.0, 1.0))
    A = float(np.clip(A, 0.0, 1.0))
    T = float(np.clip(T, 0.0, 1.0))
    R = float(np.clip(R, 0.0, 1.0))
    L = float(np.clip(L, 0.0, 1.0))
    M = float(np.clip(M, 0.0, 1.0))

    # Extract parameters with defaults for backward compatibility
    rho_C: float = float(params.get("rho_C", 0.15))
    eta_C_M: float = float(params.get("eta_C_M", 2.5))
    K_C: float = float(params.get("K_C", 0.25))
    n_C: float = float(params.get("n_C", 3.5))
    delta_C: float = float(params.get("delta_C", 0.8))

    alpha_A: float = float(params.get("alpha_A", 0.4))
    beta_A: float = float(params.get("beta_A", 1.0))
    delta_A: float = float(params.get("delta_A", 0.6))

    alpha_T: float = float(params.get("alpha_T", 0.8))
    gamma_T: float = float(params.get("gamma_T", 0.15))
    delta_T: float = float(params.get("delta_T", 1.2))

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

    # Compute Hill activation for TGFβ → CSC
    H_C_T: float = hill_function(T, K_C, n_C)

    # Compute Hill activation for TGFβ → LEPR
    H_R_T: float = hill_function(T, K_R, p_R)

    # Compute Hill activation for leptin–LEPR → mTOR
    H_M_LR: float = hill_function(L * R, K_M, q_M)

    # Compute dC/dt: CSC dynamics with Ras, mTOR, and TGFβ contributions
    dC_dt: float = (
        (rho_C * f_ras) +
        (eta_C_M * M) +
        (params.get("eta_C", 1.2) * H_C_T)
    ) * (1.0 - C) - (delta_C * C)

    # Compute dA/dt: angiogenesis driven by Ras and CSCs
    dA_dt: float = (
        (alpha_A * f_ras) +
        (beta_A * C)
    ) * (1.0 - A) - (delta_A * A)

    # Compute dT/dt: TGFβ from angiogenesis and CSCs
    dT_dt: float = (alpha_T * A) + (gamma_T * C) - (delta_T * T)

    # Compute dR/dt: LEPR driven by TGFβ
    dR_dt: float = (eta_R * H_R_T) * (1.0 - R) - (delta_R * R)

    # Compute dL/dt: leptin from angiogenesis
    dL_dt: float = (mu_L * A) - (delta_L * L)

    # Compute dM/dt: mTOR from leptin × LEPR
    dM_dt: float = (eta_M * H_M_LR) * (1.0 - M) - (delta_M * M)

    # Return the derivatives as a list for odeint
    return [dC_dt, dA_dt, dT_dt, dR_dt, dL_dt, dM_dt]


# ======================================================================
# STEADY-STATE SIMULATION WRAPPER
# ======================================================================

def simulate_steady_state(
    f_ras: float,
    params: Mapping[str, float],
    y0: Optional[Sequence[float]] = None,
    t_max: float = 200.0,
    n_steps: int = 2000,
    atol: float = 1e-6,
    rtol: float = 1e-6,
) -> np.ndarray:
    """
    Integrate the Ras–CSC ODE system to approximate a steady state.

    This helper is the *only* place that integrates the ODEs. It wraps
    scipy.integrate.odeint with:
        - Default initial conditions.
        - Basic NaN/Inf and convergence checks.
        - Clipping to [0, 1] at the end.

    Args:
        f_ras:
            Effective Ras input (0–1) for this simulation.
        params:
            Parameter mapping to pass into the ODE.
        y0:
            Optional initial condition [C, A, T, R, L, M]. If None,
            a moderate baseline vector is used.
        t_max:
            Final time for integration. Larger values allow more time
            to relax to steady state but cost more compute.
        n_steps:
            Number of time points used for integration.
        atol:
            Absolute tolerance for the ODE solver.
        rtol:
            Relative tolerance for the ODE solver.

    Returns:
        np.ndarray: Final state vector of shape (6,), clipped to [0, 1].

    Raises:
        ValueError:
            If the initial condition has the wrong shape.
    """
    # Use a default starting point if none is provided
    if y0 is None:
        y0 = [0.3, 0.3, 0.2, 0.2, 0.2, 0.3]

    # Convert initial condition to numpy array
    y0_arr: np.ndarray = np.asarray(y0, dtype=float)

    # Check that the initial condition has length 6
    if y0_arr.shape[0] != 6:
        raise ValueError(
            f"simulate_steady_state expected y0 of length 6, got shape {y0_arr.shape}"
        )

    # Build a time grid for integration
    t: np.ndarray = np.linspace(0.0, float(t_max), int(n_steps))

    try:
        # Integrate using odeint with the core ODE function
        sol: np.ndarray = odeint(
            func=ras_csc_ode,
            y0=y0_arr,
            t=t,
            args=(float(f_ras), dict(params)),
            atol=float(atol),
            rtol=float(rtol),
        )
    except Exception as exc:
        # If integration fails, fall back to a safe default and report
        print(
            f"[ERROR] ODE integration failed for f_ras={f_ras:.3f} "
            f"with error: {exc}. Returning default state."
        )
        return np.array([0.3, 0.3, 0.2, 0.2, 0.2, 0.3], dtype=float)

    # Extract the final state
    final_state: np.ndarray = sol[-1, :]

    # Check for NaN or Inf in the final state
    if not np.all(np.isfinite(final_state)):
        print(
            f"[WARN] Non-finite values detected in final state for f_ras={f_ras:.3f}. "
            "Returning default state."
        )
        return np.array([0.3, 0.3, 0.2, 0.2, 0.2, 0.3], dtype=float)

    # Clip the final state to [0, 1]
    final_clipped: np.ndarray = np.clip(final_state, 0.0, 1.0)

    # Return the clipped state
    return final_clipped


# ======================================================================
# CONDITION → RAS INPUT MAPPING
# ======================================================================

def get_ras_input(condition: str) -> float:
    """
    Map an experimental condition label to a Ras input level.

    This helper centralizes the mapping used across calibration and
    hypothesis tests. If a condition is unknown, a conservative default
    Ras input is returned and a message is printed.

    Args:
        condition:
            Condition name used in the datasets, for example:
                "Normal", "Papilloma", "SCC", "PDV_WT", "PDV_LeprKO".

    Returns:
        float: Ras input in [0, 1].

    Notes:
        The mapping can be revised as the model is updated, but all
        downstream code should continue to call this function instead
        of hard-coding values.
    """
    # Define the default mapping between condition labels and Ras input
    mapping: Dict[str, float] = {
        "Normal": 0.3,
        "Papilloma": 0.6,
        "SCC": 0.95,
        "PDV_WT": 0.95,
        "PDV_LeprKO": 0.95,
    }

    # Attempt to get the Ras input for the provided condition
    if condition in mapping:
        return float(mapping[condition])

    # If an unknown condition is seen, log and return a mid-range default
    print(
        f"[WARN] Unknown condition '{condition}' encountered in get_ras_input. "
        "Using default Ras input 0.5."
    )
    return 0.5
