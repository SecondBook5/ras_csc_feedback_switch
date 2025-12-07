#!/usr/bin/env python3
"""
Core Ras–CSC–microenvironment feedback ODE model for calibration.

This module defines:
  - RasCSCParams: dataclass of all model parameters
  - ras_csc_rhs:   ODE right-hand side dy/dt for [C, A, T, R]
  - simulate_trajectory: RK4 time integrator

State vector:
  y = [C, A, T, R]
    C = Cancer stem cell fraction / activity
    A = Angiogenesis / vessel density surrogate
    T = TGFβ level
    R = Lepr receptor level

Algebraic quasi-steady-state variables:
  L = leptin (systemic + angiogenesis-driven)
  M = mTOR activity

Equations (schematic):

1) dC/dt = [ k_C_ras * f_ras
             + k_C_TGFb * Hill_T(T)
             + k_C_M * M ] * (1 - C)  -  k_C_deg * C

2) dA/dt = k_A_ras * f_ras + k_A_C * C - k_A_deg * A

3) dT/dt = k_T_A * A + k_T_C * C - k_T_deg * T

4) dR/dt = k_R_prod * Hill_R(T) - k_R_deg * R

5) L = L_sys + k_L_A * A

6) M = k_M_max * Hill_M(S),      S = L * R

Where each Hill(.) is of the form:
  Hill(x) = x^n / (x^n + K^n)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Tuple

import numpy as np


@dataclass
class RasCSCParams:
    """
    Container for Ras–CSC feedback model parameters.

    The fields are grouped by equation but all live in a single dataclass
    so that calibration code can treat this as a plain Python object
    (convertible to/from dict via dataclasses.asdict()).

    Most time-scale / threshold / Hill exponents will typically be fixed.
    Calibration focuses on the gain parameters (k_*) that control the
    strength of Ras, TGFβ, and mTOR feedback.

    Attributes
    ----------
    f_ras:
        Effective Ras input (0–1). Condition specific (Normal vs tumor).
    k_C_ras:
        Baseline Ras-driven CSC induction strength.
    k_C_TGFb:
        Maximal TGFβ-driven CSC induction strength.
    K_C:
        TGFβ half-max for CSC Hill response.
    n_C:
        Hill exponent for TGFβ → CSC.
    k_C_M:
        mTOR feedback gain into CSC equation (C ← M).
    k_C_deg:
        CSC decay / differentiation rate.

    k_A_ras:
        Ras-driven baseline angiogenesis.
    k_A_C:
        CSC-driven angiogenesis gain.
    k_A_deg:
        Angiogenesis regression rate.

    k_T_A:
        Angiogenesis-driven TGFβ production.
    k_T_C:
        CSC-driven TGFβ production (often set to 0).
    k_T_deg:
        TGFβ clearance rate.

    k_R_prod:
        Maximal TGFβ-induced Lepr production.
    K_R:
        TGFβ half-max for Lepr Hill response.
    p_R:
        Hill exponent for TGFβ → Lepr.
    k_R_deg:
        Lepr receptor turnover rate.

    L_sys:
        Systemic leptin baseline (independent of tumor).
    k_L_A:
        Angiogenesis-driven leptin gain.

    k_M_max:
        Maximal mTOR activity.
    K_M:
        Half-max of the L*R signal for mTOR activation.
    q_M:
        Hill exponent for mTOR activation.

    clip_state:
        If True, state variables are clipped to be non-negative
        after each RK4 step.
    """

    # CSC equation
    f_ras: float
    k_C_ras: float
    k_C_TGFb: float
    K_C: float
    n_C: float
    k_C_M: float
    k_C_deg: float

    # Angiogenesis equation
    k_A_ras: float
    k_A_C: float
    k_A_deg: float

    # TGFβ equation
    k_T_A: float
    k_T_C: float
    k_T_deg: float

    # Lepr equation
    k_R_prod: float
    K_R: float
    p_R: float
    k_R_deg: float

    # Leptin algebraic
    L_sys: float
    k_L_A: float

    # mTOR algebraic
    k_M_max: float
    K_M: float
    q_M: float

    # Numerical behavior
    clip_state: bool = True


def ras_csc_rhs(t: float, y: np.ndarray, params: RasCSCParams) -> np.ndarray:
    """
    Compute the right-hand side dy/dt for the Ras–CSC model.

    Parameters
    ----------
    t:
        Time (not used explicitly but required by ODE signature).
    y:
        State vector [C, A, T, R] as a 1D numpy array.
    params:
        RasCSCParams instance containing model parameters.

    Returns
    -------
    dy_dt:
        Numpy array of derivatives [dC/dt, dA/dt, dT/dt, dR/dt].
    """
    # Unpack state
    C = float(y[0])
    A = float(y[1])
    T = float(y[2])
    R = float(y[3])

    # Optionally enforce non-negativity of state variables
    if params.clip_state:
        C = max(C, 0.0)
        A = max(A, 0.0)
        T = max(T, 0.0)
        R = max(R, 0.0)

    # -----------------------------
    # Algebraic variables: L, M
    # -----------------------------

    # Leptin: systemic baseline plus angiogenesis-driven component
    L = params.L_sys + params.k_L_A * A

    # mTOR: Hill activation of signal S = L * R
    S = L * R

    # Compute Hill(S; K_M, q_M)
    S_q = S ** params.q_M
    K_Mq = params.K_M ** params.q_M
    denom_M = S_q + K_Mq

    if denom_M > 1e-12:
        M = params.k_M_max * (S_q / denom_M)
    else:
        M = 0.0

    # -----------------------------
    # Differential equations
    # -----------------------------

    # TGFβ Hill for CSC: Hill(T; K_C, n_C)
    T_n = T ** params.n_C
    K_Cn = params.K_C ** params.n_C
    denom_C = T_n + K_Cn
    if denom_C > 1e-12:
        term_TGFb = params.k_C_TGFb * (T_n / denom_C)
    else:
        term_TGFb = 0.0

    # Total CSC activation: Ras + TGFβ + mTOR feedback
    activation_C = (
        params.k_C_ras * params.f_ras
        + term_TGFb
        + params.k_C_M * M
    )

    # CSC logistic growth with decay
    dC = activation_C * (1.0 - C) - params.k_C_deg * C

    # Angiogenesis: Ras + CSC-driven - regression
    dA = (
        params.k_A_ras * params.f_ras
        + params.k_A_C * C
        - params.k_A_deg * A
    )

    # TGFβ: angiogenesis + optional CSC contribution - clearance
    dT = (
        params.k_T_A * A
        + params.k_T_C * C
        - params.k_T_deg * T
    )

    # Lepr: TGFβ-induced Hill - receptor turnover
    T_p = T ** params.p_R
    K_Rp = params.K_R ** params.p_R
    denom_R = T_p + K_Rp
    if denom_R > 1e-12:
        term_Lepr = params.k_R_prod * (T_p / denom_R)
    else:
        term_Lepr = 0.0

    dR = term_Lepr - params.k_R_deg * R

    # Pack derivatives into numpy array
    return np.array([dC, dA, dT, dR], dtype=float)


def simulate_trajectory(
    rhs: Callable[[float, np.ndarray, RasCSCParams], np.ndarray],
    y0: np.ndarray,
    params: RasCSCParams,
    t_span: Tuple[float, float],
    dt: float = 0.1,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Integrate the Ras–CSC ODE using a simple RK4 scheme.

    Parameters
    ----------
    rhs:
        Function implementing the ODE right-hand side with signature
        f(t, y, params) -> dy/dt.
    y0:
        Initial state vector [C, A, T, R] as 1D numpy array.
    params:
        RasCSCParams instance.
    t_span:
        Tuple (t_start, t_end).
    dt:
        Time step for the fixed-step RK4 solver.

    Returns
    -------
    times:
        1D numpy array of time points.
    traj:
        2D numpy array of shape (n_steps, 4) with the trajectory.
    """
    # Unpack time span
    t_start, t_end = float(t_span[0]), float(t_span[1])

    # Create time grid (inclusive of t_start, exclusive of t_end)
    times = np.arange(t_start, t_end, dt, dtype=float)
    n_steps = times.shape[0]

    # Allocate trajectory array
    y0 = np.asarray(y0, dtype=float)
    if y0.shape != (4,):
        raise ValueError(
            f"y0 must be length-4 array [C, A, T, R], got shape {y0.shape}"
        )

    traj = np.zeros((n_steps, 4), dtype=float)
    traj[0, :] = y0

    # Current state
    y_curr = y0.copy()

    # RK4 loop
    for i in range(n_steps - 1):
        t = times[i]

        k1 = rhs(t, y_curr, params)
        k2 = rhs(t + 0.5 * dt, y_curr + 0.5 * dt * k1, params)
        k3 = rhs(t + 0.5 * dt, y_curr + 0.5 * dt * k2, params)
        k4 = rhs(t + dt, y_curr + dt * k3, params)

        y_next = y_curr + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

        if params.clip_state:
            y_next = np.maximum(y_next, 0.0)

        traj[i + 1, :] = y_next
        y_curr = y_next

    return times, traj


if __name__ == "__main__":
    # Minimal sanity check so you can run:
    #   python src/ras_csc_model_calib.py
    # to verify it doesn't crash.
    params = RasCSCParams(
        # CSC
        f_ras=1.0,
        k_C_ras=0.05,
        k_C_TGFb=0.5,
        K_C=0.5,
        n_C=4.0,
        k_C_M=0.8,
        k_C_deg=1.0,
        # Angio
        k_A_ras=0.05,
        k_A_C=0.8,
        k_A_deg=1.0,
        # TGFb
        k_T_A=1.0,
        k_T_C=0.0,
        k_T_deg=1.0,
        # Lepr
        k_R_prod=1.0,
        K_R=0.5,
        p_R=4.0,
        k_R_deg=1.0,
        # Leptin
        L_sys=0.1,
        k_L_A=1.0,
        # mTOR
        k_M_max=1.0,
        K_M=0.5,
        q_M=4.0,
        # Numerical
        clip_state=True,
    )

    y0_test = np.array([0.01, 0.01, 0.01, 0.01], dtype=float)
    t, y = simulate_trajectory(
        rhs=ras_csc_rhs,
        y0=y0_test,
        params=params,
        t_span=(0.0, 50.0),
        dt=0.1,
    )
    print("Sanity check – final state [C, A, T, R]:")
    print(y[-1, :])
