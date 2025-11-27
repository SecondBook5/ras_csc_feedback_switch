#!/usr/bin/env python3
"""
Ras–CSC–microenvironment feedback ODE model.

1. Topology: C -> A -> T -> R -> (L,R) -> M -> C
2. Dynamics:
   - C, A, T, R are differential variables (slow/medium timescale).
   - L (Leptin) and M (mTOR) are Quasi-Steady State (QSS) algebraic variables.
3. C-Equation:
   - Uses Logistic growth constraint (1 - C).
   - Driven by Ras + TGFb (Hill) + mTOR (Linear feedback).

State Vector y = [C, A, T, R]
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Tuple
import numpy as np
import sys

@dataclass
class RasCSCParams:
    """
    Parameter container matching the coefficients in the PDF equations.
    """
    # --- Equation 1: CSC Dynamics (dC/dt) ---
    # dC/dt = [ rho_C*f_ras + eta_C*Hill(T) + eta_C_M*M ] * (1 - C) - delta_C*C
    
    f_ras: float      # f_RAS: Input (0 or 1)
    k_C_ras: float    # rho_C: Baseline Ras contribution
    k_C_TGFb: float   # eta_C: Max contribution from TGFb switch
    K_C: float        # K_C:   TGFb threshold for CSC induction
    n_C: float        # n:     Hill coeff for TGFb->CSC
    k_C_M: float      # eta_C,M: Strength of mTOR feedback (THE FEEDBACK LOOP)
    k_C_deg: float    # delta_C: Decay/Differentiation rate

    # --- Equation 2: Angiogenesis (dA/dt) ---
    # dA/dt = alpha_A*f_ras + beta_A*C - delta_A*A
    
    k_A_ras: float    # alpha_A: Baseline Ras angiogenesis
    k_A_C: float      # beta_A:  CSC-driven angiogenesis (k_V_C from RNAseq)
    k_A_deg: float    # delta_A: Vessel regression rate

    # --- Equation 3: TGFb (dT/dt) ---
    # dT/dt = alpha_T*A + gamma_T*C - delta_T*T
    
    k_T_A: float      # alpha_T: Angiogenesis-driven TGFb (Stromal recruitment)
    k_T_C: float      # gamma_T: CSC-secreted TGFb (Optional, usually 0)
    k_T_deg: float    # delta_T: TGFb clearance rate

    # --- Equation 4: Lepr Dynamics (dR/dt) ---
    # dR/dt = eta_R * Hill(T) - delta_R*R
    
    k_R_prod: float   # eta_R: Max receptor production (k_R_TGFb from RNAseq)
    K_R: float        # K_R:   TGFb threshold for Lepr induction
    p_R: float        # p:     Hill coeff for TGFb->Lepr
    k_R_deg: float    # delta_R: Receptor turnover

    # --- Equation 5: Leptin QSS (L) ---
    # L = L_sys + (mu_L/delta_L) * A
    
    L_sys: float      # L_sys: Systemic baseline leptin
    k_L_A: float      # mu_L/delta_L: Angiogenesis-dependent leptin gain

    # --- Equation 6: mTOR QSS (M) ---
    # M = eta_M * Hill(L*R)
    
    k_M_max: float    # eta_M: Max mTOR activity
    K_M: float        # K_M:   Half-max for Signal S = L*R
    q_M: float        # q:     Hill coeff for mTOR activation

    # Numerical
    clip_state: bool = True


def ras_csc_rhs(t: float, y: np.ndarray, params: RasCSCParams) -> np.ndarray:
    """
    Compute dy/dt for the system [C, A, T, R].
    Matches equations 1-6 in the PDF.
    """
    # 1. Unpack State
    C = y[0]  # Cancer Stem Cells
    A = y[1]  # Angiogenesis
    T = y[2]  # TGFb
    R = y[3]  # Lepr

    if params.clip_state:
        C, A, T, R = np.maximum([C, A, T, R], 0.0)

    # 2. Calculate QSS Variables (Algebraic)
    
    # Eq 5: Leptin (L)
    # Driven by Systemic levels + Angiogenesis
    L = params.L_sys + params.k_L_A * A

    # Eq 6: mTOR (M)
    # Signal S = L * R
    # Hill function activation
    S = L * R
    S_q = S ** params.q_M
    K_Mq = params.K_M ** params.q_M
    
    # Avoid divide by zero
    denom_M = S_q + K_Mq
    if denom_M > 1e-12:
        M = params.k_M_max * (S_q / denom_M)
    else:
        M = 0.0

    # 3. Calculate Differential Equations
    
    # --- dC/dt (Eq 1) ---
    # Terms: Ras + TGFb_Hill + mTOR_Feedback
    
    # TGFb Hill Term
    T_n = T ** params.n_C
    K_Cn = params.K_C ** params.n_C
    denom_C = T_n + K_Cn
    term_TGFb = 0.0
    if denom_C > 1e-12:
        term_TGFb = params.k_C_TGFb * (T_n / denom_C)
        
    # Total Activation Rate
    activation = (params.k_C_ras * params.f_ras) + term_TGFb + (params.k_C_M * M)
    
    # Logistic Growth - Decay
    dC = activation * (1.0 - C) - params.k_C_deg * C

    # --- dA/dt (Eq 2) ---
    dA = (params.k_A_ras * params.f_ras) + (params.k_A_C * C) - (params.k_A_deg * A)

    # --- dT/dt (Eq 3) ---
    dT = (params.k_T_A * A) + (params.k_T_C * C) - (params.k_T_deg * T)

    # --- dR/dt (Eq 4) ---
    # TGFb induces Receptor
    T_p = T ** params.p_R
    K_Rp = params.K_R ** params.p_R
    denom_R = T_p + K_Rp
    
    term_Induction = 0.0
    if denom_R > 1e-12:
        term_Induction = params.k_R_prod * (T_p / denom_R)
        
    dR = term_Induction - (params.k_R_deg * R)

    return np.array([dC, dA, dT, dR])


def simulate_trajectory(rhs, y0, params, t_span, dt=0.1):
    """
    Simple RK4 solver.
    """
    t_start, t_end = t_span
    times = np.arange(t_start, t_end, dt)
    n_steps = len(times)
    
    y = np.zeros((n_steps, len(y0)))
    y[0] = y0
    
    curr_y = y0.copy()
    
    for i in range(n_steps - 1):
        t = times[i]
        
        k1 = rhs(t, curr_y, params)
        k2 = rhs(t + 0.5*dt, curr_y + 0.5*dt*k1, params)
        k3 = rhs(t + 0.5*dt, curr_y + 0.5*dt*k2, params)
        k4 = rhs(t + dt, curr_y + dt*k3, params)
        
        curr_y = curr_y + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)
        
        if params.clip_state:
            curr_y = np.maximum(curr_y, 0.0)
            
        y[i+1] = curr_y
        
    return times, y

if __name__ == "__main__":
    print("Model file loaded. Import 'RasCSCParams' and 'ras_csc_rhs' to use.")
    
    # --- Quick Demo with Placeholder Params ---
    # You should overwrite these with your RNA-seq derived values!
    
    params = RasCSCParams(
        f_ras=1.0,
        k_C_ras=0.05, k_C_TGFb=0.5, K_C=0.5, n_C=4.0, k_C_M=0.8, k_C_deg=1.0, # CSC
        k_A_ras=0.1, k_A_C=0.8, k_A_deg=1.0,                                  # Angio
        k_T_A=1.0, k_T_C=0.0, k_T_deg=1.0,                                    # TGFb
        k_R_prod=1.0, K_R=0.5, p_R=4.0, k_R_deg=1.0,                          # Lepr
        L_sys=0.1, k_L_A=1.0,                                                 # Leptin
        k_M_max=1.0, K_M=0.5, q_M=4.0                                         # mTOR
    )
    
    # Initial Condition: Benign Papilloma (Low everything)
    y0 = [0.01, 0.01, 0.01, 0.01] 
    
    t_span = (0, 100)
    t, y = simulate_trajectory(ras_csc_rhs, y0, params, t_span)
    
    print("Final State (C, A, T, R):")
    print(y[-1])