#!/usr/bin/env python3
"""
test_hypothesis.py

Tests the Alternative Hypothesis (H1):
  "Disruption of LEPR or mTOR signalling will drive the system from a
   high-C, high-A steady state toward a low-C, low-A steady state,
   even in the continued presence of oncogenic Ras."

Simulation Procedure:
  Phase 1 (Tumor Growth):  Start Benign (Ras=1). Let system reach Malignant Steady State.
  Phase 2 (Intervention):  Keep Ras=1. CUT the feedback loop (k_C_M = 0 or k_R_prod = 0).
  Check: Does the system collapse back to a Benign-like state?
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.integrate import solve_ivp
from dataclasses import dataclass, asdict

# --------------------------------------------------------------------
# 1. Model Definition
# --------------------------------------------------------------------

@dataclass
class RasCSCParams:
    f_ras: float = 1.0
    k_C_ras: float = 0.1
    k_C_TGFb: float = 0.5
    K_C: float = 0.5
    n_C: float = 4.0
    k_C_M: float = 0.5
    k_C_deg: float = 1.0
    
    k_A_ras: float = 0.1
    k_A_C: float = 0.5
    k_A_deg: float = 1.0
    
    k_T_A: float = 1.0
    k_T_C: float = 0.0
    k_T_deg: float = 1.0
    
    k_R_prod: float = 1.0
    K_R: float = 0.5
    p_R: float = 4.0
    k_R_deg: float = 1.0
    
    L_sys: float = 0.1
    k_L_A: float = 1.0
    
    k_M_max: float = 1.0
    K_M: float = 0.5
    q_M: float = 4.0
    
    clip_state: bool = True

def ras_csc_rhs(t, y, params):
    C, A, T, R = y
    if params.clip_state:
        C, A, T, R = np.maximum([C, A, T, R], 0.0)
        
    # Algebraic
    L = params.L_sys + params.k_L_A * A
    S = L * R
    denom_M = (S ** params.q_M) + (params.K_M ** params.q_M)
    M = params.k_M_max * (S ** params.q_M) / denom_M if denom_M > 1e-12 else 0.0
    
    # ODEs
    denom_C = (T ** params.n_C) + (params.K_C ** params.n_C)
    term_TGFb = params.k_C_TGFb * (T ** params.n_C) / denom_C if denom_C > 1e-12 else 0.0
    
    # dC/dt: Ras + TGFb + mTOR_Feedback
    activation = params.k_C_ras * params.f_ras + term_TGFb + params.k_C_M * M
    dC = activation * (1 - C) - params.k_C_deg * C
    
    # dA/dt
    dA = params.k_A_ras * params.f_ras + params.k_A_C * C - params.k_A_deg * A
    
    # dT/dt
    dT = params.k_T_A * A + params.k_T_C * C - params.k_T_deg * T
    
    # dR/dt
    denom_R = (T ** params.p_R) + (params.K_R ** params.p_R)
    term_R = params.k_R_prod * (T ** params.p_R) / denom_R if denom_R > 1e-12 else 0.0
    dR = term_R - params.k_R_deg * R
    
    return [dC, dA, dT, dR]

# --------------------------------------------------------------------
# 2. Simulation Logic
# --------------------------------------------------------------------

def run_intervention_test(param_dict):
    """
    Simulates:
      1. Growth: Start from small initial state with Ras ON.
      2. Intervention: At t=100, cut feedback (simulating LEPR/mTOR disruption).
    """
    # --- Phase 1: Tumor Establishment ---
    params = RasCSCParams(**param_dict)
    params.f_ras = 1.0  # Oncogenic Ras ON
    
    y0 = [0.01, 0.01, 0.01, 0.01] # Initial benign state
    t_span1 = (0, 100)
    
    sol1 = solve_ivp(lambda t, y: ras_csc_rhs(t, y, params), t_span1, y0, method='RK45', rtol=1e-6, dense_output=True)
    t1 = np.linspace(0, 100, 200)
    y1 = sol1.sol(t1).T
    
    # Check steady state at t=100 (Pre-treatment)
    csc_high = y1[-1, 0]
    angio_high = y1[-1, 1]
    
    # --- Phase 2: Intervention (Disrupt Feedback) ---
    # Scenario: Knockout LEPR (reduce R prod) OR Block mTOR (reduce k_C_M)
    # We will simulate blocking mTOR feedback directly (k_C_M -> 0)
    params_treated = RasCSCParams(**param_dict)
    params_treated.f_ras = 1.0 # Ras stays ON
    params_treated.k_C_M = 0.0 # Feedback CUT (mTOR inhibition)
    
    # Start from established tumor state
    y_initial_2 = y1[-1]
    t_span2 = (100, 200)
    
    sol2 = solve_ivp(lambda t, y: ras_csc_rhs(t, y, params_treated), t_span2, y_initial_2, method='RK45', rtol=1e-6, dense_output=True)
    t2 = np.linspace(100, 200, 200)
    y2 = sol2.sol(t2).T
    
    # Check final state (Post-treatment)
    csc_low = y2[-1, 0]
    angio_low = y2[-1, 1]
    
    # Combine for plotting
    t_full = np.concatenate([t1, t2])
    y_full = np.concatenate([y1, y2])
    
    return t_full, y_full, (csc_high, angio_high), (csc_low, angio_low)

# --------------------------------------------------------------------
# 3. Main Execution
# --------------------------------------------------------------------

PROJECT_ROOT = Path(".").resolve()
PARAMS_PATH = PROJECT_ROOT / "data" / "processed" / "model_fits" / "model_calibration_results.json"

if not PARAMS_PATH.exists():
    print(f"[ERROR] Params file not found: {PARAMS_PATH}")
    exit(1)

with open(PARAMS_PATH) as f:
    data = json.load(f)
    best_params = data["model_A_params"]
    if "description" in best_params: del best_params["description"]

print("[INFO] Loaded calibrated parameters.")

# Run Test
t, y, pre, post = run_intervention_test(best_params)

# --- Evaluation ---
# H1 requires: High State -> Low State transition
# Define "Low" roughly as < 50% of the High state, or an absolute threshold.
csc_drop = (pre[0] - post[0]) / pre[0] * 100
angio_drop = (pre[1] - post[1]) / pre[1] * 100

print(f"\n--- Hypothesis Test Results ---")
print(f"Condition: Ras stays ON (f_ras=1.0) throughout.")
print(f"Intervention: At t=100, Feedback is blocked (k_C_M -> 0).")
print(f"\n[Pre-Treatment Steady State]")
print(f"  CSC Fraction: {pre[0]:.4f}")
print(f"  Angiogenesis: {pre[1]:.4f}")
print(f"\n[Post-Treatment Steady State]")
print(f"  CSC Fraction: {post[0]:.4f}")
print(f"  Angiogenesis: {post[1]:.4f}")
print(f"\n[Effect Size]")
print(f"  CSC Reduction: {csc_drop:.1f}%")
print(f"  Angio Reduction: {angio_drop:.1f}%")

if csc_drop > 50 and angio_drop > 50:
    print("\n>>> CONCLUSION: Data SUPPORTS Alternative Hypothesis (H1).")
    print("    Disrupting feedback collapsed the malignant state despite Oncogenic Ras.")
else:
    print("\n>>> CONCLUSION: Data SUPPORTS Null Hypothesis (H0).")
    print("    Tumor persisted or did not collapse significantly.")

# --- Plotting ---
plt.figure(figsize=(10, 6))
plt.plot(t, y[:, 0], label='CSC Fraction (C)', linewidth=2.5, color='firebrick')
plt.plot(t, y[:, 1], label='Angiogenesis (A)', linewidth=2.5, color='darkblue')
plt.plot(t, y[:, 3], label='Lepr (R)', linewidth=1.5, color='green', alpha=0.6)

# Add Intervention Line
plt.axvline(x=100, color='black', linestyle='--', label='Feedback Disruption (t=100)')
plt.text(105, 0.9, 'Ras ON + No Feedback', fontsize=10, fontweight='bold')
plt.text(5, 0.9, 'Ras ON + Feedback Intact', fontsize=10, fontweight='bold')

plt.title('Hypothesis Test: Does Feedback Disruption Collapse the Tumor?')
plt.xlabel('Time (arbitrary units)')
plt.ylabel('Normalized Activity')
plt.legend()
plt.grid(True, alpha=0.3)

out_path = PROJECT_ROOT / "figures" / "hypothesis_test_result.png"
plt.savefig(out_path, dpi=300)
print(f"\n[INFO] Plot saved to {out_path}")