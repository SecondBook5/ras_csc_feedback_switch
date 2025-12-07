#!/usr/bin/env python3
"""
run_step8_system_analysis.py

Step 8: Analyze the system (parameter sensitivity, static and temporal behaviors)
Per standard systems biology workflow.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
from scipy.integrate import odeint
from scipy.optimize import fsolve
import json
from pathlib import Path
from typing import Dict, Tuple, List

sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

ROOT = Path(".")
RESULTS_DIR = ROOT / "results"
FIG_DIR = ROOT / "figures/main"

RESULTS_DIR.mkdir(exist_ok=True, parents=True)
FIG_DIR.mkdir(exist_ok=True, parents=True)


# ============================================================================
# ODE MODEL
# ============================================================================

def ras_csc_model(y: np.ndarray, t: float, f_ras: float, params: Dict) -> List[float]:
    """
    6-variable Ras-CSC ODE system.
    
    Args:
        y: State vector [C, A, T, R, L, M]
        t: Time
        f_ras: Ras activity input
        params: Parameter dictionary
        
    Returns:
        List of derivatives [dC/dt, dA/dt, dT/dt, dR/dt, dL/dt, dM/dt]
    """
    # Extract state variables
    C, A, T, R, L, M = y

    # Clip to valid range to prevent numerical issues
    C = np.clip(C, 0, 1)
    A = np.clip(A, 0, 1)
    T = np.clip(T, 0, 1)
    R = np.clip(R, 0, 1)
    L = np.clip(L, 0, 1)
    M = np.clip(M, 0, 1)

    def hill(x: float, K: float, n: float) -> float:
        """
        Hill function for cooperative binding.
        
        Args:
            x: Input concentration
            K: Half-maximal concentration
            n: Hill coefficient (cooperativity)
            
        Returns:
            Normalized Hill response
        """
        # Prevent numerical overflow/underflow
        x = np.clip(x, 0, 10)
        K = np.clip(K, 0.01, 10)
        n = np.clip(n, 1, 6)
        denom = K**n + x**n
        if denom == 0:
            return 0
        return (x**n) / denom

    # CSC dynamics: Ras seeding + mTOR feedback + TGFβ activation
    dC_dt = (params['rho_C'] * f_ras +
             params['eta_C_M'] * M +
             params['eta_C'] * hill(T, params['K_C'], params['n_C'])) * (1 - C) - params['delta_C'] * C

    # Angiogenesis: Ras activation + CSC promotion
    dA_dt = (params['alpha_A'] * f_ras + params['beta_A']
             * C) * (1 - A) - params['delta_A'] * A

    # TGFβ: Angiogenesis promotion + CSC secretion
    dT_dt = params['alpha_T'] * A + \
        params['gamma_T'] * C - params['delta_T'] * T

    # LEPR expression: TGFβ induction
    dR_dt = params['eta_R'] * hill(T, params['K_R'],
                                   params['p_R']) * (1 - R) - params['delta_R'] * R

    # Leptin import: Angiogenesis-mediated delivery
    dL_dt = params['mu_L'] * A - params['delta_L'] * L

    # mTOR activation: Leptin-LEPR signaling
    dM_dt = params['eta_M'] * hill(L * R, params['K_M'],
                                   params['q_M']) * (1 - M) - params['delta_M'] * M

    return [dC_dt, dA_dt, dT_dt, dR_dt, dL_dt, dM_dt]


def simulate_steady_state(f_ras: float, params: Dict, y0: np.ndarray = None,
                          t_max: float = 200) -> np.ndarray:
    """
    Simulate to steady state.
    
    Args:
        f_ras: Ras activity input
        params: Parameter dictionary
        y0: Initial condition (default: [0.3, 0.3, 0.2, 0.2, 0.2, 0.3])
        t_max: Maximum simulation time
        
    Returns:
        Steady state vector [C, A, T, R, L, M]
    """
    if y0 is None:
        y0 = [0.3, 0.3, 0.2, 0.2, 0.2, 0.3]

    try:
        # Simulate to steady state
        t = np.linspace(0, t_max, 2000)
        sol = odeint(ras_csc_model, y0, t, args=(f_ras, params),
                     atol=1e-8, rtol=1e-8)
        final = sol[-1, :]

        # Check for NaN/Inf
        if np.any(np.isnan(final)) or np.any(np.isinf(final)):
            print(f"  [WARN] NaN/Inf at f_ras={f_ras}, using defaults")
            return np.array([0.3, 0.3, 0.2, 0.2, 0.2, 0.3])

        return np.clip(final, 0, 1)

    except Exception as e:
        print(f"  [ERROR] Integration failed: {e}")
        return np.array([0.3, 0.3, 0.2, 0.2, 0.2, 0.3])


# ============================================================================
# SENSITIVITY ANALYSIS
# ============================================================================

def local_sensitivity_analysis(params: Dict, f_ras: float = 0.95) -> pd.DataFrame:
    """
    Local parameter sensitivity analysis.
    
    Perturb each parameter by ±10% and measure effect on steady-state CSC fraction.
    
    Args:
        params: Nominal parameter set
        f_ras: Ras input level for analysis (default: malignant state)
        
    Returns:
        DataFrame with sensitivity coefficients
    """
    print("\n[SENSITIVITY] Computing local parameter sensitivities...")

    # Get baseline steady state
    ss_baseline = simulate_steady_state(f_ras, params)
    C_baseline = ss_baseline[0]
    M_baseline = ss_baseline[5]

    # Parameters to analyze (free parameters from calibration)
    param_names = [
        'rho_C', 'eta_C_M', 'K_C', 'n_C', 'delta_C',
        'beta_A', 'delta_A', 'alpha_T', 'delta_T',
        'eta_M', 'delta_M'
    ]

    results = []

    for param_name in param_names:
        # Store original value
        original_value = params[param_name]

        # Perturb +10%
        params_plus = params.copy()
        params_plus[param_name] = original_value * 1.1
        ss_plus = simulate_steady_state(f_ras, params_plus)
        C_plus = ss_plus[0]
        M_plus = ss_plus[5]

        # Perturb -10%
        params_minus = params.copy()
        params_minus[param_name] = original_value * 0.9
        ss_minus = simulate_steady_state(f_ras, params_minus)
        C_minus = ss_minus[0]
        M_minus = ss_minus[5]

        # Compute normalized sensitivity: (dC/C) / (dp/p)
        # Using central difference approximation
        if C_baseline > 1e-6:
            sens_C = ((C_plus - C_minus) / C_baseline) / \
                0.2  # 0.2 = 20% total change
        else:
            sens_C = 0

        if M_baseline > 1e-6:
            sens_M = ((M_plus - M_minus) / M_baseline) / 0.2
        else:
            sens_M = 0

        results.append({
            'parameter': param_name,
            'nominal_value': original_value,
            'C_baseline': C_baseline,
            'C_plus10': C_plus,
            'C_minus10': C_minus,
            'sensitivity_C': sens_C,
            'M_baseline': M_baseline,
            'M_plus10': M_plus,
            'M_minus10': M_minus,
            'sensitivity_M': sens_M,
            'abs_sensitivity': abs(sens_C)
        })

    df = pd.DataFrame(results)
    df = df.sort_values('abs_sensitivity', ascending=False)

    print(f"  [COMPLETE] Analyzed {len(param_names)} parameters")
    print(f"\n  Top 5 sensitive parameters (CSC):")
    for _, row in df.head(5).iterrows():
        print(f"    {row['parameter']:10s}: S_C = {row['sensitivity_C']:7.3f}")

    return df


# ============================================================================
# TEMPORAL DYNAMICS
# ============================================================================

def temporal_dynamics_analysis(params: Dict) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Analyze temporal dynamics: step responses and approach to steady state.
    
    Tests:
    1. Step-up response (Normal → SCC)
    2. Step-down response (SCC → Normal)
    3. Oscillation check
    
    Args:
        params: Parameter dictionary
        
    Returns:
        Tuple of (step_up_df, step_down_df)
    """
    print("\n[TEMPORAL] Analyzing temporal dynamics...")

    # Time vector
    t = np.linspace(0, 200, 1000)

    # Test 1: Step-up response (Ras 0.3 → 0.9)
    print("  [1/2] Step-up response (Ras: 0.3 → 0.9)...")
    y0_benign = simulate_steady_state(0.3, params)

    # Create step function
    f_ras_step_up = np.where(t < 50, 0.3, 0.9)

    # Integrate with time-varying Ras
    sol_up = np.zeros((len(t), 6))
    sol_up[0, :] = y0_benign

    for i in range(1, len(t)):
        dt = t[i] - t[i-1]
        t_span = [0, dt]
        sol_temp = odeint(ras_csc_model, sol_up[i-1, :], t_span,
                          args=(f_ras_step_up[i], params))
        sol_up[i, :] = sol_temp[-1, :]

    df_up = pd.DataFrame({
        'time': t,
        'f_ras': f_ras_step_up,
        'C': sol_up[:, 0],
        'A': sol_up[:, 1],
        'T': sol_up[:, 2],
        'R': sol_up[:, 3],
        'L': sol_up[:, 4],
        'M': sol_up[:, 5]
    })

    # Test 2: Step-down response (Ras 0.9 → 0.3)
    print("  [2/2] Step-down response (Ras: 0.9 → 0.3)...")
    y0_malignant = simulate_steady_state(0.9, params)

    f_ras_step_down = np.where(t < 50, 0.9, 0.3)

    sol_down = np.zeros((len(t), 6))
    sol_down[0, :] = y0_malignant

    for i in range(1, len(t)):
        dt = t[i] - t[i-1]
        t_span = [0, dt]
        sol_temp = odeint(ras_csc_model, sol_down[i-1, :], t_span,
                          args=(f_ras_step_down[i], params))
        sol_down[i, :] = sol_temp[-1, :]

    df_down = pd.DataFrame({
        'time': t,
        'f_ras': f_ras_step_down,
        'C': sol_down[:, 0],
        'A': sol_down[:, 1],
        'T': sol_down[:, 2],
        'R': sol_down[:, 3],
        'L': sol_down[:, 4],
        'M': sol_down[:, 5]
    })

    # Check for oscillations
    C_final_up = df_up['C'].iloc[-50:].std()
    C_final_down = df_down['C'].iloc[-50:].std()

    if C_final_up < 0.01 and C_final_down < 0.01:
        print("  ✓ No sustained oscillations detected (steady-state convergence)")
    else:
        print(
            f"  ⚠ Possible oscillations: std(C_up)={C_final_up:.4f}, std(C_down)={C_final_down:.4f}")

    return df_up, df_down


# ============================================================================
# PHASE PLANE ANALYSIS
# ============================================================================

def compute_nullclines(params: Dict, f_ras: float = 0.9) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute C and M nullclines.
    
    Args:
        params: Parameter dictionary
        f_ras: Ras input level
        
    Returns:
        Tuple of (C_range, M_nullcline) where M_nullcline is M value that makes dM/dt=0
    """
    print("\n[PHASE PLANE] Computing nullclines...")

    C_range = np.linspace(0, 1, 100)
    M_nullcline = np.zeros_like(C_range)

    for i, C_val in enumerate(C_range):
        # For each C, find M where dM/dt = 0
        # Need to solve through full system since M depends on L, R which depend on A, T

        # Use full steady state with fixed C as approximation
        def equations(vars):
            A, T, R, L, M = vars

            # Clip values
            A = np.clip(A, 0, 1)
            T = np.clip(T, 0, 1)
            R = np.clip(R, 0, 1)
            L = np.clip(L, 0, 1)
            M = np.clip(M, 0, 1)

            def hill(x, K, n):
                x = np.clip(x, 0, 10)
                K = np.clip(K, 0.01, 10)
                n = np.clip(n, 1, 6)
                denom = K**n + x**n
                return (x**n) / denom if denom > 0 else 0

            # Equations (with C fixed)
            dA = (params['alpha_A'] * f_ras + params['beta_A']
                  * C_val) * (1 - A) - params['delta_A'] * A
            dT = params['alpha_T'] * A + params['gamma_T'] * \
                C_val - params['delta_T'] * T
            dR = params['eta_R'] * hill(T, params['K_R'],
                                        params['p_R']) * (1 - R) - params['delta_R'] * R
            dL = params['mu_L'] * A - params['delta_L'] * L
            dM = params['eta_M'] * hill(L * R, params['K_M'],
                                        params['q_M']) * (1 - M) - params['delta_M'] * M

            return [dA, dT, dR, dL, dM]

        try:
            # Solve for steady state of other variables
            sol = fsolve(equations, [0.5, 0.5, 0.5, 0.5, 0.5])
            M_nullcline[i] = np.clip(sol[4], 0, 1)
        except:
            M_nullcline[i] = 0.5

    return C_range, M_nullcline


def phase_plane_analysis(params: Dict, f_ras: float = 0.9) -> pd.DataFrame:
    """
    Generate phase plane data: vector field and trajectories.
    
    Args:
        params: Parameter dictionary
        f_ras: Ras input level
        
    Returns:
        DataFrame with phase plane data
    """
    print("\n[PHASE PLANE] Computing vector field...")

    # Create grid
    C_grid = np.linspace(0, 1, 20)
    M_grid = np.linspace(0, 1, 20)

    C_mesh, M_mesh = np.meshgrid(C_grid, M_grid)

    dC_mesh = np.zeros_like(C_mesh)
    dM_mesh = np.zeros_like(M_mesh)

    for i in range(len(M_grid)):
        for j in range(len(C_grid)):
            C_val = C_mesh[i, j]
            M_val = M_mesh[i, j]

            # Compute derivatives at this point
            # Need full system state
            y = simulate_steady_state(
                f_ras, params, y0=[C_val, 0.5, 0.5, 0.5, 0.5, M_val], t_max=50)
            dydt = ras_csc_model(y, 0, f_ras, params)

            dC_mesh[i, j] = dydt[0]
            dM_mesh[i, j] = dydt[5]

    # Store as DataFrame
    data = []
    for i in range(len(M_grid)):
        for j in range(len(C_grid)):
            data.append({
                'C': C_mesh[i, j],
                'M': M_mesh[i, j],
                'dC': dC_mesh[i, j],
                'dM': dM_mesh[i, j]
            })

    return pd.DataFrame(data)


# ============================================================================
# ENERGY LANDSCAPE
# ============================================================================

def compute_energy_landscape(params: Dict, f_ras_range: np.ndarray) -> pd.DataFrame:
    """
    Compute quasi-potential landscape for bistability visualization.
    
    Uses the approximation: U(C) ≈ -∫ dC/dt dC
    
    Args:
        params: Parameter dictionary
        f_ras_range: Range of Ras values to scan
        
    Returns:
        DataFrame with energy landscape
    """
    print("\n[LANDSCAPE] Computing energy landscape...")

    C_range = np.linspace(0, 1, 100)

    data = []

    for f_ras in f_ras_range:
        potential = np.zeros_like(C_range)

        for i, C_val in enumerate(C_range):
            # Get full steady state with this C value (approximately)
            y = simulate_steady_state(
                f_ras, params, y0=[C_val, 0.5, 0.5, 0.5, 0.5, 0.5], t_max=50)
            dydt = ras_csc_model(y, 0, f_ras, params)

            # Potential is negative integral of dC/dt
            # Approximate as cumulative sum
            if i > 0:
                dC = C_range[i] - C_range[i-1]
                potential[i] = potential[i-1] - dydt[0] * dC
            else:
                potential[i] = 0

        # Normalize potential
        potential = potential - potential.min()

        for C_val, U_val in zip(C_range, potential):
            data.append({
                'f_ras': f_ras,
                'C': C_val,
                'potential': U_val
            })

    return pd.DataFrame(data)


# ============================================================================
# MAIN ANALYSIS
# ============================================================================

def main():
    """Run complete Step 8 system analysis."""
    print("\n" + "="*70)
    print("STEP 8: SYSTEM ANALYSIS")
    print("Parameter Sensitivity | Temporal Dynamics | Phase Plane | Landscape")
    print("="*70)

    # Load parameters
    with open(RESULTS_DIR / "optimized_parameters_CORRECTED.json", 'r') as f:
        params = json.load(f)

    print("\n[INFO] Loaded optimized parameters")

    # 1. Sensitivity Analysis
    df_sens = local_sensitivity_analysis(params, f_ras=0.95)
    sens_file = RESULTS_DIR / "sensitivity_analysis.csv"
    df_sens.to_csv(sens_file, index=False)
    print(f"[SAVED] {sens_file}")

    # 2. Temporal Dynamics
    df_up, df_down = temporal_dynamics_analysis(params)
    up_file = RESULTS_DIR / "temporal_step_up.csv"
    down_file = RESULTS_DIR / "temporal_step_down.csv"
    df_up.to_csv(up_file, index=False)
    df_down.to_csv(down_file, index=False)
    print(f"[SAVED] {up_file}")
    print(f"[SAVED] {down_file}")

    # 3. Phase Plane
    df_phase = phase_plane_analysis(params, f_ras=0.95)
    phase_file = RESULTS_DIR / "phase_plane_data.csv"
    df_phase.to_csv(phase_file, index=False)
    print(f"[SAVED] {phase_file}")

    # 4. Energy Landscape
    f_ras_range = np.array([0.3, 0.5, 0.7, 0.9])
    df_landscape = compute_energy_landscape(params, f_ras_range)
    landscape_file = RESULTS_DIR / "energy_landscape.csv"
    df_landscape.to_csv(landscape_file, index=False)
    print(f"[SAVED] {landscape_file}")

    print("\n" + "="*70)
    print("STEP 8 ANALYSIS COMPLETE")
    print("="*70)
    print("\nNext: Generate Step 8 figures (Figures 5-8)")


if __name__ == "__main__":
    main()
