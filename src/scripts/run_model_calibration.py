#!/usr/bin/env python3
"""
run_model_calibration_CORRECTED.py

FIXED: Properly handles z-scored RNA-seq targets
Strategy: Transform model outputs to z-scores BEFORE comparison
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.integrate import odeint
from scipy.optimize import least_squares
import json
from pathlib import Path

sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

ROOT = Path(".")
TARGETS_CSV = ROOT / "data/processed/rnaseq/ras_csc_calibration_targets.csv"
RESULTS_DIR = ROOT / "results"
FIG_DIR = ROOT / "figures/main"

RESULTS_DIR.mkdir(exist_ok=True, parents=True)
FIG_DIR.mkdir(exist_ok=True, parents=True)


# ============================================================================
# LOAD TARGETS
# ============================================================================

def load_targets():
    """Load RNA-seq targets (z-scored)."""
    df = pd.read_csv(TARGETS_CSV)
    print(f"\n[INFO] Loaded targets from: {TARGETS_CSV}")

    targets = {}
    for _, row in df.iterrows():
        key = f"{row['dataset']}_{row['condition']}"
        targets[key] = {
            'C': row['C_target'],
            'A': row['A_target'],
            'T': row['T_target'],
            'M': row['M_target'],
            'condition': row['condition'],
            'dataset': row['dataset']
        }

    return targets


# ============================================================================
# ODE MODEL
# ============================================================================

def ras_csc_model(y, t, f_ras, params):
    """Ras-CSC ODE system."""
    C, A, T, R, L, M = y

    # Clip to valid range
    C = np.clip(C, 0, 1)
    A = np.clip(A, 0, 1)
    T = np.clip(T, 0, 1)
    R = np.clip(R, 0, 1)
    L = np.clip(L, 0, 1)
    M = np.clip(M, 0, 1)

    def hill(x, K, n):
        x = np.clip(x, 0, 10)  # Prevent overflow
        K = np.clip(K, 0.01, 10)
        n = np.clip(n, 1, 6)
        denom = K**n + x**n
        if denom == 0:
            return 0
        return (x**n) / denom

    dC_dt = (params['rho_C'] * f_ras +
             params['eta_C_M'] * M +
             params['eta_C'] * hill(T, params['K_C'], params['n_C'])) * (1 - C) - params['delta_C'] * C

    dA_dt = (params['alpha_A'] * f_ras + params['beta_A']
             * C) * (1 - A) - params['delta_A'] * A

    dT_dt = params['alpha_T'] * A + \
        params['gamma_T'] * C - params['delta_T'] * T

    dR_dt = params['eta_R'] * hill(T, params['K_R'],
                                   params['p_R']) * (1 - R) - params['delta_R'] * R

    dL_dt = params['mu_L'] * A - params['delta_L'] * L

    dM_dt = params['eta_M'] * hill(L * R, params['K_M'],
                                   params['q_M']) * (1 - M) - params['delta_M'] * M

    return [dC_dt, dA_dt, dT_dt, dR_dt, dL_dt, dM_dt]


def simulate_steady_state(f_ras, params, y0=None, t_max=200):
    """Simulate to steady state with error handling."""
    if y0 is None:
        y0 = [0.3, 0.3, 0.2, 0.2, 0.2, 0.3]

    try:
        t = np.linspace(0, t_max, 2000)
        sol = odeint(ras_csc_model, y0, t, args=(f_ras, params),
                     atol=1e-6, rtol=1e-6)
        final = sol[-1, :]

        # Check for NaN/Inf
        if np.any(np.isnan(final)) or np.any(np.isinf(final)):
            print(
                f"  [WARN] NaN/Inf detected at f_ras={f_ras}, using defaults")
            return np.array([0.3, 0.3, 0.2, 0.2, 0.2, 0.3])

        return np.clip(final, 0, 1)

    except Exception as e:
        print(f"  [ERROR] Integration failed: {e}")
        return np.array([0.3, 0.3, 0.2, 0.2, 0.2, 0.3])


# ============================================================================
# CONDITION MAPPING
# ============================================================================

def get_ras_input(condition):
    """Map condition to Ras input."""
    mapping = {
        'Normal': 0.3,
        'Papilloma': 0.6,
        'SCC': 0.95,
        'PDV_WT': 0.95,
        'PDV_LeprKO': 0.95
    }
    return mapping.get(condition, 0.5)


# ============================================================================
# CALIBRATION (WITH Z-SCORE TRANSFORM)
# ============================================================================

def objective_function(param_vec, targets, fixed_params):
    """
    Objective function with PROPER z-score handling.
    
    Strategy:
    1. Simulate model for all conditions → raw outputs
    2. Z-score transform model outputs across conditions
    3. Compare z-scored model to z-scored data
    """
    # Unpack parameters
    params = fixed_params.copy()
    params.update({
        'rho_C': param_vec[0],
        'eta_C_M': param_vec[1],
        'K_C': param_vec[2],
        'n_C': param_vec[3],
        'delta_C': param_vec[4],
        'beta_A': param_vec[5],
        'delta_A': param_vec[6],
        'alpha_T': param_vec[7],
        'delta_T': param_vec[8],
        'eta_M': param_vec[9],
        'delta_M': param_vec[10],
    })

    # Simulate all conditions
    results = {}
    for key, target in targets.items():
        f_ras = get_ras_input(target['condition'])

        # Handle LeprKO specially
        if 'LeprKO' in target['condition']:
            params_ko = params.copy()
            params_ko['mu_L'] = params['mu_L'] * 0.2  # Strong reduction
            ss = simulate_steady_state(f_ras, params_ko)
        else:
            ss = simulate_steady_state(f_ras, params)

        results[key] = {
            'C_raw': ss[0],
            'A_raw': ss[1],
            'T_raw': ss[2],
            'M_raw': ss[5]
        }

    # Z-score transform model outputs
    for module in ['C', 'A', 'T', 'M']:
        raw_key = f'{module}_raw'
        values = [results[k][raw_key] for k in results.keys()]
        mean_val = np.mean(values)
        std_val = np.std(values)

        if std_val < 1e-6:  # Avoid division by zero
            std_val = 1.0

        for key in results.keys():
            z_key = f'{module}_zscore'
            results[key][z_key] = (results[key][raw_key] - mean_val) / std_val

    # Compute residuals (z-scored model - z-scored data)
    residuals = []
    for key, target in targets.items():
        residuals.append(results[key]['C_zscore'] - target['C'])
        residuals.append(results[key]['A_zscore'] - target['A'])
        residuals.append(results[key]['T_zscore'] - target['T'])
        residuals.append(results[key]['M_zscore'] - target['M'])

    return np.array(residuals)


def fit_model(targets):
    """Fit model with proper z-score handling."""
    print("\n" + "="*70)
    print("PARAMETER FITTING (CORRECTED)")
    print("="*70)

    # Fixed parameters
    fixed_params = {
        'eta_C': 1.2,      # TGFβ → CSC (moderate)
        'alpha_A': 0.4,    # Ras → Angio
        'gamma_T': 0.15,   # CSC → TGFβ
        'eta_R': 0.6,      # TGFβ → LEPR
        'K_R': 0.35,
        'p_R': 2.5,
        'delta_R': 0.6,
        'mu_L': 0.9,       # Angio → Leptin
        'delta_L': 0.6,
        'K_M': 0.25,
        'q_M': 3.5,
    }

    # Initial guess (REVISED for stability)
    x0 = np.array([
        0.15,   # rho_C (moderate Ras seeding)
        2.5,    # eta_C_M (strong mTOR feedback)
        0.25,   # K_C
        3.5,    # n_C
        0.8,    # delta_C (moderate decay)
        1.0,    # beta_A (CSC→Angio)
        0.6,    # delta_A
        0.8,    # alpha_T (Angio→TGFβ)
        1.2,    # delta_T
        1.0,    # eta_M
        0.6,    # delta_M
    ])

    # Bounds (tighter to prevent collapse)
    bounds = (
        [0.05, 0.5, 0.1, 2.0, 0.3, 0.3, 0.3, 0.3, 0.5, 0.3, 0.3],  # lower
        [0.5, 5.0, 0.5, 5.0, 2.0, 2.0, 1.5, 2.0, 3.0, 2.0, 1.5]   # upper
    )

    print("\n[INFO] Running optimization with z-score transform...")
    result = least_squares(
        objective_function,
        x0,
        args=(targets, fixed_params),
        bounds=bounds,
        method='trf',
        verbose=2,
        max_nfev=500,  # Limit iterations
        ftol=1e-6,
        xtol=1e-6
    )

    # Extract optimized parameters
    opt_params = fixed_params.copy()
    opt_params.update({
        'rho_C': result.x[0],
        'eta_C_M': result.x[1],
        'K_C': result.x[2],
        'n_C': result.x[3],
        'delta_C': result.x[4],
        'beta_A': result.x[5],
        'delta_A': result.x[6],
        'alpha_T': result.x[7],
        'delta_T': result.x[8],
        'eta_M': result.x[9],
        'delta_M': result.x[10],
    })

    residuals = objective_function(result.x, targets, fixed_params)
    rss = np.sum(residuals**2)
    rmse = np.sqrt(rss / len(residuals))

    print(f"\n[SUCCESS] Optimization complete:")
    print(f"  RSS = {rss:.3f}")
    print(f"  RMSE = {rmse:.3f} z-score units")

    # Check biological sanity
    print(f"\n[SANITY CHECK] Key parameter ratios:")
    print(
        f"  rho_C/delta_C = {opt_params['rho_C']/opt_params['delta_C']:.3f} (should be 0.1-0.5)")
    print(
        f"  beta_A/delta_A = {opt_params['beta_A']/opt_params['delta_A']:.3f} (should be 0.5-2.0)")
    print(
        f"  eta_M/delta_M = {opt_params['eta_M']/opt_params['delta_M']:.3f} (should be 0.5-2.0)")

    return opt_params, rss, rmse


# ============================================================================
# VALIDATION
# ============================================================================

def validate_fit(params, targets):
    """Validate with z-score transform."""
    # Simulate all conditions (raw)
    results_raw = {}
    for key, target in targets.items():
        f_ras = get_ras_input(target['condition'])

        if 'LeprKO' in target['condition']:
            params_ko = params.copy()
            params_ko['mu_L'] = params['mu_L'] * 0.2
            ss = simulate_steady_state(f_ras, params_ko)
        else:
            ss = simulate_steady_state(f_ras, params)

        results_raw[key] = {
            'C': ss[0],
            'A': ss[1],
            'T': ss[2],
            'M': ss[5],
            'condition': target['condition'],
            'dataset': target['dataset']
        }

    # Z-score transform
    for module in ['C', 'A', 'T', 'M']:
        values = [results_raw[k][module] for k in results_raw.keys()]
        mean_val = np.mean(values)
        std_val = np.std(values) if np.std(values) > 1e-6 else 1.0

        for key in results_raw.keys():
            zscore_key = f'{module}_model_zscore'
            results_raw[key][zscore_key] = (
                results_raw[key][module] - mean_val) / std_val

    # Build comparison dataframe
    rows = []
    for key, target in targets.items():
        for module in ['C', 'A', 'T', 'M']:
            rows.append({
                'dataset': target['dataset'],
                'condition': target['condition'],
                'module': module,
                'data_zscore': target[module],
                'model_zscore': results_raw[key][f'{module}_model_zscore'],
                'model_raw': results_raw[key][module],
                'residual': results_raw[key][f'{module}_model_zscore'] - target[module]
            })

    return pd.DataFrame(rows)


# ============================================================================
# MAIN
# ============================================================================

def main():
    """Main calibration pipeline."""
    print("\n" + "="*70)
    print("RAS-CSC MODEL CALIBRATION (CORRECTED)")
    print("="*70)

    targets = load_targets()
    params, rss, rmse = fit_model(targets)

    # Save parameters (as JSON to avoid numpy issues)
    param_file = RESULTS_DIR / "optimized_parameters_CORRECTED.json"
    with open(param_file, 'w') as f:
        json.dump(params, f, indent=2)
    print(f"\n[SAVED] Parameters: {param_file}")

    # Validate
    df_results = validate_fit(params, targets)

    results_csv = RESULTS_DIR / "model_vs_data_CORRECTED.csv"
    df_results.to_csv(results_csv, index=False)
    print(f"[SAVED] Results: {results_csv}")

    # Quick sanity check
    print(f"\n[SANITY CHECK] Model raw outputs:")
    for cond in df_results['condition'].unique():
        df_cond = df_results[df_results['condition'] == cond]
        C_raw = df_cond[df_cond['module'] == 'C']['model_raw'].values[0]
        M_raw = df_cond[df_cond['module'] == 'M']['model_raw'].values[0]
        print(f"  {cond:15s}: C={C_raw:.3f}, M={M_raw:.3f}")

    return params, df_results


if __name__ == "__main__":
    params, df_results = main()
