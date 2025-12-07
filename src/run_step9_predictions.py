#!/usr/bin/env python3
"""
run_step9_predictions.py

Step 9: Make predictions for scenarios NOT used to construct the model
- Gene knockouts/overexpression
- Drug combinations
- Alternative perturbations
- Novel treatment protocols
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
from scipy.integrate import odeint
import json
from pathlib import Path
from typing import Dict, Tuple

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

def ras_csc_model(y: np.ndarray, t: float, f_ras: float, params: Dict) -> list:
    """6-variable Ras-CSC ODE system."""
    C, A, T, R, L, M = y

    # Clip to valid range
    C = np.clip(C, 0, 1)
    A = np.clip(A, 0, 1)
    T = np.clip(T, 0, 1)
    R = np.clip(R, 0, 1)
    L = np.clip(L, 0, 1)
    M = np.clip(M, 0, 1)

    def hill(x: float, K: float, n: float) -> float:
        """Hill function."""
        x = np.clip(x, 0, 10)
        K = np.clip(K, 0.01, 10)
        n = np.clip(n, 1, 6)
        denom = K**n + x**n
        return (x**n) / denom if denom > 0 else 0

    # ODEs
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


def simulate_trajectory(f_ras: float, params: Dict, y0: np.ndarray,
                        t_max: float = 200) -> Tuple[np.ndarray, np.ndarray]:
    """
    Simulate full trajectory.
    
    Returns:
        Tuple of (time, solution_array)
    """
    t = np.linspace(0, t_max, 1000)

    try:
        sol = odeint(ras_csc_model, y0, t, args=(f_ras, params),
                     atol=1e-8, rtol=1e-8)
        return t, sol
    except Exception as e:
        print(f"  [ERROR] Integration failed: {e}")
        return t, np.tile(y0, (len(t), 1))


# ============================================================================
# PREDICTION 1: GENE KNOCKOUTS
# ============================================================================

def predict_gene_knockouts(params: Dict) -> pd.DataFrame:
    """
    Predict effects of gene knockouts NOT tested in calibration.
    
    Calibration used: LEPR knockout (mu_L reduction)
    
    New predictions:
    1. TGFβ receptor knockout (blocks T→R)
    2. Angiogenesis knockout (blocks CSC→A)
    3. mTOR knockout (blocks M activation)
    4. Double knockouts
    """
    print("\n[PREDICTION 1] Gene knockout scenarios...")

    # Malignant baseline
    f_ras = 0.95
    y0_malignant = [0.665, 0.500, 0.300, 0.400, 0.400, 0.564]

    results = []

    # Baseline
    t, sol = simulate_trajectory(f_ras, params, y0_malignant)
    results.append({
        'scenario': 'Baseline',
        'gene': 'WT',
        'perturbation': 'None',
        'C_final': sol[-1, 0],
        'M_final': sol[-1, 5],
        'reduction_C': 0,
        'reduction_M': 0
    })
    baseline_C = sol[-1, 0]
    baseline_M = sol[-1, 5]

    # 1. TGFβ receptor knockout (eta_R = 0)
    params_tgfbr_ko = params.copy()
    params_tgfbr_ko['eta_R'] = 0
    t, sol = simulate_trajectory(f_ras, params_tgfbr_ko, y0_malignant)
    results.append({
        'scenario': 'TGFβR_KO',
        'gene': 'TGFBR',
        'perturbation': 'eta_R = 0',
        'C_final': sol[-1, 0],
        'M_final': sol[-1, 5],
        'reduction_C': baseline_C - sol[-1, 0],
        'reduction_M': baseline_M - sol[-1, 5]
    })

    # 2. Angiogenesis knockout (beta_A = 0, CSC cannot promote angio)
    params_angio_ko = params.copy()
    params_angio_ko['beta_A'] = 0
    t, sol = simulate_trajectory(f_ras, params_angio_ko, y0_malignant)
    results.append({
        'scenario': 'Angio_KO',
        'gene': 'VEGF',
        'perturbation': 'beta_A = 0',
        'C_final': sol[-1, 0],
        'M_final': sol[-1, 5],
        'reduction_C': baseline_C - sol[-1, 0],
        'reduction_M': baseline_M - sol[-1, 5]
    })

    # 3. mTOR knockout (eta_M = 0, complete loss)
    params_mtor_ko = params.copy()
    params_mtor_ko['eta_M'] = 0
    t, sol = simulate_trajectory(f_ras, params_mtor_ko, y0_malignant)
    results.append({
        'scenario': 'mTOR_KO',
        'gene': 'MTOR',
        'perturbation': 'eta_M = 0',
        'C_final': sol[-1, 0],
        'M_final': sol[-1, 5],
        'reduction_C': baseline_C - sol[-1, 0],
        'reduction_M': baseline_M - sol[-1, 5]
    })

    # 4. TGFβ secretion knockout (gamma_T = 0, CSC cannot secrete TGFβ)
    params_tgfb_ko = params.copy()
    params_tgfb_ko['gamma_T'] = 0
    t, sol = simulate_trajectory(f_ras, params_tgfb_ko, y0_malignant)
    results.append({
        'scenario': 'TGFb_secretion_KO',
        'gene': 'TGFB1',
        'perturbation': 'gamma_T = 0',
        'C_final': sol[-1, 0],
        'M_final': sol[-1, 5],
        'reduction_C': baseline_C - sol[-1, 0],
        'reduction_M': baseline_M - sol[-1, 5]
    })

    # 5. Double knockout: LEPR + mTOR
    params_double1 = params.copy()
    params_double1['mu_L'] = params['mu_L'] * 0.1
    params_double1['eta_M'] = 0
    t, sol = simulate_trajectory(f_ras, params_double1, y0_malignant)
    results.append({
        'scenario': 'LEPR_mTOR_double_KO',
        'gene': 'LEPR+MTOR',
        'perturbation': 'mu_L×0.1, eta_M=0',
        'C_final': sol[-1, 0],
        'M_final': sol[-1, 5],
        'reduction_C': baseline_C - sol[-1, 0],
        'reduction_M': baseline_M - sol[-1, 5]
    })

    # 6. Double knockout: TGFβR + Angio
    params_double2 = params.copy()
    params_double2['eta_R'] = 0
    params_double2['beta_A'] = 0
    t, sol = simulate_trajectory(f_ras, params_double2, y0_malignant)
    results.append({
        'scenario': 'TGFbR_Angio_double_KO',
        'gene': 'TGFBR+VEGF',
        'perturbation': 'eta_R=0, beta_A=0',
        'C_final': sol[-1, 0],
        'M_final': sol[-1, 5],
        'reduction_C': baseline_C - sol[-1, 0],
        'reduction_M': baseline_M - sol[-1, 5]
    })

    df = pd.DataFrame(results)
    print(f"  [COMPLETE] Predicted {len(results)-1} knockout scenarios")
    print(f"\n  Most effective knockout (CSC reduction):")
    best_idx = df['reduction_C'].idxmax()
    print(
        f"    {df.loc[best_idx, 'scenario']}: ΔC = {df.loc[best_idx, 'reduction_C']:.3f}")

    return df


# ============================================================================
# PREDICTION 2: DRUG COMBINATIONS
# ============================================================================

def predict_drug_combinations(params: Dict) -> pd.DataFrame:
    """
    Predict synergistic drug combinations NOT tested.
    
    New combinations:
    1. mTOR + TGFβ inhibitors
    2. mTOR + Angio inhibitors
    3. Triple combination
    4. Sequential vs simultaneous
    """
    print("\n[PREDICTION 2] Drug combination scenarios...")

    f_ras = 0.95
    y0_malignant = [0.665, 0.500, 0.300, 0.400, 0.400, 0.564]

    results = []

    # Helper function for dose combinations
    def apply_drugs(params_base, mtor_dose=0, tgfb_dose=0, angio_dose=0):
        p = params_base.copy()
        p['eta_M'] = p['eta_M'] * (1 - mtor_dose/100)
        p['delta_T'] = p['delta_T'] * \
            (1 + tgfb_dose/100)  # Increase TGFβ decay
        p['beta_A'] = p['beta_A'] * (1 - angio_dose/100)
        return p

    # Baseline
    baseline_C = 0.665
    baseline_M = 0.564

    # Single agents at 50% dose
    dose = 50

    params_mtor = apply_drugs(params, mtor_dose=dose)
    t, sol = simulate_trajectory(f_ras, params_mtor, y0_malignant)
    results.append({
        'scenario': f'mTOR_inh_{dose}pct',
        'drugs': 'mTOR inhibitor',
        'dose': f'{dose}%',
        'C_final': sol[-1, 0],
        'M_final': sol[-1, 5],
        'reduction_C': baseline_C - sol[-1, 0]
    })

    params_tgfb = apply_drugs(params, tgfb_dose=dose)
    t, sol = simulate_trajectory(f_ras, params_tgfb, y0_malignant)
    results.append({
        'scenario': f'TGFb_inh_{dose}pct',
        'drugs': 'TGFβ inhibitor',
        'dose': f'{dose}%',
        'C_final': sol[-1, 0],
        'M_final': sol[-1, 5],
        'reduction_C': baseline_C - sol[-1, 0]
    })

    params_angio = apply_drugs(params, angio_dose=dose)
    t, sol = simulate_trajectory(f_ras, params_angio, y0_malignant)
    results.append({
        'scenario': f'Angio_inh_{dose}pct',
        'drugs': 'Angio inhibitor',
        'dose': f'{dose}%',
        'C_final': sol[-1, 0],
        'M_final': sol[-1, 5],
        'reduction_C': baseline_C - sol[-1, 0]
    })

    # Dual combinations (50% each)
    params_mtor_tgfb = apply_drugs(params, mtor_dose=dose, tgfb_dose=dose)
    t, sol = simulate_trajectory(f_ras, params_mtor_tgfb, y0_malignant)
    results.append({
        'scenario': f'mTOR_TGFb_combo_{dose}pct',
        'drugs': 'mTOR + TGFβ',
        'dose': f'{dose}% each',
        'C_final': sol[-1, 0],
        'M_final': sol[-1, 5],
        'reduction_C': baseline_C - sol[-1, 0]
    })

    params_mtor_angio = apply_drugs(params, mtor_dose=dose, angio_dose=dose)
    t, sol = simulate_trajectory(f_ras, params_mtor_angio, y0_malignant)
    results.append({
        'scenario': f'mTOR_Angio_combo_{dose}pct',
        'drugs': 'mTOR + Angio',
        'dose': f'{dose}% each',
        'C_final': sol[-1, 0],
        'M_final': sol[-1, 5],
        'reduction_C': baseline_C - sol[-1, 0]
    })

    params_tgfb_angio = apply_drugs(params, tgfb_dose=dose, angio_dose=dose)
    t, sol = simulate_trajectory(f_ras, params_tgfb_angio, y0_malignant)
    results.append({
        'scenario': f'TGFb_Angio_combo_{dose}pct',
        'drugs': 'TGFβ + Angio',
        'dose': f'{dose}% each',
        'C_final': sol[-1, 0],
        'M_final': sol[-1, 5],
        'reduction_C': baseline_C - sol[-1, 0]
    })

    # Triple combination (50% each)
    params_triple = apply_drugs(
        params, mtor_dose=dose, tgfb_dose=dose, angio_dose=dose)
    t, sol = simulate_trajectory(f_ras, params_triple, y0_malignant)
    results.append({
        'scenario': f'Triple_combo_{dose}pct',
        'drugs': 'mTOR + TGFβ + Angio',
        'dose': f'{dose}% each',
        'C_final': sol[-1, 0],
        'M_final': sol[-1, 5],
        'reduction_C': baseline_C - sol[-1, 0]
    })

    df = pd.DataFrame(results)
    print(f"  [COMPLETE] Predicted {len(results)} combination scenarios")
    print(f"\n  Most effective combination:")
    best_idx = df['reduction_C'].idxmax()
    print(
        f"    {df.loc[best_idx, 'drugs']}: ΔC = {df.loc[best_idx, 'reduction_C']:.3f}")

    return df


# ============================================================================
# PREDICTION 3: OVEREXPRESSION SCENARIOS
# ============================================================================

def predict_overexpression(params: Dict) -> pd.DataFrame:
    """
    Predict effects of gene overexpression.
    
    NOT tested in calibration:
    1. mTOR overexpression (2× eta_M)
    2. TGFβ overexpression (2× gamma_T)
    3. LEPR overexpression (2× eta_R)
    """
    print("\n[PREDICTION 3] Overexpression scenarios...")

    # Start from benign state
    f_ras = 0.3
    y0_benign = [0.039, 0.100, 0.100, 0.100, 0.100, 0.001]

    results = []

    # Baseline benign
    t, sol = simulate_trajectory(f_ras, params, y0_benign)
    results.append({
        'scenario': 'Benign_baseline',
        'gene': 'WT',
        'fold_change': '1×',
        'C_final': sol[-1, 0],
        'M_final': sol[-1, 5],
        'transforms_to_malignant': 'No'
    })
    baseline_C = sol[-1, 0]

    # mTOR overexpression (2×)
    params_mtor_oe = params.copy()
    params_mtor_oe['eta_M'] = params['eta_M'] * 2
    t, sol = simulate_trajectory(f_ras, params_mtor_oe, y0_benign)
    transforms = 'Yes' if sol[-1, 0] > 0.3 else 'No'
    results.append({
        'scenario': 'mTOR_overexpression_2x',
        'gene': 'MTOR',
        'fold_change': '2×',
        'C_final': sol[-1, 0],
        'M_final': sol[-1, 5],
        'transforms_to_malignant': transforms
    })

    # TGFβ overexpression (2×)
    params_tgfb_oe = params.copy()
    params_tgfb_oe['gamma_T'] = params['gamma_T'] * 2
    t, sol = simulate_trajectory(f_ras, params_tgfb_oe, y0_benign)
    transforms = 'Yes' if sol[-1, 0] > 0.3 else 'No'
    results.append({
        'scenario': 'TGFb_overexpression_2x',
        'gene': 'TGFB1',
        'fold_change': '2×',
        'C_final': sol[-1, 0],
        'M_final': sol[-1, 5],
        'transforms_to_malignant': transforms
    })

    # LEPR overexpression (2×)
    params_lepr_oe = params.copy()
    params_lepr_oe['eta_R'] = params['eta_R'] * 2
    t, sol = simulate_trajectory(f_ras, params_lepr_oe, y0_benign)
    transforms = 'Yes' if sol[-1, 0] > 0.3 else 'No'
    results.append({
        'scenario': 'LEPR_overexpression_2x',
        'gene': 'LEPR',
        'fold_change': '2×',
        'C_final': sol[-1, 0],
        'M_final': sol[-1, 5],
        'transforms_to_malignant': transforms
    })

    # mTOR overexpression (3×)
    params_mtor_oe3 = params.copy()
    params_mtor_oe3['eta_M'] = params['eta_M'] * 3
    t, sol = simulate_trajectory(f_ras, params_mtor_oe3, y0_benign)
    transforms = 'Yes' if sol[-1, 0] > 0.3 else 'No'
    results.append({
        'scenario': 'mTOR_overexpression_3x',
        'gene': 'MTOR',
        'fold_change': '3×',
        'C_final': sol[-1, 0],
        'M_final': sol[-1, 5],
        'transforms_to_malignant': transforms
    })

    df = pd.DataFrame(results)
    print(f"  [COMPLETE] Predicted {len(results)-1} overexpression scenarios")

    # Check transformations
    transforms = df[df['transforms_to_malignant'] == 'Yes']
    if len(transforms) > 0:
        print(f"\n  Transformative overexpressions:")
        for _, row in transforms.iterrows():
            print(
                f"    {row['gene']} ({row['fold_change']}): C = {row['C_final']:.3f}")

    return df


# ============================================================================
# MAIN
# ============================================================================

def main():
    """Run Step 9 predictions."""
    print("\n" + "="*70)
    print("STEP 9: PREDICTIONS FOR NEW SCENARIOS")
    print("NOT used in model calibration")
    print("="*70)

    # Load parameters
    with open(RESULTS_DIR / "optimized_parameters_CORRECTED.json", 'r') as f:
        params = json.load(f)

    print("\n[INFO] Loaded optimized parameters")

    # Prediction 1: Gene knockouts
    df_ko = predict_gene_knockouts(params)
    ko_file = RESULTS_DIR / "step9_gene_knockouts.csv"
    df_ko.to_csv(ko_file, index=False)
    print(f"[SAVED] {ko_file}")

    # Prediction 2: Drug combinations
    df_combo = predict_drug_combinations(params)
    combo_file = RESULTS_DIR / "step9_drug_combinations.csv"
    df_combo.to_csv(combo_file, index=False)
    print(f"[SAVED] {combo_file}")

    # Prediction 3: Overexpression
    df_oe = predict_overexpression(params)
    oe_file = RESULTS_DIR / "step9_overexpression.csv"
    df_oe.to_csv(oe_file, index=False)
    print(f"[SAVED] {oe_file}")

    print("\n" + "="*70)
    print("STEP 9 PREDICTIONS COMPLETE")
    print("="*70)
    print("\nKey predictions:")
    print(
        f"  1. Best knockout: {df_ko.iloc[df_ko['reduction_C'].idxmax()]['scenario']}")
    print(
        f"  2. Best drug combo: {df_combo.iloc[df_combo['reduction_C'].idxmax()]['drugs']}")
    print(f"  3. Transformative OE: Check {oe_file}")
    print("\nNext: Compare these predictions to Yuan 2022 experimental data (Step 10)")


if __name__ == "__main__":
    main()
