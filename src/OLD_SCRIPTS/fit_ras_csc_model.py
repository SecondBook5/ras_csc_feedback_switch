#!/usr/bin/env python3
"""
fit_ras_csc_model.py

Calibrate Ras–CSC–microenvironment feedback model parameters to
COMBINED RNA-Seq, ATAC-Seq, and scRNA-Seq derived module scores.

This script:
  1) Loads the combined RNA/ATAC summary.
  2) Loads the single-cell scores, calculates the mean (representing the SCC attractor),
     and merges it into the dataset.
  3) Pivots data to map:
       - TGFb  ~ T
       - Angio ~ A
       - CSC   ~ C
       - mTOR  ~ M
  4) Z-score normalizes the observed data across all conditions.
  5) Uses a random search to find parameters minimizing the error between
     simulated steady states and observed omics profiles.
  6) Exports results to data/processed/model_fits/model_calibration_results.json
"""

from __future__ import annotations

import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Tuple, List

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

# --- Import Model Definitions ---
# Ensure this file exists in the same directory or python path
try:
    from ras_csc_model_calib import RasCSCParams, ras_csc_rhs
except ImportError:
    sys.path.append(str(Path(__file__).parent))
    from ras_csc_model_calib import RasCSCParams, ras_csc_rhs

# --------------------------------------------------------------------
# Paths and Constants
# --------------------------------------------------------------------

PROJECT_ROOT = Path(".").resolve()

# Input Files
COMBINED_DATA_PATH = PROJECT_ROOT / "data" / "processed" / "omics" / "module_summary_rna_atac_combined.csv"
SCRNA_DATA_PATH = PROJECT_ROOT / "data" / "processed" / "omics_summaries" / "scc_scRNA_module_scores_per_cell.csv"

# Output Files
OUTPUT_DIR = PROJECT_ROOT / "data" / "processed" / "model_fits"
OUTPUT_PATH = OUTPUT_DIR / "model_calibration_results.json"

# Calibration settings
N_SAMPLES = 5000
RANDOM_SEED = 42

# --------------------------------------------------------------------
# Helper Functions
# --------------------------------------------------------------------

def get_condition_config(dataset: str, condition: str) -> Tuple[float, float]:
    """
    Map biological conditions to model inputs:
    Returns (f_ras, lepr_knockout_fraction)
    
    f_ras: 0.1 (Normal/Low) to 1.0 (High/SCC)
    lepr_ko: 0.0 (WT) to 1.0 (Full Knockout)
    """
    d = str(dataset)
    c = str(condition).lower()

    f_ras = 1.0
    lepr_ko = 0.0

    # --- RNA-Seq Mappings ---
    if d == "Bl6":
        if "normal" in c: f_ras = 0.1
        elif "scc" in c: f_ras = 1.0
        
    elif d == "PAP_SCC":
        if "papilloma" in c: f_ras = 0.6
        elif "scc" in c: f_ras = 1.0
        
    elif d == "PDV":
        if "leprko" in c: lepr_ko = 1.0
        # PDV is inherently Ras-driven (f_ras=1.0)

    # --- ATAC-Seq Mappings ---
    elif d == "ATAC":
        if "hfsc" in c or "ife" in c:
            f_ras = 0.1 # Normal stem cell states
        elif "papilloma" in c:
            f_ras = 0.6 # Benign
        elif "scc" in c:
            f_ras = 1.0 # Malignant

    # --- scRNA-Seq Mappings ---
    elif d == "scRNA":
        # The scRNA data comes from SCC tumors
        if "scc" in c:
            f_ras = 1.0

    return f_ras, lepr_ko

def simulate_steady_state(params: RasCSCParams) -> Dict[str, float]:
    """
    Run ODE to steady state and return dictionary of variables.
    """
    y0 = [0.1, 0.1, 0.1, 0.1]
    t_span = (0, 200)
    
    fun = lambda t, y: ras_csc_rhs(t, y, params)
    
    # Use BDF for stiff systems or RK45 for general
    sol = solve_ivp(fun, t_span, y0, method='RK45', rtol=1e-6, atol=1e-8)
    
    C, A, T, R = sol.y[:, -1]
    
    # Recompute algebraic variables
    L = params.L_sys + params.k_L_A * A
    S = L * R
    denom = (S ** params.q_M) + (params.K_M ** params.q_M)
    
    if denom > 1e-12:
        M = params.k_M_max * (S ** params.q_M) / denom
    else:
        M = 0.0
        
    return {"C": C, "A": A, "T": T, "M": M}

def load_and_prep_data() -> pd.DataFrame:
    """
    Load combined summary AND scRNA-seq data, aggregate, pivoting, and z-score.
    """
    # 1. Load Combined RNA/ATAC
    if not COMBINED_DATA_PATH.exists():
        raise FileNotFoundError(f"Combined data not found: {COMBINED_DATA_PATH}")
    df_combined = pd.read_csv(COMBINED_DATA_PATH)
    
    # 2. Load and Aggregate scRNA-seq
    if SCRNA_DATA_PATH.exists():
        print(f"[INFO] Loading scRNA-seq data from {SCRNA_DATA_PATH}...")
        df_sc = pd.read_csv(SCRNA_DATA_PATH)
        
        # Columns in scRNA: TGFb_module_score, mTOR_module_score, Angio_module_score, CSC_module_score
        # Map them to standard names
        sc_map = {
            'TGFb_module_score': 'TGFb_bulk',
            'mTOR_module_score': 'mTOR_bulk',
            'Angio_module_score': 'Angio_bulk',
            'CSC_module_score': 'CSC_bulk'
        }
        
        # Create a list of summary rows
        summary_rows = []
        for sc_col, std_col in sc_map.items():
            mean_val = df_sc[sc_col].mean()
            summary_rows.append({
                'omics_type': 'scRNA',
                'dataset': 'scRNA',
                'condition': 'SCC',
                'module': std_col,
                'n_samples': len(df_sc),
                'mean_score': mean_val,
                'sd_score': df_sc[sc_col].std()
            })
            
        df_sc_summary = pd.DataFrame(summary_rows)
        
        # Merge with main dataframe
        df_combined = pd.concat([df_combined, df_sc_summary], ignore_index=True)
        print(f"[INFO] Added scRNA-seq SCC aggregate (n={len(df_sc)} cells) to dataset.")
    else:
        print(f"[WARN] scRNA-seq file not found at {SCRNA_DATA_PATH}. Proceeding without it.")

    # 3. Pivot to wide format
    # Columns: omics_type, dataset, condition, module, mean_score
    df_wide = df_combined.pivot_table(
        index=['omics_type', 'dataset', 'condition'], 
        columns='module', 
        values='mean_score'
    ).reset_index()
    
    # Rename to model variables
    rename_map = {
        'CSC_bulk': 'obs_C',
        'Angio_bulk': 'obs_A',
        'TGFb_bulk': 'obs_T',
        'mTOR_bulk': 'obs_M'
    }
    df_wide = df_wide.rename(columns=rename_map)
    
    # 4. Z-score normalization across the entire dataset
    for col in ['obs_C', 'obs_A', 'obs_T', 'obs_M']:
        if col in df_wide.columns:
            vals = df_wide[col].values
            mu = np.nanmean(vals)
            sigma = np.nanstd(vals) + 1e-9
            df_wide[col] = (vals - mu) / sigma
        else:
            df_wide[col] = np.nan # Should not happen if data is complete

    return df_wide

def evaluate_params(base_params: RasCSCParams, data: pd.DataFrame) -> float:
    """
    Calculate MSE error for a parameter set.
    """
    sim_results = []
    
    # Simulate for every row (condition) in the data
    for _, row in data.iterrows():
        f_ras, lepr_ko = get_condition_config(row['dataset'], row['condition'])
        
        # Clone parameters
        p = RasCSCParams(**asdict(base_params))
        p.f_ras = f_ras
        
        # Apply Knockout Logic
        if lepr_ko > 0:
            p.k_R_prod *= (1.0 - lepr_ko)
            
        res = simulate_steady_state(p)
        sim_results.append([res['C'], res['A'], res['T'], res['M']])
        
    sim_arr = np.array(sim_results)
    
    # Z-score simulation results to match data scale
    sim_mean = np.mean(sim_arr, axis=0)
    sim_std = np.std(sim_arr, axis=0) + 1e-9
    sim_z = (sim_arr - sim_mean) / sim_std
    
    # Compute Error (MSE) excluding NaNs
    obs_arr = data[['obs_C', 'obs_A', 'obs_T', 'obs_M']].values
    mask = ~np.isnan(obs_arr)
    
    diff = (sim_z - obs_arr)
    mse = np.sum(diff[mask]**2) / np.sum(mask)
    
    return mse

# --------------------------------------------------------------------
# Main
# --------------------------------------------------------------------

def main():
    print(f"[INFO] Loading and preparing multi-omics data...")
    df_data = load_and_prep_data()
    
    print(f"[INFO] Data loaded. Found {len(df_data)} unique conditions across RNA, ATAC, scRNA.")
    print(f"[INFO] Starting Random Search calibration ({N_SAMPLES} iterations)...")
    
    rng = np.random.default_rng(RANDOM_SEED)
    best_error = np.inf
    best_params = None
    
    for i in range(N_SAMPLES):
        # Sample parameter candidates
        # We vary the gains (k) and feedback strengths significantly
        candidate = RasCSCParams(
            f_ras=1.0, # Base, modified by condition
            
            # CSC Dynamics
            k_C_ras = rng.uniform(0.01, 0.5),
            k_C_TGFb = rng.uniform(0.1, 1.0),
            k_C_M = rng.uniform(0.1, 3.0),      # <--- The key feedback parameter
            K_C = 0.5, n_C = 4.0,               # Fixed Hills for stability
            k_C_deg = 1.0,
            
            # Angiogenesis
            k_A_ras = rng.uniform(0.01, 0.5),
            k_A_C = rng.uniform(0.1, 2.0),
            k_A_deg = 1.0,
            
            # TGFb
            k_T_A = rng.uniform(0.5, 2.0),
            k_T_C = 0.0, # Stromal source assumption
            k_T_deg = 1.0,
            
            # Lepr
            k_R_prod = rng.uniform(0.1, 1.5),
            K_R = 0.5, p_R = 4.0,
            k_R_deg = 1.0,
            
            # Leptin/mTOR
            L_sys = 0.1,
            k_L_A = 1.0,
            k_M_max = rng.uniform(0.5, 1.5),
            K_M = 0.5, q_M = 4.0
        )
        
        try:
            error = evaluate_params(candidate, df_data)
            if error < best_error:
                best_error = error
                best_params = candidate
                if i % 500 == 0:
                    print(f"  [Iter {i}] New Best MSE: {best_error:.5f}")
        except Exception:
            continue

    print("-" * 40)
    print(f"[DONE] Final Best MSE: {best_error:.5f}")
    print("Best Parameters:")
    print(asdict(best_params))
    
    # Save Results
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    params_A_dict = asdict(best_params)
    params_B_dict = params_A_dict.copy()
    params_B_dict['k_C_M'] = 0.0 # No feedback version
    
    results = {
        "model_A_params": params_A_dict,
        "model_B_params": params_B_dict,
        "meta": {
            "description": "Calibrated to RNA+ATAC+scRNA combined dataset",
            "n_samples": N_SAMPLES,
            "best_mse": best_error
        }
    }
    
    with open(OUTPUT_PATH, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"[INFO] Calibration saved to {OUTPUT_PATH}")

if __name__ == "__main__":
    main()