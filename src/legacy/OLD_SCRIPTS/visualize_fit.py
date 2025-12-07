#!/usr/bin/env python3
"""
visualize_fit.py

Visualizes the "Model vs Data" comparison for the calibrated Ras-CSC model.
It loads the parameters you just found and plots the steady-state values
alongside the Z-scored biological data.
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.integrate import solve_ivp
from dataclasses import dataclass, asdict  # <--- Fixed: Added asdict here

# --------------------------------------------------------------------
# 1. Define Model & Simulation Logic
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
    # C (Cancer Stem Cell)
    denom_C = (T ** params.n_C) + (params.K_C ** params.n_C)
    term_TGFb = params.k_C_TGFb * (T ** params.n_C) / denom_C if denom_C > 1e-12 else 0.0
    activation = params.k_C_ras * params.f_ras + term_TGFb + params.k_C_M * M
    dC = activation * (1 - C) - params.k_C_deg * C
    
    # A (Angiogenesis)
    dA = params.k_A_ras * params.f_ras + params.k_A_C * C - params.k_A_deg * A
    
    # T (TGFb)
    dT = params.k_T_A * A + params.k_T_C * C - params.k_T_deg * T
    
    # R (Lepr Receptor)
    denom_R = (T ** params.p_R) + (params.K_R ** params.p_R)
    term_R = params.k_R_prod * (T ** params.p_R) / denom_R if denom_R > 1e-12 else 0.0
    dR = term_R - params.k_R_deg * R
    
    return [dC, dA, dT, dR]

def simulate_steady_state(params):
    # Solve to steady state
    y0 = [0.1, 0.1, 0.1, 0.1]
    # Integrate for enough time to reach steady state
    sol = solve_ivp(lambda t, y: ras_csc_rhs(t, y, params), (0, 200), y0, method='RK45', rtol=1e-5)
    C, A, T, R = sol.y[:, -1]
    
    # Recompute M (algebraic variable)
    L = params.L_sys + params.k_L_A * A
    S = L * R
    denom_M = (S ** params.q_M) + (params.K_M ** params.q_M)
    M = params.k_M_max * (S ** params.q_M) / denom_M if denom_M > 1e-12 else 0.0
    
    return {"C": C, "A": A, "T": T, "M": M}

def get_config(dataset, condition):
    # Maps dataset conditions to f_ras and lepr_ko
    d, c = str(dataset), str(condition).lower()
    f_ras, lepr_ko = 1.0, 0.0

    if d == "Bl6":
        if "normal" in c: f_ras = 0.1
        elif "scc" in c: f_ras = 1.0
    elif d == "PAP_SCC":
        if "papilloma" in c: f_ras = 0.6
        elif "scc" in c: f_ras = 1.0
    elif d == "PDV":
        if "leprko" in c: lepr_ko = 1.0
    elif d == "ATAC":
        if "hfsc" in c or "ife" in c: f_ras = 0.1
        elif "papilloma" in c: f_ras = 0.6
        elif "scc" in c: f_ras = 1.0
    elif d == "scRNA":
        if "scc" in c: f_ras = 1.0

    return f_ras, lepr_ko

# --------------------------------------------------------------------
# 2. Load Data and Parameters
# --------------------------------------------------------------------

PROJECT_ROOT = Path(".").resolve()
DATA_PATH = PROJECT_ROOT / "data" / "processed" / "omics" / "module_summary_rna_atac_combined.csv"
SCRNA_PATH = PROJECT_ROOT / "data" / "processed" / "omics_summaries" / "scc_scRNA_module_scores_per_cell.csv"
PARAMS_PATH = PROJECT_ROOT / "data" / "processed" / "model_fits" / "model_calibration_results.json"

# A. Load Params
if not PARAMS_PATH.exists():
    print(f"[ERROR] Calibration results not found at {PARAMS_PATH}")
    exit(1)

with open(PARAMS_PATH) as f:
    calib_data = json.load(f)
    # Load Model A (Feedback ON) params
    best_params_dict = calib_data["model_A_params"]
    # Remove description field if present (it's not in the dataclass)
    if "description" in best_params_dict: del best_params_dict["description"]
    
    # Create param object
    fitted_params = RasCSCParams(**best_params_dict)
    print("[INFO] Loaded fitted parameters.")

# B. Load & Prep Data (Same logic as fit script)
df_combined = pd.read_csv(DATA_PATH)

# Add scRNA if available
if SCRNA_PATH.exists():
    df_sc = pd.read_csv(SCRNA_PATH)
    sc_map = {'TGFb_module_score': 'TGFb_bulk', 'mTOR_module_score': 'mTOR_bulk',
              'Angio_module_score': 'Angio_bulk', 'CSC_module_score': 'CSC_bulk'}
    summary_rows = []
    for sc_col, std_col in sc_map.items():
        summary_rows.append({
            'omics_type': 'scRNA', 'dataset': 'scRNA', 'condition': 'SCC',
            'module': std_col, 'mean_score': df_sc[sc_col].mean()
        })
    df_combined = pd.concat([df_combined, pd.DataFrame(summary_rows)], ignore_index=True)

# Pivot to wide
df_wide = df_combined.pivot_table(
    index=['omics_type', 'dataset', 'condition'], 
    columns='module', values='mean_score'
).reset_index()

# Rename columns to model variables
rename_map = {'CSC_bulk': 'C', 'Angio_bulk': 'A', 'TGFb_bulk': 'T', 'mTOR_bulk': 'M'}
df_wide = df_wide.rename(columns=rename_map)

# Z-score Normalize Data (Store mean/std to invert if needed, but we plot Z-scores)
data_z = df_wide.copy()
cols_to_plot = ['C', 'A', 'T', 'M']
stats = {}

for col in cols_to_plot:
    if col in df_wide.columns:
        mu = df_wide[col].mean()
        sigma = df_wide[col].std()
        # Avoid divide by zero
        if sigma == 0: sigma = 1.0
        data_z[col] = (df_wide[col] - mu) / sigma
        stats[col] = (mu, sigma)

# --------------------------------------------------------------------
# 3. Run Model on All Conditions
# --------------------------------------------------------------------

sim_rows = []
for _, row in df_wide.iterrows():
    f_ras, lepr_ko = get_config(row['dataset'], row['condition'])
    
    # Adjust params using the fitted values
    p = RasCSCParams(**asdict(fitted_params))
    p.f_ras = f_ras
    if lepr_ko > 0: p.k_R_prod *= (1.0 - lepr_ko)
    
    # Simulate
    res = simulate_steady_state(p)
    
    # Store raw result
    sim_entry = {
        'omics_type': row['omics_type'],
        'dataset': row['dataset'],
        'condition': row['condition'],
        'C': res['C'], 'A': res['A'], 'T': res['T'], 'M': res['M']
    }
    sim_rows.append(sim_entry)

df_sim = pd.DataFrame(sim_rows)

# Z-score Model Results (using model's own internal stats to match pattern)
df_sim_z = df_sim.copy()
for col in cols_to_plot:
    mu = df_sim[col].mean()
    sigma = df_sim[col].std()
    if sigma == 0: sigma = 1.0
    df_sim_z[col] = (df_sim[col] - mu) / sigma

# --------------------------------------------------------------------
# 4. Plotting
# --------------------------------------------------------------------

# Melt for plotting
data_melt = data_z.melt(id_vars=['dataset', 'condition'], value_vars=cols_to_plot, 
                        var_name='Variable', value_name='Z_Score')
data_melt['Source'] = 'Data'

sim_melt = df_sim_z.melt(id_vars=['dataset', 'condition'], value_vars=cols_to_plot, 
                         var_name='Variable', value_name='Z_Score')
sim_melt['Source'] = 'Model'

# Combine
plot_df = pd.concat([data_melt, sim_melt], ignore_index=True)

# Clean up condition names for plot
plot_df['Condition_Label'] = plot_df['dataset'] + "_" + plot_df['condition']

# Plot
plt.figure(figsize=(14, 8))
sns.set_style("whitegrid")

# Create barplot comparing Data vs Model side-by-side
g = sns.catplot(
    data=plot_df, 
    x='Condition_Label', 
    y='Z_Score', 
    hue='Source', 
    col='Variable', 
    col_wrap=2,
    kind='bar', 
    height=4, 
    aspect=1.5,
    palette={'Data': 'grey', 'Model': 'firebrick'},
    sharex=False 
)

g.set_xticklabels(rotation=45, ha='right')
g.fig.suptitle(f"Model Calibration Fit (MSE={calib_data.get('meta', {}).get('best_mse', 0):.3f})", y=1.02)

# Save
out_fig = PROJECT_ROOT / "figures" / "model_fit_validation.png"
out_fig.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(out_fig, bbox_inches='tight', dpi=300)

print(f"[INFO] Plot saved to: {out_fig}")


# --- NEW: Print the data table ---
print("\n[DATA] Comparison of Z-Scores (Model vs Data):")
print(plot_df.pivot_table(index=['dataset', 'condition', 'Variable'], columns='Source', values='Z_Score'))

print("[DONE] Open the image to verify the 'Low -> High' trend.")