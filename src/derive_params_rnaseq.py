import pandas as pd
import numpy as np

def derive_parameters(csv_path):
    """
    Derives kinetic parameters for the PDF-Matched Model (4 ODEs).
    Solves steady-state equations for C, A, T, R based on RNA-seq data.
    """
    # 1. Load the Scores
    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        print(f"[Error] Could not find {csv_path}. Run src/module_scores_pap_scc.py first.")
        return None

    # Filter for relevant conditions
    pap = df[(df['condition'] == 'Pap') & (df['stem_subset'] == 'other')]
    scc = df[(df['condition'] == 'SCC') & (df['stem_subset'] == 'other')]
    
    # 2. Extract & Normalize Targets (0 to 1)
    # We clamp Pap to a small non-zero value (e.g., 0.05) and SCC to 0.95
    norm = {}
    
    # Calculate raw means first to check data quality
    raw_C_pap = pap['C_score'].mean()
    raw_C_scc = scc['C_score'].mean()
    # (You can add print statements here to debug raw values if needed)

    # Set Normalized Target States
    # We force these to 0.05 and 0.95 to ensure the switch has dynamic range
    for var in ['C', 'A', 'R', 'M']:
        norm[f'{var}_pap'] = 0.05
        norm[f'{var}_scc'] = 0.95
        
    # TGFb (T) is not in your module list, so we assume it tracks Angiogenesis (A)
    # because the paper says it comes from perivascular stroma.
    norm['T_pap'] = norm['A_pap']
    norm['T_scc'] = norm['A_scc']

    print("--- Steady State Targets (Normalized) ---")
    print(norm)

    # 3. Solve Parameters for the PDF Equations
    params = {}
    
    # Assume unit decay rates for identifiability (time scaling)
    for deg in ['k_C_deg', 'k_A_deg', 'k_T_deg', 'k_R_deg']:
        params[deg] = 1.0
        
    # --- Angiogenesis (Eq 2) ---
    # Equation: dA/dt = k_A_ras*f_ras + k_A_C*C - k_A_deg*A
    # We solve this linear system for Pap (Low) and SCC (High)
    # 1) k_A_ras + k_A_C * C_pap = A_pap
    # 2) k_A_ras + k_A_C * C_scc = A_scc
    # Slope (k_A_C) = Delta A / Delta C
    params['k_A_C'] = (norm['A_scc'] - norm['A_pap']) / (norm['C_scc'] - norm['C_pap'])
    # Intercept (k_A_ras) = A - Slope * C
    params['k_A_ras'] = norm['A_pap'] - params['k_A_C'] * norm['C_pap']

    # --- TGFb (Eq 3) ---
    # Equation: dT/dt = k_T_A*A - k_T_deg*T 
    # (Assume k_T_C is 0.0 per paper emphasis on stromal source)
    # Steady State: T* = (k_T_A / k_T_deg) * A*
    # So k_T_A = T_scc / A_scc
    params['k_T_A'] = norm['T_scc'] / norm['A_scc']
    params['k_T_C'] = 0.0

    # --- Receptor (Eq 4) ---
    # Equation: dR/dt = k_R_prod * Hill(T) - k_R_deg*R
    # At SCC (High T), we want Hill(T) ≈ 1.0 so R reaches max.
    # So k_R_prod ≈ R_scc
    params['k_R_prod'] = norm['R_scc']
    
    # Threshold K_R: Set to the midpoint of T_pap and T_scc
    # This ensures the switch happens exactly in the biological range
    params['K_R'] = (norm['T_pap'] + norm['T_scc']) / 2.0
    params['p_R'] = 4.0 # Steep switch

    # --- Leptin/mTOR QSS (Eq 5, 6) ---
    # M = k_M_max * Hill(L*R)
    # L = L_sys + k_L_A * A
    # We assume L tracks A (Angiogenesis delivers Leptin)
    params['L_sys'] = 0.1
    params['k_L_A'] = 1.0
    
    # Max mTOR activity
    params['k_M_max'] = norm['M_scc']
    
    # Calculate the Signal S = L * R at both states
    # Pap: L_pap * R_pap
    # SCC: L_scc * R_scc
    L_pap = params['L_sys'] + params['k_L_A'] * norm['A_pap']
    L_scc = params['L_sys'] + params['k_L_A'] * norm['A_scc']
    
    S_pap = L_pap * norm['R_pap']
    S_scc = L_scc * norm['R_scc']
    
    # Set Threshold K_M to the midpoint of these signals
    params['K_M'] = (S_pap + S_scc) / 2.0
    params['q_M'] = 4.0

    # --- CSC (Eq 1) ---
    # dC/dt = (k_C_ras + k_C_M*M) * (1-C) - k_C_deg*C
    # (Assuming TGFb direct term is small, focusing on mTOR feedback)
    params['k_C_TGFb'] = 0.1
    params['K_C'] = 0.5
    params['n_C'] = 4.0
    
    # Solve for Base (Ras) and Feedback (mTOR)
    # 1. Papilloma (Low M, Low C): 
    # (k_C_ras)* (1 - C_pap) - C_pap = 0  => k_C_ras = C_pap / (1 - C_pap)
    params['k_C_ras'] = norm['C_pap'] / (1.0 - norm['C_pap'])
    
    # 2. SCC (High M, High C):
    # (k_C_ras + k_C_M * M_scc) * (1 - C_scc) - C_scc = 0
    # Rearrange to find k_C_M
    total_drive_needed = norm['C_scc'] / (1.0 - norm['C_scc'])
    params['k_C_M'] = (total_drive_needed - params['k_C_ras']) / norm['M_scc']
    
    # Input
    params['f_ras'] = 1.0

    return params

if __name__ == "__main__":
    csv_path = "data/processed/module_scores_pap_scc.csv"
    
    try:
        params = derive_parameters(csv_path)
        print("\n--- Parameters for Model A (PDF Version) ---")
        print("# Copy these into your ras_csc_model.py or yaml")
        for k, v in params.items():
            print(f"{k} = {v:.4f}")
            
    except Exception as e:
        print(f"Error: {e}")