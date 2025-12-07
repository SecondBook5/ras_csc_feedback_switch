#!/usr/bin/env python3
"""
run_simulation.py

This script:
1. Loads the parameters from config/model_params.yaml.
2. Sets up two scenarios:
   - Model A: Full Feedback (Malignant/SCC)
   - Model B: Broken Feedback (Inhibited/Benign)
3. Simulates the ODEs for 60 days.
4. Plots the CSC Fraction (C) and Tumor Volume to visualize the switch collapse.
"""

import os
import yaml
import numpy as np
import matplotlib.pyplot as plt
from ras_csc_model import RasCSCParams, ras_csc_rhs, simulate_trajectory

def load_config(yaml_path):
    """Load the YAML config file."""
    if not os.path.exists(yaml_path):
        raise FileNotFoundError(f"Config not found at {yaml_path}")
    with open(yaml_path, 'r') as f:
        return yaml.safe_load(f)

def calculate_tumor_volume(time, c_state, r_benign=0.03, r_malignant=0.08):
    """
    Estimates tumor volume based on growth rates derived from Figure 3/4.
    Volume grows at 'r_benign' when C=0, and 'r_malignant' when C=1.
    """
    vol = [50.0] # Start palpable at 50mm3
    for i in range(1, len(time)):
        dt = time[i] - time[i-1]
        current_c = c_state[i]
        # Interpolate rate
        r = r_benign + current_c * (r_malignant - r_benign)
        # Exponential growth step
        new_vol = vol[-1] * np.exp(r * dt)
        vol.append(new_vol)
    return np.array(vol)

def main():
    # 1. Setup Paths
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    config_path = os.path.join(base_dir, "config", "model_params.yaml")
    
    print(f"Loading configuration from: {config_path}")
    config = load_config(config_path)

    # 2. Run Model A (Full Feedback - Malignant)
    print("Running Model A (Malignant Feedback Loop)...")
    params_dict_A = config['models']['model_A_feedback']['params']
    params_A = RasCSCParams(**params_dict_A)
    
    # Initial Condition: Start with a small "seed" of CSCs (Papilloma transition)
    # y = [C, A, T, R]
    y0 = [0.60, 0.60, 0.60, 0.60] 
    t_span = (0, 60) # 60 Days
    t, y_A = simulate_trajectory(ras_csc_rhs, y0, params_A, t_span)
    vol_A = calculate_tumor_volume(t, y_A[:, 0])

    # 3. Run Model B (No Feedback - Therapeutic)
    print("Running Model B (Feedback Broken/Inhibited)...")
    params_dict_B = config['models']['model_B_no_feedback']['params']
    params_B = RasCSCParams(**params_dict_B)
    
    # Same initial condition
    _, y_B = simulate_trajectory(ras_csc_rhs, y0, params_B, t_span)
    vol_B = calculate_tumor_volume(t, y_B[:, 0])

    # 4. Visualize the Hypothesis
    print("Plotting results...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Plot 1: The Internal Switch (CSC State)
    ax1.plot(t, y_A[:, 0], 'r-', linewidth=3, label='Full Feedback (SCC)')
    ax1.plot(t, y_B[:, 0], 'b--', linewidth=2, label='Feedback Broken (Inhibited)')
    ax1.set_title('Hypothesis Test: Stability of CSC State', fontsize=14)
    ax1.set_ylabel('Cancer Stem Cell Fraction (C)', fontsize=12)
    ax1.set_xlabel('Time (Days)', fontsize=12)
    ax1.set_ylim(-0.05, 1.05)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)

    # Plot 2: The Phenotypic Outcome (Tumor Volume)
    ax2.plot(t, vol_A, 'r-', linewidth=3, label='Malignant Growth')
    ax2.plot(t, vol_B, 'b--', linewidth=2, label='Benign Growth')
    ax2.set_title('Clinical Outcome: Tumor Volume', fontsize=14)
    ax2.set_ylabel('Volume (mm3)', fontsize=12)
    ax2.set_xlabel('Time (Days)', fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)

    plt.tight_layout()
    
    # Save plot
    out_file = os.path.join(base_dir, "hypothesis_test_result.png")
    plt.savefig(out_file)
    print(f"Graph saved to: {out_file}")
    plt.show()

if __name__ == "__main__":
    main()