# Ras–CSC Feedback Switch Model

This repository implements a small, self-contained Ras–CSC–microenvironment model calibrated to the cutaneous SCC system of:

> Yuan, S., Stewart, K.S., Yang, Y. et al. **Ras drives malignancy through stem cell crosstalk with the microenvironment.** *Nature* 612, 555–563 (2022). https://doi.org/10.1038/s41586-022-05475-6

The goal is to test a specific hypothesis: that a Ras-driven positive feedback loop between CSCs and their microenvironment creates a bistable “malignant switch” and that perturbing angiogenesis or LEPR/mTOR can collapse the high-CSC state even when Ras stays high.

---

## Model in one paragraph

The dynamical system tracks six variables:

- **C** – CSC fraction  
- **A** – Angiogenesis / vascular support  
- **T** – TGFβ signaling  
- **R** – LEPR (leptin receptor)  
- **L** – Leptin ligand  
- **M** – mTOR activity  

with the interaction chain:

> Ras → A → T → R → L → M → C,  
> and C → A closes the loop.

All ODEs and Hill functions live in `src/ras_csc_core.py`, and every downstream script imports this core so there is exactly one mathematical model.

---

## Repository layout (only the parts that matter)

```text
config/
  model_params.yaml          # Priors, bounds, Ras mapping
  gene_sets_rnaseq.yaml      # Bulk module gene sets (C/A/T/M)

pipelines/data_processing/
  01_process_bulk_rnaseq.R   # Bulk module scores → calibration targets
  02_process_scrna.R         # CSC calling, module scores, pseudotime
  03_process_atacseq.R       # ATAC module scores + rough checks

src/
  ras_csc_core.py            # Hill fn, ODE RHS, steady-state solver, Ras mapping
  scripts/
    04_run_model_calibration.py
    05_run_hypothesis_tests.py
    06_run_sensitivity_analysis.py
    07_run_phase_plane_plots.py
    08_run_temporal_analysis.py
    09_run_predictions.py

data/
  raw/                       # Source data (not tracked here)
  interim/                   # Cleaned matrices
  processed/                 # Module scores, calibration targets, validation logs

results/
  calibration/               # Fitted params + residuals
  model/                     # Hysteresis, time courses, derivatives
  predictions/               # Therapy scenarios, bimodality summary
  sensitivity/               # Parameter sensitivity and perturbation grids
  temporal/                  # scRNA pseudotime summaries

figures/
  main/                      # “Manuscript” figures
  phase_plane/
  temporal/
  sensitivity/
````

---

## How to reproduce the full pipeline

After setting up R and Python environments and placing raw data under `data/raw/` in the expected format:

```bash
# 1) Process omics data
Rscript pipelines/data_processing/01_process_bulk_rnaseq.R
Rscript pipelines/data_processing/02_process_scrna.R
Rscript pipelines/data_processing/03_process_atacseq.R

# 2) Calibrate model to bulk RNA-seq targets
python src/scripts/04_run_model_calibration.py

# 3) Core analyses
python src/scripts/05_run_hypothesis_tests.py        # Bistability + LEPR/mTOR tests
python src/scripts/06_run_sensitivity_analysis.py    # Parameter + perturbation sensitivity
python src/scripts/07_run_phase_plane_plots.py       # C–M phase plane + Jacobian spectra
python src/scripts/08_run_temporal_analysis.py       # Time course + scRNA pseudotime
python src/scripts/09_run_predictions.py             # Ras inhibition, anti-angiogenic therapy
```

All figures used in the writeup are generated under `figures/`, and all numbers cited in the text (calibration targets, CSC gaps, ΔC under perturbations, bimodality BICs, etc.) come from CSV/JSON files in `results/`.

---

## Core files to read first

If you are trying to understand the model rather than the plumbing:

1. `src/ras_csc_core.py` – the actual ODE system.
2. `config/model_params.yaml` – priors, bounds, and Ras mapping.
3. `results/calibration/optimized_parameters_CORRECTED.json` – fitted parameters.
4. `figures/main/ras_hysteresis_CORRECTED.png` – bistability.
5. `figures/main/lepr_mtor_perturbations_CORRECTED.png` – LEPR/mTOR vs CSC.

---

## Citation

If you reuse this model or pipeline, please cite both:

* Yuan et al., *Nature* 2022 (above) for the biological system and data.
* This repository as the source of the Ras–CSC feedback model implementation and calibration strategy.
