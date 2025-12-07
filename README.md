# Ras–CSC Feedback Switch Model

By AJ Book and Shrinithi Natarajan

This repository implements a compact Ras–CSC–microenvironment model calibrated to the cutaneous SCC system of:

> Yuan, S., Stewart, K.S., Yang, Y. et al. **Ras drives malignancy through stem cell crosstalk with the microenvironment.**  
> *Nature* 612, 555–563 (2022). https://doi.org/10.1038/s41586-022-05475-6

The central hypothesis is that a Ras-driven positive feedback loop between CSCs and their microenvironment produces a **bistable malignant switch**.  
The model predicts that perturbing angiogenesis or the LEPR–mTOR axis can collapse the high-CSC state even if Ras remains high.

---

## Model in one paragraph

The dynamical system tracks:

- **C** – CSC fraction  
- **A** – Angiogenesis  
- **T** – TGFβ signaling  
- **R** – LEPR  
- **L** – Leptin  
- **M** – mTOR  

with the feedforward chain:

> Ras → A → T → R → L → M → C, with CSCs feeding back to angiogenesis.

The ODEs and Hill activations live in `src/ras_csc_core.py`, and every analysis script imports this single source to avoid drift.

---

## Software requirements

You need both R and Python:

- **R** ≥ 4.x with CRAN + Bioconductor packages installed via:

  ```bash
  Rscript R/install_packages.R
````

* **Python** ≥ 3.10 with the packages listed in `requirements.txt`:

  ```bash
  python -m venv .venv_ras-csc-switch-mm
  source .venv_ras-csc-switch-mm/bin/activate

  pip install -r requirements.txt
  ```

`requirements.txt` includes:

```text
numpy
pandas
scipy
scikit-learn
matplotlib
seaborn
pyyaml
```

---

## Repository layout

High-level structure:

```text
R/
  install_packages.R          # Installs required CRAN/Bioconductor packages

config/
  model_params.yaml           # Priors, bounds, Ras mapping
  gene_sets_rnaseq.yaml       # Gene sets for bulk module scoring
  scrna_module_genesets.csv   # scRNA module gene sets

data/
  raw/                        # GEO data + Yuan MOESM Excel files (as downloaded)
    GSE190411/                # Bulk RNA-seq (Normal, Papilloma, SCC, PDV, LEPR KO)
    GSE190414/                # ATAC-seq narrowPeak files
    GSE207975/                # SCC scRNA-seq (counts + TPM)
    yuan2022_moesm/           # Supplementary tables from Yuan et al.
  interim/
    rnaseq/                   # Parsed / harmonized bulk expression tables
    atac/                     # Gene-level accessibility + ATAC sample metadata
    omics_qc/                 # QC tables for omics processing
  processed/
    rnaseq/                   # Calibration targets + validation reports
    scrna/                    # CSC calls, module scores, mechanistic correlations
    atac/                     # ATAC module scores + cross-omics comparisons
    model_fits/               # Model calibration outputs (JSON/CSV)

pipelines/
  data_processing/
    01_process_bulk_rnaseq.R  # Bulk module scores → calibration targets
    02_process_scrna.R        # CSC calling, module scores, pseudotime, correlations
    03_process_atacseq.R      # ATAC module scores + cross-omics with RNA
  analysis/
    check_rnaseq_module_scores.R  # Optional QC of RNA-seq module scores

src/
  __init__.py
  ras_csc_core.py             # Core ODE system + steady-state solver
  scripts/
    04_run_model_calibration.py
    05_run_hypothesis_tests.py
    06_run_sensitivity_analysis.py
    07_run_phase_plane_plots.py
    08_run_temporal_analysis.py
    09_run_predictions.py

results/
  calibration/                # Fitted parameters, residuals, model vs data
  model/                      # Hysteresis, flux field, time courses, steady states
  predictions/                # Therapy predictions + bimodality summary
  perturbations/              # Full perturbation grid (LEPR/mTOR/angiogenesis)
  hypothesis_tests/           # Summary table of hypothesis test outcomes
  sensitivity/                # RSS splits + sensitivity to parameters and perturbations
  temporal/                   # Temporal alignment with scRNA pseudotime

figures/
  main/                       # Primary “manuscript-style” figures (Figures 1,5–7)
  model/                      # Additional model diagnostics (derivatives, seeds)
  phase_plane/                # Phase planes + Jacobian spectra
  sensitivity/                # Calibration / hysteresis / perturbation sensitivity
  temporal/                   # Temporal model vs scRNA pseudotime panels

docs/
  equations.md                # Full equation set and parameter definitions
  hypothesis.md               # Verbal and mathematical statement of the hypothesis
  outline.md                  # Paper outline
  paper_figure1.md            # Text for main figures (1–5)
  paper_figure2.md
  paper_figure3.md
  paper_figure4.md
  paper_figure5.md

external/
  README_external.md          # Notes on external scripts/notebooks
  yuan2022/
    Technical_noise_test_Leo_*.R
    Yuanetal2022_Nature_scRNAseq_jupyternotebook.ipynb

notebooks/
  scRNA_SCC_RasCSC_state_summary.ipynb  # Exploratory scRNA summaries (not used by pipeline)

run_ras_csc_feedforward_hypothesis.sh   # One-click pipeline script (R + Python)
requirements.txt                         # Python dependencies
README.md
```

---

## Quick start: run the entire hypothesis pipeline

After installing R and Python dependencies:

```bash
# Activate the Python environment
source .venv_ras-csc-switch-mm/bin/activate

# Run the full pipeline
./run_ras_csc_feedforward_hypothesis.sh
```

This script performs, in order:

1. Bulk RNA-seq processing and module scoring
2. scRNA-seq processing, CSC calling, and mechanistic correlations
3. ATAC-seq processing and cross-omics comparison
4. Model calibration to the bulk calibration targets
5. Hypothesis tests (bistability, Ras hysteresis, perturbations)
6. Sensitivity analyses
7. Phase-plane and Jacobian diagnostics
8. Temporal alignment between model trajectories and scRNA pseudotime
9. Therapy predictions, including Ras inhibition and anti-angiogenic scenarios

All numeric outputs go to `results/`, and manuscript-style figures go to `figures/main/`.

---

## Manual reproduction of each step

If you prefer to run steps individually:

```bash
# 1–3: Omics → calibration targets
Rscript pipelines/data_processing/01_process_bulk_rnaseq.R
Rscript pipelines/data_processing/02_process_scrna.R
Rscript pipelines/data_processing/03_process_atacseq.R

# (Optional) Check RNA-seq module scores
Rscript pipelines/analysis/check_rnaseq_module_scores.R

# 4–9: Model → results
python src/scripts/04_run_model_calibration.py
python src/scripts/05_run_hypothesis_tests.py
python src/scripts/06_run_sensitivity_analysis.py
python src/scripts/07_run_phase_plane_plots.py
python src/scripts/08_run_temporal_analysis.py
python src/scripts/09_run_predictions.py
```

Key outputs:

* `results/calibration/optimized_parameters_CORRECTED.json`
* `results/model/ras_hysteresis_CORRECTED.csv`
* `results/predictions/prediction_bimodality_summary_CORRECTED.json`
* `results/temporal/*.csv`

Figures used in the manuscript/live write-up are under `figures/main/`.

---

## Key files to understand first

1. `docs/hypothesis.md` and `docs/equations.md` – verbal and mathematical statement of the model.
2. `src/ras_csc_core.py` – the actual ODE model and steady-state solvers.
3. `config/model_params.yaml` – priors, bounds, Ras mapping, and parameter configuration.
4. `results/calibration/optimized_parameters_CORRECTED.json` – fitted parameter set.
5. `figures/main/ras_hysteresis_CORRECTED.png` – bistability and Ras hysteresis.
6. `figures/main/lepr_mtor_perturbations_CORRECTED.png` – LEPR/mTOR/angiogenesis perturbation outcomes.

---

## Data sources

All datasets are public and come directly from the Yuan et al. system:

* **GSE190411** – bulk RNA-seq (Bl6 Normal, Papilloma, SCC, PDV, LEPR KO).
* **GSE207975** – SCC scRNA-seq.
* **GSE190414** – HFSC/IFE/PAP/SCC ATAC-seq.
* **MOESM Excel files** – supplementary tables from the Yuan et al. paper.

Raw files are stored under `data/raw/` exactly as downloaded; all intermediate and processed derivatives live under `data/interim/` and `data/processed/`.

---

## Citation

If this model or workflow is used in your work, please cite:

**Biological system + data:**
Yuan et al., *Nature* 2022.

**Model + pipeline implementation:**
This repository: *Ras–CSC Feedback Switch Model*.


