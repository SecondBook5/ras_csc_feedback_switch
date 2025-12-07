#!/usr/bin/env bash
# run_ras_csc_feedforward_hypothesis.sh
#
# Run the full Ras–CSC feedforward loop hypothesis pipeline:
#   01_process_bulk_rnaseq.R
#   02_process_scrna.R
#   03_process_atacseq.R
#   04_run_model_calibration.py
#   05_run_hypothesis_tests.py
#   06_run_sensitivity_analysis.py
#   07_run_phase_plane_plots.py
#   08_run_temporal_analysis.py
#   09_run_predictions.py
#
# Assumes:
#   - You run this from the repository root
#   - Rscript is on PATH
#   - python points to the environment where dependencies are installed

set -euo pipefail

# Simple error trap to make failures obvious
trap 'echo "[ERROR] Pipeline failed at line $LINENO" >&2' ERR

echo "======================================================================"
echo "Ras–CSC Feedforward Loop Hypothesis Pipeline"
echo "======================================================================"
echo "Working directory: $(pwd)"
echo

# Optional: activate a venv if you want (uncomment and edit if needed)
# source .venv_ras-csc-switch-mm/bin/activate

# ----------------------------------------------------------------------
# STEP 01–03: Omics processing (R)
# ----------------------------------------------------------------------

echo "[STEP 01/09] Bulk RNA-seq processing (01_process_bulk_rnaseq.R)"
Rscript pipelines/data_processing/01_process_bulk_rnaseq.R
echo "[DONE] Step 01 complete."
echo

echo "[STEP 02/09] scRNA-seq processing (02_process_scrna.R)"
Rscript pipelines/data_processing/02_process_scrna.R
echo "[DONE] Step 02 complete."
echo

echo "[STEP 03/09] ATAC-seq processing (03_process_atacseq.R)"
Rscript pipelines/data_processing/03_process_atacseq.R
echo "[DONE] Step 03 complete."
echo

# ----------------------------------------------------------------------
# STEP 04–09: Model calibration, tests, and predictions (Python)
# ----------------------------------------------------------------------

echo "[STEP 04/09] Model calibration (04_run_model_calibration.py)"
python src/scripts/04_run_model_calibration.py
echo "[DONE] Step 04 complete."
echo

echo "[STEP 05/09] Hypothesis tests (05_run_hypothesis_tests.py)"
python src/scripts/05_run_hypothesis_tests.py
echo "[DONE] Step 05 complete."
echo

echo "[STEP 06/09] Sensitivity analysis (06_run_sensitivity_analysis.py)"
python src/scripts/06_run_sensitivity_analysis.py
echo "[DONE] Step 06 complete."
echo

echo "[STEP 07/09] Phase-plane analysis (07_run_phase_plane_plots.py)"
python src/scripts/07_run_phase_plane_plots.py
echo "[DONE] Step 07 complete."
echo

echo "[STEP 08/09] Temporal analysis (08_run_temporal_analysis.py)"
python src/scripts/08_run_temporal_analysis.py
echo "[DONE] Step 08 complete."
echo

echo "[STEP 09/09] Predictions + falsification (09_run_predictions.py)"
python src/scripts/09_run_predictions.py
echo "[DONE] Step 09 complete."
echo

echo "======================================================================"
echo "Pipeline complete: all 09 steps ran without error."
echo "Check ./results and ./figures/main for outputs."
echo "======================================================================"
