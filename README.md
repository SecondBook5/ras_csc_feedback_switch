ras_csc_feedback_switch
=======================

Minimal Ras–CSC feedback model and calibration pipeline (Yuan et al. 2022 inspired).

Quick start
-----------

1. Create and activate a Python virtual environment and install dependencies:

```bash
python -m venv .venv
.venv/bin/pip install --upgrade pip
.venv/bin/pip install numpy pandas
```

2. Compute Pap/SCC module scores (requires raw TPM in `data/raw`):

```bash
python src/module_scores_pap_scc.py
```

3. Calibrate model parameters (writes `data/processed/model_calibration_results.json`):

```bash
python src/calibrate_params.py
```

4. Run the model demo (requires calibrated params):

```bash
python src/ras_csc_model.py
```

Notes for contributors
----------------------
- Scripts assume they live under `src/` and use `get_repo_root()` to find `data/` and `config/`.
- `data/raw/` contains large raw files and is git-ignored by default; keep those files out of the repo.
- See `.github/copilot-instructions.md` for guidance aimed at AI coding agents.

How to publish to GitHub
------------------------
Replace `USERNAME` and `REPO` with your GitHub username and repository name:

```bash
git branch -M main
git remote add origin git@github.com:USERNAME/REPO.git
git push -u origin main
```
