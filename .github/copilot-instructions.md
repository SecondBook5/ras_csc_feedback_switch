# Copilot / AI agent instructions for ras_csc_feedback_switch

Purpose
- Help an AI coding agent become productive quickly in this repo by documenting the "why", key files, and common workflows.

Quick architecture & intent
- Core model: `src/ras_csc_model.py` implements a minimal 5-variable ODE model (C, V, L, R, M) and a small RK4 simulator. It exposes:
  - `RasCSCParams` dataclass (all kinetic parameters)
  - `ras_csc_rhs(t, y, params)` — ODE RHS; state vector order: `[C, V, L, R, M]`
  - `simulate_trajectory(...)` — fixed-step RK4 (exploratory use)
  - `load_calibrated_params()` — loads calibrated parameter sets from `data/processed/model_calibration_results.json`

- Calibration pipeline:
  1. `src/module_scores_pap_scc.py` reads raw TPM (`data/raw/GSE190411_*.txt.gz`) and writes `data/processed/module_scores_pap_scc.csv` (group-level Pap vs SCC mpos scores).
  2. `src/calibrate_params.py` reads `module_scores_pap_scc.csv` and writes `data/processed/model_calibration_results.json` (Model A / Model B parameters + targets).
  3. `src/ras_csc_model.py` expects that JSON to exist before running its demo or being imported.

Key files to inspect
- `src/ras_csc_model.py` — model implementation, parameter dataclass, CLI demo.
- `src/calibrate_params.py` — builds parameter dicts from module scores.
- `src/module_scores_pap_scc.py` — computes `module_scores_pap_scc.csv` from raw TPM.
- `config/model_params.yaml` — placeholder for manual parameter sets; currently empty.
- `data/processed/` — runtime artifacts: `module_scores_pap_scc.csv`, `model_calibration_results.json`.
- `docs/equations.md` — human-readable model equations (partial).

Important conventions & patterns
- Path resolution: scripts use a `get_repo_root()` helper that assumes each script lives under `src/` and computes paths relative to that file. This means scripts can be invoked from any CWD and will still find `data/` and other repo files.
- Data-first workflow: processing raw data (`data/raw`) → compute module scores (`data/processed/module_scores_pap_scc.csv`) → calibrate (`data/processed/model_calibration_results.json`) → simulate/analysis. Many scripts will raise helpful errors if upstream files are missing.
- Parameters are portable dataclasses/dicts: `calibrate_params.py` writes plain dicts; `ras_csc_model.py` converts them into `RasCSCParams` with float casting and an optional `clip_state` boolean flag (calibration may store `clip_state` as 0/1).
- State ordering and units: model is dimensionless; maps to data via linear offsets/scales (`Pap`->0, `SCC`->1) produced by `calibrate_params.py`.

Dependencies & environment
- Minimal Python packages required (discoverable from imports): `numpy`, `pandas`.
- There is no `requirements.txt` in repo root. Typical quick setup:

```bash
python -m venv .venv
.venv/bin/pip install --upgrade pip
.venv/bin/pip install numpy pandas
```

Common commands / workflows
- Compute module scores (produces `data/processed/module_scores_pap_scc.csv`):

```bash
python src/module_scores_pap_scc.py
```

- Calibrate parameters (produces `data/processed/model_calibration_results.json`):

```bash
python src/calibrate_params.py
```

- Run the model demo (requires calibrated JSON):

```bash
python src/ras_csc_model.py
```

- Programmatic usage example (from repo root):

```python
from src.ras_csc_model import load_calibrated_params, ras_csc_rhs, simulate_trajectory
params_A, params_B = load_calibrated_params()
# call simulate_trajectory(ras_csc_rhs, ...) using params_A
```

Project-specific gotchas & notes for agents (must-read)
- JSON key mismatch detected: `ras_csc_model.load_calibrated_params()` expects the calibration JSON to contain keys named `model_A_params` and `model_B_params`, but `src/calibrate_params.py` currently writes keys `model_A` and `model_B`. This will cause `load_calibrated_params()` to raise a `KeyError` unless the JSON is adjusted or one of the scripts is patched. Two remediation options:
  - Update `src/calibrate_params.py` to write `model_A_params` / `model_B_params` instead of `model_A` / `model_B` (small, local change).
  - Or update `src/ras_csc_model.py`'s `load_calibrated_params()` to look for `model_A` / `model_B` (or both variants) for backward compatibility.

- Upstream data requirement: `src/module_scores_pap_scc.py` expects `data/raw/GSE190411_Yuan2021_PAP_SCC_RNAseq_TPM.txt.gz`. If that file is not present, the script exits with a clear message instructing to download the GEO file.

- Numeric coercions: `ras_csc_model` aggressively casts JSON numeric values to `float` and tolerates `clip_state` as 0/1 or boolean — follow that format when creating or editing parameter JSON.

Examples to reference in code
- State vector order in `ras_csc_model.py`:
  - `y[0] = C`, `y[1] = V`, `y[2] = L`, `y[3] = R`, `y[4] = M`.
- Calibration choices in `calibrate_params.py`: many decay rates are set to `1.0` to fix time scales; `u_TGFb=1.0` is used for SCC.

When making edits
- Preserve the repo-relative path pattern used by `get_repo_root()` (scripts assume `src/` layout).
- If touching the JSON contract between calibrator and model, update both `calibrate_params.py` and `ras_csc_model.py` together or make `load_calibrated_params()` tolerant of both key names.
- Add a `requirements.txt` if you add more dependencies; keep installs minimal.

What I could not discover automatically
- CI/test commands (no tests or CI config present). If you need unit tests or CI, add a minimal `pytest` harness and `github/workflows` entries.
- Exact gene sets intended for modules — `module_scores_pap_scc.py` contains placeholder gene lists; confirm with domain experts or Yuan et al. supplementary material.

Ask me to refine
- Would you like me to (a) patch `calibrate_params.py` to write `model_A_params` / `model_B_params`, (b) make `ras_csc_model.load_calibrated_params()` accept both keys, or (c) add a `requirements.txt` and a short `Makefile` with the common commands above? Provide preference and I'll implement it.
