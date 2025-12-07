#!/usr/bin/env python3
"""
Compute module scores (CSC program, angiogenesis, LEPR, mTOR) from GSE190411 TPM data.

This script reads the processed bulk RNA-seq TPM matrix from Yuan et al. (Nature 2022, GSE190411),
assigns samples to biological conditions (Normal skin vs Tumor in this Bl6 subset),
and computes simple gene-set module scores as mean log2(TPM + 1) for each module.

The goal here is to obtain approximate "low" (normal) and "high" (tumor) levels for the
CSC, angiogenesis, LEPR, and mTOR modules, which will be used to:
    1) Sanity-check that the chosen gene sets behave in the expected direction.
    2) Provide early anchoring for the Ras–CSC feedback ODE model.
    3) Demonstrate that the model parameters can be tied to real expression differences.

Papilloma vs SCC module scores will be derived later from the separate
'GSE190411_Yuan2021_PAP_SCC_RNAseq_TPM.txt.gz' file; this script is focused
on the Bl6 Norm/Tumor matrix.

The script is deliberately defensive:
    - It checks that the data file exists.
    - It parses gene symbols from the 'transcript_id' field.
    - It drops junk columns and aggregates duplicate symbols.
    - It prints a summary of conditions and module scores to the console.

Args:
    None, all paths are resolved relative to the repo root.

Outputs:
    data/processed/module_scores_bl6_norm_tumor.csv
        Columns: condition, n_samples, C_score, V_score, R_score, M_score
"""

from __future__ import annotations

import os
import sys
from typing import Dict, List

import numpy as np
import pandas as pd


def get_repo_root() -> str:
    """
    Determine the repository root directory based on this script's location.

    The script is expected to live in the 'src/' subdirectory of the project.
    By taking the parent directory of the directory containing this file,
    we obtain the project root in a way that does not depend on the current
    working directory. This improves reproducibility when running the script
    from different shells or IDEs.

    Returns:
        str: Absolute path to the repository root directory.

    Raises:
        RuntimeError: If the computed root directory does not exist.
    """
    # Get absolute path of this script
    script_path: str = os.path.abspath(__file__)
    # Parent directory should be 'src/'
    src_dir: str = os.path.dirname(script_path)
    # One level up is the repo root
    repo_root: str = os.path.dirname(src_dir)

    # Verify that the directory actually exists
    if not os.path.isdir(repo_root):
        raise RuntimeError(
            f"[ERROR] Computed repo root does not exist: {repo_root}. "
            f"Check that the script is located in 'src/' under the project root."
        )

    return repo_root


def load_tpm_matrix() -> pd.DataFrame:
    """
    Load the Bl6 Norm/Tumor TPM matrix for GSE190411 from data/raw.

    The input file 'GSE190411_Yuan2021_Bl6_Norm_Pap_SCC_RNAseq_TPM.csv.gz'
    stores transcript identifiers in a 'transcript_id' column. Each entry
    looks like 'ENSMUSG00000000001_Gnai3'. This function:

        1) Reads the file without presetting an index.
        2) Extracts gene symbols by splitting 'transcript_id' on '_' and
           taking the last field (e.g. 'Gnai3' -> 'GNAI3').
        3) Uppercases the gene symbols and uses them as the DataFrame index.
        4) Drops non-sample columns (including 'transcript_id' and any
           'Unnamed: x' junk columns).
        5) Keeps only numeric columns as TPM samples.
        6) Aggregates duplicate gene symbols by taking the mean across
           transcripts.

    Returns:
        pd.DataFrame: TPM matrix with genes as index (uppercase symbols)
                      and samples as columns.

    Raises:
        FileNotFoundError: If the expected TPM file does not exist.
        RuntimeError: If parsing fails or if there are too few sample columns.
    """
    # Resolve the repository root
    repo_root: str = get_repo_root()

    # Build the full path to the TPM file
    tpm_path: str = os.path.join(
        repo_root,
        "data",
        "raw",
        "GSE190411_Yuan2021_Bl6_Norm_Pap_SCC_RNAseq_TPM.csv.gz",
    )

    # Check for existence of the TPM file
    if not os.path.exists(tpm_path):
        raise FileNotFoundError(
            f"[ERROR] TPM file not found at: {tpm_path}\n"
            f"Download 'GSE190411_Yuan2021_Bl6_Norm_Pap_SCC_RNAseq_TPM.csv.gz' "
            f"from GEO (GSE190411) and place it in data/raw/."
        )

    # Read the raw file without specifying an index
    try:
        raw_df: pd.DataFrame = pd.read_csv(tpm_path)
    except Exception as exc:
        raise RuntimeError(
            f"[ERROR] Failed to read TPM matrix from {tpm_path}: {exc}"
        ) from exc

    # Ensure the 'transcript_id' column is present
    if "transcript_id" not in raw_df.columns:
        raise RuntimeError(
            "[ERROR] Column 'transcript_id' not found in TPM file. "
            "Inspect the file header to confirm the column names."
        )

    # Extract transcript IDs as strings
    transcript_ids = raw_df["transcript_id"].astype(str)

    # Parse gene symbols by splitting on '_' and taking the last segment
    gene_symbols = transcript_ids.str.split("_").str[-1]

    # Uppercase gene symbols to standardise and enable case-insensitive matching
    gene_symbols = gene_symbols.str.upper()

    # Set parsed gene symbols as the DataFrame index
    raw_df.index = gene_symbols

    # Drop the original 'transcript_id' column
    raw_df = raw_df.drop(columns=["transcript_id"])

    # Drop any junk columns that start with 'Unnamed'
    cols_to_drop: List[str] = [c for c in raw_df.columns if c.startswith("Unnamed")]
    if len(cols_to_drop) > 0:
        raw_df = raw_df.drop(columns=cols_to_drop)

    # Keep only numeric columns as TPM samples
    tpm_numeric = raw_df.select_dtypes(include=[np.number])

    # Sanity check: require at least a few sample columns
    if tpm_numeric.shape[1] < 3:
        raise RuntimeError(
            f"[ERROR] After selecting numeric sample columns, only "
            f"{tpm_numeric.shape[1]} columns remain. Please inspect the file."
        )

    # Aggregate duplicate gene symbols by mean across transcripts
    tpm_numeric = tpm_numeric.groupby(tpm_numeric.index).mean()

    # Print summary of the cleaned TPM matrix
    print(
        f"[INFO] TPM matrix after parsing gene symbols: "
        f"{tpm_numeric.shape[0]} genes x {tpm_numeric.shape[1]} samples.",
        file=sys.stderr,
    )
    print(
        f"[INFO] Sample columns: {list(tpm_numeric.columns)}",
        file=sys.stderr,
    )

    return tpm_numeric


def assign_conditions(sample_names: List[str]) -> Dict[str, str]:
    """
    Assign biological condition labels to samples using name patterns.

    For the Bl6 Norm/Tumor TPM file, the sample names look like:
        - '102F_TGFbRIIhet_TetOHras_normal'
        - '100F_TGFbRIIcKO_TetOHras_normal'
        - '47_M_Tumor'
        - '85F_TGFbRIIflfl_TetOHras'
        - '46_M_Tumor'
        - '83F_TGFbRIIfl_TetOHras'

    At this stage we distinguish only between:
        - 'Norm'  : if the name contains 'norm' or 'normal'
        - 'Tumor' : if the name contains 'tumor'
        - 'Other' : all remaining samples

    You can refine these labels later (e.g., by using metadata) or
    handle Pap vs SCC separately using the PAP/SCC TPM file.

    Args:
        sample_names (List[str]): List of sample identifiers (column names).

    Returns:
        Dict[str, str]: Mapping from sample name to condition label.
    """
    condition_map: Dict[str, str] = {}

    for s in sample_names:
        label: str = "Other"
        s_lower: str = s.lower()

        if "norm" in s_lower or "normal" in s_lower:
            label = "Norm"
        elif "tumor" in s_lower:
            label = "Tumor"

        condition_map[s] = label

    # Summarise counts per condition
    counts: Dict[str, int] = {}
    for _, cond in condition_map.items():
        counts[cond] = counts.get(cond, 0) + 1

    print("[INFO] Condition assignment summary:", file=sys.stderr)
    for cond, n in counts.items():
        print(f"  {cond}: {n} samples", file=sys.stderr)

    return condition_map


def _mean_log2_tpm(
    tpm: pd.DataFrame,
    gene_set: List[str],
    sample: str,
) -> float:
    """
    Compute mean log2(TPM + 1) for a given gene set in a single sample.

    This helper extracts TPM values for a specified list of genes in one sample,
    applies log2(TPM + 1) to avoid issues at zero, and returns the arithmetic mean.
    Genes not present in the TPM matrix are ignored (with a warning if all are
    missing).

    Args:
        tpm (pd.DataFrame): TPM matrix with genes as index and samples as columns.
        gene_set (List[str]): List of gene symbols to include in the module.
        sample (str): Column name of the sample in the TPM matrix.

    Returns:
        float: Mean log2(TPM + 1) value over the present genes in the module.

    Raises:
        KeyError: If the specified sample is not a column in the TPM matrix.
    """
    if sample not in tpm.columns:
        raise KeyError(
            f"[ERROR] Sample '{sample}' not found in TPM matrix columns. "
            f"Available columns (first few): {list(tpm.columns)[:5]}"
        )

    # Filter to genes actually present in the TPM index
    present_genes: List[str] = [g for g in gene_set if g in tpm.index]

    if len(present_genes) == 0:
        print(
            f"[WARN] None of the genes in the gene set are present for sample "
            f"{sample}. Returning NaN for this module.",
            file=sys.stderr,
        )
        return float("nan")

    values = tpm.loc[present_genes, sample].astype(float)
    log2_vals = np.log2(values + 1.0)

    return float(log2_vals.mean())


def compute_module_scores(
    tpm: pd.DataFrame,
    condition_map: Dict[str, str],
) -> pd.DataFrame:
    """
    Compute module scores for each sample and aggregate by condition.

    For each sample, this function calculates:
        - C_score: CSC program module
        - V_score: Angiogenesis / vasculature module
        - R_score: LEPR module
        - M_score: mTOR activity module

    Each score is defined as the mean log2(TPM + 1) over a predefined gene set.
    After computing these scores for every sample, they are aggregated by
    biological condition to obtain condition-level averages.

    Args:
        tpm (pd.DataFrame): TPM matrix (genes as rows, samples as columns).
        condition_map (Dict[str, str]): Mapping from sample name to condition label.

    Returns:
        pd.DataFrame: Per-condition module scores with columns:
            condition, n_samples, C_score, V_score, R_score, M_score
    """
    # NOTE: These gene sets are placeholders and should be refined using
    # the Yuan et al. paper and pathway annotations. They are uppercase to
    # match the parsed gene index.
    c_genes: List[str] = [
        "SOX2",
        "KLF4",
        "MYC",
    ]

    v_genes: List[str] = [
        "VEGFA",
        "KDR",
        "ANGPT2",
    ]

    r_genes: List[str] = [
        "LEPR",
    ]

    m_genes: List[str] = [
        "RPS6KB1",
        "EIF4EBP1",
    ]

    rows: List[Dict[str, float]] = []

    for sample in tpm.columns:
        cond: str = condition_map.get(sample, "Other")

        c_score: float = _mean_log2_tpm(tpm, c_genes, sample)
        v_score: float = _mean_log2_tpm(tpm, v_genes, sample)
        r_score: float = _mean_log2_tpm(tpm, r_genes, sample)
        m_score: float = _mean_log2_tpm(tpm, m_genes, sample)

        rows.append(
            {
                "sample": sample,
                "condition": cond,
                "C_score": c_score,
                "V_score": v_score,
                "R_score": r_score,
                "M_score": m_score,
            }
        )

    sample_scores = pd.DataFrame(rows)

    grouped = (
        sample_scores
        .groupby("condition", dropna=False)
        .agg(
            n_samples=("sample", "count"),
            C_score=("C_score", "mean"),
            V_score=("V_score", "mean"),
            R_score=("R_score", "mean"),
            M_score=("M_score", "mean"),
        )
        .reset_index()
    )

    print("[INFO] Per-condition module scores:", file=sys.stderr)
    with pd.option_context("display.max_columns", None):
        print(grouped, file=sys.stderr)

    return grouped


def main() -> None:
    """
    Orchestrate the module score computation for GSE190411 Bl6 Norm/Tumor.

    Steps:
        1. Load the TPM matrix from data/raw and parse gene symbols.
        2. Assign each sample to a biological condition (Norm, Tumor, Other).
        3. Compute module scores for each sample and aggregate by condition.
        4. Write the per-condition scores to
           data/processed/module_scores_bl6_norm_tumor.csv.
        5. Print a brief summary of the results to the console.

    The output CSV will be used to define Normal vs Tumor targets and to
    normalise ODE model variables for the Ras–CSC feedback switch analysis.
    """
    tpm: pd.DataFrame = load_tpm_matrix()
    condition_map: Dict[str, str] = assign_conditions(list(tpm.columns))
    module_scores: pd.DataFrame = compute_module_scores(tpm, condition_map)

    repo_root: str = get_repo_root()
    out_dir: str = os.path.join(repo_root, "data", "processed")
    os.makedirs(out_dir, exist_ok=True)

    out_path: str = os.path.join(out_dir, "module_scores_bl6_norm_tumor.csv")

    try:
        module_scores.to_csv(out_path, index=False)
    except Exception as exc:
        raise RuntimeError(
            f"[ERROR] Failed to write module scores to {out_path}: {exc}"
        ) from exc

    print(f"[OK] Wrote per-condition module scores to: {out_path}")
    print(module_scores.to_string(index=False))


if __name__ == "__main__":
    main()
