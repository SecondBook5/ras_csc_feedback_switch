#!/usr/bin/env python3
"""
Compute module scores (CSC program, angiogenesis, LEPR, mTOR) from
GSE190411 PAP/SCC TPM data.

This script reads the processed bulk RNA-seq TPM matrix
'GSE190411_Yuan2021_PAP_SCC_RNAseq_TPM.txt.gz' from Yuan et al.
(Nature 2022, GSE190411), assigns samples to biological conditions
(Papilloma vs SCC, and mneg vs mpos subsets), and computes simple
gene-set module scores as mean log2(TPM + 1) for each module.

The goals are to:
    1) Obtain approximate "low" (papilloma) and "high" (SCC) levels
       for angiogenesis / vasculature (V), LEPR (R), and mTOR (M).
    2) Optionally distinguish mneg vs mpos subsets to inform the CSC
       program variable C(t).
    3) Use these condition-level means as anchors to scale variables
       and steady states in the Ras–CSC feedback ODE model.

The script is defensive:
    - Checks that the data file exists.
    - Verifies the expected 'gene' column is present.
    - Drops any non-numeric columns beyond the gene symbol column.
    - Aggregates duplicate gene symbols (if any) by mean.
    - Prints a summary of conditions and module scores to stderr.

Output:
    data/processed/module_scores_pap_scc.csv

Columns:
    condition      : 'Pap' or 'SCC'
    stem_subset    : 'mneg', 'mpos', or 'other'
    n_samples      : number of samples in that group
    C_score        : CSC program module mean log2(TPM + 1)
    V_score        : Angiogenesis/vasculature module mean log2(TPM + 1)
    R_score        : LEPR module mean log2(TPM + 1)
    M_score        : mTOR module mean log2(TPM + 1)
"""

from __future__ import annotations

# Import standard libraries
import os
import sys
from typing import Dict, List

# Import third-party libraries
import numpy as np
import pandas as pd


def get_repo_root() -> str:
    """
    Determine the repository root directory based on this script's location.

    The script is expected to live in the 'src/' subdirectory of the project.
    By taking the parent directory of the directory containing this file,
    we obtain the project root in a way that does not depend on the current
    working directory.

    Returns:
        str: Absolute path to the repository root directory.

    Raises:
        RuntimeError: If the computed root directory does not exist.
    """
    # Get absolute path to this script
    script_path: str = os.path.abspath(__file__)
    # Get the directory that contains this script (expected 'src/')
    src_dir: str = os.path.dirname(script_path)
    # One level up should be the repo root
    repo_root: str = os.path.dirname(src_dir)
    # Check that directory actually exists
    if not os.path.isdir(repo_root):
        # Raise error if layout is not as expected
        raise RuntimeError(
            f"[ERROR] Computed repo root does not exist: {repo_root}. "
            f"Check that the script is located in 'src/' under the project root."
        )
    # Return verified root
    return repo_root


def load_tpm_matrix_pap_scc() -> pd.DataFrame:
    """
    Load the PAP/SCC TPM matrix for GSE190411 from data/raw.

    The input file 'GSE190411_Yuan2021_PAP_SCC_RNAseq_TPM.txt.gz'
    has the following structure (tab-separated):

        gene    PAPmneg1  PAPmneg2  PAPmpos1  PAPmpos2
                SCCmneg1 SCCmneg2 SCCmpos1 SCCmpos2

    where:
        - 'gene' holds gene symbols (e.g. 'SOX2', 'VEGFA', etc.).
        - Remaining columns are TPM values for each sample.

    This function:
        1) Checks that the file exists.
        2) Reads the table as tab-separated.
        3) Sets the 'gene' column as index after uppercasing symbols.
        4) Keeps only numeric columns as samples.
        5) Aggregates duplicate gene symbols by taking the mean.

    Returns:
        pd.DataFrame: TPM matrix with gene symbols (uppercase) as index
                      and samples as columns.

    Raises:
        FileNotFoundError: If the file is missing.
        RuntimeError: If parsing fails or the structure is not as expected.
    """
    # Resolve repo root
    repo_root: str = get_repo_root()
    # Build path to TPM file
    tpm_path: str = os.path.join(
        repo_root,
        "data",
        "raw",
        "GSE190411_Yuan2021_PAP_SCC_RNAseq_TPM.txt.gz",
    )
    # Check existence
    if not os.path.exists(tpm_path):
        # Raise if TPM file missing
        raise FileNotFoundError(
            f"[ERROR] PAP/SCC TPM file not found at: {tpm_path}\n"
            f"Download 'GSE190411_Yuan2021_PAP_SCC_RNAseq_TPM.txt.gz' "
            f"from GEO (GSE190411) and place it in data/raw/."
        )
    # Attempt to read file as tab-separated
    try:
        raw_df: pd.DataFrame = pd.read_csv(tpm_path, sep="\t")
    except Exception as exc:
        # Wrap any parsing exception
        raise RuntimeError(
            f"[ERROR] Failed to read PAP/SCC TPM matrix from {tpm_path}: {exc}"
        ) from exc
    # Ensure expected 'gene' column exists
    if "gene" not in raw_df.columns:
        # Raise if 'gene' column not present
        raise RuntimeError(
            "[ERROR] Column 'gene' not found in PAP/SCC TPM file. "
            "Inspect the file header to confirm the column names."
        )
    # Extract gene symbols as strings
    gene_symbols = raw_df["gene"].astype(str)
    # Uppercase gene symbols for consistent matching
    gene_symbols = gene_symbols.str.upper()
    # Set uppercase gene symbols as DataFrame index
    raw_df.index = gene_symbols
    # Drop the original 'gene' column
    raw_df = raw_df.drop(columns=["gene"])
    # Keep only numeric columns as TPM samples
    tpm_numeric = raw_df.select_dtypes(include=[np.number])
    # Sanity check: require at least a few sample columns
    if tpm_numeric.shape[1] < 4:
        # Raise if too few samples
        raise RuntimeError(
            f"[ERROR] PAP/SCC TPM matrix has only {tpm_numeric.shape[1]} "
            f"sample columns. Inspect the file structure."
        )
    # Group duplicates (if any) by gene symbol and average TPM
    tpm_numeric = tpm_numeric.groupby(tpm_numeric.index).mean()
    # Print summary to stderr
    print(
        f"[INFO] PAP/SCC TPM matrix: {tpm_numeric.shape[0]} genes x "
        f"{tpm_numeric.shape[1]} samples.",
        file=sys.stderr,
    )
    print(
        f"[INFO] Sample columns: {list(tpm_numeric.columns)}",
        file=sys.stderr,
    )
    # Return cleaned TPM matrix
    return tpm_numeric


def assign_conditions_pap_scc(
    sample_names: List[str],
) -> pd.DataFrame:
    """
    Assign biological conditions (Pap vs SCC) and stem subsets (mneg vs mpos).

    The PAP/SCC TPM file uses column names like:
        - 'PAPmneg1', 'PAPmneg2'
        - 'PAPmpos1', 'PAPmpos2'
        - 'SCCmneg1', 'SCCmneg2'
        - 'SCCmpos1', 'SCCmpos2'

    This function parses each sample name and assigns:
        - condition   : 'Pap' if name starts with 'PAP',
                        'SCC' if name starts with 'SCC',
                        'Other' otherwise.
        - stem_subset : 'mneg' if 'mneg' in name,
                        'mpos' if 'mpos' in name,
                        'other' otherwise.

    It returns a small DataFrame mapping each sample to its condition
    and stem subset, and prints a summary to stderr.

    Args:
        sample_names (List[str]): List of sample column names.

    Returns:
        pd.DataFrame: DataFrame with columns
            'sample', 'condition', 'stem_subset'.
    """
    # Prepare containers for each attribute
    records: List[Dict[str, str]] = []
    # Iterate over sample names
    for s in sample_names:
        # Default condition
        condition: str = "Other"
        # Default subset
        stem_subset: str = "other"
        # Lowercase name for pattern matching
        s_lower: str = s.lower()
        # Determine condition
        if s_lower.startswith("pap"):
            condition = "Pap"
        elif s_lower.startswith("scc"):
            condition = "SCC"
        # Determine stem subset
        if "mneg" in s_lower:
            stem_subset = "mneg"
        elif "mpos" in s_lower:
            stem_subset = "mpos"
        # Append record
        records.append(
            {
                "sample": s,
                "condition": condition,
                "stem_subset": stem_subset,
            }
        )
    # Convert to DataFrame
    cond_df = pd.DataFrame.from_records(records)
    # Summarise counts for quick inspection
    summary = (
        cond_df.groupby(["condition", "stem_subset"], dropna=False)
        .size()
        .reset_index(name="n_samples")
    )
    # Print summary to stderr
    print("[INFO] PAP/SCC condition assignment summary:", file=sys.stderr)
    with pd.option_context("display.max_rows", None, "display.max_columns", None):
        print(summary, file=sys.stderr)
    # Return mapping DataFrame
    return cond_df


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
        gene_set (List[str]): List of gene symbols to include in the module
                              (must be uppercase to match the index).
        sample (str): Column name of the sample in the TPM matrix.

    Returns:
        float: Mean log2(TPM + 1) value over the present genes in the module.

    Raises:
        KeyError: If the specified sample is not a column in the TPM matrix.
    """
    # Check that sample exists
    if sample not in tpm.columns:
        # Raise helpful error if sample missing
        raise KeyError(
            f"[ERROR] Sample '{sample}' not found in PAP/SCC TPM matrix. "
            f"Available columns (first few): {list(tpm.columns)[:5]}"
        )
    # Filter gene set by those present in index
    present_genes: List[str] = [g for g in gene_set if g in tpm.index]
    # If none present, warn and return NaN
    if len(present_genes) == 0:
        print(
            f"[WARN] None of the genes in the gene set are present for sample "
            f"{sample}. Returning NaN for this module.",
            file=sys.stderr,
        )
        return float("nan")
    # Extract TPM values for present genes in chosen sample
    values = tpm.loc[present_genes, sample].astype(float)
    # Compute log2(TPM + 1)
    log2_vals = np.log2(values + 1.0)
    # Return mean of transformed values
    return float(log2_vals.mean())


def compute_module_scores_pap_scc(
    tpm: pd.DataFrame,
    cond_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Compute module scores for each PAP/SCC sample and aggregate by condition.

    For each sample, this function calculates:
        - C_score: CSC program module
        - V_score: Angiogenesis / vasculature module
        - R_score: LEPR module
        - M_score: mTOR activity module

    Each score is defined as the mean log2(TPM + 1) over a predefined gene set.
    After computing these scores for every sample, they are aggregated by
    biological condition ('Pap' vs 'SCC') and stem subset ('mneg' vs 'mpos')
    to obtain group-level averages.

    Args:
        tpm (pd.DataFrame): TPM matrix (genes as rows, samples as columns).
        cond_df (pd.DataFrame): DataFrame mapping samples to condition and
                                stem subset. Must have columns:
                                'sample', 'condition', 'stem_subset'.

    Returns:
        pd.DataFrame: Per-group module scores with columns:
            condition, stem_subset, n_samples,
            C_score, V_score, R_score, M_score
    """
    # Define gene sets (uppercase to match TPM index)
    # NOTE: These are placeholders and should be refined using Yuan et al.
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
    # Prepare list for per-sample scores
    rows: List[Dict[str, float]] = []
    # Iterate over samples
    for _, row in cond_df.iterrows():
        # Extract sample ID
        sample: str = str(row["sample"])
        # Extract condition
        condition: str = str(row["condition"])
        # Extract stem subset
        stem_subset: str = str(row["stem_subset"])
        # Compute module scores for this sample
        c_score: float = _mean_log2_tpm(tpm, c_genes, sample)
        v_score: float = _mean_log2_tpm(tpm, v_genes, sample)
        r_score: float = _mean_log2_tpm(tpm, r_genes, sample)
        m_score: float = _mean_log2_tpm(tpm, m_genes, sample)
        # Append record
        rows.append(
            {
                "sample": sample,
                "condition": condition,
                "stem_subset": stem_subset,
                "C_score": c_score,
                "V_score": v_score,
                "R_score": r_score,
                "M_score": m_score,
            }
        )
    # Convert to DataFrame
    sample_scores = pd.DataFrame.from_records(rows)
    # Group by condition and stem subset, compute means and counts
    grouped = (
        sample_scores
        .groupby(["condition", "stem_subset"], dropna=False)
        .agg(
            n_samples=("sample", "count"),
            C_score=("C_score", "mean"),
            V_score=("V_score", "mean"),
            R_score=("R_score", "mean"),
            M_score=("M_score", "mean"),
        )
        .reset_index()
    )
    # Print summary to stderr
    print("[INFO] PAP/SCC per-group module scores:", file=sys.stderr)
    with pd.option_context("display.max_columns", None):
        print(grouped, file=sys.stderr)
    # Return grouped scores
    return grouped


def main() -> None:
    """
    Main entry point for PAP/SCC module score computation.

    Steps:
        1. Load PAP/SCC TPM matrix from data/raw and parse gene symbols.
        2. Assign each sample to a biological condition (Pap vs SCC) and
           stem subset (mneg vs mpos).
        3. Compute module scores for each sample and aggregate by group.
        4. Write per-group scores to
           data/processed/module_scores_pap_scc.csv.
        5. Print a brief summary to stdout.
    """
    # Load PAP/SCC TPM matrix
    tpm: pd.DataFrame = load_tpm_matrix_pap_scc()
    # Assign conditions and stem subsets
    cond_df: pd.DataFrame = assign_conditions_pap_scc(list(tpm.columns))
    # Compute per-group module scores
    module_scores: pd.DataFrame = compute_module_scores_pap_scc(tpm, cond_df)
    # Resolve output directory
    repo_root: str = get_repo_root()
    out_dir: str = os.path.join(repo_root, "data", "processed")
    # Ensure directory exists
    os.makedirs(out_dir, exist_ok=True)
    # Build output path
    out_path: str = os.path.join(out_dir, "module_scores_pap_scc.csv")
    # Attempt to write CSV
    try:
        module_scores.to_csv(out_path, index=False)
    except Exception as exc:
        # Wrap write error
        raise RuntimeError(
            f"[ERROR] Failed to write PAP/SCC module scores to {out_path}: {exc}"
        ) from exc
    # Print success message and table to stdout
    print(f"[OK] Wrote PAP/SCC module scores to: {out_path}")
    print(module_scores.to_string(index=False))


if __name__ == "__main__":
    main()
