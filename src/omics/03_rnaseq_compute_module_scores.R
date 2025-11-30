#!/usr/bin/env Rscript

# ===========================================================
# 03_rnaseq_compute_module_scores.R
#
# Pipeline step: MAIN STEP 2 (module scores)
#
# Purpose:
#   1) Load log2-transformed expression matrices from step 02:
#        - data/interim/rnaseq/bl6_expression.tsv
#        - data/interim/rnaseq/pap_scc_expression.tsv
#        - data/interim/rnaseq/pdv_expression.tsv
#   2) Load sample metadata from step 02:
#        - data/interim/rnaseq/sample_metadata_GSE190411.csv
#   3) Load gene sets:
#        - config/gene_sets_rnaseq.yaml
#   4) For each dataset (Bl6, PAP_SCC, PDV):
#        - Z-score each gene across samples.
#        - Compute module score per sample = mean z-score of genes in set.
#   5) Write:
#        - data/processed/rnaseq/module_scores_by_sample.csv
#
# Notes:
#   - Output feeds directly into model calibration (Python side).
# ===========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(yaml)
})

#------------------------- paths ---------------------------

expr_dir   <- "data/interim/rnaseq"
meta_path  <- file.path(expr_dir, "sample_metadata_GSE190411.csv")
geneset_yml <- "config/gene_sets_rnaseq.yaml"
out_dir    <- "data/processed/rnaseq"

if (!dir.exists(expr_dir)) {
  stop("[ERROR] Expression directory not found: ", expr_dir)
}
if (!file.exists(meta_path)) {
  stop("[ERROR] Metadata file not found: ", meta_path)
}
if (!file.exists(geneset_yml)) {
  stop("[ERROR] Gene set YAML not found: ", geneset_yml)
}
if (!dir.exists(out_dir)) {
  message("[INFO] Creating output directory: ", out_dir)
  dir.create(out_dir, recursive = TRUE)
}

#---------------------- load data --------------------------

message("[STEP] Loading expression matrices")

bl6_expr <- read_tsv(file.path(expr_dir, "bl6_expression.tsv"), show_col_types = FALSE)
pap_expr <- read_tsv(file.path(expr_dir, "pap_scc_expression.tsv"), show_col_types = FALSE)
pdv_expr <- read_tsv(file.path(expr_dir, "pdv_expression.tsv"), show_col_types = FALSE)

meta <- read_csv(meta_path, show_col_types = FALSE)

# Ensure gene_symbol exists
stopifnot("gene_symbol" %in% colnames(bl6_expr))
stopifnot("gene_symbol" %in% colnames(pap_expr))
stopifnot("gene_symbol" %in% colnames(pdv_expr))

#---------------------- load gene sets ---------------------

message("[STEP] Loading gene sets from YAML: ", geneset_yml)
gene_sets <- yaml::read_yaml(geneset_yml)

if (length(gene_sets) == 0) {
  stop("[ERROR] No gene sets found in YAML.")
}

module_names <- names(gene_sets)
message("[INFO] Found gene sets: ", paste(module_names, collapse = ", "))

#---------------- helper: compute module scores ------------

compute_module_scores_dataset <- function(expr_df, dataset_name, gene_sets) {
  # expr_df: tibble with gene_symbol + samples as columns
  # dataset_name: "Bl6", "PAP_SCC", or "PDV"

  message("\n[DATASET] Computing module scores for: ", dataset_name)

  # Long format
  expr_long <- expr_df %>%
    pivot_longer(
      cols = -gene_symbol,
      names_to = "sample_id",
      values_to = "expr"
    )

  # Check that all sample_ids exist in metadata
  unknown_samples <- setdiff(unique(expr_long$sample_id), meta$sample_id)
  if (length(unknown_samples) > 0) {
    stop("[ERROR] Some ", dataset_name, " samples missing in metadata: ",
         paste(unknown_samples, collapse = ", "))
  }

  # Restrict metadata to this dataset
  meta_sub <- meta %>%
    filter(dataset == dataset_name)

  # Z-score per gene across samples in THIS dataset
  expr_long <- expr_long %>%
    group_by(gene_symbol) %>%
    mutate(
      expr_z = (expr - mean(expr, na.rm = TRUE)) /
               sd(expr, na.rm = TRUE)
    ) %>%
    ungroup()

  # If sd is zero, expr_z becomes NA; set back to 0
  expr_long <- expr_long %>%
    mutate(
      expr_z = ifelse(is.na(expr_z), 0, expr_z)
    )

  # Module scores
  module_scores_list <- list()

  for (mod in names(gene_sets)) {
    genes <- gene_sets[[mod]]

    # Restrict to genes present in this dataset
    present_genes <- intersect(genes, unique(expr_long$gene_symbol))

    if (length(present_genes) == 0) {
      warning("[WARN] Module ", mod, " has no genes present in dataset ", dataset_name)
      next
    }

    # Filter and compute mean z per sample
    mod_scores <- expr_long %>%
      filter(gene_symbol %in% present_genes) %>%
      group_by(sample_id) %>%
      summarise(
        score = mean(expr_z, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        module = mod
      )

    module_scores_list[[mod]] <- mod_scores
  }

  if (length(module_scores_list) == 0) {
    stop("[ERROR] No module scores computed for dataset: ", dataset_name)
  }

  scores_all <- bind_rows(module_scores_list)

  # Join metadata for dataset / condition
  scores_all <- scores_all %>%
    left_join(meta_sub, by = "sample_id")

  # Pivot modules wide
  scores_wide <- scores_all %>%
    select(sample_id, dataset, condition, module, score) %>%
    pivot_wider(
      names_from = module,
      values_from = score
    )

  return(scores_wide)
}

#---------------- compute for all datasets -----------------

scores_bl6 <- compute_module_scores_dataset(bl6_expr, "Bl6", gene_sets)
scores_pap <- compute_module_scores_dataset(pap_expr, "PAP_SCC", gene_sets)
scores_pdv <- compute_module_scores_dataset(pdv_expr, "PDV", gene_sets)

module_scores_all <- bind_rows(scores_bl6, scores_pap, scores_pdv)

out_scores <- file.path(out_dir, "module_scores_by_sample.csv")
write_csv(module_scores_all, out_scores)

message("\n[OK] Wrote module scores to: ", out_scores)
message("[DONE] RNA-Seq module scoring completed.\n")
