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
  library(yaml)
})

#------------------------- paths ---------------------------

expr_dir <- "data/interim/rnaseq"
meta_path <- file.path(expr_dir, "sample_metadata_GSE190411.csv")
geneset_yml <- "config/gene_sets_rnaseq.yaml"
out_dir <- "data/processed/rnaseq"

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

# basic sanity
stopifnot("gene_symbol" %in% colnames(bl6_expr))
stopifnot("gene_symbol" %in% colnames(pap_expr))
stopifnot("gene_symbol" %in% colnames(pdv_expr))
stopifnot(all(c("sample_id", "dataset", "condition") %in% colnames(meta)))

#---------------------- load gene sets ---------------------

message("[STEP] Loading gene sets from YAML: ", geneset_yml)
gene_sets <- yaml::read_yaml(geneset_yml)

if (!is.list(gene_sets) || length(gene_sets) == 0L) {
  stop("[ERROR] No gene sets found in YAML or YAML did not parse as list.")
}

module_names <- names(gene_sets)
message("[INFO] Found gene sets: ", paste(module_names, collapse = ", "))

#---------------- helper: compute module scores ------------

compute_module_scores_dataset <- function(expr_df, dataset_name, gene_sets) {
  message("\n[DATASET] Computing module scores for: ", dataset_name)

  # long format
  expr_long <- expr_df %>%
    pivot_longer(
      cols = -gene_symbol,
      names_to = "sample_id",
      values_to = "expr"
    )

  # metadata restricted to this dataset
  meta_sub <- meta %>%
    filter(dataset == dataset_name)

  if (nrow(meta_sub) == 0L) {
    stop("[ERROR] No metadata rows for dataset ", dataset_name, ".")
  }

  # enforce 1:1 mapping between expression sample IDs and dataset-specific metadata
  missing_in_meta <- setdiff(unique(expr_long$sample_id), meta_sub$sample_id)
  if (length(missing_in_meta) > 0L) {
    stop(
      "[ERROR] Dataset ", dataset_name,
      " has expression samples with no metadata rows: ",
      paste(missing_in_meta, collapse = ", ")
    )
  }

  # also check for meta samples with no expression (not fatal but warn)
  missing_in_expr <- setdiff(meta_sub$sample_id, unique(expr_long$sample_id))
  if (length(missing_in_expr) > 0L) {
    warning(
      "[WARN] Dataset ", dataset_name,
      " has metadata samples with no expression columns: ",
      paste(missing_in_expr, collapse = ", ")
    )
    # drop those from meta_sub so they don't appear later with all-NA modules
    meta_sub <- meta_sub %>%
      filter(sample_id %in% unique(expr_long$sample_id))
  }

  # z-score per gene across samples within this dataset
  expr_long <- expr_long %>%
    group_by(gene_symbol) %>%
    mutate(
      expr_z = (expr - mean(expr, na.rm = TRUE)) /
        sd(expr, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      expr_z = ifelse(is.na(expr_z), 0, expr_z)
    )

  module_scores_list <- list()

  for (mod in names(gene_sets)) {
    genes <- gene_sets[[mod]]

    # genes present in this dataset
    present_genes <- intersect(genes, unique(expr_long$gene_symbol))

    if (length(present_genes) == 0L) {
      warning("[WARN] Module ", mod, " has no genes present in dataset ", dataset_name)
      next
    }

    mod_scores <- expr_long %>%
      filter(gene_symbol %in% present_genes) %>%
      group_by(sample_id) %>%
      summarise(
        score = mean(expr_z, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(module = mod)

    module_scores_list[[mod]] <- mod_scores
  }

  if (length(module_scores_list) == 0L) {
    stop("[ERROR] No module scores computed for dataset: ", dataset_name)
  }

  scores_all <- bind_rows(module_scores_list)

  # attach dataset/condition
  scores_all <- scores_all %>%
    left_join(meta_sub, by = "sample_id")

  if (any(is.na(scores_all$dataset)) || any(is.na(scores_all$condition))) {
    stop(
      "[ERROR] After join, some rows in dataset ", dataset_name,
      " have NA dataset/condition. Check metadata."
    )
  }

  # wide modules
  scores_wide <- scores_all %>%
    select(sample_id, dataset, condition, module, score) %>%
    pivot_wider(
      names_from  = module,
      values_from = score
    )

  scores_wide
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
