#!/usr/bin/env Rscript

# ===========================================================
# 08_compare_rna_atac_modules.R
#
# Purpose:
#   1) Read RNA-Seq module summaries:
#        - data/processed/rnaseq/module_summary_by_condition.csv
#   2) Read ATAC-seq module summaries:
#        - data/processed/atac/atac_module_summary_by_condition.csv
#   3) Build a joint long-format table with:
#        - omics_type (RNA / ATAC)
#        - condition
#        - module
#        - n_samples
#        - mean_score
#        - sd_score
#   4) Write:
#        - data/processed/omics/module_summary_rna_atac_combined.csv
#
# Notes:
#   - ATAC n_samples = 1, so sd is NA; that’s fine as long as you
#     treat it as qualitative.
# ===========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

#------------------------- paths ---------------------------

rna_path  <- "data/processed/rnaseq/module_summary_by_condition.csv"
atac_path <- "data/processed/atac/atac_module_summary_by_condition.csv"
out_dir   <- "data/processed/omics"
out_path  <- file.path(out_dir, "module_summary_rna_atac_combined.csv")

if (!file.exists(rna_path)) {
  stop("[ERROR] RNA summary not found: ", rna_path)
}
if (!file.exists(atac_path)) {
  stop("[ERROR] ATAC summary not found: ", atac_path)
}
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

#--------------------- load RNA summary --------------------

message("[STEP] Reading RNA-Seq module summary: ", rna_path)
rna <- read_csv(rna_path, show_col_types = FALSE)

req_rna <- c("dataset", "condition", "module", "n_samples", "mean_score", "sd_score")
missing_rna <- setdiff(req_rna, colnames(rna))
if (length(missing_rna) > 0L) {
  stop("[ERROR] RNA summary missing columns: ",
       paste(missing_rna, collapse = ", "))
}

rna_clean <- rna %>%
  mutate(
    omics_type = "RNA"
  ) %>%
  select(omics_type, dataset, condition, module, n_samples, mean_score, sd_score)

#--------------------- load ATAC summary -------------------

message("[STEP] Reading ATAC module summary: ", atac_path)
atac <- read_csv(atac_path, show_col_types = FALSE)

req_atac <- c("condition", "module", "n_samples", "mean_score", "sd_score")
missing_atac <- setdiff(req_atac, colnames(atac))
if (length(missing_atac) > 0L) {
  stop("[ERROR] ATAC summary missing columns: ",
       paste(missing_atac, collapse = ", "))
}

atac_clean <- atac %>%
  mutate(
    omics_type = "ATAC",
    dataset    = "ATAC"
  ) %>%
  select(omics_type, dataset, condition, module, n_samples, mean_score, sd_score)

#----------------------- combine ---------------------------

combined <- bind_rows(rna_clean, atac_clean) %>%
  arrange(module, omics_type, condition)

if (nrow(combined) == 0L) {
  stop("[ERROR] Combined RNA+ATAC table is empty.")
}

message("[STEP] Writing combined RNA+ATAC summary to: ", out_path)
write_csv(combined, out_path)

message("[DONE] Built combined RNA+ATAC module summary.")
