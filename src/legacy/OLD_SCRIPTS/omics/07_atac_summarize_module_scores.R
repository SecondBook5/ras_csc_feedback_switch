#!/usr/bin/env Rscript

# ===========================================================
# 07_atac_summarize_module_scores.R
#
# Purpose:
#   1) Read ATAC module scores:
#        - data/processed/atac/atac_module_scores_by_sample.csv
#   2) Summarize by condition:
#        - n_samples, mean, sd for each module
#   3) Write:
#        - data/processed/atac/atac_module_summary_by_condition.csv
#
# This mirrors the RNA-Seq summary so you can compare RNA vs ATAC
# for TGFb_bulk, mTOR_bulk, Angio_bulk, CSC_bulk.
# ===========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

#------------------------- paths ---------------------------

in_scores  <- "data/processed/atac/atac_module_scores_by_sample.csv"
out_dir    <- "data/processed/atac"
out_summary <- file.path(out_dir, "atac_module_summary_by_condition.csv")

if (!file.exists(in_scores)) {
  stop("[ERROR] ATAC module scores not found: ", in_scores,
       "\nRun 06_atac_from_narrowPeak_to_module_scores.R first.")
}
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

#---------------------- load data --------------------------

message("[STEP] Loading ATAC module scores from: ", in_scores)
scores <- read_csv(in_scores, show_col_types = FALSE)

# Expect columns: sample_id, dataset, condition, TGFb_bulk, mTOR_bulk, Angio_bulk, CSC_bulk
required_cols <- c("sample_id", "dataset", "condition")
missing_req <- setdiff(required_cols, colnames(scores))
if (length(missing_req) > 0) {
  stop("[ERROR] Missing required columns in ATAC scores: ",
       paste(missing_req, collapse = ", "))
}

# Identify module columns automatically (anything not id cols)
id_cols <- c("sample_id", "dataset", "condition")
module_cols <- setdiff(colnames(scores), id_cols)

if (length(module_cols) == 0L) {
  stop("[ERROR] No module columns found in ATAC scores.")
}

message("[INFO] Detected module columns: ", paste(module_cols, collapse = ", "))

#------------------ summarize by condition -----------------

message("[STEP] Summarizing ATAC module scores by condition")

summary_long <- scores %>%
  tidyr::pivot_longer(
    cols = all_of(module_cols),
    names_to = "module",
    values_to = "score"
  ) %>%
  group_by(condition, module) %>%
  summarise(
    n_samples = n(),
    mean_score = mean(score, na.rm = TRUE),
    sd_score   = sd(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(condition, module)

readr::write_csv(summary_long, out_summary)
message("[OK] Wrote ATAC module summary to: ", out_summary)

message("[DONE] ATAC module summary by condition completed.")
