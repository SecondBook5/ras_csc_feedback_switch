#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

in_path  <- "data/processed/rnaseq/module_summary_by_condition.csv"
out_path <- "data/processed/rnaseq/ras_csc_calibration_targets.csv"

if (!file.exists(in_path)) {
  stop("[ERROR] module_summary_by_condition.csv not found at: ", in_path)
}

targets <- readr::read_csv(in_path, show_col_types = FALSE)

keep_rows <- targets %>%
  dplyr::filter(
    (dataset == "Bl6" & condition %in% c("Normal", "SCC")) |
    (dataset == "PAP_SCC" & condition %in% c("Papilloma", "SCC")) |
    (dataset == "PDV" & condition %in% c("PDV_WT", "PDV_LeprKO"))
  )

readr::write_csv(keep_rows, out_path)
message("[OK] Wrote Ras–CSC calibration targets to: ", out_path)
