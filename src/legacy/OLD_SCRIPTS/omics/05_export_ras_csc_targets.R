#!/usr/bin/env Rscript

#-----------------------------------------------------------
# 05_export_ras_csc_targets.R
#
# Purpose:
#   Take the condition-level RNA module summary
#     data/processed/rnaseq/module_summary_by_condition.csv
#   and construct a *wide* calibration target table for the
#   Ras–CSC ODE model:
#
#     dataset, condition, C_target, A_target, T_target, M_target
#
#   where:
#     C_target = CSC_bulk mean score
#     A_target = Angio_bulk mean score
#     T_target = TGFb_bulk mean score
#     M_target = mTOR_bulk mean score
#
# Output:
#   data/processed/rnaseq/ras_csc_calibration_targets.csv
#
# This is what evaluate_ras_csc_fit.py expects.
#-----------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

#---------- paths ----------

root_dir <- "."
in_path <- file.path(
  root_dir,
  "data", "processed", "rnaseq",
  "module_summary_by_condition.csv"
)
out_path <- file.path(
  root_dir,
  "data", "processed", "rnaseq",
  "ras_csc_calibration_targets.csv"
)

if (!file.exists(in_path)) {
  stop(
    "[ERROR] Input summary not found: ", in_path,
    "\nRun 04_rnaseq_qc_and_summary.R first."
  )
}

message("[STEP] Reading RNA module summary from: ", in_path)
summary_df <- readr::read_csv(in_path, show_col_types = FALSE)

required_cols <- c("dataset", "condition", "module", "mean_score")
missing_cols <- setdiff(required_cols, colnames(summary_df))

if (length(missing_cols) > 0) {
  stop(
    "[ERROR] module_summary_by_condition.csv is missing columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# We only care about the four Ras–CSC modules
needed_modules <- c("TGFb_bulk", "mTOR_bulk", "Angio_bulk", "CSC_bulk")

present_modules <- unique(summary_df$module)
if (!all(needed_modules %in% present_modules)) {
  missing_mods <- setdiff(needed_modules, present_modules)
  stop(
    "[ERROR] Some required modules are missing from summary: ",
    paste(missing_mods, collapse = ", ")
  )
}

# Filter to the relevant modules only
summary_filt <- summary_df %>%
  dplyr::filter(module %in% needed_modules)

# Pivot to wide format:
#   TGFb_bulk, mTOR_bulk, Angio_bulk, CSC_bulk columns
calib_wide <- summary_filt %>%
  dplyr::select(dataset, condition, module, mean_score) %>%
  tidyr::pivot_wider(
    names_from  = module,
    values_from = mean_score
  )

# Sanity check: no NAs in the module columns
for (mod in needed_modules) {
  if (any(is.na(calib_wide[[mod]]))) {
    bad_rows <- calib_wide %>%
      dplyr::filter(is.na(.data[[mod]]))
    msg <- paste0(
      "[ERROR] NA detected in module '", mod, "' for some rows.\n",
      "Offending rows:\n",
      paste(
        paste0(
          "  dataset=", bad_rows$dataset,
          ", condition=", bad_rows$condition
        ),
        collapse = "\n"
      )
    )
    stop(msg)
  }
}

# Map module means to model targets
# C ↔ CSC, A ↔ Angio, T ↔ TGFb, M ↔ mTOR
calib_targets <- calib_wide %>%
  dplyr::mutate(
    C_target = CSC_bulk,
    A_target = Angio_bulk,
    T_target = TGFb_bulk,
    M_target = mTOR_bulk
  ) %>%
  # Keep both the model targets and the original module columns
  dplyr::select(
    dataset, condition,
    C_target, A_target, T_target, M_target,
    TGFb_bulk, mTOR_bulk, Angio_bulk, CSC_bulk
  )

message("[STEP] Writing calibration targets to: ", out_path)
readr::write_csv(calib_targets, out_path)

message("[DONE] Exported Ras–CSC calibration targets.")
