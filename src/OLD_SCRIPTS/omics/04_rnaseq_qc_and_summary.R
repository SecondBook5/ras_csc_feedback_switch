#!/usr/bin/env Rscript

#-----------------------------------------------------------
# 04_rnaseq_qc_and_summary.R
#
# Goal:
#   1) QC: visualize module scores by dataset x condition.
#   2) Summarize to condition-level means/SD:
#        data/processed/rnaseq/module_summary_by_condition.csv
#
# This is the bridge between RNA-seq and the Ras–CSC model:
#   - Model sees per-condition module targets, not raw samples.
#-----------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

in_path  <- "data/processed/rnaseq/module_scores_by_sample.csv"
out_dir  <- "data/processed/rnaseq"

if (!file.exists(in_path)) {
  stop("[ERROR] module_scores_by_sample.csv not found at: ", in_path)
}
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

message("[STEP] Reading module scores: ", in_path)
scores <- readr::read_csv(in_path, show_col_types = FALSE)

required_cols <- c("sample_id", "dataset", "condition",
                   "TGFb_bulk", "mTOR_bulk", "Angio_bulk", "CSC_bulk")
missing_cols <- setdiff(required_cols, colnames(scores))
if (length(missing_cols) > 0) {
  stop("[ERROR] Missing columns in module_scores_by_sample.csv: ",
       paste(missing_cols, collapse = ", "))
}

#-----------------------------------------------------------
# 1) QC plots
#-----------------------------------------------------------

message("[STEP] Building QC plots for module scores")

scores_long <- scores %>%
  tidyr::pivot_longer(
    cols = c(TGFb_bulk, mTOR_bulk, Angio_bulk, CSC_bulk),
    names_to = "module",
    values_to = "score"
  )

# Order conditions roughly along benign -> invasive axis
scores_long <- scores_long %>%
  dplyr::mutate(
    condition = factor(
      condition,
      levels = c("Normal", "Papilloma", "SCC", "PDV_WT", "PDV_LeprKO")
    ),
    module = factor(
      module,
      levels = c("TGFb_bulk", "Angio_bulk", "CSC_bulk", "mTOR_bulk")
    )
  )

qc_plot_path <- file.path(out_dir, "qc_module_scores_by_condition.png")

p_qc <- ggplot(scores_long,
               aes(x = condition, y = score, fill = condition)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.7) +
  facet_grid(module ~ dataset, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Module scores by dataset and condition",
    x = "Condition",
    y = "Module z-score (within dataset)"
  )

ggsave(qc_plot_path, p_qc, width = 10, height = 7, dpi = 300)
message("[OK] Wrote QC plot to: ", qc_plot_path)

#-----------------------------------------------------------
# 2) Condition-level summary for model calibration
#-----------------------------------------------------------

message("[STEP] Computing condition-level module summaries")

summary_df <- scores_long %>%
  group_by(dataset, condition, module) %>%
  summarize(
    n_samples = dplyr::n(),
    mean_score = mean(score, na.rm = TRUE),
    sd_score   = sd(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(dataset, condition, module)

summary_path <- file.path(out_dir, "module_summary_by_condition.csv")
readr::write_csv(summary_df, summary_path)
message("[OK] Wrote module summary to: ", summary_path)

message("[DONE] QC + summary of RNA-seq module scores complete.")
