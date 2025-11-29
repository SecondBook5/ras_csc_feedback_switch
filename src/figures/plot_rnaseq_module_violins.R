#!/usr/bin/env Rscript

#===========================================================
# plot_rnaseq_module_violins.R
#
# Violin + box + jitter panels for module scores by
# dataset and condition.
#
# Input:
#   data/processed/rnaseq/module_scores_by_sample.csv
#
# Output:
#   figures/main/fig2_rnaseq_module_violins.pdf
#   figures/main/fig2_rnaseq_module_violins.png
#===========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

scores_path <- "data/processed/rnaseq/module_scores_by_sample.csv"
out_dir     <- "figures/main"

if (!file.exists(scores_path)) {
  stop("[ERROR] module_scores_by_sample.csv not found at: ", scores_path)
}
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

scores <- readr::read_csv(scores_path, show_col_types = FALSE)

required_cols <- c("sample_id", "dataset", "condition")
module_cols   <- setdiff(colnames(scores), required_cols)

if (!all(required_cols %in% colnames(scores))) {
  stop(
    "[ERROR] module_scores_by_sample.csv missing: ",
    paste(setdiff(required_cols, colnames(scores)), collapse = ", ")
  )
}
if (length(module_cols) == 0L) {
  stop("[ERROR] No module columns in module_scores_by_sample.csv.")
}

scores <- scores %>%
  mutate(
    condition = factor(
      condition,
      levels = c("Normal", "Papilloma", "SCC", "PDV_WT", "PDV_LeprKO")
    ),
    dataset = factor(
      dataset,
      levels = c("Bl6", "PAP_SCC", "PDV")
    )
  ) %>%
  filter(!is.na(condition))

# module labels
module_labels <- c(
  TGFb_bulk  = "TGFb module",
  mTOR_bulk  = "mTOR module",
  Angio_bulk = "Angiogenesis module",
  CSC_bulk   = "CSC module"
)

condition_palette <- c(
  "Normal"     = "#1f78b4",
  "Papilloma"  = "#33a02c",
  "SCC"        = "#e31a1c",
  "PDV_WT"     = "#6a3d9a",
  "PDV_LeprKO" = "#b15928"
)

theme_ras_csc <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(linewidth = 0.2, colour = "grey88"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(fill = NA, colour = "grey70", linewidth = 0.4),
      axis.title       = ggplot2::element_text(face = "bold"),
      axis.text        = ggplot2::element_text(colour = "black"),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      strip.text       = ggplot2::element_text(face = "bold"),
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle    = ggplot2::element_text(hjust = 0)
    )
}

scores_long <- scores %>%
  tidyr::pivot_longer(
    cols      = all_of(module_cols),
    names_to  = "module",
    values_to = "score"
  ) %>%
  mutate(
    module = factor(module, levels = names(module_labels), labels = module_labels)
  )

p_violin <- ggplot(
  scores_long,
  aes(x = condition, y = score, fill = condition)
) +
  geom_violin(
    alpha = 0.7,
    trim  = TRUE
  ) +
  geom_boxplot(
    width         = 0.18,
    alpha         = 0.95,
    outlier.shape = NA,
    color         = "black"
  ) +
  geom_jitter(
    width = 0.12,
    size  = 1.0,
    alpha = 0.55,
    color = "black"
  ) +
  scale_fill_manual(values = condition_palette, guide = "none") +
  facet_grid(
    rows   = vars(dataset),
    cols   = vars(module),
    scales = "free_y"
  ) +
  labs(
    title    = "Module activity per condition",
    subtitle = "Violin + box + points for TGFb, mTOR, Angiogenesis, CSC modules",
    x        = "Condition",
    y        = "Module z-score"
  ) +
  theme_ras_csc(base_size = 12) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(0.5, "lines")
  )

ggsave(file.path(out_dir, "fig2_rnaseq_module_violins.pdf"),
       p_violin, width = 8.5, height = 4.8, units = "in")
ggsave(file.path(out_dir, "fig2_rnaseq_module_violins.png"),
       p_violin, width = 8.5, height = 4.8, units = "in", dpi = 300)

message("[DONE] Violin panels written to figures/main.")
