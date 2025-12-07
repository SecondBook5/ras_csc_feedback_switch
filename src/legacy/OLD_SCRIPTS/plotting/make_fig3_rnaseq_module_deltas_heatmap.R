#!/usr/bin/env Rscript

# ===========================================================
# make_fig3_rnaseq_module_deltas_heatmap.R
#
# Figure 3B: Tile heatmap of delta module scores (high - low
# condition) for TGFb, mTOR, Angio, CSC across Bl6, PAP_SCC,
# and PDV contrasts.
#
# Input:
#   data/processed/rnaseq/module_scores_by_sample.csv
#
# Output:
#   figures/main/fig3_rnaseq_module_deltas_heatmap.pdf
#   figures/main/fig3_rnaseq_module_deltas_heatmap.png
# ===========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

scores_path <- "data/processed/rnaseq/module_scores_by_sample.csv"
out_dir <- "figures/main"

if (!file.exists(scores_path)) {
  stop("[ERROR] module_scores_by_sample.csv not found at: ", scores_path)
}
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

#------------------- load data -----------------------------

scores <- readr::read_csv(scores_path, show_col_types = FALSE)

required_cols <- c("sample_id", "dataset", "condition")
module_cols <- setdiff(colnames(scores), required_cols)

if (!all(required_cols %in% colnames(scores))) {
  stop(
    "[ERROR] module_scores_by_sample.csv missing columns: ",
    paste(setdiff(required_cols, colnames(scores)), collapse = ", ")
  )
}
if (length(module_cols) == 0L) {
  stop("[ERROR] No module score columns found in module_scores_by_sample.csv.")
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

# readable module labels
module_labels <- c(
  TGFb_bulk  = "TGFb module",
  mTOR_bulk  = "mTOR module",
  Angio_bulk = "Angiogenesis module",
  CSC_bulk   = "CSC module"
)

# simple theme
theme_ras_csc <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      axis.title       = ggplot2::element_text(face = "bold"),
      axis.text        = ggplot2::element_text(colour = "black"),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      strip.text       = ggplot2::element_text(face = "bold"),
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle    = ggplot2::element_text(hjust = 0)
    )
}

#------------------- summarize means -----------------------

summary_tbl <- scores %>%
  group_by(dataset, condition) %>%
  summarise(
    across(
      all_of(module_cols),
      ~ mean(.x, na.rm = TRUE),
      .names = "{.col}_mean"
    ),
    .groups = "drop"
  )

get_row <- function(ds, cond) {
  summary_tbl %>%
    filter(dataset == ds, condition == cond)
}

build_contrast <- function(ds, cond_high, cond_low, label) {
  hi <- get_row(ds, cond_high)
  lo <- get_row(ds, cond_low)

  if (nrow(hi) != 1L || nrow(lo) != 1L) {
    warning(
      "[WARN] Contrast ", label, " for dataset ", ds,
      " could not be formed (missing rows). Skipping."
    )
    return(NULL)
  }

  res_list <- list()

  for (m in module_cols) {
    mean_hi <- hi[[paste0(m, "_mean")]]
    mean_lo <- lo[[paste0(m, "_mean")]]

    if (any(is.na(c(mean_hi, mean_lo)))) next

    delta <- mean_hi - mean_lo

    res_list[[m]] <- data.frame(
      dataset    = as.character(ds),
      contrast   = label,
      module_raw = m,
      delta      = delta
    )
  }

  if (length(res_list) == 0L) {
    return(NULL)
  }

  do.call(rbind, res_list)
}

delta_bl6 <- build_contrast(
  ds         = "Bl6",
  cond_high  = "SCC",
  cond_low   = "Normal",
  label      = "Bl6: SCC - Normal"
)

delta_pap <- build_contrast(
  ds         = "PAP_SCC",
  cond_high  = "SCC",
  cond_low   = "Papilloma",
  label      = "PAP/SCC: SCC - Papilloma"
)

delta_pdv <- build_contrast(
  ds         = "PDV",
  cond_high  = "PDV_LeprKO",
  cond_low   = "PDV_WT",
  label      = "PDV: LeprKO - WT"
)

delta_all <- bind_rows(delta_bl6, delta_pap, delta_pdv)

if (nrow(delta_all) == 0L) {
  stop("[ERROR] No contrasts could be computed; check dataset/condition names.")
}

delta_all <- delta_all %>%
  mutate(
    module = factor(
      module_raw,
      levels = names(module_labels),
      labels = module_labels
    ),
    contrast = factor(
      contrast,
      levels = c(
        "Bl6: SCC - Normal",
        "PAP/SCC: SCC - Papilloma",
        "PDV: LeprKO - WT"
      )
    )
  )

# symmetric color scale so zero is neutral
max_abs <- max(abs(delta_all$delta), na.rm = TRUE)

p_delta_heat <- ggplot(delta_all, aes(x = module, y = contrast, fill = delta)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = sprintf("%.2f", delta)),
    color = "black",
    size = 3
  ) +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = c(-max_abs, max_abs),
    name     = "Delta mean\nz-score"
  ) +
  labs(
    title    = "Contrasts in module activity between key conditions",
    subtitle = "Tiles show change in mean module z-score (high condition minus low condition)",
    x        = "Module",
    y        = "Contrast"
  ) +
  theme_ras_csc(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(out_dir, "fig3_rnaseq_module_deltas_heatmap.pdf"),
  p_delta_heat,
  width = 6.8, height = 3.4, units = "in"
)
ggsave(
  file.path(out_dir, "fig3_rnaseq_module_deltas_heatmap.png"),
  p_delta_heat,
  width = 6.8, height = 3.4, units = "in", dpi = 300
)

message("[DONE] Delta heatmap written to figures/main.")
