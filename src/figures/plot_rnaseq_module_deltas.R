#!/usr/bin/env Rscript

#===========================================================
# plot_rnaseq_module_deltas.R
#
# Fancier delta barplot:
#   - Computes delta = mean(module) in condition_high - condition_low
#   - For 3 contrasts:
#       * Bl6:     SCC - Normal
#       * PAP_SCC: SCC - Papilloma
#       * PDV:     PDV_LeprKO - PDV_WT
#   - Shows bars with SE error bars and clean facets.
#
# Input:
#   data/processed/rnaseq/module_scores_by_sample.csv
#
# Output:
#   figures/main/fig3_rnaseq_module_deltas.pdf
#   figures/main/fig3_rnaseq_module_deltas.png
#===========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

scores_path <- "data/processed/rnaseq/module_scores_by_sample.csv"
out_dir     <- "figures/main"

if (!file.exists(scores_path)) {
  stop("[ERROR] module_scores_by_sample.csv not found at: ", scores_path)
}
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

#------------------- load data -----------------------------

scores <- readr::read_csv(scores_path, show_col_types = FALSE)

required_cols <- c("sample_id", "dataset", "condition")
module_cols   <- setdiff(colnames(scores), required_cols)

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

# palette for modules
module_palette <- c(
  "TGFb module"         = "#1f78b4",
  "mTOR module"         = "#33a02c",
  "Angiogenesis module" = "#fb9a99",
  "CSC module"          = "#6a3d9a"
)

#------------------- theme helper --------------------------

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

#------------------- summarize means -----------------------

summary_tbl <- scores %>%
  group_by(dataset, condition) %>%
  summarise(
    across(
      all_of(module_cols),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd   = ~ sd(.x, na.rm = TRUE),
        n    = ~ sum(!is.na(.x))
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# convenience accessor: pick one dataset/condition row
get_row <- function(ds, cond) {
  summary_tbl %>%
    filter(dataset == ds, condition == cond)
}

#------------------- build contrasts -----------------------

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

  # for each module, columns are e.g. TGFb_bulk_mean, TGFb_bulk_sd, TGFb_bulk_n
  res_list <- list()

  for (m in module_cols) {
    mean_hi <- hi[[paste0(m, "_mean")]]
    mean_lo <- lo[[paste0(m, "_mean")]]

    sd_hi   <- hi[[paste0(m, "_sd")]]
    sd_lo   <- lo[[paste0(m, "_sd")]]

    n_hi    <- hi[[paste0(m, "_n")]]
    n_lo    <- lo[[paste0(m, "_n")]]

    if (any(is.na(c(mean_hi, mean_lo, sd_hi, sd_lo, n_hi, n_lo))) ||
        n_hi == 0 || n_lo == 0) {
      next
    }

    delta <- mean_hi - mean_lo

    # standard error of difference of means
    se_delta <- sqrt((sd_hi^2 / n_hi) + (sd_lo^2 / n_lo))

    res_list[[m]] <- data.frame(
      dataset    = as.character(ds),
      contrast   = label,
      module_raw = m,
      delta      = delta,
      se_delta   = se_delta
    )
  }

  if (length(res_list) == 0L) return(NULL)

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
      levels = c("Bl6: SCC - Normal",
                 "PAP/SCC: SCC - Papilloma",
                 "PDV: LeprKO - WT")
    )
  )

#------------------- plot -----------------------

p_delta <- ggplot(delta_all, aes(x = module, y = delta, fill = module)) +
  # light band around zero to emphasize sign
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = -0.05, ymax = 0.05,
    fill = "grey92", alpha = 0.7
  ) +
  geom_col(width = 0.55, color = "black", linewidth = 0.25) +
  geom_errorbar(
    aes(ymin = delta - se_delta, ymax = delta + se_delta),
    width     = 0.15,
    linewidth = 0.3
  ) +
  geom_hline(
    yintercept = 0,
    linetype   = "dashed",
    color      = "#6B7280",
    linewidth  = 0.35
  ) +
  scale_fill_manual(values = module_palette, guide = "none") +
  facet_wrap(
    ~ contrast,
    nrow   = 1,
    scales = "free_y"
  ) +
  labs(
    title    = "Contrasts in module activity between key conditions",
    subtitle = "Bars show change in mean module z-score (high condition minus low condition)\nError bars show standard error of the difference of means",
    x        = "Module",
    y        = "Change in module z-score"
  ) +
  theme_ras_csc(base_size = 12) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(0.6, "lines")
  )

ggsave(
  file.path(out_dir, "fig3_rnaseq_module_deltas.pdf"),
  p_delta, width = 8.2, height = 2.9, units = "in"
)
ggsave(
  file.path(out_dir, "fig3_rnaseq_module_deltas.png"),
  p_delta, width = 8.2, height = 2.9, units = "in", dpi = 300
)

message("[DONE] Delta module barplot written to figures/main.")
