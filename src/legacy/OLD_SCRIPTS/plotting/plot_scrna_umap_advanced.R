#!/usr/bin/env Rscript

#===========================================================
# plot_scrna_umap_advanced.R
#
# Uses per-cell metadata (with UMAP coords, cluster labels,
# Lepr expression, CSC flag) to make publication-grade:
#   1) Cluster UMAP with refined palette and layout
#   2) Lepr expression UMAP with continuous heat overlay
#   3) CSC island UMAP with ghost background and density
#
# Input:
#   data/processed/omics_summaries/
#       scc_scRNA_K14pos_metadata_with_CSC_labels.csv
#
# Output:
#   figures/main/fig_scRNA_umap_clusters_advanced.[pdf|png]
#   figures/main/fig_scRNA_umap_lepr_advanced.[pdf|png]
#   figures/main/fig_scRNA_umap_csc_island_advanced.[pdf|png]
#===========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

meta_path <- "data/processed/omics_summaries/scc_scRNA_K14pos_metadata_with_CSC_labels.csv"
out_dir   <- "figures/main"

if (!file.exists(meta_path)) {
  stop("[ERROR] Metadata CSV not found: ", meta_path,
       "\nRun scrna_build_seurat_and_module_scores.R first.")
}
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

#------------------------ load metadata ---------------------

meta <- readr::read_csv(meta_path, show_col_types = FALSE)

needed_cols <- c(
  "UMAP_1", "UMAP_2",
  "louvain_cluster",
  "is_CSC",
  "Lepr_expr"
)

missing_cols <- setdiff(needed_cols, colnames(meta))
if (length(missing_cols) > 0L) {
  stop("[ERROR] Metadata is missing columns: ",
       paste(missing_cols, collapse = ", "))
}

meta <- meta %>%
  mutate(
    louvain_cluster = factor(louvain_cluster),
    is_CSC          = as.logical(is_CSC)
  )

# simple house theme
theme_ras_csc <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid       = element_blank(),
      axis.title       = element_text(face = "bold"),
      axis.text        = element_text(colour = "black"),
      plot.title       = element_text(face = "bold", hjust = 0),
      plot.subtitle    = element_text(hjust = 0),
      legend.title     = element_text(face = "bold")
    )
}

# Define a discrete palette for clusters (use many but gentle colors)
cluster_cols <- c(
  "#F8766D", "#C49A00", "#7CAE00", "#00BFC4",
  "#619CFF", "#AA00AA", "#FFB6C1", "#00A08A"
)

#------------------------ 1) cluster UMAP -------------------

p_clusters <- ggplot(meta, aes(x = UMAP_1, y = UMAP_2)) +
  # light density underlay for shape
  stat_density_2d(
    aes(fill = after_stat(level)),
    geom = "polygon",
    colour = NA,
    alpha  = 0.18,
    contour_var = "ndensity"
  ) +
  scale_fill_gradient(low = "grey92", high = "grey70", guide = "none") +
  geom_point(
    aes(colour = louvain_cluster),
    size = 1.1,
    alpha = 0.9
  ) +
  scale_colour_manual(
    values = cluster_cols[seq_along(levels(meta$louvain_cluster))],
    name   = "Cluster"
  ) +
  labs(
    title    = "SCC K14+ tumour cells",
    subtitle = "UMAP coloured by Louvain clusters",
    x        = "UMAP 1",
    y        = "UMAP 2"
  ) +
  coord_equal() +
  theme_ras_csc(base_size = 13) +
  theme(
    legend.position = c(0.98, 0.15),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", colour = "grey80")
  )

ggsave(
  file.path(out_dir, "fig_scRNA_umap_clusters_advanced.pdf"),
  p_clusters, width = 6.2, height = 4.8, units = "in"
)
ggsave(
  file.path(out_dir, "fig_scRNA_umap_clusters_advanced.png"),
  p_clusters, width = 6.2, height = 4.8, units = "in", dpi = 300
)

#------------------------ 2) Lepr UMAP ----------------------

# cap Lepr for colour scale sanity
lepr_cap <- quantile(meta$Lepr_expr, 0.99, na.rm = TRUE)
meta <- meta %>%
  mutate(Lepr_expr_capped = pmin(Lepr_expr, lepr_cap))

p_lepr <- ggplot(meta, aes(x = UMAP_1, y = UMAP_2)) +
  stat_density_2d(
    aes(fill = after_stat(ndensity)),
    geom = "raster",
    contour = FALSE,
    alpha = 0.20
  ) +
  scale_fill_gradient(low = "black", high = "grey70", guide = "none") +
  geom_point(
    aes(colour = Lepr_expr_capped),
    size = 1.0,
    alpha = 0.95
  ) +
  scale_colour_gradientn(
    colours = c("#440154", "#3B528B", "#21918C", "#5DC963", "#FDE725"),
    name    = "Lepr\n(log-normalized)",
    limits  = c(min(meta$Lepr_expr_capped, na.rm = TRUE), lepr_cap),
    oob     = squish
  ) +
  labs(
    title    = "Lepr signalling in SCC K14+ tumour cells",
    subtitle = "UMAP coloured by Lepr expression",
    x        = "UMAP 1",
    y        = "UMAP 2"
  ) +
  coord_equal() +
  theme_ras_csc(base_size = 13) +
  theme(
    legend.position = c(0.98, 0.2),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", colour = "grey80")
  )

ggsave(
  file.path(out_dir, "fig_scRNA_umap_lepr_advanced.pdf"),
  p_lepr, width = 6.2, height = 4.8, units = "in"
)
ggsave(
  file.path(out_dir, "fig_scRNA_umap_lepr_advanced.png"),
  p_lepr, width = 6.2, height = 4.8, units = "in", dpi = 300
)

#------------------------ 3) CSC island UMAP ----------------

meta <- meta %>%
  mutate(
    CSC_state = ifelse(is_CSC, "CSC", "non-CSC"),
    CSC_state = factor(CSC_state, levels = c("non-CSC", "CSC"))
  )

# restrict density overlay to CSC cells only
meta_csc <- meta %>% filter(is_CSC)

p_csc <- ggplot() +
  # background non-CSC silhouette
  geom_point(
    data = meta %>% filter(!is_CSC),
    aes(x = UMAP_1, y = UMAP_2),
    colour = "grey80",
    size   = 0.8,
    alpha  = 0.4
  ) +
  # density of CSC island
  stat_density_2d(
    data = meta_csc,
    aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level)),
    geom        = "polygon",
    colour      = "grey40",
    alpha       = 0.55,
    contour_var = "ndensity"
  ) +
  scale_fill_gradientn(
    colours = c("#4C1D4B", "#AE305C", "#F5764C", "#FDC527"),
    name    = "CSC density"
  ) +
  # CSC points on top
  geom_point(
    data = meta_csc,
    aes(x = UMAP_1, y = UMAP_2),
    colour = "#B2182B",
    size   = 1.3,
    alpha  = 0.9
  ) +
  labs(
    title    = "CSC island within K14+ SCC state space",
    subtitle = "UMAP with CSC cells and density overlay",
    x        = "UMAP 1",
    y        = "UMAP 2"
  ) +
  coord_equal() +
  theme_ras_csc(base_size = 13) +
  theme(
    legend.position = c(0.98, 0.18),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", colour = "grey80")
  )

ggsave(
  file.path(out_dir, "fig_scRNA_umap_csc_island_advanced.pdf"),
  p_csc, width = 6.2, height = 4.8, units = "in"
)
ggsave(
  file.path(out_dir, "fig_scRNA_umap_csc_island_advanced.png"),
  p_csc, width = 6.2, height = 4.8, units = "in", dpi = 300
)

message("[DONE] Advanced UMAP figures written to ", out_dir)
