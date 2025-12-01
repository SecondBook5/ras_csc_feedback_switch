#!/usr/bin/env Rscript

# ===========================================================
# scrna_build_seurat_and_module_scores.R
#
# Purpose:
#   - Take K14+ SCC scRNA-Seq (prepared by scrna_prepare_k14_and_hvgs.R)
#   - Build a Seurat object from raw counts
#   - Use Yuan-style ERCC HVGs as VariableFeatures
#   - Run PCA -> neighbors -> clustering -> UMAP
#   - Identify a CSC cluster as the Lepr-high cluster
#   - Compute CSC fractions and marker summaries
#   - Save per-cell metadata, a compact JSON summary, and the Seurat RDS
#   - Generate polished UMAP figures
#
# Inputs:
#   data/interim/omics_qc/scc_scRNA_K14pos_expression_logged_filtered.csv
#   data/interim/omics_qc/scc_scRNA_K14pos_counts_filtered.csv
#   data/interim/omics_qc/all_tumor_K14morethan7_Counts_geneCollapsed_filtered_variable.csv
#
# Outputs:
#   data/processed/omics_summaries/scc_scRNA_K14pos_metadata_with_CSC_labels.csv
#   data/processed/omics_summaries/scc_scRNA_CSC_summary.json
#   data/processed/omics_summaries/scc_scRNA_seurat.rds
#   figures/main/scc_scRNA_umap_louvain_clusters.png
#   figures/main/scc_scRNA_umap_Lepr_expr.png
#   figures/main/scc_scRNA_umap_CSC_vs_nonCSC.png
# ===========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(Seurat)
  library(ggplot2)
  library(jsonlite)
})

message("[STEP] Setting paths")

expr_log_path <- "data/interim/omics_qc/scc_scRNA_K14pos_expression_logged_filtered.csv"
counts_path <- "data/interim/omics_qc/scc_scRNA_K14pos_counts_filtered.csv"
hvg_path <- "data/interim/omics_qc/all_tumor_K14morethan7_Counts_geneCollapsed_filtered_variable.csv"

out_dir_proc <- "data/processed/omics_summaries"
out_dir_fig <- "figures/main"

meta_out_csv <- file.path(out_dir_proc, "scc_scRNA_K14pos_metadata_with_CSC_labels.csv")
summary_out_json <- file.path(out_dir_proc, "scc_scRNA_CSC_summary.json")
seurat_rds_path <- file.path(out_dir_proc, "scc_scRNA_seurat.rds")

if (!dir.exists(out_dir_proc)) dir.create(out_dir_proc, recursive = TRUE)
if (!dir.exists(out_dir_fig)) dir.create(out_dir_fig, recursive = TRUE)

# Marker Ensembl IDs (consistent with the original notebook)
lepr_id <- "ENSMUSG00000059316.2"
krt14_id <- "ENSMUSG00000045545.8"

# ----------------------------------------------------------
# 1. Load K14+ log expression and counts
# ----------------------------------------------------------

if (!file.exists(expr_log_path)) {
  stop("[ERROR] K14+ log expression file not found at: ", expr_log_path)
}
if (!file.exists(counts_path)) {
  stop("[ERROR] K14+ counts file not found at: ", counts_path)
}
if (!file.exists(hvg_path)) {
  stop("[ERROR] HVG file not found at: ", hvg_path)
}

message("[STEP] Reading K14+ log2(TPM+1) expression from: ", expr_log_path)
expr_df <- readr::read_csv(expr_log_path, show_col_types = FALSE)

if (!"gene_id" %in% colnames(expr_df)) {
  stop("[ERROR] Expression file is missing 'gene_id' column.")
}
gene_ids_expr <- expr_df$gene_id
expr_mat_log <- as.matrix(expr_df[, -1, drop = FALSE])
rownames(expr_mat_log) <- gene_ids_expr

message("[STEP] Reading K14+ counts from: ", counts_path)
counts_df <- readr::read_csv(counts_path, show_col_types = FALSE)

if (!"gene_id" %in% colnames(counts_df)) {
  stop("[ERROR] Counts file is missing 'gene_id' column.")
}
gene_ids_counts <- counts_df$gene_id
counts_mat <- as.matrix(counts_df[, -1, drop = FALSE])
rownames(counts_mat) <- gene_ids_counts

# Sanity checks
if (!identical(rownames(expr_mat_log), rownames(counts_mat))) {
  stop("[ERROR] Gene IDs differ between log expression and counts matrices.")
}
if (!identical(colnames(expr_mat_log), colnames(counts_mat))) {
  stop("[ERROR] Cell columns differ between log expression and counts matrices.")
}

message("[INFO] K14+ counts dimensions: ", paste(dim(counts_mat), collapse = " x "))

# ----------------------------------------------------------
# 2. Load HVGs from Yuan technical-noise model
# ----------------------------------------------------------

message("[STEP] Reading HVG flags from: ", hvg_path)
hvg_df <- readr::read_csv(hvg_path, show_col_types = FALSE)

if (!all(c("gene_id", "sig") %in% colnames(hvg_df))) {
  stop("[ERROR] HVG file must have columns 'gene_id' and 'sig'.")
}

hvg_ids <- hvg_df %>%
  filter(sig) %>%
  pull(gene_id)

hvg_ids_in_obj <- intersect(hvg_ids, rownames(counts_mat))

if (length(hvg_ids_in_obj) < 200L) {
  warning(
    "[WARN] Fewer than 200 HVGs intersect with expression matrix (n = ",
    length(hvg_ids_in_obj),
    "). Clustering may be unstable."
  )
}

message("[INFO] HVGs from Yuan model present in object: ", length(hvg_ids_in_obj))

# ----------------------------------------------------------
# 3. Build Seurat object and run standard workflow
# ----------------------------------------------------------

message("[STEP] Creating Seurat object")

scc_obj <- CreateSeuratObject(
  counts       = counts_mat,
  project      = "SCC_K14",
  min.cells    = 0,
  min.features = 0
)

# Sample ID from cell barcode (prefix before first underscore)
scc_obj$sample_id <- sub("_.*", "", colnames(scc_obj))

# Use Yuan HVGs as VariableFeatures
VariableFeatures(scc_obj) <- hvg_ids_in_obj

message("[STEP] Normalizing, scaling, PCA, neighbors, clustering, UMAP")

scc_obj <- NormalizeData(
  scc_obj,
  normalization.method = "LogNormalize",
  scale.factor         = 10000,
  verbose              = FALSE
)

scc_obj <- ScaleData(
  scc_obj,
  features = VariableFeatures(scc_obj),
  verbose  = FALSE
)

scc_obj <- RunPCA(
  scc_obj,
  features = VariableFeatures(scc_obj),
  npcs     = 30,
  verbose  = FALSE
)

scc_obj <- FindNeighbors(
  scc_obj,
  dims    = 1:20,
  verbose = FALSE
)

scc_obj <- FindClusters(
  scc_obj,
  resolution = 0.4,
  verbose    = FALSE
)

scc_obj <- RunUMAP(
  scc_obj,
  dims    = 1:20,
  verbose = FALSE
)

# Copy UMAP coordinates into metadata columns
umap_emb <- Embeddings(scc_obj, reduction = "umap")
scc_obj$UMAP_1 <- umap_emb[, 1]
scc_obj$UMAP_2 <- umap_emb[, 2]

# ----------------------------------------------------------
# 4. Add Lepr and Krt14 expression to metadata
# ----------------------------------------------------------

message("[STEP] Adding Lepr and Krt14 expression to metadata")

missing_markers <- setdiff(c(lepr_id, krt14_id), rownames(scc_obj))
if (length(missing_markers) > 0L) {
  stop(
    "[ERROR] Marker gene IDs not found in Seurat object rownames: ",
    paste(missing_markers, collapse = ", ")
  )
}

# Use normalized data ("data" slot) for markers
norm_data <- GetAssayData(scc_obj, slot = "data")

scc_obj$Lepr_expr <- as.numeric(norm_data[lepr_id, colnames(scc_obj), drop = TRUE])
scc_obj$Krt14_expr <- as.numeric(norm_data[krt14_id, colnames(scc_obj), drop = TRUE])

# ----------------------------------------------------------
# 5. Identify CSC cluster as Lepr-high cluster
# ----------------------------------------------------------

message("[STEP] Summarizing markers by cluster and defining CSC cluster")

meta_df <- scc_obj@meta.data %>%
  as.data.frame() %>%
  rownames_to_column(var = "cell_id")

cluster_summary <- meta_df %>%
  group_by(seurat_clusters) %>%
  summarize(
    n_cells      = n(),
    Lepr_mean    = mean(Lepr_expr, na.rm = TRUE),
    Lepr_median  = median(Lepr_expr, na.rm = TRUE),
    Krt14_mean   = mean(Krt14_expr, na.rm = TRUE),
    Krt14_median = median(Krt14_expr, na.rm = TRUE),
    .groups      = "drop"
  )

if (nrow(cluster_summary) == 0L) {
  stop("[ERROR] No clusters found in seurat_clusters.")
}

csc_cluster_id <- cluster_summary %>%
  filter(Lepr_mean == max(Lepr_mean, na.rm = TRUE)) %>%
  slice(1) %>%
  pull(seurat_clusters)

message("[INFO] Selected CSC cluster ID (highest mean Lepr_expr): ", csc_cluster_id)

scc_obj$is_CSC <- scc_obj$seurat_clusters == csc_cluster_id

# ----------------------------------------------------------
# 6. Global and per-sample CSC summaries
# ----------------------------------------------------------

meta_df <- scc_obj@meta.data %>%
  as.data.frame() %>%
  rownames_to_column(var = "cell_id")

global_csc_fraction <- mean(meta_df$is_CSC)

sample_stats <- meta_df %>%
  group_by(sample_id) %>%
  summarize(
    n_cells      = n(),
    n_CSC        = sum(is_CSC),
    CSC_fraction = ifelse(n_cells > 0, n_CSC / n_cells, NA_real_),
    .groups      = "drop"
  )

cluster_csc_stats <- meta_df %>%
  group_by(seurat_clusters) %>%
  summarize(
    n_cells      = n(),
    n_CSC        = sum(is_CSC),
    CSC_fraction = ifelse(n_cells > 0, n_CSC / n_cells, NA_real_),
    .groups      = "drop"
  )

message("[INFO] Global CSC fraction among K14+ cells: ", sprintf("%.3f", global_csc_fraction))

# ----------------------------------------------------------
# 7. Save per-cell metadata and JSON summary
# ----------------------------------------------------------

message("[STEP] Writing per-cell metadata and CSC summary")

meta_for_export <- meta_df %>%
  select(
    cell_id,
    sample_id,
    seurat_clusters,
    is_CSC,
    UMAP_1,
    UMAP_2,
    Lepr_expr,
    Krt14_expr,
    nCount_RNA,
    nFeature_RNA
  )

readr::write_csv(meta_for_export, meta_out_csv)
message("[OK] Per-cell metadata with CSC labels written to: ", meta_out_csv)

summary_list <- list(
  n_cells_total = nrow(meta_df),
  n_CSC = sum(meta_df$is_CSC),
  n_non_CSC = sum(!meta_df$is_CSC),
  CSC_fraction_K14 = global_csc_fraction,
  CSC_cluster_id = as.character(csc_cluster_id),
  cluster_marker_summary = cluster_summary,
  per_sample_CSC = sample_stats,
  per_cluster_CSC = cluster_csc_stats
)

writeLines(
  jsonlite::toJSON(summary_list, pretty = TRUE, auto_unbox = TRUE),
  con = summary_out_json
)

message("[OK] CSC summary JSON written to: ", summary_out_json)

# ----------------------------------------------------------
# 8. Save Seurat object for downstream module scoring
# ----------------------------------------------------------

saveRDS(scc_obj, seurat_rds_path)
message("[OK] Saved Seurat object with CSC labels to: ", seurat_rds_path)

# ----------------------------------------------------------
# 9. UMAP figures (polished)
# ----------------------------------------------------------

theme_sc_advanced <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_rect(fill = "white", colour = NA),
    axis.text          = element_text(color = "black"),
    axis.title         = element_text(face = "bold"),
    axis.ticks         = element_line(colour = "grey70", linewidth = 0.3),
    axis.line          = element_line(colour = "grey70", linewidth = 0.3),
    legend.title       = element_text(face = "bold"),
    legend.background  = element_rect(fill = "white", colour = NA),
    legend.key         = element_rect(fill = "white", colour = NA),
    plot.title         = element_text(face = "bold", hjust = 0, size = 14),
    plot.subtitle      = element_text(hjust = 0, size = 11)
  )

umap_df <- meta_df %>%
  select(
    cell_id, UMAP_1, UMAP_2,
    seurat_clusters, is_CSC, Lepr_expr, sample_id
  ) %>%
  mutate(
    seurat_clusters = factor(seurat_clusters),
    sample_id       = factor(sample_id)
  )

# UMAP by cluster
p_umap_cluster <- ggplot(
  umap_df,
  aes(x = UMAP_1, y = UMAP_2, colour = seurat_clusters)
) +
  geom_point(size = 0.5, alpha = 0.9) +
  scale_color_hue(name = "Cluster") +
  labs(
    title    = "SCC K14+ tumour cells",
    subtitle = "UMAP coloured by Louvain clusters",
    x        = "UMAP 1",
    y        = "UMAP 2"
  ) +
  coord_equal() +
  theme_sc_advanced +
  guides(
    colour = guide_legend(
      override.aes = list(size = 2, alpha = 1),
      ncol = 2
    )
  )

ggsave(
  filename = file.path(out_dir_fig, "scc_scRNA_umap_louvain_clusters.png"),
  plot = p_umap_cluster,
  width = 5.8, height = 4.8, units = "in", dpi = 400
)

# UMAP by Lepr expression
p_umap_lepr <- ggplot(
  umap_df,
  aes(x = UMAP_1, y = UMAP_2, colour = Lepr_expr)
) +
  geom_point(size = 0.55, alpha = 0.95) +
  scale_colour_viridis_c(
    name   = "Lepr\n(log-normalized)",
    option = "C"
  ) +
  labs(
    title    = "Lepr signalling in SCC K14+ tumour cells",
    subtitle = "UMAP coloured by Lepr expression",
    x        = "UMAP 1",
    y        = "UMAP 2"
  ) +
  coord_equal() +
  theme_sc_advanced +
  theme(
    legend.position = "right"
  )

ggsave(
  filename = file.path(out_dir_fig, "scc_scRNA_umap_Lepr_expr.png"),
  plot = p_umap_lepr,
  width = 5.8, height = 4.8, units = "in", dpi = 400
)

# UMAP CSC vs non-CSC with density overlay
umap_bg <- umap_df %>% filter(!is.na(is_CSC))
umap_csc <- umap_bg %>% filter(is_CSC)

p_umap_csc <- ggplot() +
  geom_point(
    data = umap_bg %>% filter(!is_CSC),
    mapping = aes(x = UMAP_1, y = UMAP_2),
    colour = "grey85",
    size = 0.4,
    alpha = 0.7
  ) +
  geom_point(
    data = umap_csc,
    mapping = aes(x = UMAP_1, y = UMAP_2),
    colour = "red3",
    size = 0.6,
    alpha = 0.95
  ) +
  stat_density_2d(
    data = umap_csc,
    aes(x = UMAP_1, y = UMAP_2, fill = after_stat(level)),
    geom = "polygon",
    alpha = 0.25,
    colour = NA
  ) +
  scale_fill_viridis_c(
    name   = "CSC density",
    option = "B",
    guide  = "none"
  ) +
  labs(
    title    = "CSC island within K14+ SCC state space",
    subtitle = "UMAP with CSC cells and density overlay",
    x        = "UMAP 1",
    y        = "UMAP 2"
  ) +
  coord_equal() +
  theme_sc_advanced

ggsave(
  filename = file.path(out_dir_fig, "scc_scRNA_umap_CSC_vs_nonCSC.png"),
  plot = p_umap_csc,
  width = 5.8, height = 4.8, units = "in", dpi = 400
)

message("[DONE] scRNA Seurat CSC summary completed.")
