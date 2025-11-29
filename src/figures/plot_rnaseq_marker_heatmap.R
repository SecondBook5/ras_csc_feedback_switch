#!/usr/bin/env Rscript

#===========================================================
# plot_rnaseq_marker_heatmap.R
#
# Build a marker-gene pheatmap for TGFb, mTOR, Angio, CSC
# across all bulk RNA-Seq samples.
#
# Outputs:
#   figures/main/fig1_rnaseq_marker_heatmap.pdf
#   figures/main/fig1_rnaseq_marker_heatmap.png
#===========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(yaml)
  library(pheatmap)
})

#---------------- paths ------------------------------------

expr_paths <- list(
  Bl6     = "data/interim/rnaseq/bl6_expression.tsv",
  PAP_SCC = "data/interim/rnaseq/pap_scc_expression.tsv",
  PDV     = "data/interim/rnaseq/pdv_expression.tsv"
)

sample_meta_path <- "data/interim/rnaseq/sample_metadata_GSE190411.csv"
geneset_path     <- "config/gene_sets_rnaseq.yaml"
out_dir          <- "figures/main"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

#---------------- helpers ----------------------------------

detect_gene_col <- function(df, stage = "detect_gene_col") {
  cand <- intersect(c("gene_symbol", "GeneSymbol", "gene", "Gene"), colnames(df))
  if (length(cand) == 0L) {
    stop(
      "[", stage, "] Could not find a gene column among: gene_symbol, GeneSymbol, gene, Gene"
    )
  }
  cand[[1]]
}

#---------------- load sample metadata ---------------------

if (!file.exists(sample_meta_path)) {
  stop("[ERROR] Sample metadata file not found at: ", sample_meta_path)
}

meta <- readr::read_csv(sample_meta_path, show_col_types = FALSE)

required_meta <- c("sample_id", "dataset", "condition")
if (!all(required_meta %in% colnames(meta))) {
  stop(
    "[ERROR] Sample metadata missing columns: ",
    paste(setdiff(required_meta, colnames(meta)), collapse = ", ")
  )
}

meta <- meta %>%
  mutate(
    dataset   = factor(dataset,   levels = c("Bl6", "PAP_SCC", "PDV")),
    condition = factor(condition, levels = c("Normal", "Papilloma", "SCC",
                                             "PDV_WT", "PDV_LeprKO"))
  )

#---------------- load expression and stack ----------------

expr_long_all <- list()

for (ds in names(expr_paths)) {
  path <- expr_paths[[ds]]
  if (!file.exists(path)) {
    stop("[ERROR] Expression file for dataset ", ds, " not found at: ", path)
  }

  message("[INFO] Reading expression for dataset ", ds, " from: ", path)
  df <- readr::read_tsv(path, show_col_types = FALSE)

  gene_col <- detect_gene_col(df, stage = paste0("expr_", ds))

  # unify gene column name
  df <- df %>%
    dplyr::rename(gene_symbol = !!gene_col)

  # long format: gene, sample_id, expr
  df_long <- df %>%
    tidyr::pivot_longer(
      cols      = -gene_symbol,
      names_to  = "sample_id",
      values_to = "expr"
    ) %>%
    mutate(
      dataset = ds
    )

  expr_long_all[[ds]] <- df_long
}

expr_long <- dplyr::bind_rows(expr_long_all)

# keep only samples that appear in metadata
expr_long <- expr_long %>%
  dplyr::semi_join(meta, by = "sample_id")

# summarise in case of duplicates
expr_long <- expr_long %>%
  dplyr::group_by(gene_symbol, sample_id) %>%
  dplyr::summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop")

# wide gene x sample matrix
expr_wide <- expr_long %>%
  tidyr::pivot_wider(
    id_cols     = gene_symbol,
    names_from  = sample_id,
    values_from = expr
  )

if (nrow(expr_wide) == 0L) {
  stop("[ERROR] Expression matrix ended up empty after joins.")
}

#---------------- load gene sets ---------------------------

if (!file.exists(geneset_path)) {
  stop("[ERROR] Gene set YAML file not found at: ", geneset_path)
}

genesets <- yaml::read_yaml(geneset_path)

modules_of_interest <- c("TGFb_bulk", "mTOR_bulk", "Angio_bulk", "CSC_bulk")

for (m in modules_of_interest) {
  if (is.null(genesets[[m]])) {
    stop("[ERROR] Gene set '", m, "' not found in ", geneset_path)
  }
}

# uppercase for intersection robustness
expr_wide <- expr_wide %>%
  mutate(gene_upper = toupper(gene_symbol))

geneset_upper <- lapply(genesets[modules_of_interest], function(v) {
  unique(toupper(unlist(v)))
})

#---------------- select marker genes ----------------------

top_n_per_module <- 10L

marker_rows <- list()
row_module  <- list()

for (m in modules_of_interest) {
  gs <- geneset_upper[[m]]

  idx <- which(expr_wide$gene_upper %in% gs)
  if (length(idx) == 0L) {
    warning("[WARN] No genes from module ", m, " found in expression matrix.")
    next
  }

  sub_mat <- as.matrix(expr_wide[idx, -c(1, ncol(expr_wide)), drop = FALSE])
  # row variance
  row_var <- apply(sub_mat, 1L, stats::var, na.rm = TRUE)

  # order by variance and take top N
  ord <- order(row_var, decreasing = TRUE)
  use <- ord[seq_len(min(top_n_per_module, length(ord)))]

  marker_rows[[m]] <- idx[use]
  row_module[[m]]  <- rep(m, length(use))
}

if (length(marker_rows) == 0L) {
  stop("[ERROR] No marker genes could be selected for any module.")
}

all_idx <- unlist(marker_rows, use.names = FALSE)
all_mod <- unlist(row_module, use.names = FALSE)

# subset expression
expr_markers <- expr_wide[all_idx, , drop = FALSE]

# expression numeric matrix (drop gene_symbol, gene_upper)
sample_cols <- setdiff(colnames(expr_markers), c("gene_symbol", "gene_upper"))

expr_mat <- as.matrix(expr_markers[, sample_cols, drop = FALSE])

# rownames as unique gene symbols
rownames(expr_mat) <- make.unique(expr_markers$gene_symbol)

#---------------- annotations ------------------------------

# column annotation from metadata
ann_col <- meta %>%
  dplyr::filter(sample_id %in% colnames(expr_mat)) %>%
  dplyr::select(sample_id, dataset, condition) %>%
  tibble::column_to_rownames("sample_id")

# ensure column order match
expr_mat <- expr_mat[, rownames(ann_col), drop = FALSE]

# row annotation: module label
ann_row <- data.frame(
  module = factor(
    all_mod,
    levels = modules_of_interest,
    labels = c("TGFb", "mTOR", "Angio", "CSC")
  )
)
rownames(ann_row) <- rownames(expr_mat)

#---------------- plotting (pheatmap) ----------------------

heat_colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)

pdf(file.path(out_dir, "fig1_rnaseq_marker_heatmap.pdf"),
    width = 7, height = 6)
pheatmap::pheatmap(
  mat                      = expr_mat,
  scale                    = "row",
  color                    = heat_colors,
  cluster_rows             = TRUE,
  cluster_cols             = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete",
  annotation_col           = ann_col,
  annotation_row           = ann_row,
  show_rownames            = TRUE,
  show_colnames            = FALSE,
  fontsize                 = 11,
  fontsize_row             = 6.5,
  fontsize_col             = 9,
  border_color             = NA,
  main                     = "Marker genes for TGFb, mTOR, Angio and CSC modules"
)
dev.off()

png(file.path(out_dir, "fig1_rnaseq_marker_heatmap.png"),
    width = 7, height = 6, units = "in", res = 300)
pheatmap::pheatmap(
  mat                      = expr_mat,
  scale                    = "row",
  color                    = heat_colors,
  cluster_rows             = TRUE,
  cluster_cols             = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete",
  annotation_col           = ann_col,
  annotation_row           = ann_row,
  show_rownames            = TRUE,
  show_colnames            = FALSE,
  fontsize                 = 11,
  fontsize_row             = 6.5,
  fontsize_col             = 9,
  border_color             = NA,
  main                     = "Marker genes for TGFb, mTOR, Angio and CSC modules"
)
dev.off()

message("[DONE] Marker heatmap written to figures/main.")
