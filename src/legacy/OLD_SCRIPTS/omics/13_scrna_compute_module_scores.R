#!/usr/bin/env Rscript

# ===========================================================
# scrna_compute_module_scores.R
#
# Purpose:
#   - Load SCC K14+ Seurat object created by
#       scrna_build_seurat_and_module_scores.R
#   - Read pathway gene sets (TGFb, mTOR, Angiogenesis, CSC, etc.)
#     from a CSV with columns: module_name, gene_id
#   - For each module, intersect genes with Seurat rownames,
#     compute AddModuleScore, and store per-cell scores
#   - Export per-cell module scores to a CSV for ODE / network model
#
# Inputs:
#   data/processed/omics_summaries/scc_scRNA_seurat.rds
#   config/scrna_module_genesets.csv
#
# Output:
#   data/processed/omics_summaries/scc_scRNA_module_scores_per_cell.csv
# ===========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
})

message("[STEP] Setting paths")

seurat_rds_path <- "data/processed/omics_summaries/scc_scRNA_seurat.rds"
geneset_path    <- "config/scrna_module_genesets.csv"
out_scores_csv  <- "data/processed/omics_summaries/scc_scRNA_module_scores_per_cell.csv"

#-----------------------------------------------------------
# 1. Load Seurat object
#-----------------------------------------------------------

if (!file.exists(seurat_rds_path)) {
  stop("[ERROR] Seurat RDS not found at: ", seurat_rds_path,
       "\nRun scrna_build_seurat_and_module_scores.R first.")
}

message("[STEP] Loading Seurat object from: ", seurat_rds_path)
scc_obj <- readRDS(seurat_rds_path)

#-----------------------------------------------------------
# 2. Load curated gene sets
#-----------------------------------------------------------

if (!file.exists(geneset_path)) {
  stop("[ERROR] Gene set CSV not found at: ", geneset_path,
       "\nExpected a file with columns: module_name, gene_id")
}

message("[STEP] Reading gene sets from: ", geneset_path)
gs_df <- readr::read_csv(geneset_path, show_col_types = FALSE)

if (!all(c("module_name", "gene_id") %in% colnames(gs_df))) {
  stop("[ERROR] Gene set file must have columns 'module_name' and 'gene_id'.")
}

# Clean up any whitespace
gs_df <- gs_df %>%
  mutate(
    module_name = trimws(module_name),
    gene_id     = trimws(gene_id)
  ) %>%
  filter(module_name != "", gene_id != "")

if (nrow(gs_df) == 0L) {
  stop("[ERROR] Gene set table is empty after trimming.")
}

#-----------------------------------------------------------
# 3. Build feature lists for AddModuleScore
#-----------------------------------------------------------

all_genes <- rownames(scc_obj)
if (is.null(all_genes) || length(all_genes) == 0L) {
  stop("[ERROR] Seurat object has no gene rownames.")
}

message("[STEP] Intersecting gene sets with Seurat feature space")

module_names <- unique(gs_df$module_name)

feature_lists <- list()
effective_names <- character()

for (m in module_names) {
  genes_m <- gs_df %>%
    filter(module_name == m) %>%
    pull(gene_id) %>%
    unique()

  intersect_m <- intersect(genes_m, all_genes)

  message(
    "[INFO] Module '", m, "': ",
    length(genes_m), " genes in file; ",
    length(intersect_m), " present in Seurat object."
  )

  if (length(intersect_m) < 5L) {
    warning(
      "[WARN] Module '", m, "' has fewer than 5 overlapping genes; ",
      "skipping this module for AddModuleScore."
    )
    next
  }

  feature_lists[[length(feature_lists) + 1L]] <- intersect_m
  effective_names[length(effective_names) + 1L] <- m
}

if (length(feature_lists) == 0L) {
  stop(
    "[ERROR] No modules had >=5 overlapping genes with the Seurat object.\n",
    "Check that your gene IDs match rownames(scc_obj) (e.g. Ensembl+version)."
  )
}

message("[INFO] Number of modules to score: ", length(feature_lists))
message("[INFO] Modules: ", paste(effective_names, collapse = ", "))

#-----------------------------------------------------------
# 4. Compute module scores using AddModuleScore
#-----------------------------------------------------------

message("[STEP] Computing module scores via AddModuleScore")

# Seurat will name columns like 'TGFb1', 'mTOR2', etc. depending on length
# We call it once with all features lists to keep consistent scaling.
scc_obj <- AddModuleScore(
  object   = scc_obj,
  features = feature_lists,
  name     = effective_names,
  assay    = DefaultAssay(scc_obj),
  nbin     = 24,
  ctrl     = 100,
  seed     = 1,
  search   = TRUE
)

# After AddModuleScore, each module will add one column per feature list.
# For simplicity we assume one list per module, so we keep the first column
# that starts with 'moduleName'.

meta_df <- scc_obj@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell_id")

module_cols <- list()

for (m in effective_names) {
  # Seurat pattern: paste0(name, seq_along(feature_lists)) => e.g. "TGFb1"
  pattern_m <- paste0("^", m, "[0-9]+$")
  cols_m <- grep(pattern_m, colnames(meta_df), value = TRUE)

  if (length(cols_m) == 0L) {
    warning("[WARN] No AddModuleScore columns found for module '", m, "'.")
    next
  }

  # take the first (usually just one)
  chosen <- cols_m[[1]]
  module_cols[[m]] <- chosen
  message("[INFO] Module '", m, "' -> using column: ", chosen)
}

if (length(module_cols) == 0L) {
  stop(
    "[ERROR] AddModuleScore did not produce any detectable module score columns.\n",
    "Check naming / effective_names and rerun."
  )
}

# Build export table: cell_id + one column per module
scores_df <- meta_df %>%
  select(cell_id)

for (m in names(module_cols)) {
  col_src <- module_cols[[m]]
  new_name <- paste0(m, "_module_score")
  scores_df[[new_name]] <- meta_df[[col_src]]
}

#-----------------------------------------------------------
# 5. Write per-cell module scores to CSV
#-----------------------------------------------------------

readr::write_csv(scores_df, out_scores_csv)
message("[OK] Per-cell module scores written to: ", out_scores_csv)

message("[DONE] scrna_compute_module_scores.R completed.")
