#!/usr/bin/env Rscript

# ===========================================================
# build_scrna_module_genesets_from_yaml.R
#
# Purpose:
#   - Read pathway gene sets (TGFb_bulk, mTOR_bulk, Angio_bulk, CSC_bulk)
#     from config/gene_sets_rnaseq.yaml (symbols).
#   - Map symbols to Ensembl IDs using biomaRt.
#   - Harmonize with scRNA Seurat object's feature names
#     (Ensembl+version, e.g. ENSMUSG00000045545.8).
#   - Write config/scrna_module_genesets.csv with columns:
#       module_name, gene_id
#     which scrna_compute_module_scores.R already expects.
#
# Inputs:
#   config/gene_sets_rnaseq.yaml
#   data/processed/omics_summaries/scc_scRNA_seurat.rds
#
# Output:
#   config/scrna_module_genesets.csv
# ===========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
  library(biomaRt)
  library(yaml)
  library(purrr)
})

message("[STEP] Setting paths")

yaml_path <- "config/gene_sets_rnaseq.yaml"
seurat_rds <- "data/processed/omics_summaries/scc_scRNA_seurat.rds"
out_geneset_csv <- "config/scrna_module_genesets.csv"

#-----------------------------------------------------------
# 1. Sanity checks and load Seurat object
#-----------------------------------------------------------

if (!file.exists(seurat_rds)) {
  stop(
    "[ERROR] Seurat RDS not found at: ", seurat_rds, "\n",
    "Run scrna_build_seurat_and_module_scores.R first."
  )
}

if (!file.exists(yaml_path)) {
  stop(
    "[ERROR] YAML gene set config not found at: ", yaml_path, "\n",
    "Expected your TGFb_bulk / mTOR_bulk / Angio_bulk / CSC_bulk here."
  )
}

message("[STEP] Loading Seurat object from: ", seurat_rds)
scc_obj <- readRDS(seurat_rds)

genes_sc <- rownames(scc_obj)
if (is.null(genes_sc) || length(genes_sc) == 0L) {
  stop("[ERROR] Seurat object has no feature rownames.")
}

message("[INFO] Seurat feature count: ", length(genes_sc))

# Pre-compute bare Ensembl IDs (without version) from Seurat rownames
lookup_features <- tibble(
  gene_id_seurat = genes_sc,
  bare_id        = sub("\\..*$", "", genes_sc)
)

#-----------------------------------------------------------
# 2. Read YAML gene sets
#-----------------------------------------------------------

message("[STEP] Reading gene sets from YAML: ", yaml_path)
gs_yml <- yaml::read_yaml(yaml_path)

if (length(gs_yml) == 0L) {
  stop("[ERROR] gene_sets_rnaseq.yaml is empty or not a list.")
}

# Flatten into a data frame: module_name_raw, symbol
gene_symbol_df <- purrr::imap_dfr(
  gs_yml,
  ~ tibble(
    module_name_raw = .y,
    symbol          = as.character(.x)
  )
)

# Drop empties / NAs
gene_symbol_df <- gene_symbol_df %>%
  dplyr::filter(!is.na(symbol), symbol != "")

if (nrow(gene_symbol_df) == 0L) {
  stop("[ERROR] No non-empty gene symbols found in YAML.")
}

# Clean module names: drop _bulk suffix for scoring if desired
gene_symbol_df <- gene_symbol_df %>%
  dplyr::mutate(
    module_name = gsub("_bulk$", "", module_name_raw)
  )

message(
  "[INFO] Modules found in YAML (raw): ",
  paste(unique(gene_symbol_df$module_name_raw), collapse = ", ")
)
message(
  "[INFO] Cleaned module names: ",
  paste(unique(gene_symbol_df$module_name), collapse = ", ")
)

#-----------------------------------------------------------
# 3. SYMBOL -> Ensembl mapping via biomaRt
#-----------------------------------------------------------

message("[STEP] Mapping symbols -> Ensembl via biomaRt (mmusculus_gene_ensembl)")

mart <- biomaRt::useEnsembl(
  biomart = "ensembl",
  dataset = "mmusculus_gene_ensembl"
)

map_df <- biomaRt::getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters    = "external_gene_name",
  values     = unique(gene_symbol_df$symbol),
  mart       = mart
)

if (nrow(map_df) == 0L) {
  stop(
    "[ERROR] biomaRt returned no mappings.\n",
    "Check internet access and that the symbols in YAML are valid mouse gene names."
  )
}

# Join mappings back to module table
gene_symbol_df <- gene_symbol_df %>%
  dplyr::left_join(
    map_df,
    by = c("symbol" = "external_gene_name")
  )

# Remove unmapped symbols
unmapped <- gene_symbol_df %>%
  dplyr::filter(is.na(ensembl_gene_id)) %>%
  dplyr::distinct(symbol)

if (nrow(unmapped) > 0L) {
  warning(
    "[WARN] Some symbols had no Ensembl mapping and will be dropped: ",
    paste(unmapped$symbol, collapse = ", ")
  )
}

gene_symbol_df <- gene_symbol_df %>%
  dplyr::filter(!is.na(ensembl_gene_id))

if (nrow(gene_symbol_df) == 0L) {
  stop("[ERROR] After mapping, no genes had valid Ensembl IDs.")
}

#-----------------------------------------------------------
# 4. Harmonize Ensembl IDs with Seurat (add versioned IDs)
#-----------------------------------------------------------

message("[STEP] Intersecting Ensembl IDs with Seurat features")

gene_symbol_df <- gene_symbol_df %>%
  dplyr::left_join(
    lookup_features,
    by = c("ensembl_gene_id" = "bare_id")
  )

# Drop any that don't appear in the Seurat object
gene_symbol_df <- gene_symbol_df %>%
  dplyr::filter(!is.na(gene_id_seurat))

if (nrow(gene_symbol_df) == 0L) {
  stop(
    "[ERROR] After intersection, none of the YAML gene sets overlap Seurat features.\n",
    "Check that Seurat rownames are Ensembl+version and that we are stripping versions correctly."
  )
}

# De-duplicate module_name x gene_id_seurat
gene_symbol_df <- gene_symbol_df %>%
  dplyr::distinct(module_name, gene_id = gene_id_seurat)

# Print counts per module
module_counts <- gene_symbol_df %>%
  dplyr::count(module_name)

message("[INFO] Overlapping genes per module:")
print(module_counts)

small_modules <- module_counts %>% dplyr::filter(n < 5)
if (nrow(small_modules) > 0L) {
  warning(
    "[WARN] Some modules have < 5 overlapping genes with the scRNA object:\n",
    paste(
      paste0("  - ", small_modules$module_name, ": ", small_modules$n, " genes"),
      collapse = "\n"
    ),
    "\nAdd more genes if you want stable AddModuleScore behaviour."
  )
}

#-----------------------------------------------------------
# 5. Write config/scrna_module_genesets.csv
#-----------------------------------------------------------

if (!dir.exists("config")) {
  dir.create("config", recursive = TRUE)
}

out_df <- gene_symbol_df %>%
  dplyr::select(module_name, gene_id)

readr::write_csv(out_df, out_geneset_csv)

message("[OK] Wrote gene set CSV for scRNA module scores to: ", out_geneset_csv)
message("[DONE] build_scrna_module_genesets_from_yaml.R completed.")
