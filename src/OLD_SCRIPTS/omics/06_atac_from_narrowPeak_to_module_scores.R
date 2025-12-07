#!/usr/bin/env Rscript

# ===========================================================
# 06_atac_from_narrowPeak_to_module_scores.R
#
# Purpose (real ATAC-seq pipeline for this project):
#   1) Read GSE190414 narrowPeak files:
#        - HFSC, IFE, PAP_mneg, PAP_mpos, SCC_mneg, SCC_mpos
#   2) Annotate peaks to nearest genes using mm10 TxDb + org.Mm.eg.db.
#   3) Collapse peaks to gene-level accessibility per sample
#        (log2-transformed signalValue).
#   4) Z-score each gene across ATAC samples.
#   5) Compute module scores for:
#        - TGFb_bulk, mTOR_bulk, Angio_bulk, CSC_bulk
#      using config/gene_sets_rnaseq.yaml
#   6) Write:
#        - data/interim/atac/atac_gene_accessibility.tsv
#        - data/interim/atac/sample_metadata_ATAC.csv
#        - data/processed/atac/atac_module_scores_by_sample.csv
#
# This uses Bioconductor:
#   GenomicRanges, rtracklayer, TxDb.Mmusculus.UCSC.mm10.knownGene,
#   org.Mm.eg.db, ChIPseeker
#
# Install (once) in R:
#   if (!require("BiocManager")) install.packages("BiocManager")
#   BiocManager::install(c(
#     "GenomicRanges",
#     "rtracklayer",
#     "TxDb.Mmusculus.UCSC.mm10.knownGene",
#     "org.Mm.eg.db",
#     "ChIPseeker"
#   ))
# ===========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(yaml)

  library(GenomicRanges)
  library(rtracklayer)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(ChIPseeker)
})

#------------------------- paths ---------------------------

root_atac <- "data/raw/GSE190414"
expr_out_dir <- "data/interim/atac"
proc_out_dir <- "data/processed/atac"
geneset_yml  <- "config/gene_sets_rnaseq.yaml"

if (!dir.exists(root_atac)) {
  stop("[ERROR] ATAC raw dir not found: ", root_atac)
}
if (!dir.exists(expr_out_dir)) {
  message("[INFO] Creating interim ATAC dir: ", expr_out_dir)
  dir.create(expr_out_dir, recursive = TRUE)
}
if (!dir.exists(proc_out_dir)) {
  message("[INFO] Creating processed ATAC dir: ", proc_out_dir)
  dir.create(proc_out_dir, recursive = TRUE)
}
if (!file.exists(geneset_yml)) {
  stop("[ERROR] Gene set YAML not found: ", geneset_yml)
}

#----------------- sample / file mapping -------------------

# These must match exactly what is in data/raw/GSE190414
samples <- tibble::tribble(
  ~sample_id,              ~npeak_file,                                         ~condition,
  "HFSC_ATAC_Rep1",        "GSE190414_HFSC_ATAC_Rep1.narrowPeak.gz",           "HFSC",
  "IFE_ATAC_Rep1",         "GSE190414_IFE_ATAC_Rep1.narrowPeak.gz",            "IFE",
  "PAP_mneg_ATAC_POOL",    "GSE190414_PAP_mneg_ATAC_Reps_POOL.narrowPeak.gz", "Papilloma_mneg",
  "PAP_mpos_ATAC_POOL",    "GSE190414_PAP_mpos_ATAC_Reps_POOL.narrowPeak.gz", "Papilloma_mpos",
  "SCC_mneg_ATAC_POOL",    "GSE190414_SCC_mneg_ATAC_Reps_POOL.narrowPeak.gz", "SCC_mneg",
  "SCC_mpos_ATAC_POOL",    "GSE190414_SCC_mpos_ATAC_Reps_POOL.narrowPeak.gz", "SCC_mpos"
) %>%
  mutate(
    path = file.path(root_atac, npeak_file),
    dataset = "ATAC"
  )

for (i in seq_len(nrow(samples))) {
  if (!file.exists(samples$path[i])) {
    stop("[ERROR] Missing narrowPeak file: ", samples$path[i])
  }
}

#--------------------- TxDb / annotation -------------------

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# promoter definition: -2kb to +2kb around TSS
promoter_upstream  <- 2000
promoter_downstream <- 2000

#---------------------- gene sets -------------------------

message("[STEP] Loading gene sets from YAML: ", geneset_yml)
gene_sets <- yaml::read_yaml(geneset_yml)

if (length(gene_sets) == 0L) {
  stop("[ERROR] No gene sets found in YAML.")
}
msg_modules <- paste(names(gene_sets), collapse = ", ")
message("[INFO] Gene sets: ", msg_modules)

#---------------- helper: narrowPeak -> gene table --------

annotate_sample_peaks <- function(sample_row, txdb) {
  sample_id <- sample_row$sample_id
  peak_file <- sample_row$path

  message("[SAMPLE] Importing peaks for: ", sample_id)

  # Read narrowPeak
  gr <- rtracklayer::import(peak_file, format = "narrowPeak")

  if (length(gr) == 0L) {
    warning("[WARN] No peaks found in file: ", peak_file)
    return(NULL)
  }

  # Annotate peaks to nearest genes (mm10)
  peak_anno <- ChIPseeker::annotatePeak(
    gr,
    TxDb = txdb,
    tssRegion = c(-promoter_upstream, promoter_downstream),
    annoDb = "org.Mm.eg.db",
    verbose = FALSE
  )

  anno_df <- as.data.frame(peak_anno)

  # Expect SYMBOL column for gene symbols
  if (!"SYMBOL" %in% colnames(anno_df)) {
    stop("[ERROR] Annotation result missing SYMBOL column for sample: ", sample_id)
  }

  # Use signalValue as accessibility measure
  # If missing, fall back to score
  if ("signalValue" %in% colnames(anno_df)) {
    signal_col <- "signalValue"
  } else if ("score" %in% colnames(anno_df)) {
    signal_col <- "score"
  } else {
    stop("[ERROR] narrowPeak lacks 'signalValue' or 'score' columns for sample: ", sample_id)
  }

  gene_access <- anno_df %>%
    dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(
      access_raw = mean(.data[[signal_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      gene_symbol = SYMBOL,
      access_log2 = log2(access_raw + 1),
      sample_id = sample_id
    ) %>%
    dplyr::select(sample_id, gene_symbol, access_log2)

  gene_access
}

#----------------- build gene x sample matrix --------------

message("[STEP] Annotating all ATAC samples and building gene x sample matrix")

access_list <- list()

for (i in seq_len(nrow(samples))) {
  res <- annotate_sample_peaks(samples[i, ], txdb)
  if (!is.null(res)) {
    access_list[[samples$sample_id[i]]] <- res
  }
}

if (length(access_list) == 0L) {
  stop("[ERROR] No annotated ATAC peaks found for any sample.")
}

access_long <- bind_rows(access_list)

# Wide: genes x samples
access_wide <- access_long %>%
  tidyr::pivot_wider(
    id_cols = gene_symbol,
    names_from = sample_id,
    values_from = access_log2
  )

# Replace NAs (genes lacking peaks in a given sample) with 0
access_wide[is.na(access_wide)] <- 0

out_access <- file.path(expr_out_dir, "atac_gene_accessibility.tsv")
message("[OK] Writing gene-level ATAC accessibility to: ", out_access)
readr::write_tsv(access_wide, out_access)

#----------------- write ATAC metadata ---------------------

meta_atac <- samples %>%
  dplyr::select(sample_id, dataset, condition)

meta_out <- file.path(expr_out_dir, "sample_metadata_ATAC.csv")
readr::write_csv(meta_atac, meta_out)

message("[OK] Wrote ATAC sample metadata to: ", meta_out)

#----------------- compute ATAC module scores --------------

message("[STEP] Computing ATAC module scores from accessibility")

expr_df <- access_wide
meta_df <- meta_atac

# Long format for z-scoring
expr_long <- expr_df %>%
  tidyr::pivot_longer(
    cols = -gene_symbol,
    names_to = "sample_id",
    values_to = "access"
  )

# Ensure metadata coverage
unknown_samples <- setdiff(unique(expr_long$sample_id), meta_df$sample_id)
if (length(unknown_samples) > 0L) {
  stop(
    "[ERROR] Some ATAC samples not found in metadata: ",
    paste(unknown_samples, collapse = ", ")
  )
}

# Z-score per gene across ATAC samples
expr_long <- expr_long %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::mutate(
    access_z = (access - mean(access, na.rm = TRUE)) /
               sd(access, na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    access_z = ifelse(is.na(access_z), 0, access_z)
  )

# Module scores
module_scores_list <- list()

for (mod in names(gene_sets)) {
  genes_mod <- gene_sets[[mod]]

  present_genes <- intersect(genes_mod, unique(expr_long$gene_symbol))

  if (length(present_genes) == 0L) {
    warning("[WARN] Module ", mod, " has no genes present in ATAC accessibility.")
    next
  }

  mod_scores <- expr_long %>%
    dplyr::filter(gene_symbol %in% present_genes) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(
      score = mean(access_z, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      module = mod
    )

  module_scores_list[[mod]] <- mod_scores
}

if (length(module_scores_list) == 0L) {
  stop("[ERROR] No ATAC module scores computed (no overlap with gene sets).")
}

scores_all <- bind_rows(module_scores_list)

# Join metadata
scores_all <- scores_all %>%
  dplyr::left_join(meta_df, by = "sample_id")

# Pivot modules wide
scores_wide <- scores_all %>%
  dplyr::select(sample_id, dataset, condition, module, score) %>%
  tidyr::pivot_wider(
    names_from = module,
    values_from = score
  )

out_scores <- file.path(proc_out_dir, "atac_module_scores_by_sample.csv")
readr::write_csv(scores_wide, out_scores)

message("[OK] Wrote ATAC module scores to: ", out_scores)
message("[DONE] ATAC pipeline (narrowPeak → gene → modules) completed.\n")
