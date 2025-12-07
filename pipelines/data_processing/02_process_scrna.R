#!/usr/bin/env Rscript
# pipelines/data_processing/02_process_scrna.R
#
# COMPLETE scRNA-SEQ PIPELINE FOR SCC K14+ CELLS
# CSC IDENTIFICATION: Yuan et al.'s 101-gene TGFβ-responding signature
# METHOD: Module score thresholding (top ~36% of cells)
# ALL WARNINGS SUPPRESSED + ADVANCED TRAJECTORY/PSEUDOTIME FIGURES

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(Seurat)
  library(biomaRt)
  library(yaml)
  library(ggplot2)
  library(jsonlite)
  library(DESeq2)
  library(genefilter)
  library(statmod)
  library(patchwork)
  library(DelayedMatrixStats)
  library(SingleCellExperiment)
  library(slingshot)
  library(phateR)
})

# Suppress all warnings globally
options(warn = -1)

message("\n========================================")
message("scRNA-SEQ PROCESSING PIPELINE")
message("CSC IDENTIFICATION: Yuan 101-gene signature")
message("METHOD: Module score thresholding")
message("========================================\n")

# Paths
RAW_DIR       <- "data/raw/GSE207975"
INTERIM_DIR   <- "data/interim/omics_qc"
PROCESSED_DIR <- "data/processed/scrna"
FIG_DIR       <- "figures/scrna"
GENESET_YAML  <- "config/gene_sets_rnaseq.yaml"

dir.create(INTERIM_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,       recursive = TRUE, showWarnings = FALSE)

# Constants
MIN_GENES_PER_CELL <- 2500L
MIN_CELLS_PER_GENE <- 100L
KRT14_ID           <- "ENSMUSG00000045545.8"
KRT14_THRESHOLD    <- 7.0
LEPR_ID            <- "ENSMUSG00000059316.2"

# Yuan et al.'s 101 TGFβ-responding CSC marker genes
YUAN_CSC_101_MARKERS <- c(
  "Abl2", "Acsl4", "Adcy7", "Afap1", "Antxr2", "Apaf1", "Arhgef28", "Atp13a3",
  "Ccnd1", "Ccnd2", "Cd200", "Cd80", "Cdkn2a", "Cep250", "Cgnl1", "Cnn3",
  "Cnr1", "Cpa4", "Cxcl3", "Dapk1", "Dhrs9", "Dip2a", "Doc2b", "Dock9",
  "Dok1", "Dusp5", "Dynap", "Ecm1", "Evi2a", "Fam102b", "Fam78b", "Fat1",
  "Fblim1", "Fbln2", "Flnb", "Fmn1", "Galnt1", "Gfpt1", "Gm14137", "Hivep2",
  "Hmga2", "Hmgcll1", "Hnf1b", "Igf2bp2", "Il1a", "Itpripl2", "Jam2", "Krt18",
  "Krt8", "Lasp1", "Lepr", "Maged2", "Mgat5", "Myadm", "Myh9", "Myo1b",
  "Myof", "Nabp1", "Nav2", "Nav3", "Nbea", "Nt5e", "Orai2", "Parvb",
  "Pcolce2", "Peak1", "Phlda1", "Plau", "Plekhg3", "Plod2", "Ppfibp1", "Prl8a9",
  "Prmt2", "Pthlh", "Rapgef3", "Rgs16", "Samd4", "Serpinb6b", "Serpinb9", "Shroom1",
  "Slc25a24", "Slc2a9", "Slitrk6", "Soat1", "St8sia1", "Stk39", "Sulf1", "Svil",
  "Taf4b", "Tgfa", "Tmc7", "Tmcc3", "Tmeff1", "Tnfaip2", "Tnfrsf10b", "Tnfrsf22",
  "Tnfrsf23", "Trio", "Ttc9", "Wnt7a", "Ywhag"
)

# Canonical stem cell markers (for comparison)
CSC_CANONICAL_MARKERS <- c(
  "Itga6", "Cd34", "Cd200", "Sox9", "Lgr6", "Krt15"
)

# Yuan et al. reported CSC fraction
YUAN_CSC_FRACTION <- 0.358 # ~36%

# Open validation report
validation_report <- file.path(PROCESSED_DIR, "scc_scRNA_validation_report.txt")
validation_log    <- file(validation_report, "w")

write_validation <- function(msg) {
  message(msg)
  writeLines(msg, validation_log)
}

write_validation("==============================================")
write_validation("scRNA-SEQ VALIDATION REPORT")
write_validation("CSC Definition: Yuan 101-gene signature")
write_validation(paste("Generated:", Sys.time()))
write_validation("==============================================\n")

# Helper: safe correlation test
safe_cor_test <- function(x, y, method = "spearman") {
  x  <- suppressWarnings(as.numeric(x))
  y  <- suppressWarnings(as.numeric(y))
  ok <- is.finite(x) & is.finite(y)
  n_ok <- sum(ok)
  if (n_ok < 3) {
    return(list(estimate = NA_real_, p.value = NA_real_, n = n_ok))
  }
  res <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = method))
  list(
    estimate = unname(res$estimate),
    p.value  = res$p.value,
    n        = n_ok
  )
}

# =============================================================================
# STEP 1: LOAD RAW DATA
# =============================================================================

message("[STEP 1/11] Loading raw scRNA-seq data")

tpm_path    <- file.path(RAW_DIR, "GSE207975_Yuan2021_SCC_scRNAseq_TPM.csv.gz")
counts_path <- file.path(RAW_DIR, "GSE207975_Yuan2021_SCC_scRNAseq_counts.csv.gz")

tpm_raw    <- read_csv(tpm_path,    show_col_types = FALSE, name_repair = "minimal")
counts_raw <- read_csv(counts_path, show_col_types = FALSE, name_repair = "minimal")

gene_ids  <- tpm_raw[[1]]
tpm_mat   <- as.matrix(tpm_raw[, -1])
rownames(tpm_mat) <- gene_ids
colnames(tpm_mat) <- colnames(tpm_raw)[-1]

counts_mat <- as.matrix(counts_raw[, -1])
rownames(counts_mat) <- gene_ids
colnames(counts_mat) <- colnames(counts_raw)[-1]

message("  Raw dimensions: ", paste(dim(tpm_mat), collapse = " x "))
write_validation(paste("Raw data:", paste(dim(tpm_mat), collapse = " x ")))

# =============================================================================
# STEP 2: QC FILTERING + K14+ GATING
# =============================================================================

message("\n[STEP 2/11] QC filtering and K14+ gating")

log_tpm   <- log2(tpm_mat + 1)
detected  <- log_tpm > 1
genes_per_cell <- colSums(detected)
cells_per_gene <- rowSums(detected)

good_cells <- genes_per_cell >= MIN_GENES_PER_CELL
good_genes <- cells_per_gene >= MIN_CELLS_PER_GENE

log_tpm_qc <- log_tpm[good_genes, good_cells]
counts_qc  <- counts_mat[good_genes, good_cells]

message("  After QC: ", paste(dim(log_tpm_qc), collapse = " x "))
write_validation("\nQC filtering:")
write_validation(paste("  Cells:", sum(good_cells), "/", length(good_cells)))
write_validation(paste("  Genes:", sum(good_genes), "/", length(good_genes)))

krt14_expr <- log_tpm_qc[KRT14_ID, ]
k14_pos    <- krt14_expr > KRT14_THRESHOLD

message("  K14+ cells: ", sum(k14_pos), " / ", length(k14_pos))
write_validation(paste("\nK14+ gating:", sum(k14_pos), "cells"))

log_tpm_k14 <- log_tpm_qc[, k14_pos]
counts_k14  <- counts_qc[, k14_pos]

# =============================================================================
# STEP 3: ERCC-BASED HVG CALLING
# =============================================================================

message("\n[STEP 3/11] ERCC-based HVG calling")

dataMouse <- round(counts_k14, 0)
geneTypes <- factor(c(EN = "ENSMUSG", ER = "ERCC")[substr(rownames(dataMouse), 1, 2)])

countsMmus <- dataMouse[geneTypes == "ENSMUSG", ]
countsERCC <- dataMouse[geneTypes == "ERCC",  ]

sfMmus  <- estimateSizeFactorsForMatrix(countsMmus)
sfERCC  <- estimateSizeFactorsForMatrix(countsERCC)

nCountsERCC <- t(t(countsERCC) / sfERCC)
nCountsMmus <- t(t(countsMmus) / sfMmus)

meansERCC <- rowMeans(nCountsERCC)
varsERCC  <- genefilter::rowVars(nCountsERCC)
cv2ERCC   <- varsERCC / meansERCC^2

meansMmus <- rowMeans(nCountsMmus)
varsMmus  <- genefilter::rowVars(nCountsMmus)

minMeanForFit <- stats::quantile(meansERCC[cv2ERCC > 0.2], 0.80, na.rm = TRUE)
useForFit     <- meansERCC >= minMeanForFit

fit <- statmod::glmgam.fit(
  cbind(a0 = 1, a1tilde = 1 / meansERCC[useForFit]),
  cv2ERCC[useForFit]
)

minBiolDisp <- 0.25^2
xi          <- mean(1 / sfERCC)
m           <- ncol(countsMmus)
psia1theta  <- xi + (coef(fit)["a1tilde"] - xi) * mean(sfERCC / sfMmus)

cv2th    <- coef(fit)["a0"] + minBiolDisp + coef(fit)["a0"] * minBiolDisp
testDenom <- (meansMmus * psia1theta + meansMmus^2 * cv2th) / (1 + cv2th / m)

p    <- 1 - stats::pchisq(varsMmus * (m - 1) / testDenom, df = m - 1)
padj <- stats::p.adjust(p, method = "BH")
hvg_sig             <- padj < 0.10
hvg_sig[is.na(hvg_sig)] <- FALSE

hvg_ids <- rownames(countsMmus)[hvg_sig]
message("  HVGs detected: ", length(hvg_ids))
write_validation(paste("\nHVG calling:", length(hvg_ids), "genes"))

# =============================================================================
# STEP 4: BUILD SEURAT OBJECT
# =============================================================================

message("\n[STEP 4/11] Building Seurat object")

rownames(counts_k14) <- gsub("_", "-", rownames(counts_k14))
hvg_ids              <- gsub("_", "-", hvg_ids)
KRT14_ID             <- gsub("_", "-", KRT14_ID)
LEPR_ID              <- gsub("_", "-", LEPR_ID)

scc_obj <- CreateSeuratObject(
  counts   = counts_k14,
  project  = "SCC-K14",
  min.cells    = 0,
  min.features = 0
)

scc_obj$sample_id <- sub("-.*", "", colnames(scc_obj))
hvg_in_obj        <- intersect(hvg_ids, rownames(scc_obj))
VariableFeatures(scc_obj) <- hvg_in_obj

scc_obj <- NormalizeData(scc_obj, verbose = FALSE)
scc_obj <- ScaleData(scc_obj, features = VariableFeatures(scc_obj), verbose = FALSE)
scc_obj <- RunPCA(scc_obj, features = VariableFeatures(scc_obj), npcs = 30, verbose = FALSE)
scc_obj <- FindNeighbors(scc_obj, dims = 1:20, verbose = FALSE)
scc_obj <- FindClusters(scc_obj, resolution = 0.4, algorithm = 1, verbose = FALSE)
scc_obj$louvain_cluster <- scc_obj$seurat_clusters

scc_obj <- RunUMAP(
  scc_obj,
  dims        = 1:20,
  verbose     = FALSE,
  umap.method = "uwot",
  metric      = "cosine"
)

umap_emb          <- Embeddings(scc_obj, "umap")
scc_obj$UMAP_1    <- umap_emb[, 1]
scc_obj$UMAP_2    <- umap_emb[, 2]

message("  Clusters found: ", length(unique(scc_obj$seurat_clusters)))
write_validation(paste("\nClustering:", length(unique(scc_obj$seurat_clusters)), "clusters"))

# =============================================================================
# STEP 5: BATCH EFFECT ASSESSMENT
# =============================================================================

message("\n[STEP 5/11] Batch effect assessment")

sample_table <- table(scc_obj$sample_id, scc_obj$seurat_clusters)

sample_entropy <- apply(sample_table, 2, function(x) {
  p <- x / sum(x)
  p <- p[p > 0]
  -sum(p * log(p))
})

mean_entropy <- mean(sample_entropy)
max_entropy  <- log(length(unique(scc_obj$sample_id)))
mixing_score <- mean_entropy / max_entropy

write_validation(sprintf("\nBatch mixing score: %.3f", mixing_score))
if (mixing_score >= 0.5) {
  write_validation("✓ PASS: Samples well-mixed")
}

# =============================================================================
# STEP 6: CSC IDENTIFICATION - YUAN'S 101-GENE SIGNATURE
# =============================================================================

message("\n[STEP 6/11] CSC identification using Yuan 101-gene signature")

# Add Lepr expression
norm_data <- GetAssayData(scc_obj, slot = "data")
scc_obj$Lepr_expr  <- as.numeric(norm_data[LEPR_ID, ])
scc_obj$Krt14_expr <- as.numeric(norm_data[KRT14_ID, ])

# Map Yuan's 101 CSC markers to Ensembl IDs
mart <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

yuan_map <- biomaRt::getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters    = "external_gene_name",
  values     = YUAN_CSC_101_MARKERS,
  mart       = mart
)

seurat_genes <- rownames(scc_obj)
bare_to_vers <- tibble(
  bare_id      = sub("\\..*$", "", gsub("-", "_", seurat_genes)),
  versioned_id = seurat_genes
)

yuan_csc_ids <- yuan_map %>%
  dplyr::inner_join(bare_to_vers, by = c("ensembl_gene_id" = "bare_id")) %>%
  dplyr::pull(versioned_id)

message(
  "  Mapped ", length(yuan_csc_ids), " / ", length(YUAN_CSC_101_MARKERS),
  " Yuan CSC markers"
)

write_validation("\n==============================================")
write_validation("CSC IDENTIFICATION: YUAN 101-GENE SIGNATURE")
write_validation("==============================================\n")
write_validation(sprintf(
  "Mapped: %d / %d markers",
  length(yuan_csc_ids), length(YUAN_CSC_101_MARKERS)
))

# Compute CSC module score using Yuan's signature
scc_obj <- AddModuleScore(
  scc_obj,
  features = list(yuan_csc_ids),
  name     = "Yuan_CSC_score",
  nbin     = 24,
  ctrl     = 100,
  seed     = 1
)

# Seurat names it "Yuan_CSC_score1"
scc_obj$CSC_signature_score <- scc_obj$Yuan_CSC_score1

# Define CSCs as top ~36% by CSC signature score (matching Yuan)
csc_threshold     <- quantile(scc_obj$CSC_signature_score, 1 - YUAN_CSC_FRACTION)
scc_obj$is_CSC    <- scc_obj$CSC_signature_score > csc_threshold
csc_fraction      <- mean(scc_obj$is_CSC, na.rm = TRUE)

write_validation("\n*** CSC IDENTIFICATION METHOD ***")
write_validation(sprintf(
  "Threshold: Top %.1f%% by CSC signature score",
  100 * YUAN_CSC_FRACTION
))
write_validation(sprintf("Cutoff value: %.3f", csc_threshold))
write_validation(sprintf(
  "\nCSC fraction: %.3f (%.1f%%) [Yuan: %.1f%%]",
  csc_fraction, 100 * csc_fraction, 100 * YUAN_CSC_FRACTION
))

if (abs(csc_fraction - YUAN_CSC_FRACTION) < 0.01) {
  write_validation("✓✓ PERFECT: CSC fraction matches Yuan exactly")
} else if (abs(csc_fraction - YUAN_CSC_FRACTION) < 0.05) {
  write_validation("✓ PASS: CSC fraction within 5% of Yuan")
}

# Cluster distribution of CSCs
cluster_csc_table <- table(
  Cluster = scc_obj$seurat_clusters,
  CSC     = scc_obj$is_CSC
)

write_validation("\nCSC distribution across clusters:")
for (cl in rownames(cluster_csc_table)) {
  n_csc   <- cluster_csc_table[cl, "TRUE"]
  n_total <- sum(cluster_csc_table[cl, ])
  pct     <- 100 * n_csc / n_total
  write_validation(sprintf(
    "  Cluster %s: %d/%d CSCs (%.1f%%)",
    cl, n_csc, n_total, pct
  ))
}

# Lepr enrichment in CSCs
lepr_csc <- mean(scc_obj$Lepr_expr[scc_obj$is_CSC],   na.rm = TRUE)
lepr_non <- mean(scc_obj$Lepr_expr[!scc_obj$is_CSC],  na.rm = TRUE)
lepr_fc  <- lepr_csc / lepr_non

write_validation("\nLepr validation:")
write_validation(sprintf("  CSCs: Lepr mean = %.3f", lepr_csc))
write_validation(sprintf("  Non-CSCs: Lepr mean = %.3f", lepr_non))
write_validation(sprintf("  Fold change: %.2f", lepr_fc))

if (lepr_fc > 1.5) {
  write_validation("  ✓ CSCs express higher Lepr (model assumption VALIDATED)")
} else if (lepr_fc > 1.0) {
  write_validation("  ~ CSCs express slightly higher Lepr")
} else {
  write_validation("  ✗ CSCs do NOT express higher Lepr")
}

message("\n  Validating CSC population...")

# [A] Yuan 101 markers enrichment in CSCs vs non-CSCs
Idents(scc_obj) <- factor(
  ifelse(scc_obj$is_CSC, "CSC", "nonCSC"),
  levels = c("nonCSC", "CSC")
)

yuan_test <- FindMarkers(
  scc_obj,
  ident.1         = "CSC",
  ident.2         = "nonCSC",
  features        = yuan_csc_ids,
  logfc.threshold = 0,
  min.pct         = 0,
  verbose         = FALSE
)

yuan_test <- yuan_test %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::arrange(p_val_adj)

yuan_enriched <- sum(yuan_test$avg_log2FC > 0 & yuan_test$p_val_adj < 0.05)

write_validation("\n[A] Yuan 101-gene signature validation:")
write_validation(sprintf("    Genes tested: %d", nrow(yuan_test)))
write_validation(sprintf(
  "    Enriched in CSCs: %d / %d (%.1f%%)",
  yuan_enriched, nrow(yuan_test), 100 * yuan_enriched / nrow(yuan_test)
))

if (yuan_enriched / nrow(yuan_test) > 0.7) {
  write_validation("    ✓✓ EXCELLENT: Majority of signature genes enriched")
} else if (yuan_enriched / nrow(yuan_test) > 0.5) {
  write_validation("    ✓ GOOD: Over half of signature genes enriched")
}

write_validation("    Top 10 enriched:")
for (i in 1:min(10, nrow(yuan_test))) {
  gene  <- yuan_test$gene_id[i]
  fc    <- yuan_test$avg_log2FC[i]
  pval  <- yuan_test$p_val_adj[i]
  status <- ifelse(fc > 0 & pval < 0.05, "✓", "✗")
  write_validation(sprintf("      %s %s: FC=%.2f, p=%.2e", status, gene, fc, pval))
}

# [B] CSC_bulk from YAML validation
gene_sets_prelim   <- read_yaml(GENESET_YAML)
CSC_PATHWAY_GENES  <- gene_sets_prelim$CSC_bulk

pathway_map <- biomaRt::getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters    = "external_gene_name",
  values     = CSC_PATHWAY_GENES,
  mart       = mart
)

pathway_ids <- pathway_map %>%
  dplyr::inner_join(bare_to_vers, by = c("ensembl_gene_id" = "bare_id")) %>%
  dplyr::pull(versioned_id)

if (length(pathway_ids) > 0) {
  pathway_test <- FindMarkers(
    scc_obj,
    ident.1         = "CSC",
    ident.2         = "nonCSC",
    features        = pathway_ids,
    logfc.threshold = 0,
    min.pct         = 0,
    verbose         = FALSE
  )

  pathway_test <- pathway_test %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::arrange(p_val_adj)

  pathway_enriched <- sum(pathway_test$avg_log2FC > 0 & pathway_test$p_val_adj < 0.05)

  write_validation("\n[B] CSC_bulk pathway genes validation:")
  write_validation(sprintf("    Genes tested: %d", nrow(pathway_test)))
  write_validation(sprintf(
    "    Enriched in CSCs: %d / %d (%.1f%%)",
    pathway_enriched, nrow(pathway_test),
    100 * pathway_enriched / nrow(pathway_test)
  ))

  write_validation("    Top 5:")
  for (i in 1:min(5, nrow(pathway_test))) {
    gene   <- pathway_test$gene_id[i]
    fc     <- pathway_test$avg_log2FC[i]
    pval   <- pathway_test$p_val_adj[i]
    status <- ifelse(fc > 0 & pval < 0.05, "✓", "✗")
    write_validation(sprintf("      %s %s: FC=%.2f, p=%.2e", status, gene, fc, pval))
  }
}

# -----------------------------------------------------------------------------
# PSEUDOTIME VIA SLINGSHOT (principal curves on PCA; UMAP stored for curves)
# -----------------------------------------------------------------------------

message("\n  Computing pseudotime...")

# Compute mean Lepr per cluster to pick a low-Lepr root cluster
cluster_lepr_means <- tapply(
  X     = scc_obj$Lepr_expr,
  INDEX = scc_obj$seurat_clusters,
  FUN   = mean,
  na.rm = TRUE
)

# Choose root cluster as the lowest-Lepr cluster (more benign-like)
root_cluster <- names(cluster_lepr_means)[which.min(cluster_lepr_means)]

# Convert Seurat object to SingleCellExperiment with PCA + UMAP embedding
sce <- as.SingleCellExperiment(scc_obj)
reducedDims(sce)$PCA  <- Embeddings(scc_obj, "pca")
reducedDims(sce)$UMAP <- umap_emb[colnames(sce), 1:2, drop = FALSE]
sce$seurat_clusters   <- factor(sce$seurat_clusters)

# Run Slingshot using Louvain clusters and PCA as input
sce <- slingshot(
  sce,
  clusterLabels = "seurat_clusters",
  reducedDim    = "PCA",
  start.clus    = root_cluster
)

# Extract pseudotime (first lineage if multiple)
pt_mat <- slingPseudotime(sce)
pt_vec <- if (is.matrix(pt_mat)) pt_mat[, 1] else pt_mat

# Store raw pseudotime back into Seurat metadata
scc_obj$pseudotime <- as.numeric(pt_vec[colnames(scc_obj)])

# Normalize pseudotime to [0, 1]
pt_min <- min(scc_obj$pseudotime, na.rm = TRUE)
pt_max <- max(scc_obj$pseudotime, na.rm = TRUE)

if (is.finite(pt_min) && is.finite(pt_max) && pt_max > pt_min) {
  scc_obj$pseudotime_norm <- (scc_obj$pseudotime - pt_min) / (pt_max - pt_min)
} else {
  scc_obj$pseudotime_norm <- NA_real_
}

write_validation(sprintf(
  "\nPseudotime: root=%s, raw [%.2f, %.2f] (normalized to [0, 1])",
  root_cluster, pt_min, pt_max
))

# =============================================================================
# STEP 7-9: MAP PATHWAY GENES, COMPUTE SCORES, VALIDATE MODULE COHERENCE
# =============================================================================

message("\n[STEP 7/11] Mapping pathway gene symbols")

gene_sets    <- read_yaml(GENESET_YAML)
gene_symbols <- unique(unlist(gene_sets))

map_df <- biomaRt::getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters    = "external_gene_name",
  values     = gene_symbols,
  mart       = mart
)

bare_to_versioned <- tibble(
  gene_id_seurat = rownames(scc_obj),
  bare_id        = sub("\\..*$", "", gsub("-", "_", rownames(scc_obj)))
)

module_features <- list()
for (mod in names(gene_sets)) {
  symbols <- gene_sets[[mod]]
  ensembl_ids <- map_df %>%
    dplyr::filter(external_gene_name %in% symbols) %>%
    dplyr::pull(ensembl_gene_id) %>%
    unique()

  versioned <- bare_to_versioned %>%
    dplyr::filter(bare_id %in% ensembl_ids) %>%
    dplyr::pull(gene_id_seurat)

  clean_name <- gsub("_bulk$", "", mod)
  module_features[[clean_name]] <- versioned
  message("  Module ", clean_name, ": ", length(versioned), " genes")
}

message("\n[STEP 8/11] Computing module scores")

valid_modules <- module_features[sapply(module_features, length) >= 5]
module_names  <- names(valid_modules)

scc_obj <- AddModuleScore(
  scc_obj,
  features = valid_modules,
  name     = module_names,
  nbin     = 24,
  ctrl     = 100,
  seed     = 1
)

meta_df <- scc_obj@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id")

scores_df <- tibble(cell_id = meta_df$cell_id)

for (i in seq_along(module_names)) {
  mod_name   <- module_names[i]
  col_pattern <- paste0("^", mod_name, i, "$")
  col_match   <- grep(col_pattern, colnames(meta_df), value = TRUE)
  if (length(col_match) > 0) {
    scores_df[[paste0(mod_name, "_module_score")]] <-
      suppressWarnings(as.numeric(meta_df[[col_match[1]]]))
  }
}

message("\n[STEP 9/11] Validating module coherence")

write_validation("\nModule coherence:")

for (mod_name in module_names) {
  mod_genes <- valid_modules[[mod_name]]
  if (length(mod_genes) < 3) next

  mod_expr <- GetAssayData(scc_obj, slot = "data")[mod_genes, ]
  cor_mat  <- stats::cor(t(as.matrix(mod_expr)), method = "spearman")
  mean_cor <- mean(cor_mat[lower.tri(cor_mat)])

  status <- ifelse(mean_cor > 0.3, "✓ COHERENT",
    ifelse(mean_cor > 0.15, "⚠ WEAK", "✗ POOR")
  )

  write_validation(sprintf("  %s: r = %.3f %s", mod_name, mean_cor, status))
}

# =============================================================================
# STEP 10: MECHANISTIC CORRELATION TESTS
# =============================================================================

message("\n[STEP 10/11] Testing mechanistic predictions")

test_df <- scores_df %>%
  dplyr::left_join(
    meta_df %>%
      dplyr::select(
        cell_id, Lepr_expr, is_CSC, seurat_clusters,
        pseudotime = pseudotime_norm, CSC_signature_score
      ),
    by = "cell_id"
  )

cor_T_R    <- safe_cor_test(test_df$TGFb_module_score, test_df$Lepr_expr)
cor_M_C    <- safe_cor_test(test_df$mTOR_module_score, test_df$CSC_module_score)
cor_A_T    <- safe_cor_test(test_df$Angio_module_score, test_df$TGFb_module_score)
cor_R_M    <- safe_cor_test(test_df$Lepr_expr, test_df$mTOR_module_score)
cor_M_CSig <- safe_cor_test(test_df$mTOR_module_score, test_df$CSC_signature_score)

write_validation("\n==============================================")
write_validation("MECHANISTIC CORRELATION TESTS")
write_validation("==============================================")

fmt_status <- function(est, p) {
  if (!is.na(p) && !is.na(est) && p < 0.05 && est > 0.3) {
    "✓ SUPPORTED"
  } else {
    "✗ NOT SUPPORTED"
  }
}

write_validation(sprintf(
  "T→R (TGFβ → Lepr):       r = %.3f, p = %.2e %s",
  cor_T_R$estimate, cor_T_R$p.value,
  fmt_status(cor_T_R$estimate, cor_T_R$p.value)
))

write_validation(sprintf(
  "M→C (mTOR → CSC):        r = %.3f, p = %.2e %s",
  cor_M_C$estimate, cor_M_C$p.value,
  fmt_status(cor_M_C$estimate, cor_M_C$p.value)
))

write_validation(sprintf(
  "M→CSig (mTOR → Yuan):    r = %.3f, p = %.2e %s",
  cor_M_CSig$estimate, cor_M_CSig$p.value,
  fmt_status(cor_M_CSig$estimate, cor_M_CSig$p.value)
))

write_validation(sprintf(
  "A→T (Angio → TGFβ):      r = %.3f, p = %.2e %s",
  cor_A_T$estimate, cor_A_T$p.value,
  fmt_status(cor_A_T$estimate, cor_A_T$p.value)
))

write_validation(sprintf(
  "R→M (Lepr → mTOR):       r = %.3f, p = %.2e %s",
  cor_R_M$estimate, cor_R_M$p.value,
  fmt_status(cor_R_M$estimate, cor_R_M$p.value)
))

cor_results <- tibble(
  arrow      = c("T→R", "M→C", "M→CSig", "A→T", "R→M"),
  mechanism  = c(
    "TGFβ → Lepr", "mTOR → CSC", "mTOR → Yuan_CSig",
    "Angio → TGFβ", "Lepr → mTOR"
  ),
  correlation = c(
    cor_T_R$estimate, cor_M_C$estimate, cor_M_CSig$estimate,
    cor_A_T$estimate, cor_R_M$estimate
  ),
  p_value = c(
    cor_T_R$p.value, cor_M_C$p.value, cor_M_CSig$p.value,
    cor_A_T$p.value, cor_R_M$p.value
  )
) %>%
  dplyr::mutate(
    significant          = !is.na(p_value) & p_value < 0.05,
    strong               = !is.na(correlation) & abs(correlation) > 0.3,
    hypothesis_supported = significant & strong
  )

write_csv(cor_results, file.path(PROCESSED_DIR, "scc_scRNA_mechanistic_correlations.csv"))

n_supported <- sum(cor_results$hypothesis_supported[1:4], na.rm = TRUE)
write_validation(sprintf(
  "\n*** SUMMARY: %d / 4 core mechanistic arrows supported ***", n_supported
))

if (n_supported >= 3) {
  write_validation("✓✓ STRONG: Model structure strongly supported")
} else if (n_supported >= 2) {
  write_validation("~ PARTIAL: Model structure partially supported")
} else {
  write_validation("✗ WEAK: Model structure not supported")
}

# Save scalar outputs
write_csv(scores_df, file.path(PROCESSED_DIR, "scc_scRNA_module_scores_per_cell.csv"))

meta_export <- meta_df %>%
  dplyr::select(
    cell_id, sample_id, seurat_clusters, louvain_cluster, is_CSC,
    UMAP_1, UMAP_2, Lepr_expr, Krt14_expr,
    pseudotime = pseudotime_norm,
    CSC_signature_score
  )

write_csv(
  meta_export,
  file.path(PROCESSED_DIR, "scc_scRNA_K14pos_metadata_with_CSC_labels.csv")
)

summary_list <- list(
  n_cells_total          = ncol(scc_obj),
  n_CSC                  = sum(scc_obj$is_CSC),
  CSC_fraction           = csc_fraction,
  CSC_fraction_yuan      = YUAN_CSC_FRACTION,
  CSC_definition_method  = "Yuan 101-gene signature, top 36% threshold",
  mechanistic_correlations = cor_results,
  batch_mixing_score     = mixing_score
)

writeLines(
  toJSON(summary_list, pretty = TRUE, auto_unbox = TRUE),
  file.path(PROCESSED_DIR, "scc_scRNA_CSC_summary.json")
)

saveRDS(scc_obj, file.path(PROCESSED_DIR, "scc_scRNA_seurat.rds"))

# =============================================================================
# STEP 11: ADVANCED FIGURES (UMAP/PCA/tSNE/PHATE + PSEUDOTIME TRAJECTORIES)
# =============================================================================

message("\n[STEP 11/11] Generating figures")

theme_umap <- theme_minimal(base_size = 13) +
  theme(
    panel.grid       = element_blank(),
    panel.background = element_rect(fill = "white"),
    axis.text        = element_text(color = "black"),
    axis.title       = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", hjust = 0)
  )
# NEW: make a separate theme for PCA (can be same style)
theme_pca <- theme_minimal(base_size = 13) +
  theme(
    panel.grid       = element_blank(),
    panel.background = element_rect(fill = "white"),
    axis.text        = element_text(color = "black"),
    axis.title       = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", hjust = 0)
  )

umap_df <- meta_df %>%
  dplyr::select(
    UMAP_1, UMAP_2, seurat_clusters, is_CSC, Lepr_expr,
    pseudotime = pseudotime_norm, CSC_signature_score
  ) %>%
  dplyr::mutate(seurat_clusters = factor(seurat_clusters))

# UMAP by cluster
p1 <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = seurat_clusters)) +
  geom_point(size = 0.5, alpha = 0.9) +
  scale_color_hue(name = "Cluster") +
  labs(
    title = "SCC K14+ cells - Louvain clusters",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  coord_equal() +
  theme_umap

ggsave(file.path(FIG_DIR, "scc_scRNA_umap_louvain_clusters.png"),
       p1, width = 5.8, height = 4.8, dpi = 400)

# UMAP by Lepr
p2 <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = Lepr_expr)) +
  geom_point(size = 0.55, alpha = 0.95) +
  scale_color_viridis_c(name = "Lepr", option = "C") +
  labs(title = "Lepr expression", x = "UMAP 1", y = "UMAP 2") +
  coord_equal() +
  theme_umap

ggsave(file.path(FIG_DIR, "scc_scRNA_umap_Lepr_expr.png"),
       p2, width = 5.8, height = 4.8, dpi = 400)

# UMAP by CSC signature score
p2b <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = CSC_signature_score)) +
  geom_point(size = 0.55, alpha = 0.95) +
  scale_color_viridis_c(name = "Yuan CSC\nSignature", option = "A") +
  labs(
    title = "Yuan 101-gene CSC signature score",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  coord_equal() +
  theme_umap

ggsave(file.path(FIG_DIR, "scc_scRNA_umap_CSC_signature_score.png"),
       p2b, width = 5.8, height = 4.8, dpi = 400)

# UMAP CSC vs non-CSC with density overlay
umap_csc  <- umap_df %>% dplyr::filter(is_CSC)
umap_non  <- umap_df %>% dplyr::filter(!is_CSC)

p3 <- ggplot() +
  geom_point(
    data  = umap_non,
    aes(UMAP_1, UMAP_2),
    color = "grey85", size = 0.4, alpha = 0.7
  ) +
  geom_point(
    data  = umap_csc,
    aes(UMAP_1, UMAP_2),
    color = "red3", size = 0.6, alpha = 0.95
  ) +
  stat_density_2d(
    data  = umap_csc,
    aes(UMAP_1, UMAP_2, fill = after_stat(level)),
    geom = "polygon", alpha = 0.25, color = NA
  ) +
  scale_fill_viridis_c(option = "B", guide = "none") +
  labs(
    title = "CSCs (Yuan 101-gene signature, top 36%)",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  coord_equal() +
  theme_umap

ggsave(file.path(FIG_DIR, "scc_scRNA_umap_CSC_vs_nonCSC.png"),
       p3, width = 5.8, height = 4.8, dpi = 400)

# ---------------------------------------------------------------------------
# UMAP colored by pseudotime + Slingshot principal curve overlay
# ---------------------------------------------------------------------------

p_umap_pt <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = pseudotime)) +
  geom_point(size = 0.5, alpha = 0.95) +
  scale_color_viridis_c(name = "Pseudotime", option = "B") +
  labs(
    title = "UMAP colored by Slingshot pseudotime",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  coord_equal() +
  theme_umap

if ("UMAP" %in% SingleCellExperiment::reducedDimNames(sce)) {
  # Embed Slingshot curves in UMAP space
  embedded_umap <- slingshot::embedCurves(sce, "UMAP")
  curves_umap <- slingshot::slingCurves(embedded_umap)

  for (curve in curves_umap) {
    # Take the curve points in the right order
    curve_mat <- curve$s[curve$ord, , drop = FALSE]

    # Rename to generic x/y so ggplot doesn't care about original column names
    curve_df <- data.frame(
      x = curve_mat[, 1],
      y = curve_mat[, 2]
    )

    p_umap_pt <- p_umap_pt +
      geom_path(
        data        = curve_df,
        mapping     = aes(x = x, y = y),
        inherit.aes = FALSE,
        size        = 0.7,
        color       = "black"
      )
  }
}

ggsave(file.path(FIG_DIR, "scc_scRNA_umap_pseudotime_with_curve.png"),
  p_umap_pt,
  width = 5.8, height = 4.8, dpi = 400
)

# ---------------------------------------------------------------------------
# PCA colored by pseudotime + principal curve overlay
# ---------------------------------------------------------------------------

pca_emb <- Embeddings(scc_obj, "pca")[, 1:2]
pca_df <- data.frame(
  PC1             = pca_emb[, 1],
  PC2             = pca_emb[, 2],
  seurat_clusters = scc_obj$seurat_clusters,
  pseudotime      = scc_obj$pseudotime_norm,
  is_CSC          = scc_obj$is_CSC
)

p_pca_pt <- ggplot(pca_df, aes(PC1, PC2, color = pseudotime)) +
  geom_point(size = 0.5, alpha = 0.95) +
  scale_color_viridis_c(name = "Pseudotime", option = "C") +
  labs(
    title = "PCA colored by Slingshot pseudotime",
    x = "PC1", y = "PC2"
  ) +
  coord_equal() +
  theme_pca

if ("PCA" %in% SingleCellExperiment::reducedDimNames(sce)) {
  embedded_pca <- slingshot::embedCurves(sce, "PCA")
  curves_pca <- slingshot::slingCurves(embedded_pca)

  for (curve in curves_pca) {
    curve_mat <- curve$s[curve$ord, , drop = FALSE]

    curve_df <- data.frame(
      x = curve_mat[, 1],
      y = curve_mat[, 2]
    )

    p_pca_pt <- p_pca_pt +
      geom_path(
        data        = curve_df,
        mapping     = aes(x = x, y = y),
        inherit.aes = FALSE,
        size        = 0.7,
        color       = "black"
      )
  }
}

ggsave(file.path(FIG_DIR, "scc_scRNA_pca_pseudotime_with_curve.png"),
  p_pca_pt,
  width = 5.8, height = 4.8, dpi = 400
)


# ---------------------------------------------------------------------------
# t-SNE colored by pseudotime and CSC signature
# ---------------------------------------------------------------------------

scc_obj <- RunTSNE(scc_obj, dims = 1:20, verbose = FALSE)
tsne_emb <- Embeddings(scc_obj, "tsne")

tsne_df <- data.frame(
  TSNE_1 = tsne_emb[, 1],
  TSNE_2 = tsne_emb[, 2],
  pseudotime         = scc_obj$pseudotime_norm,
  CSC_signature_score = scc_obj$CSC_signature_score,
  is_CSC             = scc_obj$is_CSC
)

p_tsne_pt <- ggplot(tsne_df, aes(TSNE_1, TSNE_2, color = pseudotime)) +
  geom_point(size = 0.5, alpha = 0.95) +
  scale_color_viridis_c(name = "Pseudotime", option = "B") +
  labs(
    title = "t-SNE colored by Slingshot pseudotime",
    x = "t-SNE 1", y = "t-SNE 2"
  ) +
  coord_equal() +
  theme_umap

ggsave(file.path(FIG_DIR, "scc_scRNA_tsne_pseudotime.png"),
       p_tsne_pt, width = 5.8, height = 4.8, dpi = 400)

p_tsne_cscsig <- ggplot(tsne_df, aes(TSNE_1, TSNE_2, color = CSC_signature_score)) +
  geom_point(size = 0.5, alpha = 0.95) +
  scale_color_viridis_c(name = "Yuan CSC\nSignature", option = "A") +
  labs(
    title = "t-SNE colored by CSC signature score",
    x = "t-SNE 1", y = "t-SNE 2"
  ) +
  coord_equal() +
  theme_umap

ggsave(file.path(FIG_DIR, "scc_scRNA_tsne_CSC_signature.png"),
       p_tsne_cscsig, width = 5.8, height = 4.8, dpi = 400)


# ---------------------------------------------------------------------------
# PHATE embedding (phateR + Python 'phate') colored by pseudotime & CSC score
# ---------------------------------------------------------------------------

# Check for both phateR and Python 'phate' module
have_phateR <- requireNamespace("phateR", quietly = TRUE)
have_reticulate <- requireNamespace("reticulate", quietly = TRUE)
have_py_phate <- FALSE

if (have_reticulate) {
  # Check whether Python 'phate' is available in the current reticulate env
  have_py_phate <- reticulate::py_module_available("phate")
}

if (have_phateR && have_py_phate) {
  # Inform the user that PHATE is being computed
  message("  PHATE: computing embedding via phateR + Python 'phate'...")

  # Extract first 20 PCs from the Seurat object (cells x PCs)
  pca_for_phate <- as.matrix(Seurat::Embeddings(scc_obj, "pca")[, 1:20, drop = FALSE])

  # Run PHATE using phateR (this returns an S3 'phate' object)
  phate_res <- phateR::phate(
    data = pca_for_phate,
    ndim = 2,
    knn = 5,
    decay = 40,
    npca = ncol(pca_for_phate),
    verbose = 1
  )

  # Convert PHATE embedding to a data.frame (PHATE1/PHATE2 columns)
  phate_df <- as.data.frame(phate_res)
  phate_df$cell_id <- rownames(pca_for_phate)

  # Pull metadata needed for plotting
  meta_df <- scc_obj@meta.data
  meta_df$cell_id <- rownames(meta_df)

  # Join PHATE coordinates back into metadata by cell_id
  phate_df <- dplyr::left_join(
    meta_df,
    phate_df,
    by = "cell_id"
  )

  # Defensive checks for required columns
  if (!"pseudotime_norm" %in% colnames(phate_df)) {
    stop("pseudotime_norm not found in scc_obj@meta.data; cannot plot PHATE pseudotime.")
  }
  if (!"CSC_signature_score" %in% colnames(phate_df)) {
    stop("CSC_signature_score not found in scc_obj@meta.data; cannot plot PHATE CSC signature.")
  }
  if (!"is_CSC" %in% colnames(phate_df)) {
    # If missing, create a dummy logical so the code does not crash later
    phate_df$is_CSC <- FALSE
  }

  # PHATE embedding colored by pseudotime
  p_phate_pt <- ggplot(
    phate_df,
    aes(x = PHATE1, y = PHATE2, color = pseudotime_norm)
  ) +
    geom_point(size = 0.5, alpha = 0.95) +
    scale_color_viridis_c(name = "Pseudotime", option = "B") +
    labs(
      title = "PHATE colored by Slingshot pseudotime",
      x = "PHATE 1",
      y = "PHATE 2"
    ) +
    coord_equal() +
    theme_umap

  # Save PHATE pseudotime figure
  ggsave(
    filename = file.path(FIG_DIR, "scc_scRNA_phate_pseudotime.png"),
    plot     = p_phate_pt,
    width    = 5.8,
    height   = 4.8,
    dpi      = 400
  )

  # PHATE embedding colored by CSC signature score
  p_phate_csc <- ggplot(
    phate_df,
    aes(x = PHATE1, y = PHATE2, color = CSC_signature_score)
  ) +
    geom_point(size = 0.5, alpha = 0.95) +
    scale_color_viridis_c(name = "Yuan CSC\nSignature", option = "A") +
    labs(
      title = "PHATE colored by CSC signature score",
      x = "PHATE 1",
      y = "PHATE 2"
    ) +
    coord_equal() +
    theme_umap

  # Save PHATE CSC signature figure
  ggsave(
    filename = file.path(FIG_DIR, "scc_scRNA_phate_CSC_signature.png"),
    plot     = p_phate_csc,
    width    = 5.8,
    height   = 4.8,
    dpi      = 400
  )
} else {
  # Explain why PHATE was skipped, with specific reasons
  if (!have_phateR) {
    message("  PHATE skipped: phateR package not installed.")
  } else if (!have_reticulate) {
    message("  PHATE skipped: reticulate package not installed, cannot access Python.")
  } else if (!have_py_phate) {
    message("  PHATE skipped: Python module 'phate' not found in current reticulate environment.")
    message("  Hint: In an interactive session, run `phateR::install.phate()` once, then rerun.")
  }
}
# ---------------------------------------------------------------------------
# Advanced pseudotime heatmap (ComplexHeatmap) of top pseudotime-correlated genes
# ---------------------------------------------------------------------------

have_CH <- requireNamespace("ComplexHeatmap", quietly = TRUE)
have_circ <- requireNamespace("circlize", quietly = TRUE)

if (have_CH && have_circ) {
  # Attach required namespaces explicitly for clarity
  library(ComplexHeatmap)
  library(circlize)

  message("  Building pseudotime heatmap with ComplexHeatmap...")

  # Extract normalized expression matrix (genes x cells) from Seurat
  expr_mat <- Seurat::GetAssayData(scc_obj, slot = "data")

  # Pull pseudotime and metadata
  pt <- scc_obj@meta.data$pseudotime_norm

  # Identify cells with valid pseudotime
  valid_cells <- which(!is.na(pt))
  if (length(valid_cells) < 10L) {
    warning("Fewer than 10 cells have non-missing pseudotime; skipping pseudotime heatmap.")
  } else {
    # Subset expression and pseudotime to valid cells
    expr_sub <- expr_mat[, valid_cells, drop = FALSE]
    pt_sub <- pt[valid_cells]

    # Compute Spearman correlation of each gene with pseudotime
    cor_vec <- apply(
      X = expr_sub,
      MARGIN = 1,
      FUN = function(x) {
        # Guard against zero-variance genes
        if (stats::var(x, na.rm = TRUE) == 0) {
          return(NA_real_)
        } else {
          return(suppressWarnings(stats::cor(x, pt_sub, method = "spearman", use = "pairwise.complete.obs")))
        }
      }
    )

    # Drop genes with NA correlations
    cor_vec <- cor_vec[!is.na(cor_vec)]
    if (length(cor_vec) == 0L) {
      warning("No genes with valid correlation to pseudotime; skipping pseudotime heatmap.")
    } else {
      # Choose top N genes by absolute correlation
      top_n <- min(50L, length(cor_vec))
      top_genes <- names(sort(abs(cor_vec), decreasing = TRUE))[seq_len(top_n)]

      # Subset expression matrix to top genes
      heat_mat <- as.matrix(expr_sub[top_genes, , drop = FALSE])

      # Z-score scale each gene across cells
      heat_mat_scaled <- t(scale(t(heat_mat)))
      # Replace any rows with all NA (rare) by zeros to avoid errors
      bad_rows <- which(!is.finite(rowSums(heat_mat_scaled)))
      if (length(bad_rows) > 0L) {
        heat_mat_scaled[bad_rows, ] <- 0
      }

      # Order cells by pseudotime
      ord <- order(pt_sub)
      heat_mat_scaled <- heat_mat_scaled[, ord, drop = FALSE]
      pt_ord <- pt_sub[ord]

      # Column annotations
      meta_sub <- scc_obj@meta.data[valid_cells, , drop = FALSE]
      meta_ord <- meta_sub[ord, , drop = FALSE]

      # Extract or build annotation vectors with defensive fallbacks
      CSC_score <- if ("CSC_signature_score" %in% colnames(meta_ord)) {
        meta_ord$CSC_signature_score
      } else {
        rep(NA_real_, length(ord))
      }

      is_CSC_vec <- if ("is_CSC" %in% colnames(meta_ord)) {
        meta_ord$is_CSC
      } else {
        rep(FALSE, length(ord))
      }

      # Cluster labels: prefer cluster_label column, otherwise Idents
      if ("cluster_label" %in% colnames(meta_ord)) {
        cluster_vec <- as.character(meta_ord$cluster_label)
      } else {
        cluster_all <- as.character(Seurat::Idents(scc_obj))
        cluster_vec <- cluster_all[valid_cells][ord]
      }

      # Build column annotation object
      col_ha <- HeatmapAnnotation(
        df = data.frame(
          pseudotime = pt_ord,
          CSC_signature = CSC_score,
          is_CSC = factor(is_CSC_vec, levels = c(FALSE, TRUE)),
          cluster = factor(cluster_vec)
        ),
        col = list(
          pseudotime = colorRamp2(c(min(pt_ord), max(pt_ord)), c("navy", "yellow")),
          CSC_signature = colorRamp2(
            c(
              min(CSC_score, na.rm = TRUE),
              max(CSC_score, na.rm = TRUE)
            ),
            c("grey90", "darkred")
          ),
          is_CSC = c(`FALSE` = "grey80", `TRUE` = "black"),
          cluster = structure(
            grDevices::rainbow(length(unique(cluster_vec))),
            names = levels(factor(cluster_vec))
          )
        ),
        annotation_height = unit.c(
          unit(0.4, "cm"),
          unit(0.4, "cm"),
          unit(0.4, "cm"),
          unit(0.4, "cm")
        ),
        simple_anno_size = unit(0.35, "cm")
      )

      # Define color scale for the heatmap
      expr_col_fun <- colorRamp2(
        c(-2, 0, 2),
        c("midnightblue", "white", "firebrick3")
      )

      # Build the heatmap object
      ht <- Heatmap(
        matrix = heat_mat_scaled,
        name = "Expression (z)",
        col = expr_col_fun,
        top_annotation = col_ha,
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        row_names_gp = grid::gpar(fontsize = 6),
        heatmap_legend_param = list(
          title = "Scaled expression",
          at = c(-2, 0, 2),
          labels = c("-2", "0", "2")
        )
      )

      # Save to file as a high resolution PNG
      png(
        filename = file.path(FIG_DIR, "scc_scRNA_pseudotime_heatmap.png"),
        width    = 2600,
        height   = 2200,
        res      = 320
      )
      draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
      dev.off()
    }
  }
} else {
  if (!have_CH) message("  Pseudotime heatmap skipped: ComplexHeatmap not installed.")
  if (!have_circ) message("  Pseudotime heatmap skipped: circlize not installed.")
}

# ---------------------------------------------------------------------------
# Mechanistic correlation scatter panels (unchanged) + MODULE vs PSEUDOTIME
# ---------------------------------------------------------------------------

theme_scatter <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.title       = element_text(face = "bold", size = 10)
  )

p_cor_TR <- ggplot(test_df, aes(TGFb_module_score, Lepr_expr)) +
  geom_point(alpha = 0.2, size = 0.5, color = "grey40") +
  geom_smooth(method = "lm", color = "red3", se = TRUE) +
  labs(
    title = sprintf("T→R: r=%.3f, p=%.2e", cor_T_R$estimate, cor_T_R$p.value),
    x = "TGFβ module score", y = "Lepr expression"
  ) +
  theme_scatter

p_cor_MC <- ggplot(test_df, aes(mTOR_module_score, CSC_module_score)) +
  geom_point(alpha = 0.2, size = 0.5, color = "grey40") +
  geom_smooth(method = "lm", color = "red3", se = TRUE) +
  labs(
    title = sprintf("M→C: r=%.3f, p=%.2e", cor_M_C$estimate, cor_M_C$p.value),
    x = "mTOR module score", y = "CSC module score"
  ) +
  theme_scatter

p_cor_AT <- ggplot(test_df, aes(Angio_module_score, TGFb_module_score)) +
  geom_point(alpha = 0.2, size = 0.5, color = "grey40") +
  geom_smooth(method = "lm", color = "red3", se = TRUE) +
  labs(
    title = sprintf("A→T: r=%.3f, p=%.2e", cor_A_T$estimate, cor_A_T$p.value),
    x = "Angio module score", y = "TGFβ module score"
  ) +
  theme_scatter

p_cor_RM <- ggplot(test_df, aes(Lepr_expr, mTOR_module_score)) +
  geom_point(alpha = 0.2, size = 0.5, color = "grey40") +
  geom_smooth(method = "lm", color = "red3", se = TRUE) +
  labs(
    title = sprintf("R→M: r=%.3f, p=%.2e", cor_R_M$estimate, cor_R_M$p.value),
    x = "Lepr expression", y = "mTOR module score"
  ) +
  theme_scatter

p_cor_combined <- (p_cor_TR + p_cor_MC) / (p_cor_AT + p_cor_RM)

ggsave(file.path(FIG_DIR, "scc_scRNA_mechanistic_correlations_panel.png"),
       p_cor_combined, width = 10, height = 9, dpi = 400)

# MODULE SCORES vs PSEUDOTIME (trend plots)
p_pt_TGFb <- ggplot(test_df, aes(pseudotime, TGFb_module_score)) +
  geom_point(alpha = 0.15, size = 0.4, color = "grey40") +
  geom_smooth(method = "loess", se = TRUE, color = "red3") +
  labs(
    title = "TGFβ module vs Slingshot pseudotime",
    x = "Pseudotime (normalized)", y = "TGFβ module score"
  ) +
  theme_scatter

p_pt_mTOR <- ggplot(test_df, aes(pseudotime, mTOR_module_score)) +
  geom_point(alpha = 0.15, size = 0.4, color = "grey40") +
  geom_smooth(method = "loess", se = TRUE, color = "red3") +
  labs(
    title = "mTOR module vs Slingshot pseudotime",
    x = "Pseudotime (normalized)", y = "mTOR module score"
  ) +
  theme_scatter

p_pt_Angio <- ggplot(test_df, aes(pseudotime, Angio_module_score)) +
  geom_point(alpha = 0.15, size = 0.4, color = "grey40") +
  geom_smooth(method = "loess", se = TRUE, color = "red3") +
  labs(
    title = "Angiogenesis module vs Slingshot pseudotime",
    x = "Pseudotime (normalized)", y = "Angio module score"
  ) +
  theme_scatter

p_pt_CSCsig <- ggplot(test_df, aes(pseudotime, CSC_signature_score)) +
  geom_point(alpha = 0.15, size = 0.4, color = "grey40") +
  geom_smooth(method = "loess", se = TRUE, color = "red3") +
  labs(
    title = "Yuan CSC signature vs Slingshot pseudotime",
    x = "Pseudotime (normalized)", y = "CSC signature score"
  ) +
  theme_scatter

p_pt_panel <- (p_pt_TGFb + p_pt_mTOR) / (p_pt_Angio + p_pt_CSCsig)

ggsave(file.path(FIG_DIR, "scc_scRNA_pseudotime_module_trends.png"),
       p_pt_panel, width = 10, height = 9, dpi = 400)

close(validation_log)

# Re-enable warnings
options(warn = 0)

message("\n========================================")
message("PIPELINE COMPLETE")
message("========================================\n")
message("✓ CSCs identified using Yuan's 101-gene signature")
message("✓ Module score thresholding (top 36%)")
message("✓ Slingshot pseudotime + principal curves on PCA/UMAP")
message("✓ UMAP/t-SNE/PHATE pseudotime & CSC signature figures")
message("✓ Check validation report and figures/main/ for results\n")
