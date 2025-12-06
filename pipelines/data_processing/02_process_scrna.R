#!/usr/bin/env Rscript
# pipelines/data_processing/02_process_scrna.R
#
# COMPLETE scRNA-SEQ PIPELINE FOR SCC K14+ CELLS WITH RIGOROUS VALIDATION
# Consolidates scripts 10-13 + adds hypothesis testing validation:
#   1. Load raw TPM/counts from GSE207975
#   2. QC filtering + K14+ gating
#   3. ERCC-based HVG calling (Yuan method)
#   4. Build Seurat object, clustering, UMAP
#   5. Batch effect assessment
#   6. Identify CSC cluster (Lepr-high) with dual marker validation
#   7. Map gene symbols → Ensembl IDs
#   8. Compute module scores per cell
#   9. Module coherence validation
#  10. Mechanistic correlation tests (CRITICAL FOR HYPOTHESIS)
#  11. Generate publication-quality figures
#
# Outputs:
#   - data/processed/omics_summaries/scc_scRNA_seurat.rds
#   - data/processed/omics_summaries/scc_scRNA_module_scores_per_cell.csv
#   - data/processed/omics_summaries/scc_scRNA_CSC_summary.json
#   - data/processed/omics_summaries/scc_scRNA_mechanistic_correlations.csv
#   - data/processed/omics_summaries/scc_scRNA_validation_report.txt
#   - figures/main/scc_scRNA_*.png (multiple UMAP and correlation plots)

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
})

message("\n========================================")
message("scRNA-SEQ PROCESSING PIPELINE")
message("WITH RIGOROUS VALIDATION")
message("========================================\n")

# Paths
RAW_DIR <- "data/raw/GSE207975"
INTERIM_DIR <- "data/interim/omics_qc"
PROCESSED_DIR <- "data/processed/omics_summaries"
FIG_DIR <- "figures/main"
GENESET_YAML <- "config/gene_sets_rnaseq.yaml"

dir.create(INTERIM_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Constants
MIN_GENES_PER_CELL <- 2500L
MIN_CELLS_PER_GENE <- 100L
KRT14_ID <- "ENSMUSG00000045545.8"
KRT14_THRESHOLD <- 7.0
LEPR_ID <- "ENSMUSG00000059316.2"

# Canonical stem cell markers (independent validation)
CSC_CANONICAL_MARKERS <- c(
  "Itga6",   # CD49f - THE canonical CSC marker
  "Cd34",    # CD34 - hair follicle stem cell
  "Cd200",   # CD200 - Yuan identified
  "Sox9",    # SOX9 - stem cell TF
  "Lgr6",    # LGR6 - stem cell marker
  "Krt15"    # Keratin 15 - basal marker
)

# Yuan et al. reported values for validation
YUAN_CSC_FRACTION <- 0.358  # Fig 1B, C2 cluster
YUAN_LEPR_FRACTION <- 0.748 # Fig 1D, Lepr+ cells

# Open validation report file
validation_report <- file.path(PROCESSED_DIR, "scc_scRNA_validation_report.txt")
validation_log <- file(validation_report, "w")

write_validation <- function(msg) {
  message(msg)
  writeLines(msg, validation_log)
}

write_validation("==============================================")
write_validation("scRNA-SEQ VALIDATION REPORT")
write_validation(paste("Generated:", Sys.time()))
write_validation("==============================================\n")

#=============================================================================
# STEP 1: LOAD RAW DATA
#=============================================================================

message("[STEP 1/11] Loading raw scRNA-seq data")

tpm_path <- file.path(RAW_DIR, "GSE207975_Yuan2021_SCC_scRNAseq_TPM.csv.gz")
counts_path <- file.path(RAW_DIR, "GSE207975_Yuan2021_SCC_scRNAseq_counts.csv.gz")

tpm_raw <- read_csv(tpm_path, show_col_types = FALSE)
counts_raw <- read_csv(counts_path, show_col_types = FALSE)

gene_ids <- tpm_raw[[1]]
tpm_mat <- as.matrix(tpm_raw[, -1])
rownames(tpm_mat) <- gene_ids

counts_mat <- as.matrix(counts_raw[, -1])
rownames(counts_mat) <- gene_ids

message("  Raw dimensions: ", paste(dim(tpm_mat), collapse = " x "))
write_validation(paste("Raw data dimensions:", paste(dim(tpm_mat), collapse = " x ")))

#=============================================================================
# STEP 2: QC FILTERING + K14+ GATING
#=============================================================================

message("\n[STEP 2/11] QC filtering and K14+ gating")

# Log transform
log_tpm <- log2(tpm_mat + 1)

# QC metrics
detected <- log_tpm > 1
genes_per_cell <- colSums(detected)
cells_per_gene <- rowSums(detected)

message("  Genes per cell: ", paste(range(genes_per_cell), collapse = " - "))
message("  Cells per gene: ", paste(range(cells_per_gene), collapse = " - "))

# Filter
good_cells <- genes_per_cell >= MIN_GENES_PER_CELL
good_genes <- cells_per_gene >= MIN_CELLS_PER_GENE

log_tpm_qc <- log_tpm[good_genes, good_cells]
counts_qc <- counts_mat[good_genes, good_cells]

message("  After QC: ", paste(dim(log_tpm_qc), collapse = " x "))
write_validation(paste("\nQC filtering:"))
write_validation(paste("  Cells passing QC:", sum(good_cells), "/", length(good_cells)))
write_validation(paste("  Genes passing QC:", sum(good_genes), "/", length(good_genes)))

# K14+ gating
if (!KRT14_ID %in% rownames(log_tpm_qc)) {
  stop("[ERROR] Krt14 gene ID not found: ", KRT14_ID)
}

krt14_expr <- log_tpm_qc[KRT14_ID, ]
k14_pos <- krt14_expr > KRT14_THRESHOLD

message("  K14+ cells: ", sum(k14_pos), " / ", length(k14_pos))
write_validation(paste("\nK14+ gating:"))
write_validation(paste("  K14+ cells:", sum(k14_pos), "/", length(k14_pos)))
write_validation(paste("  K14+ fraction:", sprintf("%.1f%%", 100*mean(k14_pos))))

log_tpm_k14 <- log_tpm_qc[, k14_pos]
counts_k14 <- counts_qc[, k14_pos]

#=============================================================================
# STEP 3: ERCC-BASED HVG CALLING
#=============================================================================

message("\n[STEP 3/11] ERCC-based highly variable gene calling")

dataMouse <- round(counts_k14, 0)

# Split mouse genes vs ERCC
geneTypes <- factor(c(EN = "ENSMUSG", ER = "ERCC")[substr(rownames(dataMouse), 1, 2)])

countsMmus <- dataMouse[geneTypes == "ENSMUSG", ]
countsERCC <- dataMouse[geneTypes == "ERCC", ]

message("  Mouse genes: ", nrow(countsMmus))
message("  ERCC spikes: ", nrow(countsERCC))

# Size factor normalization
sfMmus <- estimateSizeFactorsForMatrix(countsMmus)
sfERCC <- estimateSizeFactorsForMatrix(countsERCC)

nCountsERCC <- t(t(countsERCC) / sfERCC)
nCountsMmus <- t(t(countsMmus) / sfMmus)

# Moments
meansERCC <- rowMeans(nCountsERCC)
varsERCC <- genefilter::rowVars(nCountsERCC)
cv2ERCC <- varsERCC / meansERCC^2

meansMmus <- rowMeans(nCountsMmus)
varsMmus <- genefilter::rowVars(nCountsMmus)
cv2Mmus <- varsMmus / meansMmus^2

# Fit technical noise on ERCC
minMeanForFit <- quantile(meansERCC[cv2ERCC > 0.2], 0.80, na.rm = TRUE)
useForFit <- meansERCC >= minMeanForFit

fit <- statmod::glmgam.fit(
  cbind(a0 = 1, a1tilde = 1 / meansERCC[useForFit]),
  cv2ERCC[useForFit]
)

# Chi-square test for HVGs
minBiolDisp <- 0.25^2
xi <- mean(1 / sfERCC)
m <- ncol(countsMmus)
psia1theta <- xi + (coef(fit)["a1tilde"] - xi) * mean(sfERCC / sfMmus)

cv2th <- coef(fit)["a0"] + minBiolDisp + coef(fit)["a0"] * minBiolDisp
testDenom <- (meansMmus * psia1theta + meansMmus^2 * cv2th) / (1 + cv2th / m)

p <- 1 - pchisq(varsMmus * (m - 1) / testDenom, df = m - 1)
padj <- p.adjust(p, method = "BH")
hvg_sig <- padj < 0.10
hvg_sig[is.na(hvg_sig)] <- FALSE

hvg_ids <- rownames(countsMmus)[hvg_sig]
message("  HVGs detected: ", length(hvg_ids))
write_validation(paste("\nHVG calling (ERCC-based):"))
write_validation(paste("  HVGs detected:", length(hvg_ids)))

#=============================================================================
# STEP 4: BUILD SEURAT OBJECT
#=============================================================================

message("\n[STEP 4/11] Building Seurat object with clustering")

scc_obj <- CreateSeuratObject(
  counts = counts_k14,
  project = "SCC_K14",
  min.cells = 0,
  min.features = 0
)

scc_obj$sample_id <- sub("_.*", "", colnames(scc_obj))

# Use Yuan HVGs
hvg_in_obj <- intersect(hvg_ids, rownames(scc_obj))
message("  Using ", length(hvg_in_obj), " HVGs for dimensionality reduction")

VariableFeatures(scc_obj) <- hvg_in_obj

# Standard workflow
scc_obj <- NormalizeData(scc_obj, verbose = FALSE)
scc_obj <- ScaleData(scc_obj, features = VariableFeatures(scc_obj), verbose = FALSE)
scc_obj <- RunPCA(scc_obj, features = VariableFeatures(scc_obj), npcs = 30, verbose = FALSE)
scc_obj <- FindNeighbors(scc_obj, dims = 1:20, verbose = FALSE)
scc_obj <- FindClusters(scc_obj, resolution = 0.4, verbose = FALSE)
scc_obj <- RunUMAP(scc_obj, dims = 1:20, verbose = FALSE)

# Add UMAP to metadata
umap_emb <- Embeddings(scc_obj, "umap")
scc_obj$UMAP_1 <- umap_emb[, 1]
scc_obj$UMAP_2 <- umap_emb[, 2]

n_clusters <- length(unique(scc_obj$seurat_clusters))
message("  Clusters found: ", n_clusters)
write_validation(paste("\nClustering:"))
write_validation(paste("  Number of clusters:", n_clusters))
write_validation(paste("  Using all ", ncol(scc_obj), " K14+ cells (doublet detection skipped)"))

#=============================================================================
# STEP 5: BATCH EFFECT ASSESSMENT
#=============================================================================

message("\n[STEP 5/11] Assessing batch effects")

# Check if samples cluster separately (potential batch effect)
sample_table <- table(scc_obj$sample_id, scc_obj$seurat_clusters)
message("  Sample x Cluster contingency table:")
print(sample_table)

# Visualize batch effect
p_batch <- DimPlot(scc_obj, group.by = "sample_id", reduction = "umap") +
  labs(title = "UMAP colored by sample ID (batch check)") +
  theme_minimal()

ggsave(file.path(FIG_DIR, "scc_scRNA_umap_batch_check.png"),
       p_batch, width = 6, height = 4.5, dpi = 300)

# Calculate entropy to assess mixing
# Low entropy = samples segregate (batch effect present)
sample_entropy <- apply(sample_table, 2, function(x) {
  p <- x / sum(x)
  p <- p[p > 0]
  -sum(p * log(p))
})

mean_entropy <- mean(sample_entropy)
max_entropy <- log(length(unique(scc_obj$sample_id)))
mixing_score <- mean_entropy / max_entropy

write_validation(paste("\nBatch effect assessment:"))
write_validation(sprintf("  Sample mixing score: %.3f (0=segregated, 1=perfect mixing)", mixing_score))

if (mixing_score < 0.5) {
  write_validation("  ⚠ WARNING: Potential batch effect detected (mixing score < 0.5)")
  write_validation("  Consider running Harmony integration if this affects conclusions")
} else {
  write_validation("  ✓ PASS: Samples are well-mixed across clusters")
}

#=============================================================================
# STEP 6: IDENTIFY CSC CLUSTER WITH DUAL VALIDATION
#=============================================================================

message("\n[STEP 6/11] Identifying and validating CSC cluster")

norm_data <- GetAssayData(scc_obj, slot = "data")
scc_obj$Lepr_expr <- as.numeric(norm_data[LEPR_ID, ])
scc_obj$Krt14_expr <- as.numeric(norm_data[KRT14_ID, ])

# Find Lepr-high cluster
cluster_summary <- scc_obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarize(
    n_cells = n(),
    Lepr_mean = mean(Lepr_expr, na.rm = TRUE),
    Lepr_median = median(Lepr_expr, na.rm = TRUE),
    Krt14_mean = mean(Krt14_expr, na.rm = TRUE),
    .groups = "drop"
  )

csc_cluster <- cluster_summary %>%
  dplyr::arrange(dplyr::desc(Lepr_mean)) %>%
  dplyr::slice(1) %>%
  dplyr::pull(seurat_clusters)


scc_obj$is_CSC <- scc_obj$seurat_clusters == csc_cluster

csc_fraction <- mean(scc_obj$is_CSC)
message("  CSC cluster: ", csc_cluster)
message("  CSC fraction: ", sprintf("%.3f", csc_fraction))

write_validation(paste("\nCSC cluster identification:"))
write_validation(paste("  CSC cluster ID:", csc_cluster))
write_validation(sprintf("  CSC fraction (ours): %.3f (%.1f%%)", csc_fraction, 100*csc_fraction))
write_validation(sprintf("  CSC fraction (Yuan): %.3f (%.1f%%)", YUAN_CSC_FRACTION, 100*YUAN_CSC_FRACTION))
write_validation(sprintf("  Difference: %.1f%%", 100*abs(csc_fraction - YUAN_CSC_FRACTION)))

if (abs(csc_fraction - YUAN_CSC_FRACTION) > 0.10) {
  write_validation("  ⚠ WARNING: CSC fraction differs by >10% from Yuan et al.")
} else {
  write_validation("  ✓ PASS: CSC fraction within 10% of published value")
}

# DUAL VALIDATION: Test both pathway genes and canonical markers
message("\n  Validating CSC cluster with dual marker sets...")

# Load CSC pathway genes from YAML
gene_sets_prelim <- read_yaml(GENESET_YAML)
CSC_PATHWAY_GENES <- gene_sets_prelim$CSC_bulk

# biomaRt setup
mart_prelim <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

# Get versioned IDs
seurat_genes <- rownames(scc_obj)
bare_to_vers <- tibble(
  bare_id = sub("\\..*$", "", seurat_genes),
  versioned_id = seurat_genes
)

#---------------------------------------------------------------------------
# VALIDATION A: Model pathway genes (CSC_bulk from YAML)
#---------------------------------------------------------------------------

message("  [A] Testing model pathway genes (CSC_bulk)...")

pathway_map <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters = "external_gene_name",
  values = CSC_PATHWAY_GENES,
  mart = mart_prelim
)

pathway_marker_ids <- pathway_map %>%
  inner_join(bare_to_vers, by = c("ensembl_gene_id" = "bare_id")) %>%
  pull(versioned_id)

if (length(pathway_marker_ids) > 0) {
  pathway_test <- FindMarkers(
    scc_obj,
    ident.1 = csc_cluster,
    features = pathway_marker_ids,
    logfc.threshold = 0,
    min.pct = 0,
    verbose = FALSE
  )
  
  pathway_test <- pathway_test %>%
    rownames_to_column("gene_id") %>%
    arrange(p_val_adj)
  
  # Count significant enrichments
  pathway_enriched <- sum(pathway_test$avg_log2FC > 0 & pathway_test$p_val_adj < 0.05)
  pathway_total <- nrow(pathway_test)
  
  write_validation("\n  [A] Model pathway genes (CSC_bulk):")
  write_validation(sprintf("      Genes tested: %d", pathway_total))
  write_validation(sprintf("      Enriched in CSC cluster: %d / %d (%.1f%%)",
                          pathway_enriched, pathway_total,
                          100 * pathway_enriched / pathway_total))
  
  write_validation("      Top 5 pathway genes:")
  for (i in 1:min(5, nrow(pathway_test))) {
    gene <- pathway_test$gene_id[i]
    fc <- pathway_test$avg_log2FC[i]
    pval <- pathway_test$p_val_adj[i]
    status <- ifelse(fc > 0 & pval < 0.05, "✓", "✗")
    write_validation(sprintf("        %s %s: FC=%.2f, p=%.2e", status, gene, fc, pval))
  }
  
  if (pathway_enriched / pathway_total > 0.5) {
    write_validation("      ✓ PASS: Majority of pathway genes enriched")
  } else {
    write_validation("      ⚠ WEAK: Less than half of pathway genes enriched")
  }
}

#---------------------------------------------------------------------------
# VALIDATION B: Canonical stem cell markers (independent)
#---------------------------------------------------------------------------

message("  [B] Testing canonical stem cell markers...")

canonical_map <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters = "external_gene_name",
  values = CSC_CANONICAL_MARKERS,
  mart = mart_prelim
)

canonical_marker_ids <- canonical_map %>%
  inner_join(bare_to_vers, by = c("ensembl_gene_id" = "bare_id")) %>%
  pull(versioned_id)

if (length(canonical_marker_ids) > 0) {
  canonical_test <- FindMarkers(
    scc_obj,
    ident.1 = csc_cluster,
    features = canonical_marker_ids,
    logfc.threshold = 0,
    min.pct = 0,
    verbose = FALSE
  )
  
  canonical_test <- canonical_test %>%
    rownames_to_column("gene_id") %>%
    arrange(p_val_adj)
  
  # Count significant enrichments
  canonical_enriched <- sum(canonical_test$avg_log2FC > 0 & canonical_test$p_val_adj < 0.05)
  canonical_total <- nrow(canonical_test)
  
  write_validation("\n  [B] Canonical stem cell markers (independent validation):")
  write_validation(sprintf("      Genes tested: %d", canonical_total))
  write_validation(sprintf("      Enriched in CSC cluster: %d / %d (%.1f%%)",
                          canonical_enriched, canonical_total,
                          100 * canonical_enriched / canonical_total))
  
  write_validation("      Individual markers:")
  for (i in 1:nrow(canonical_test)) {
    gene <- canonical_test$gene_id[i]
    fc <- canonical_test$avg_log2FC[i]
    pval <- canonical_test$p_val_adj[i]
    status <- ifelse(fc > 0 & pval < 0.05, "✓", "✗")
    write_validation(sprintf("        %s %s: FC=%.2f, p=%.2e", status, gene, fc, pval))
  }
  
  if (canonical_enriched / canonical_total >= 0.5) {
    write_validation("      ✓ PASS: CSC cluster expresses canonical stem markers")
  } else {
    write_validation("      ⚠ WEAK: Few canonical markers enriched")
  }
}

write_validation("\n  Overall CSC validation:")
if (exists("pathway_enriched") && exists("canonical_enriched")) {
  if ((pathway_enriched / pathway_total > 0.5) && (canonical_enriched / canonical_total >= 0.5)) {
    write_validation("  ✓ STRONG: Cluster enriched for BOTH pathway genes AND canonical markers")
  } else if ((pathway_enriched / pathway_total > 0.3) || (canonical_enriched / canonical_total >= 0.3)) {
    write_validation("  ⚠ MODERATE: Partial enrichment detected")
  } else {
    write_validation("  ✗ WEAK: Limited enrichment for CSC markers")
  }
}

#=============================================================================
# STEP 7: MAP GENE SYMBOLS → ENSEMBL IDs
#=============================================================================

message("\n[STEP 7/11] Mapping gene symbols to Ensembl IDs")

gene_sets <- read_yaml(GENESET_YAML)

# Flatten to symbols
gene_symbols <- unlist(gene_sets) %>% unique()

# biomaRt mapping
mart <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
map_df <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters = "external_gene_name",
  values = gene_symbols,
  mart = mart
)

message("  Mapped ", nrow(map_df), " symbols to Ensembl")

# Build feature lists for each module
seurat_genes <- rownames(scc_obj)
bare_to_versioned <- tibble(
  gene_id_seurat = seurat_genes,
  bare_id = sub("\\..*$", "", seurat_genes)
)

module_features <- list()
for (mod in names(gene_sets)) {
  symbols <- gene_sets[[mod]]
  
  # Symbol → Ensembl → Versioned Ensembl
  ensembl_ids <- map_df %>%
    filter(external_gene_name %in% symbols) %>%
    pull(ensembl_gene_id) %>%
    unique()
  
  versioned <- bare_to_versioned %>%
    filter(bare_id %in% ensembl_ids) %>%
    pull(gene_id_seurat)
  
  clean_name <- gsub("_bulk$", "", mod)
  module_features[[clean_name]] <- versioned
  
  message("  Module ", clean_name, ": ", length(versioned), " genes")
}

#=============================================================================
# STEP 8: COMPUTE MODULE SCORES
#=============================================================================

message("\n[STEP 8/11] Computing module scores per cell")

# Filter modules with >=5 genes
valid_modules <- module_features[sapply(module_features, length) >= 5]
module_names <- names(valid_modules)

scc_obj <- AddModuleScore(
  scc_obj,
  features = valid_modules,
  name = module_names,
  nbin = 24,
  ctrl = 100,
  seed = 1
)

# Extract scores
meta_df <- scc_obj@meta.data %>%
  rownames_to_column("cell_id")

# Build scores table
scores_df <- tibble(cell_id = meta_df$cell_id)

for (i in seq_along(module_names)) {
  mod_name <- module_names[i]
  col_pattern <- paste0("^", mod_name, i, "$")
  col_match <- grep(col_pattern, colnames(meta_df), value = TRUE)
  
  if (length(col_match) > 0) {
    scores_df[[paste0(mod_name, "_module_score")]] <- meta_df[[col_match[1]]]
  }
}

#=============================================================================
# STEP 9: MODULE COHERENCE VALIDATION
#=============================================================================

message("\n[STEP 9/11] Validating module gene coherence")

write_validation("\nModule coherence (gene co-expression):")

for (mod_name in module_names) {
  mod_genes <- valid_modules[[mod_name]]
  
  if (length(mod_genes) < 3) {
    write_validation(sprintf("  %s: Too few genes (%d)", mod_name, length(mod_genes)))
    next
  }
  
  # Get expression matrix for module genes
  mod_expr <- GetAssayData(scc_obj, slot = "data")[mod_genes, ]
  
  # Compute pairwise correlations
  cor_mat <- cor(t(as.matrix(mod_expr)), method = "spearman")
  
  # Mean pairwise correlation (excluding diagonal)
  mean_cor <- mean(cor_mat[lower.tri(cor_mat)])
  
  status <- ifelse(mean_cor > 0.3, "✓ COHERENT", 
                   ifelse(mean_cor > 0.15, "⚠ WEAK", "✗ POOR"))
  
  write_validation(sprintf("  %s: mean r = %.3f %s", mod_name, mean_cor, status))
}

#=============================================================================
# STEP 10: MECHANISTIC CORRELATION TESTS (CRITICAL)
#=============================================================================

message("\n[STEP 10/11] Testing mechanistic predictions (HYPOTHESIS VALIDATION)")

# Join module scores with metadata
test_df <- scores_df %>%
  dplyr::left_join(
    meta_df %>%
      dplyr::select(cell_id, Lepr_expr, is_CSC, seurat_clusters),
    by = "cell_id"
  )


# Test T→R: TGFβ induces Lepr
cor_T_R <- cor.test(
  test_df$TGFb_module_score,
  test_df$Lepr_expr,
  method = "spearman"
)

# Test M→C: mTOR drives CSC state
cor_M_C <- cor.test(
  test_df$mTOR_module_score,
  test_df$CSC_module_score,
  method = "spearman"
)

# Test A→T: Angiogenesis leads to TGFβ
cor_A_T <- cor.test(
  test_df$Angio_module_score,
  test_df$TGFb_module_score,
  method = "spearman"
)

# Test R→M: Lepr leads to mTOR
cor_R_M <- cor.test(
  test_df$Lepr_expr,
  test_df$mTOR_module_score,
  method = "spearman"
)

write_validation("\n==============================================")
write_validation("MECHANISTIC CORRELATION TESTS (HYPOTHESIS)")
write_validation("==============================================")

write_validation(sprintf(
  "\nT→R (TGFβ → Lepr): r = %.3f, p = %.2e %s",
  cor_T_R$estimate,
  cor_T_R$p.value,
  ifelse(cor_T_R$p.value < 0.05 & cor_T_R$estimate > 0.3, "✓ SUPPORTED", "✗ NOT SUPPORTED")
))

write_validation(sprintf(
  "M→C (mTOR → CSC): r = %.3f, p = %.2e %s",
  cor_M_C$estimate,
  cor_M_C$p.value,
  ifelse(cor_M_C$p.value < 0.05 & cor_M_C$estimate > 0.3, "✓ SUPPORTED", "✗ NOT SUPPORTED")
))

write_validation(sprintf(
  "A→T (Angio → TGFβ): r = %.3f, p = %.2e %s",
  cor_A_T$estimate,
  cor_A_T$p.value,
  ifelse(cor_A_T$p.value < 0.05 & cor_A_T$estimate > 0.3, "✓ SUPPORTED", "✗ NOT SUPPORTED")
))

write_validation(sprintf(
  "R→M (Lepr → mTOR): r = %.3f, p = %.2e %s",
  cor_R_M$estimate,
  cor_R_M$p.value,
  ifelse(cor_R_M$p.value < 0.05 & cor_R_M$estimate > 0.3, "✓ SUPPORTED", "✗ NOT SUPPORTED")
))

# Save correlation results
cor_results <- tibble(
  arrow = c("T→R", "M→C", "A→T", "R→M"),
  mechanism = c(
    "TGFβ → Lepr",
    "mTOR → CSC",
    "Angiogenesis → TGFβ",
    "Lepr → mTOR"
  ),
  correlation = c(
    cor_T_R$estimate,
    cor_M_C$estimate,
    cor_A_T$estimate,
    cor_R_M$estimate
  ),
  p_value = c(
    cor_T_R$p.value,
    cor_M_C$p.value,
    cor_A_T$p.value,
    cor_R_M$p.value
  ),
  significant = p_value < 0.05,
  strong = abs(correlation) > 0.3,
  hypothesis_supported = significant & strong
)

write_csv(
  cor_results,
  file.path(PROCESSED_DIR, "scc_scRNA_mechanistic_correlations.csv")
)

n_supported <- sum(cor_results$hypothesis_supported)
write_validation(sprintf("\n✓ SUMMARY: %d / 4 mechanistic arrows supported", n_supported))

if (n_supported >= 3) {
  write_validation("✓ PASS: Model structure strongly supported by scRNA data")
} else if (n_supported >= 2) {
  write_validation("⚠ PARTIAL: Model structure partially supported, consider refinement")
} else {
  write_validation("✗ FAIL: Model structure not supported by scRNA data")
}

#=============================================================================
# SAVE PROCESSED DATA
#=============================================================================

# Save scores
write_csv(
  scores_df,
  file.path(PROCESSED_DIR, "scc_scRNA_module_scores_per_cell.csv")
)

# Save metadata with CSC labels
meta_export <- meta_df %>%
  dplyr::select(cell_id, sample_id, seurat_clusters, is_CSC,
         UMAP_1, UMAP_2, Lepr_expr, Krt14_expr)

write_csv(
  meta_export,
  file.path(PROCESSED_DIR, "scc_scRNA_K14pos_metadata_with_CSC_labels.csv")
)

# Save summary JSON
summary_list <- list(
  n_cells_total = ncol(scc_obj),
  n_CSC = sum(scc_obj$is_CSC),
  CSC_fraction = csc_fraction,
  CSC_fraction_yuan = YUAN_CSC_FRACTION,
  CSC_cluster_id = as.character(csc_cluster),
  cluster_summary = cluster_summary,
  mechanistic_correlations = cor_results,
  batch_mixing_score = mixing_score
)

writeLines(
  toJSON(summary_list, pretty = TRUE, auto_unbox = TRUE),
  file.path(PROCESSED_DIR, "scc_scRNA_CSC_summary.json")
)

# Save Seurat object
saveRDS(scc_obj, file.path(PROCESSED_DIR, "scc_scRNA_seurat.rds"))

message("  ✓ Module scores written")
message("  ✓ CSC summary written")
message("  ✓ Mechanistic correlations written")
message("  ✓ Seurat object saved")

#=============================================================================
# STEP 11: GENERATE PUBLICATION FIGURES
#=============================================================================

message("\n[STEP 11/11] Generating publication-quality figures")

theme_umap <- theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    axis.ticks = element_line(color = "grey70"),
    axis.line = element_line(color = "grey70"),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0)
  )

umap_df <- meta_df %>%
  dplyr::select(UMAP_1, UMAP_2, seurat_clusters, is_CSC, Lepr_expr) %>%
  mutate(seurat_clusters = factor(seurat_clusters))

# Figure 1: UMAP by cluster
p1 <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = seurat_clusters)) +
  geom_point(size = 0.5, alpha = 0.9) +
  scale_color_hue(name = "Cluster") +
  labs(title = "SCC K14+ cells - Louvain clusters",
       x = "UMAP 1", y = "UMAP 2") +
  coord_equal() +
  theme_umap

ggsave(file.path(FIG_DIR, "scc_scRNA_umap_louvain_clusters.png"),
       p1, width = 5.8, height = 4.8, dpi = 400)

# Figure 2: UMAP by Lepr expression
p2 <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = Lepr_expr)) +
  geom_point(size = 0.55, alpha = 0.95) +
  scale_color_viridis_c(name = "Lepr\n(log-norm)", option = "C") +
  labs(title = "Lepr expression in SCC K14+ cells",
       x = "UMAP 1", y = "UMAP 2") +
  coord_equal() +
  theme_umap

ggsave(file.path(FIG_DIR, "scc_scRNA_umap_Lepr_expr.png"),
       p2, width = 5.8, height = 4.8, dpi = 400)

# Figure 3: UMAP CSC vs non-CSC
umap_csc <- umap_df %>% filter(is_CSC)
umap_non <- umap_df %>% filter(!is_CSC)

p3 <- ggplot() +
  geom_point(data = umap_non, aes(UMAP_1, UMAP_2),
             color = "grey85", size = 0.4, alpha = 0.7) +
  geom_point(data = umap_csc, aes(UMAP_1, UMAP_2),
             color = "red3", size = 0.6, alpha = 0.95) +
  stat_density_2d(data = umap_csc, aes(UMAP_1, UMAP_2, fill = after_stat(level)),
                  geom = "polygon", alpha = 0.25, color = NA) +
  scale_fill_viridis_c(option = "B", guide = "none") +
  labs(title = "CSC island in SCC K14+ state space",
       x = "UMAP 1", y = "UMAP 2") +
  coord_equal() +
  theme_umap

ggsave(file.path(FIG_DIR, "scc_scRNA_umap_CSC_vs_nonCSC.png"),
       p3, width = 5.8, height = 4.8, dpi = 400)

# Figure 4: Mechanistic correlation scatter plots
theme_scatter <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.title = element_text(face = "bold", size = 10)
  )

p_cor_TR <- ggplot(test_df, aes(TGFb_module_score, Lepr_expr)) +
  geom_point(alpha = 0.2, size = 0.5, color = "grey40") +
  geom_smooth(method = "lm", color = "red3", se = TRUE) +
  labs(
    title = sprintf("T→R: r=%.3f, p=%.2e", cor_T_R$estimate, cor_T_R$p.value),
    x = "TGFβ module score",
    y = "Lepr expression"
  ) +
  theme_scatter

p_cor_MC <- ggplot(test_df, aes(mTOR_module_score, CSC_module_score)) +
  geom_point(alpha = 0.2, size = 0.5, color = "grey40") +
  geom_smooth(method = "lm", color = "red3", se = TRUE) +
  labs(
    title = sprintf("M→C: r=%.3f, p=%.2e", cor_M_C$estimate, cor_M_C$p.value),
    x = "mTOR module score",
    y = "CSC module score"
  ) +
  theme_scatter

p_cor_AT <- ggplot(test_df, aes(Angio_module_score, TGFb_module_score)) +
  geom_point(alpha = 0.2, size = 0.5, color = "grey40") +
  geom_smooth(method = "lm", color = "red3", se = TRUE) +
  labs(
    title = sprintf("A→T: r=%.3f, p=%.2e", cor_A_T$estimate, cor_A_T$p.value),
    x = "Angio module score",
    y = "TGFβ module score"
  ) +
  theme_scatter

p_cor_RM <- ggplot(test_df, aes(Lepr_expr, mTOR_module_score)) +
  geom_point(alpha = 0.2, size = 0.5, color = "grey40") +
  geom_smooth(method = "lm", color = "red3", se = TRUE) +
  labs(
    title = sprintf("R→M: r=%.3f, p=%.2e", cor_R_M$estimate, cor_R_M$p.value),
    x = "Lepr expression",
    y = "mTOR module score"
  ) +
  theme_scatter

# Combined correlation panel
p_cor_combined <- (p_cor_TR + p_cor_MC) / (p_cor_AT + p_cor_RM)

ggsave(
  file.path(FIG_DIR, "scc_scRNA_mechanistic_correlations_panel.png"),
  p_cor_combined, width = 10, height = 9, dpi = 400
)

# Save individual correlation plots too
ggsave(file.path(FIG_DIR, "scc_scRNA_correlation_TGFb_Lepr.png"),
       p_cor_TR, width = 5, height = 4, dpi = 300)
ggsave(file.path(FIG_DIR, "scc_scRNA_correlation_mTOR_CSC.png"),
       p_cor_MC, width = 5, height = 4, dpi = 300)
ggsave(file.path(FIG_DIR, "scc_scRNA_correlation_Angio_TGFb.png"),
       p_cor_AT, width = 5, height = 4, dpi = 300)
ggsave(file.path(FIG_DIR, "scc_scRNA_correlation_Lepr_mTOR.png"),
       p_cor_RM, width = 5, height = 4, dpi = 300)

message("  ✓ All figures saved to ", FIG_DIR)

# Close validation report
close(validation_log)

message("\n========================================")
message("scRNA PIPELINE COMPLETE")
message("========================================\n")
message("Outputs:")
message("  - ", file.path(PROCESSED_DIR, "scc_scRNA_seurat.rds"))
message("  - ", file.path(PROCESSED_DIR, "scc_scRNA_module_scores_per_cell.csv"))
message("  - ", file.path(PROCESSED_DIR, "scc_scRNA_CSC_summary.json"))
message("  - ", file.path(PROCESSED_DIR, "scc_scRNA_mechanistic_correlations.csv"))
message("  - ", file.path(PROCESSED_DIR, "scc_scRNA_validation_report.txt"))
message("  - ", file.path(FIG_DIR, "scc_scRNA_*.png"))
message("\n✓ READ VALIDATION REPORT FOR HYPOTHESIS TEST RESULTS\n")