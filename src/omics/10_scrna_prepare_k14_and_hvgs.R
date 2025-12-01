#!/usr/bin/env Rscript

#===========================================================
# scrna_prepare_k14_and_hvgs.R
#
# Purpose:
#   - Start from raw scRNA TPM and counts for GSE207975
#   - Do simple QC (genes per cell, cells per gene)
#   - Gate Krt14+ tumour basal cells (log2(TPM+1) > 7)
#   - Save K14+ log2 TPM and counts matrices
#   - Run ERCC-based technical noise model (Yuan et al.)
#     to call highly variable genes
#   - Save a "variable" CSV with a logical 'sig' column
#
# Inputs:
#   data/raw/GSE207975/GSE207975_Yuan2021_SCC_scRNAseq_TPM.csv.gz
#   data/raw/GSE207975/GSE207975_Yuan2021_SCC_scRNAseq_counts.csv.gz
#
# Outputs:
#   data/interim/omics_qc/scc_scRNA_K14pos_expression_logged_filtered.csv
#   data/interim/omics_qc/scc_scRNA_K14pos_counts_filtered.csv
#   data/interim/omics_qc/all_tumor_K14morethan7_Counts_geneCollapsed_filtered_variable.csv
#
#===========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(DESeq2)
  library(genefilter)
  library(statmod)
})

message("[STEP] Setting paths")

tpm_path    <- "data/raw/GSE207975/GSE207975_Yuan2021_SCC_scRNAseq_TPM.csv.gz"
counts_path <- "data/raw/GSE207975/GSE207975_Yuan2021_SCC_scRNAseq_counts.csv.gz"

out_dir_interim <- "data/interim/omics_qc"
if (!dir.exists(out_dir_interim)) {
  dir.create(out_dir_interim, recursive = TRUE)
}

k14_expr_out   <- file.path(out_dir_interim, "scc_scRNA_K14pos_expression_logged_filtered.csv")
k14_counts_out <- file.path(out_dir_interim, "scc_scRNA_K14pos_counts_filtered.csv")
hvg_out        <- file.path(out_dir_interim, "all_tumor_K14morethan7_Counts_geneCollapsed_filtered_variable.csv")

# Constants
min_genes_per_cell <- 2500L
min_cells_per_gene <- 100L
krt14_id           <- "ENSMUSG00000045545.8"
krt14_threshold    <- 7.0  # log2(TPM+1) > 7 as in Yuan

#-----------------------------------------------------------
# 1. Load TPM and counts matrices
#-----------------------------------------------------------

if (!file.exists(tpm_path)) {
  stop("[ERROR] TPM file not found at: ", tpm_path)
}
if (!file.exists(counts_path)) {
  stop("[ERROR] Counts file not found at: ", counts_path)
}

message("[STEP] Reading TPM from: ", tpm_path)
tpm_raw <- readr::read_csv(tpm_path, show_col_types = FALSE)

message("[STEP] Reading counts from: ", counts_path)
counts_raw <- readr::read_csv(counts_path, show_col_types = FALSE)

# First column is gene ID, remaining columns are cells
if (ncol(tpm_raw) < 2L || ncol(counts_raw) < 2L) {
  stop("[ERROR] TPM or counts matrix has fewer than 2 columns. Check input files.")
}

gene_id_col_tpm    <- colnames(tpm_raw)[1]
gene_id_col_counts <- colnames(counts_raw)[1]

if (gene_id_col_tpm != gene_id_col_counts) {
  warning(
    "[WARN] First column names differ between TPM and counts (",
    gene_id_col_tpm, " vs ", gene_id_col_counts,
    "). Proceeding but assuming both are gene IDs."
  )
}

message("[INFO] Using '", gene_id_col_tpm, "' as gene ID column.")

gene_ids_tpm    <- tpm_raw[[gene_id_col_tpm]]
gene_ids_counts <- counts_raw[[gene_id_col_counts]]

if (!identical(as.character(gene_ids_tpm), as.character(gene_ids_counts))) {
  stop("[ERROR] Gene ID columns for TPM and counts do not match in content or order.")
}

# Build matrices with genes as rows, cells as columns
tpm_mat <- as.matrix(tpm_raw[, -1, drop = FALSE])
rownames(tpm_mat) <- as.character(gene_ids_tpm)

counts_mat <- as.matrix(counts_raw[, -1, drop = FALSE])
rownames(counts_mat) <- as.character(gene_ids_counts)

# Basic checks on cell names
if (!identical(colnames(tpm_mat), colnames(counts_mat))) {
  stop("[ERROR] Cell columns of TPM and counts do not match in content or order.")
}

message("[INFO] TPM matrix dimensions:    ", paste(dim(tpm_mat), collapse = " x "))
message("[INFO] Counts matrix dimensions: ", paste(dim(counts_mat), collapse = " x "))

#-----------------------------------------------------------
# 2. Compute log2(TPM+1) and basic QC
#-----------------------------------------------------------

message("[STEP] Computing log2(TPM + 1)")

if (any(tpm_mat < 0, na.rm = TRUE)) {
  stop("[ERROR] TPM matrix contains negative values. Check input.")
}

log_tpm <- log2(tpm_mat + 1)

# Detection threshold: log2(TPM+1) > 1
message("[STEP] Computing QC metrics (genes per cell and cells per gene)")
detected <- log_tpm > 1

num_genes_per_cell <- colSums(detected)
num_cells_per_gene <- rowSums(detected)

message("[INFO] Summary of detected genes per cell:")
print(summary(num_genes_per_cell))

message("[INFO] Summary of detected cells per gene:")
print(summary(num_cells_per_gene))

message("[STEP] Applying QC thresholds: genes_per_cell >= ",
        min_genes_per_cell,
        ", cells_per_gene >= ",
        min_cells_per_gene)

good_cells <- num_genes_per_cell >= min_genes_per_cell
good_genes <- num_cells_per_gene >= min_cells_per_gene

if (!any(good_cells)) {
  stop("[ERROR] No cells pass the genes-per-cell threshold. Threshold may be too strict.")
}
if (!any(good_genes)) {
  stop("[ERROR] No genes pass the cells-per-gene threshold. Threshold may be too strict.")
}

log_tpm_qc  <- log_tpm[good_genes, good_cells, drop = FALSE]
counts_qc   <- counts_mat[good_genes, good_cells, drop = FALSE]

message("[INFO] After QC - log_tpm_qc:  ", paste(dim(log_tpm_qc), collapse = " x "))
message("[INFO] After QC - counts_qc:   ", paste(dim(counts_qc), collapse = " x "))

#-----------------------------------------------------------
# 3. Krt14+ gating (K14+ tumour basal cells)
#-----------------------------------------------------------

message("[STEP] Applying Krt14+ gating with gene ID: ", krt14_id)

if (!(krt14_id %in% rownames(log_tpm_qc))) {
  stop(
    "[ERROR] Krt14 gene ID '", krt14_id,
    "' not found in filtered log_tpm_qc. Check gene IDs or adjust krt14_id."
  )
}

krt14_expr <- log_tpm_qc[krt14_id, , drop = TRUE]

krt14_pos_cells <- krt14_expr > krt14_threshold

message("[INFO] Number of Krt14+ cells (log2(TPM+1) > ", krt14_threshold, "): ",
        sum(krt14_pos_cells))

if (sum(krt14_pos_cells) == 0L) {
  stop("[ERROR] No Krt14+ cells found with the current threshold. Check data and threshold.")
}

log_tpm_k14 <- log_tpm_qc[, krt14_pos_cells, drop = FALSE]
counts_k14  <- counts_qc[, krt14_pos_cells,  drop = FALSE]

message("[INFO] K14+ log_tpm dimensions: ", paste(dim(log_tpm_k14), collapse = " x "))
message("[INFO] K14+ counts dimensions:  ", paste(dim(counts_k14), collapse = " x "))

# Save K14+ matrices for reuse in downstream steps
message("[STEP] Writing K14+ log_tpm and counts matrices to ", out_dir_interim)

# Expression is already log2(TPM+1)
readr::write_csv(
  as.data.frame(log_tpm_k14) %>%
    tibble::rownames_to_column(var = "gene_id"),
  k14_expr_out
)

readr::write_csv(
  as.data.frame(counts_k14) %>%
    tibble::rownames_to_column(var = "gene_id"),
  k14_counts_out
)

message("[OK] Saved K14+ log_tpm to:   ", k14_expr_out)
message("[OK] Saved K14+ counts to:    ", k14_counts_out)

#-----------------------------------------------------------
# 4. ERCC-based technical noise and HVG calling (Yuan style)
#-----------------------------------------------------------

message("[STEP] Running ERCC-based technical noise model to call HVGs")

# Round counts as in the original script
dataMouse <- round(counts_k14, digits = 0)

if (any(dataMouse < 0, na.rm = TRUE)) {
  stop("[ERROR] Negative counts detected after rounding. Check input counts.")
}

# Split mouse genes vs ERCC spike-ins based on prefix
# In Yuan's script they used substr(rownames, 1, 2) mapped to "ENSMUSG" or "ERCC"
geneTypes <- factor(
  c(EN = "ENSMUSG", ER = "ERCC")[substr(rownames(dataMouse), 1, 2)]
)

if (all(is.na(geneTypes))) {
  stop("[ERROR] geneTypes is all NA. Rowname prefixes do not match expected 'EN' or 'ER'.")
}

countsMmus <- dataMouse[geneTypes == "ENSMUSG", , drop = FALSE]
countsERCC <- dataMouse[geneTypes == "ERCC",    , drop = FALSE]

if (nrow(countsMmus) == 0L) {
  stop("[ERROR] No mouse genes detected by prefix. Check rownames.")
}
if (nrow(countsERCC) == 0L) {
  stop("[ERROR] No ERCC spike-in rows detected. Technical noise model cannot be applied.")
}

message("[INFO] Mouse gene rows: ", nrow(countsMmus))
message("[INFO] ERCC rows:      ", nrow(countsERCC))

# Size factor normalisation (DESeq)
sfMmus <- estimateSizeFactorsForMatrix(countsMmus)
sfERCC <- estimateSizeFactorsForMatrix(countsERCC)

# Normalise counts
nCountsERCC <- t(t(countsERCC) / sfERCC)
nCountsMmus <- t(t(countsMmus) / sfMmus)

# Sample moments
meansERCC <- rowMeans(nCountsERCC)
varsERCC  <- genefilter::rowVars(nCountsERCC)
cv2ERCC   <- varsERCC / meansERCC^2

meansMmus <- rowMeans(nCountsMmus)
varsMmus  <- genefilter::rowVars(nCountsMmus)
cv2Mmus   <- varsMmus / meansMmus^2

# Fit technical noise model on ERCC
minMeanForFit <- unname(
  stats::quantile(meansERCC[which(cv2ERCC > 0.2)], 0.80, na.rm = TRUE)
)

useForFit <- meansERCC >= minMeanForFit

if (!any(useForFit)) {
  stop("[ERROR] No ERCC spikes pass the minimum mean threshold for fit. Check data.")
}

fit <- statmod::glmgam.fit(
  cbind(a0 = 1, a1tilde = 1 / meansERCC[useForFit]),
  cv2ERCC[useForFit]
)

# Parameters for high variance test
minBiolDisp <- 0.25^2

xi         <- mean(1 / sfERCC)
m          <- ncol(countsMmus)
psia1theta <- mean(1 / sfERCC) +
  (coef(fit)["a1tilde"] - xi) * mean(sfERCC / sfMmus)

cv2th    <- coef(fit)["a0"] + minBiolDisp + coef(fit)["a0"] * minBiolDisp
testDenom <- (meansMmus * psia1theta + meansMmus^2 * cv2th) / (1 + cv2th / m)

# Chi-square test for each mouse gene
p <- 1 - stats::pchisq(varsMmus * (m - 1) / testDenom, df = m - 1)

padj <- stats::p.adjust(p, method = "BH")
sig  <- padj < 0.10
sig[is.na(sig)] <- FALSE

message("[INFO] Number of highly variable genes (padj < 0.10): ", sum(sig))

# Save logical vector with rownames = mouse gene IDs
sig_df <- data.frame(sig = sig)
rownames(sig_df) <- rownames(countsMmus)

readr::write_csv(
  sig_df %>%
    tibble::rownames_to_column(var = "gene_id"),
  hvg_out
)

message("[OK] HVG flags written to: ", hvg_out)
message("[DONE] scrna_prepare_k14_and_hvgs.R finished successfully.")
