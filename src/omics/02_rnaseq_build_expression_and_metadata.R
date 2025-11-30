#!/usr/bin/env Rscript

# ===========================================================
# 02_rnaseq_build_expression_and_metadata.R
#
# Pipeline step: MAIN STEP 1 (expression + canonical metadata)
#
# Purpose:
#   1) Read three GSE190411 RNA-Seq matrices:
#        - Bl6 Norm/Pap/SCC (ENSEMBLID_Symbol)
#        - PAP/SCC (gene symbols)
#        - PDV WT / LeprKO (gene symbols)
#   2) Convert each to gene_symbol x sample matrix.
#   3) Log2-transform values (log2(x + 1)).
#   4) Write:
#        - data/interim/rnaseq/bl6_expression.tsv
#        - data/interim/rnaseq/pap_scc_expression.tsv
#        - data/interim/rnaseq/pdv_expression.tsv
#   5) Build the *canonical* sample metadata used downstream:
#        - data/interim/rnaseq/sample_metadata_GSE190411.csv
#
# Notes:
#   - This is the upstream for 03_rnaseq_compute_module_scores.R.
#   - If 01_... is run, its metadata is overwritten here to keep the
#     module-scoring logic consistent.
# ===========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

#------------------------- paths ---------------------------

root_rnaseq <- "data/raw/GSE190411"
out_dir     <- "data/interim/rnaseq"

if (!dir.exists(root_rnaseq)) {
  stop("[ERROR] Directory not found: ", root_rnaseq,
       ". Run this from project root.")
}

if (!dir.exists(out_dir)) {
  message("[INFO] Creating output directory: ", out_dir)
  dir.create(out_dir, recursive = TRUE)
}

# Input files
f_bl6 <- file.path(root_rnaseq,
                   "GSE190411_Yuan2021_Bl6_Norm_Pap_SCC_RNAseq_counts.csv.gz")
f_pap <- file.path(root_rnaseq,
                   "GSE190411_Yuan2021_PAP_SCC_RNAseq_counts.txt.gz")
f_pdv <- file.path(root_rnaseq,
                   "GSE190411_Yuan2021_PDV_WT_leprKO_RNAseq_counts.txt.gz")

#--------------------- helper functions --------------------

log2p1 <- function(x) {
  log2(x + 1)
}

#----------------- 1) BL6 Norm/Pap/SCC ---------------------

message("\n[STEP] Processing Bl6 Norm/Pap/SCC expression")

if (!file.exists(f_bl6)) {
  stop("[ERROR] Missing Bl6 file: ", f_bl6)
}

bl6_raw <- read_csv(f_bl6, show_col_types = FALSE)

# Expect first column "transcript_id" with ENSEMBLID_Symbol
if (!"transcript_id" %in% colnames(bl6_raw)) {
  stop("[ERROR] Expected column 'transcript_id' in Bl6 file, found: ",
       paste(colnames(bl6_raw)[1:5], collapse = ", "))
}

bl6_sym <- bl6_raw %>%
  mutate(
    # transcript_id like "ENSMUSG00000000001_Gnai3" -> "Gnai3"
    gene_symbol = str_replace(transcript_id, "^[^_]+_", "")
  ) %>%
  select(gene_symbol, everything(), -transcript_id) %>%
  group_by(gene_symbol) %>%
  summarise(
    across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    across(where(is.numeric), log2p1)
  )

out_bl6 <- file.path(out_dir, "bl6_expression.tsv")
write_tsv(bl6_sym, out_bl6)
message("[OK] Wrote Bl6 expression to: ", out_bl6)

#----------------- 2) PAP / SCC expression -----------------

message("\n[STEP] Processing PAP/SCC expression")

if (!file.exists(f_pap)) {
  stop("[ERROR] Missing PAP/SCC file: ", f_pap)
}

pap_raw <- read_tsv(f_pap, show_col_types = FALSE)

if (!"gene" %in% colnames(pap_raw)) {
  stop("[ERROR] Expected column 'gene' in PAP/SCC file, found: ",
       paste(colnames(pap_raw)[1:5], collapse = ", "))
}

pap_sym <- pap_raw %>%
  rename(gene_symbol = gene) %>%
  group_by(gene_symbol) %>%
  summarise(
    across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    across(where(is.numeric), log2p1)
  )

out_pap <- file.path(out_dir, "pap_scc_expression.tsv")
write_tsv(pap_sym, out_pap)
message("[OK] Wrote PAP/SCC expression to: ", out_pap)

#----------------- 3) PDV WT / LeprKO ----------------------

message("\n[STEP] Processing PDV WT / LeprKO expression")

if (!file.exists(f_pdv)) {
  stop("[ERROR] Missing PDV file: ", f_pdv)
}

pdv_raw <- read_tsv(f_pdv, show_col_types = FALSE)

if (!"gene" %in% colnames(pdv_raw)) {
  stop("[ERROR] Expected column 'gene' in PDV file, found: ",
       paste(colnames(pdv_raw)[1:5], collapse = ", "))
}

pdv_sym <- pdv_raw %>%
  rename(gene_symbol = gene) %>%
  group_by(gene_symbol) %>%
  summarise(
    across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    across(where(is.numeric), log2p1)
  )

out_pdv <- file.path(out_dir, "pdv_expression.tsv")
write_tsv(pdv_sym, out_pdv)
message("[OK] Wrote PDV expression to: ", out_pdv)

#----------------- 4) Sample metadata ----------------------

message("\n[STEP] Building sample metadata from column names")

# Reload summaries just to grab column names cleanly
bl6_expr <- read_tsv(out_bl6, show_col_types = FALSE)
pap_expr <- read_tsv(out_pap, show_col_types = FALSE)
pdv_expr <- read_tsv(out_pdv, show_col_types = FALSE)

ids_bl6 <- colnames(bl6_expr)[colnames(bl6_expr) != "gene_symbol"]
ids_pap <- colnames(pap_expr)[colnames(pap_expr) != "gene_symbol"]
ids_pdv <- colnames(pdv_expr)[colnames(pdv_expr) != "gene_symbol"]

meta <- tibble(
  sample_id = c(ids_bl6, ids_pap, ids_pdv),
  dataset   = c(
    rep("Bl6",     length(ids_bl6)),
    rep("PAP_SCC", length(ids_pap)),
    rep("PDV",     length(ids_pdv))
  )
) %>%
  mutate(
    sample_id_lower = stringr::str_to_lower(sample_id),
    condition = case_when(
      # Bl6 normals vs tumors
      dataset == "Bl6" &
        stringr::str_detect(sample_id_lower, "normal") ~ "Normal",

      dataset == "Bl6" &
        stringr::str_detect(sample_id_lower, "tumor") ~ "SCC",

      # PAP/SCC dataset
      dataset == "PAP_SCC" &
        stringr::str_detect(sample_id_lower, "pap") ~ "Papilloma",

      dataset == "PAP_SCC" &
        stringr::str_detect(sample_id_lower, "scc") ~ "SCC",

      # PDV WT / KO
      dataset == "PDV" &
        stringr::str_detect(sample_id_lower, "^wt_") ~ "PDV_WT",

      dataset == "PDV" &
        stringr::str_detect(sample_id_lower, "^ko_") ~ "PDV_LeprKO",

      TRUE ~ "UNKNOWN"
    ),
    replicate = 1L
  ) %>%
  select(sample_id, dataset, condition, replicate)

out_meta <- file.path(out_dir, "sample_metadata_GSE190411.csv")
write_csv(meta, out_meta)

message("[OK] Wrote sample metadata to: ", out_meta)
message("\n[SUMMARY] Dataset x condition table:\n")
print(table(meta$dataset, meta$condition))

unknown <- meta$sample_id[meta$condition == "UNKNOWN"]
if (length(unknown) > 0) {
  message("\n[WARNING] Some samples are labeled UNKNOWN. ",
          "Open ", out_meta, " and fix 'condition' for these rows if needed:\n  ",
          paste(unknown, collapse = ", "))
} else {
  message("\n[INFO] No UNKNOWN conditions detected.")
}

message("\n[DONE] RNA-Seq expression + metadata build completed.\n")
