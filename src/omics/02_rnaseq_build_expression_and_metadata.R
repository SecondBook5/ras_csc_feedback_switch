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
#   - Overwrites metadata from 01_build_rnaseq_counts_and_metadata.R
#   - Condition rules must match 01 and yield ZERO UNKNOWNs.
# ===========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

#------------------------- paths ---------------------------

root_rnaseq <- "data/raw/GSE190411"
out_dir <- "data/interim/rnaseq"

if (!dir.exists(root_rnaseq)) {
  stop(
    "[ERROR] Directory not found: ", root_rnaseq,
    ". Run this from project root."
  )
}

if (!dir.exists(out_dir)) {
  message("[INFO] Creating output directory: ", out_dir)
  dir.create(out_dir, recursive = TRUE)
}

f_bl6 <- file.path(
  root_rnaseq,
  "GSE190411_Yuan2021_Bl6_Norm_Pap_SCC_RNAseq_counts.csv.gz"
)
f_pap <- file.path(
  root_rnaseq,
  "GSE190411_Yuan2021_PAP_SCC_RNAseq_counts.txt.gz"
)
f_pdv <- file.path(
  root_rnaseq,
  "GSE190411_Yuan2021_PDV_WT_leprKO_RNAseq_counts.txt.gz"
)

#--------------------- helper functions --------------------

log2p1 <- function(x) {
  log2(x + 1)
}

#----------------- 1) BL6 Norm/Pap/SCC ---------------------

message("\n[STEP] Processing Bl6 Norm/Pap/SCC expression")

if (!file.exists(f_bl6)) {
  stop("[ERROR] Missing Bl6 file: ", f_bl6)
}

bl6_raw <- readr::read_csv(f_bl6, show_col_types = FALSE)

if (!"transcript_id" %in% colnames(bl6_raw)) {
  stop(
    "[ERROR] Expected column 'transcript_id' in Bl6 file, found: ",
    paste(colnames(bl6_raw)[1:min(5, ncol(bl6_raw))], collapse = ", ")
  )
}

bl6_sym <- bl6_raw %>%
  mutate(
    gene_symbol = stringr::str_replace(transcript_id, "^[^_]+_", "")
  ) %>%
  select(gene_symbol, dplyr::everything(), -transcript_id) %>%
  group_by(gene_symbol) %>%
  summarise(
    dplyr::across(
      where(is.numeric),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  mutate(
    dplyr::across(where(is.numeric), log2p1)
  )

out_bl6 <- file.path(out_dir, "bl6_expression.tsv")
readr::write_tsv(bl6_sym, out_bl6)
message("[OK] Wrote Bl6 expression to: ", out_bl6)

#----------------- 2) PAP / SCC expression -----------------

message("\n[STEP] Processing PAP/SCC expression")

if (!file.exists(f_pap)) {
  stop("[ERROR] Missing PAP/SCC file: ", f_pap)
}

pap_raw <- readr::read_tsv(f_pap, show_col_types = FALSE)

if (!"gene" %in% colnames(pap_raw)) {
  stop(
    "[ERROR] Expected column 'gene' in PAP/SCC file, found: ",
    paste(colnames(pap_raw)[1:min(5, ncol(pap_raw))], collapse = ", ")
  )
}

pap_sym <- pap_raw %>%
  dplyr::rename(gene_symbol = gene) %>%
  group_by(gene_symbol) %>%
  summarise(
    dplyr::across(
      where(is.numeric),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  mutate(
    dplyr::across(where(is.numeric), log2p1)
  )

out_pap <- file.path(out_dir, "pap_scc_expression.tsv")
readr::write_tsv(pap_sym, out_pap)
message("[OK] Wrote PAP/SCC expression to: ", out_pap)

#----------------- 3) PDV WT / LeprKO ----------------------

message("\n[STEP] Processing PDV WT / LeprKO expression")

if (!file.exists(f_pdv)) {
  stop("[ERROR] Missing PDV file: ", f_pdv)
}

pdv_raw <- readr::read_tsv(f_pdv, show_col_types = FALSE)

if (!"gene" %in% colnames(pdv_raw)) {
  stop(
    "[ERROR] Expected column 'gene' in PDV file, found: ",
    paste(colnames(pdv_raw)[1:min(5, ncol(pdv_raw))], collapse = ", ")
  )
}

pdv_sym <- pdv_raw %>%
  dplyr::rename(gene_symbol = gene) %>%
  group_by(gene_symbol) %>%
  summarise(
    dplyr::across(
      where(is.numeric),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  mutate(
    dplyr::across(where(is.numeric), log2p1)
  )

out_pdv <- file.path(out_dir, "pdv_expression.tsv")
readr::write_tsv(pdv_sym, out_pdv)
message("[OK] Wrote PDV expression to: ", out_pdv)

#----------------- 4) Sample metadata ----------------------

message("\n[STEP] Building sample metadata from column names")

bl6_expr <- readr::read_tsv(out_bl6, show_col_types = FALSE)
pap_expr <- readr::read_tsv(out_pap, show_col_types = FALSE)
pdv_expr <- readr::read_tsv(out_pdv, show_col_types = FALSE)

ids_bl6 <- colnames(bl6_expr)[colnames(bl6_expr) != "gene_symbol"]
ids_pap <- colnames(pap_expr)[colnames(pap_expr) != "gene_symbol"]
ids_pdv <- colnames(pdv_expr)[colnames(pdv_expr) != "gene_symbol"]

meta <- tibble::tibble(
  sample_id = c(ids_bl6, ids_pap, ids_pdv),
  dataset = c(
    rep("Bl6", length(ids_bl6)),
    rep("PAP_SCC", length(ids_pap)),
    rep("PDV", length(ids_pdv))
  )
) %>%
  mutate(
    sample_id_lower = stringr::str_to_lower(sample_id),
    condition = dplyr::case_when(
      # 1) Normal first so Ras+ normals remain 'Normal'
      stringr::str_detect(sample_id_lower, "norm") ~ "Normal",

      # 2) PAP / papilloma
      stringr::str_detect(sample_id_lower, "pap") ~ "Papilloma",

      # 3) PDV LeprKO (explicit or KO_rep)
      stringr::str_detect(sample_id_lower, "leprko") ~ "PDV_LeprKO",
      stringr::str_detect(sample_id_lower, "^ko_rep") ~ "PDV_LeprKO",

      # 4) PDV WT (explicit or WT_rep)
      stringr::str_detect(sample_id_lower, "pdv") ~ "PDV_WT",
      stringr::str_detect(sample_id_lower, "^wt_rep") ~ "PDV_WT",

      # 5) Everything else that looks like tumor/SCC/Ras/TGFbR
      stringr::str_detect(sample_id_lower, "tumor") ~ "SCC",
      stringr::str_detect(sample_id_lower, "scc") ~ "SCC",
      stringr::str_detect(sample_id_lower, "tgfbr") ~ "SCC",
      stringr::str_detect(sample_id_lower, "tetohras") ~ "SCC",

      # 6) Anything unexpected => force failure so we notice
      TRUE ~ "UNKNOWN"
    ),
    replicate = 1L
  ) %>%
  dplyr::select(sample_id, dataset, condition, replicate)

out_meta <- file.path(out_dir, "sample_metadata_GSE190411.csv")
readr::write_csv(meta, out_meta)

message("[OK] Wrote sample metadata to: ", out_meta)
message("\n[SUMMARY] Dataset x condition table:\n")
print(table(meta$dataset, meta$condition))

unknown <- meta$sample_id[meta$condition == "UNKNOWN"]

if (length(unknown) > 0L) {
  stop(
    "\n[ERROR] Some samples are still labeled UNKNOWN. ",
    "Update the case_when rules to handle these IDs explicitly:\n  ",
    paste(unknown, collapse = ", ")
  )
} else {
  message("\n[INFO] No UNKNOWN conditions detected.")
}

message("\n[DONE] RNA-Seq expression + metadata build completed.\n")
