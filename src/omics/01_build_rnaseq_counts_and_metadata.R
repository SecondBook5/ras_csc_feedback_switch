#!/usr/bin/env Rscript

#-----------------------------------------------------------
# 01_build_rnaseq_counts_and_metadata.R
#
# Pipeline step: OPTIONAL (counts-only)
#
# Goal:
#   1) Read all GSE190411 count matrices
#   2) Merge into a single genes x samples matrix
#   3) Auto-generate a *provisional* metadata table
#   4) Write:
#        - data/interim/rnaseq/GSE190411_all_counts.tsv
#        - data/interim/rnaseq/sample_metadata_GSE190411.csv  (will be
#          overwritten by step 02; use this only for quick QC)
#
# Notes:
#   - This step is NOT required for the Ras–CSC model calibration.
#   - It exists to give you a unified counts matrix for DESeq2 or other
#     count-based analyses.
#-----------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
})

#---------- paths ----------

root <- "data/raw/GSE190411"
out_dir <- "data/interim/rnaseq"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

count_files <- c(
  file.path(root, "GSE190411_Yuan2021_Bl6_Norm_Pap_SCC_RNAseq_counts.csv.gz"),
  file.path(root, "GSE190411_Yuan2021_PAP_SCC_RNAseq_counts.txt.gz"),
  file.path(root, "GSE190411_Yuan2021_PDV_WT_leprKO_RNAseq_counts.txt.gz")
)

#---------- helper: read counts, handle csv vs tsv ----------

read_counts <- function(path) {
  message("[INFO] Reading counts from: ", path)

  if (!file.exists(path)) {
    stop("[ERROR] File not found: ", path)
  }

  # Choose delimiter based on extension
  df <- if (grepl("\\.csv", path, ignore.case = TRUE)) {
    readr::read_csv(path, show_col_types = FALSE)
  } else {
    readr::read_tsv(path, show_col_types = FALSE)
  }

  if (ncol(df) < 2) {
    stop("[ERROR] Counts file appears to have <2 columns: ", path)
  }

  gene_col <- colnames(df)[1]

  # Convert to data.frame with gene IDs as rownames
  df <- as.data.frame(df)
  rownames(df) <- df[[gene_col]]
  df[[gene_col]] <- NULL

  return(df)
}

#---------- read and merge all counts ----------

counts_list <- lapply(count_files, read_counts)

# merge on rownames (gene IDs)
all_counts <- Reduce(function(x, y) {
  # bring rownames into a column for join
  x_tbl <- tibble(gene_id = rownames(x)) %>%
    bind_cols(as.data.frame(x))
  y_tbl <- tibble(gene_id = rownames(y)) %>%
    bind_cols(as.data.frame(y))

  full_join(x_tbl, y_tbl, by = "gene_id")
}, counts_list)

if (anyDuplicated(all_counts$gene_id)) {
  stop("[ERROR] Duplicate gene_ids after merge. Inspect GSE190411 counts files.")
}

# convert gene_id back to rownames
rownames(all_counts) <- all_counts$gene_id
all_counts$gene_id <- NULL

# sanity: check for duplicated sample columns
dup_cols <- colnames(all_counts)[duplicated(colnames(all_counts))]
if (length(dup_cols) > 0) {
  stop(
    "[ERROR] Duplicate sample IDs in merged counts matrix: ",
    paste(unique(dup_cols), collapse = ", ")
  )
}

# write unified counts
counts_out <- file.path(out_dir, "GSE190411_all_counts.tsv")
message("[INFO] Writing merged counts to: ", counts_out)

all_counts_out <- tibble(gene_id = rownames(all_counts)) %>%
  bind_cols(as.data.frame(all_counts))

write_tsv(all_counts_out, counts_out)

#---------- auto-generate metadata ----------

sample_ids <- colnames(all_counts)

metadata <- tibble(sample_id = sample_ids) %>%
  mutate(
    condition = case_when(
      str_detect(sample_id, regex("Norm", ignore_case = TRUE)) ~ "Normal",
      str_detect(sample_id, regex("Pap", ignore_case = TRUE)) ~ "Papilloma",
      str_detect(sample_id, regex("SCC", ignore_case = TRUE)) ~ "SCC",
      str_detect(sample_id, regex("LeprKO|KO", ignore_case = TRUE)) ~ "PDV_LeprKO",
      str_detect(sample_id, regex("PDV", ignore_case = TRUE)) ~ "PDV_WT",
      TRUE ~ "UNKNOWN"
    ),
    group = case_when(
      condition %in% c("Normal", "Papilloma", "SCC") ~ "Bl6_Pap_SCC",
      condition %in% c("PDV_WT", "PDV_LeprKO") ~ "PDV_Lepr",
      TRUE ~ "UNKNOWN"
    ),
    # crude replicate index: will be refined if needed
    replicate = 1L
  )

meta_out <- file.path(out_dir, "sample_metadata_GSE190411.csv")
message("[INFO] Writing auto-generated metadata to: ", meta_out)
write_csv(metadata, meta_out)

# Quick summary to help you spot issues immediately
message("\n[SUMMARY] Auto-generated metadata condition counts:")
print(table(metadata$condition))

unknown_ids <- metadata$sample_id[metadata$condition == "UNKNOWN"]
if (length(unknown_ids) > 0) {
  message("\n[WARNING] Some samples were labeled UNKNOWN. ",
          "You should open ", meta_out, " and fix these rows manually:\n  ",
          paste(unknown_ids, collapse = ", "))
} else {
  message("\n[INFO] No UNKNOWN conditions detected.")
}

message("[DONE] GSE190411 counts + metadata build complete.")
