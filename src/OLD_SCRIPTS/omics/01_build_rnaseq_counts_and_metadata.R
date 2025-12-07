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
#        - data/interim/rnaseq/sample_metadata_GSE190411.csv
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

  # Treat first column as gene ID
  df <- as.data.frame(df)
  colnames(df)[1] <- "gene_id"
  df$gene_id <- as.character(df$gene_id)

  # Check for duplicate gene IDs (incl. Excel-mangled ones like "1-Mar")
  dup_ids <- df$gene_id[duplicated(df$gene_id)]

  if (length(dup_ids) > 0L) {
    message(
      "[WARN] Found ", length(unique(dup_ids)),
      " duplicated gene IDs in ", basename(path), "."
    )
    message(
      "[WARN] Examples of duplicated IDs: ",
      paste(head(unique(dup_ids), 10L), collapse = ", ")
    )

    # Collapse duplicates by summing counts across rows
    df <- df %>%
      dplyr::group_by(gene_id) %>%
      dplyr::summarise(
        dplyr::across(
          where(is.numeric),
          ~ sum(.x, na.rm = TRUE)
        ),
        .groups = "drop"
      )
  }

  df
}

#---------- read and merge all counts ----------

counts_list <- lapply(count_files, read_counts)

# merge on gene_id column
all_counts <- Reduce(function(x, y) {
  dplyr::full_join(x, y, by = "gene_id")
}, counts_list)

if (!"gene_id" %in% colnames(all_counts)) {
  stop("[ERROR] gene_id column missing after merge; something is wrong upstream.")
}

if (anyDuplicated(all_counts$gene_id)) {
  stop("[ERROR] Duplicate gene_ids after merge. Inspect GSE190411 counts files.")
}

# convert gene_id to rownames
all_counts_df <- as.data.frame(all_counts)
rownames(all_counts_df) <- all_counts_df$gene_id
all_counts_df$gene_id <- NULL

# sanity: check for duplicated sample columns
dup_cols <- colnames(all_counts_df)[duplicated(colnames(all_counts_df))]
if (length(dup_cols) > 0) {
  stop(
    "[ERROR] Duplicate sample IDs in merged counts matrix: ",
    paste(unique(dup_cols), collapse = ", ")
  )
}

#---------- write unified counts ----------

counts_out <- file.path(out_dir, "GSE190411_all_counts.tsv")
message("[INFO] Writing merged counts to: ", counts_out)

all_counts_out <- tibble(gene_id = rownames(all_counts_df)) %>%
  bind_cols(as.data.frame(all_counts_df))

write_tsv(all_counts_out, counts_out)

#---------- auto-generate metadata ----------

sample_ids <- colnames(all_counts_df)

metadata <- tibble(sample_id = sample_ids) %>%
  mutate(
    sample_id_lower = stringr::str_to_lower(sample_id),
    condition = case_when(
      # 1) Normals must win over everything else
      str_detect(sample_id_lower, "normal") ~ "Normal",

      # 2) Papillomas
      str_detect(sample_id_lower, "^papmneg|^papmpos|papilloma") ~ "Papilloma",

      # 3) PDV LeprKO (KO reps or explicit LeprKO)
      str_detect(sample_id_lower, "^ko_rep") ~ "PDV_LeprKO",
      str_detect(sample_id_lower, "leprko") ~ "PDV_LeprKO",

      # 4) PDV WT (WT reps or explicit pdv)
      str_detect(sample_id_lower, "^wt_rep") ~ "PDV_WT",
      str_detect(sample_id_lower, "pdv") ~ "PDV_WT",

      # 5) Remaining tumors / SCC:
      #    - generic 'tumor'
      #    - explicit SCCmneg/SCCmpos
      #    - Ras/TGFbRII tumors (TetOHras, Tgfbr)
      str_detect(sample_id_lower, "tumor") ~ "SCC",
      str_detect(sample_id_lower, "^sccmneg|^sccmpos") ~ "SCC",
      str_detect(sample_id_lower, "tetohras") ~ "SCC",
      str_detect(sample_id_lower, "tgfbr") ~ "SCC",
      TRUE ~ "UNKNOWN"
    ),
    group = case_when(
      condition %in% c("Normal", "Papilloma", "SCC") ~ "Bl6_Pap_SCC",
      condition %in% c("PDV_WT", "PDV_LeprKO") ~ "PDV_Lepr",
      TRUE ~ "UNKNOWN"
    ),
    replicate = 1L
  ) %>%
  select(sample_id, condition, group, replicate)

meta_out <- file.path(out_dir, "sample_metadata_GSE190411.csv")
message("[INFO] Writing auto-generated metadata to: ", meta_out)
write_csv(metadata, meta_out)

#---------- sanity: enforce zero UNKNOWN -------------------

message("\n[SUMMARY] Auto-generated metadata condition counts:")
print(table(metadata$condition))

unknown_ids <- metadata$sample_id[metadata$condition == "UNKNOWN"]

if (length(unknown_ids) > 0) {
  stop(
    "[ERROR] Some samples are still UNKNOWN. ",
    "Update the case_when rules to handle these IDs explicitly:\n  ",
    paste(unknown_ids, collapse = ", ")
  )
} else {
  message("\n[INFO] No UNKNOWN conditions detected.")
}

message("[DONE] GSE190411 counts + metadata build complete.")
