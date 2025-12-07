#!/usr/bin/env Rscript

# ===========================================================
# 09_correlate_rna_atac_modules.R
#
# Purpose:
#   1) Read combined RNA+ATAC module summary:
#        - data/processed/omics/module_summary_rna_atac_combined.csv
#   2) For each module, compute correlation between
#      RNA and ATAC mean scores across matching biological states.
#   3) Write:
#        - data/processed/omics/rna_atac_module_correlations.csv
#
# Notes:
#   - ATAC has n = 1 per condition; this is *not* a statistical test,
#     it is a descriptive cross-omics concordance check.
# ===========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

in_path  <- "data/processed/omics/module_summary_rna_atac_combined.csv"
out_dir  <- "data/processed/omics"
out_path <- file.path(out_dir, "rna_atac_module_correlations.csv")

if (!file.exists(in_path)) {
  stop("[ERROR] Combined RNA+ATAC summary not found: ", in_path)
}
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

message("[STEP] Reading combined RNA+ATAC summary: ", in_path)
df <- read_csv(in_path, show_col_types = FALSE)

req_cols <- c("omics_type", "dataset", "condition", "module",
              "n_samples", "mean_score", "sd_score")
missing <- setdiff(req_cols, colnames(df))
if (length(missing) > 0L) {
  stop("[ERROR] Combined table missing columns: ",
       paste(missing, collapse = ", "))
}

#------------------------------------------
# 1) Define rough mapping of conditions
#    (RNA condition -> ATAC condition)
#    This is necessarily approximate.
#------------------------------------------
mapping <- tibble::tribble(
  ~rna_condition,   ~atac_condition,
  "Normal",         "HFSC",          # quiescent / baseline
  "Papilloma",      "Papilloma_mneg",# benign papilloma
  "SCC",            "SCC_mneg"       # aggressive SCC
)

#------------------------------------------
# 2) Prepare RNA and ATAC subsets
#------------------------------------------
rna <- df %>%
  filter(omics_type == "RNA") %>%
  select(module, condition, mean_score) %>%
  rename(rna_condition = condition,
         rna_mean      = mean_score)

atac <- df %>%
  filter(omics_type == "ATAC") %>%
  select(module, condition, mean_score) %>%
  rename(atac_condition = condition,
         atac_mean      = mean_score)

#------------------------------------------
# 3) Join RNA and ATAC by mapped conditions
#------------------------------------------
joined <- mapping %>%
  inner_join(rna,  by = c("rna_condition")) %>%
  inner_join(atac, by = c("atac_condition", "module"))

if (nrow(joined) == 0L) {
  stop("[ERROR] No overlapping module-condition pairs between RNA and ATAC.")
}

#------------------------------------------
# 4) Compute per-module correlations
#------------------------------------------
cor_df <- joined %>%
  group_by(module) %>%
  summarise(
    n_pairs           = n(),
    cor_pearson       = ifelse(n() >= 2,
                               cor(rna_mean, atac_mean, method = "pearson"),
                               NA_real_),
    cor_spearman      = ifelse(n() >= 2,
                               cor(rna_mean, atac_mean, method = "spearman"),
                               NA_real_),
    .groups = "drop"
  )

message("[STEP] Writing RNA–ATAC module correlations to: ", out_path)
write_csv(cor_df, out_path)

message("[DONE] RNA–ATAC module concordance computed.")
