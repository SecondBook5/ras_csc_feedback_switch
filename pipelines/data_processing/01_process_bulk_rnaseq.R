#!/usr/bin/env Rscript
# pipelines/data_processing/01_process_bulk_rnaseq.R
#
# COMPLETE BULK RNA-SEQ PIPELINE
# Runs all 5 steps in sequence:
#   1. Build counts + metadata
#   2. Build expression matrices
#   3. Compute module scores
#   4. QC + summarize by condition
#   5. Export calibration targets
#
# Output: data/processed/rnaseq/ras_csc_calibration_targets.csv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(yaml)
  library(ggplot2)
})

ROOT <- "."
RAW_DIR <- file.path(ROOT, "data/raw/GSE190411")
INTERIM_DIR <- file.path(ROOT, "data/interim/rnaseq")
PROCESSED_DIR <- file.path(ROOT, "data/processed/rnaseq")
GENESET_YAML <- file.path(ROOT, "config/gene_sets_rnaseq.yaml")

# Create directories
dir.create(INTERIM_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)

message("\n========================================")
message("BULK RNA-SEQ PROCESSING PIPELINE")
message("========================================\n")

#=============================================================================
# STEP 1: BUILD EXPRESSION MATRICES
#=============================================================================

message("[STEP 1/5] Building log2-transformed expression matrices")

# Helper function
log2p1 <- function(x) log2(x + 1)

# 1A: Bl6 Normal/Papilloma/SCC
message("  Processing Bl6 dataset...")
bl6_raw <- read_csv(
  file.path(RAW_DIR, "GSE190411_Yuan2021_Bl6_Norm_Pap_SCC_RNAseq_counts.csv.gz"),
  show_col_types = FALSE
)
bl6_expr <- bl6_raw %>%
  mutate(gene_symbol = str_replace(transcript_id, "^[^_]+_", "")) %>%
  select(gene_symbol, everything(), -transcript_id) %>%
  group_by(gene_symbol) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(across(where(is.numeric), log2p1))
write_tsv(bl6_expr, file.path(INTERIM_DIR, "bl6_expression.tsv"))

# 1B: PAP/SCC
message("  Processing PAP/SCC dataset...")
pap_raw <- read_tsv(
  file.path(RAW_DIR, "GSE190411_Yuan2021_PAP_SCC_RNAseq_counts.txt.gz"),
  show_col_types = FALSE
)
pap_expr <- pap_raw %>%
  rename(gene_symbol = gene) %>%
  group_by(gene_symbol) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(across(where(is.numeric), log2p1))
write_tsv(pap_expr, file.path(INTERIM_DIR, "pap_scc_expression.tsv"))

# 1C: PDV WT / LeprKO
message("  Processing PDV dataset...")
pdv_raw <- read_tsv(
  file.path(RAW_DIR, "GSE190411_Yuan2021_PDV_WT_leprKO_RNAseq_counts.txt.gz"),
  show_col_types = FALSE
)
pdv_expr <- pdv_raw %>%
  rename(gene_symbol = gene) %>%
  group_by(gene_symbol) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(across(where(is.numeric), log2p1))
write_tsv(pdv_expr, file.path(INTERIM_DIR, "pdv_expression.tsv"))

message("  ✓ Expression matrices written to: ", INTERIM_DIR)

#=============================================================================
# STEP 2: BUILD SAMPLE METADATA
#=============================================================================

message("\n[STEP 2/5] Building sample metadata")

ids_bl6 <- setdiff(colnames(bl6_expr), "gene_symbol")
ids_pap <- setdiff(colnames(pap_expr), "gene_symbol")
ids_pdv <- setdiff(colnames(pdv_expr), "gene_symbol")

metadata <- tibble(
  sample_id = c(ids_bl6, ids_pap, ids_pdv),
  dataset = c(rep("Bl6", length(ids_bl6)),
              rep("PAP_SCC", length(ids_pap)),
              rep("PDV", length(ids_pdv)))
) %>%
  mutate(
    sample_id_lower = str_to_lower(sample_id),
    condition = case_when(
      str_detect(sample_id_lower, "norm") ~ "Normal",
      str_detect(sample_id_lower, "pap") ~ "Papilloma",
      str_detect(sample_id_lower, "leprko|^ko_rep") ~ "PDV_LeprKO",
      str_detect(sample_id_lower, "pdv|^wt_rep") ~ "PDV_WT",
      str_detect(sample_id_lower, "tumor|scc|tgfbr|tetohras") ~ "SCC",
      TRUE ~ "UNKNOWN"
    )
  ) %>%
  select(sample_id, dataset, condition)

# Enforce no unknowns
unknowns <- metadata %>% filter(condition == "UNKNOWN")
if (nrow(unknowns) > 0) {
  stop("[ERROR] Unknown conditions found: ", paste(unknowns$sample_id, collapse = ", "))
}

write_csv(metadata, file.path(INTERIM_DIR, "sample_metadata_GSE190411.csv"))
message("  ✓ Metadata written. Conditions found:")
print(table(metadata$condition))

#=============================================================================
# STEP 3: COMPUTE MODULE SCORES
#=============================================================================

message("\n[STEP 3/5] Computing module scores")

# Load gene sets
gene_sets <- read_yaml(GENESET_YAML)
message("  Gene sets loaded: ", paste(names(gene_sets), collapse = ", "))

# Function to compute z-scored module means per dataset
compute_modules <- function(expr_df, dataset_name, gene_sets, metadata) {
  message("  Computing modules for: ", dataset_name)
  
  # Long format
  expr_long <- expr_df %>%
    pivot_longer(-gene_symbol, names_to = "sample_id", values_to = "expr")
  
  # Z-score within dataset
  expr_long <- expr_long %>%
    group_by(gene_symbol) %>%
    mutate(expr_z = (expr - mean(expr, na.rm = TRUE)) / sd(expr, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(expr_z = if_else(is.na(expr_z), 0, expr_z))
  
  # Compute each module
  module_list <- lapply(names(gene_sets), function(mod) {
    genes <- gene_sets[[mod]]
    present <- intersect(genes, unique(expr_long$gene_symbol))
    
    if (length(present) == 0) {
      warning("  [WARN] Module ", mod, " has no genes in ", dataset_name)
      return(NULL)
    }
    
    expr_long %>%
      filter(gene_symbol %in% present) %>%
      group_by(sample_id) %>%
      summarise(score = mean(expr_z, na.rm = TRUE), .groups = "drop") %>%
      mutate(module = mod)
  })
  
  scores <- bind_rows(module_list) %>%
    left_join(metadata %>% filter(dataset == dataset_name), by = "sample_id") %>%
    pivot_wider(names_from = module, values_from = score)
  
  scores
}

# Compute for all datasets
scores_bl6 <- compute_modules(bl6_expr, "Bl6", gene_sets, metadata)
scores_pap <- compute_modules(pap_expr, "PAP_SCC", gene_sets, metadata)
scores_pdv <- compute_modules(pdv_expr, "PDV", gene_sets, metadata)

all_scores <- bind_rows(scores_bl6, scores_pap, scores_pdv)
write_csv(all_scores, file.path(PROCESSED_DIR, "module_scores_by_sample.csv"))
message("  ✓ Module scores written")

#=============================================================================
# STEP 4: SUMMARIZE BY CONDITION
#=============================================================================

message("\n[STEP 4/5] Summarizing to condition-level means")

scores_long <- all_scores %>%
  pivot_longer(
    cols = c(TGFb_bulk, mTOR_bulk, Angio_bulk, CSC_bulk),
    names_to = "module",
    values_to = "score"
  )

condition_summary <- scores_long %>%
  group_by(dataset, condition, module) %>%
  summarise(
    n_samples = n(),
    mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(condition_summary, file.path(PROCESSED_DIR, "module_summary_by_condition.csv"))
message("  ✓ Condition summary written")

# QC plot
scores_long <- scores_long %>%
  mutate(
    condition = factor(condition, levels = c("Normal", "Papilloma", "SCC", "PDV_WT", "PDV_LeprKO")),
    module = factor(module, levels = c("TGFb_bulk", "Angio_bulk", "CSC_bulk", "mTOR_bulk"))
  )

p <- ggplot(scores_long, aes(x = condition, y = score, fill = condition)) +
  geom_boxplot(alpha = 0.7) +
  facet_grid(module ~ dataset, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "Module scores by condition", x = "Condition", y = "Z-score")

ggsave(file.path(PROCESSED_DIR, "qc_module_scores.png"), p, width = 10, height = 7, dpi = 300)

#=============================================================================
# STEP 5: EXPORT CALIBRATION TARGETS
#=============================================================================

message("\n[STEP 5/5] Exporting model calibration targets")

targets <- condition_summary %>%
  select(dataset, condition, module, mean_score) %>%
  pivot_wider(names_from = module, values_from = mean_score) %>%
  mutate(
    C_target = CSC_bulk,
    A_target = Angio_bulk,
    T_target = TGFb_bulk,
    M_target = mTOR_bulk
  ) %>%
  select(dataset, condition, C_target, A_target, T_target, M_target, everything())

write_csv(targets, file.path(PROCESSED_DIR, "ras_csc_calibration_targets.csv"))

message("\n========================================")
message("PIPELINE COMPLETE")
message("========================================")
message("\nCalibration targets written to:")
message("  ", file.path(PROCESSED_DIR, "ras_csc_calibration_targets.csv"))
message("\nConditions available for model fitting:")
print(table(targets$condition))
message("\n")
