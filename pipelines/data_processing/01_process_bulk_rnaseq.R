#!/usr/bin/env Rscript
# pipelines/data_processing/01_process_bulk_rnaseq.R
#
# COMPLETE BULK RNA-SEQ PIPELINE WITH RIGOROUS VALIDATION
# Runs all steps with statistical testing:
#   1. Build expression matrices
#   2. Build sample metadata
#   3. Compute module scores with z-scoring
#   4. Statistical validation (ANOVA, effect sizes, correlations)
#   5. PCA and clustering for batch effects
#   6. Cross-dataset validation
#   7. Generate publication figures
#   8. Export calibration targets with confidence intervals
#
# Output: data/processed/rnaseq/ras_csc_calibration_targets.csv
#         data/processed/rnaseq/rnaseq_validation_report.txt

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(yaml)
  library(ggplot2)
  library(ggpubr)
  library(pheatmap)
  library(RColorBrewer)
  library(broom)
})

ROOT <- "."
RAW_DIR <- file.path(ROOT, "data/raw/GSE190411")
INTERIM_DIR <- file.path(ROOT, "data/interim/rnaseq")
PROCESSED_DIR <- file.path(ROOT, "data/processed/rnaseq")
FIG_DIR <- file.path(ROOT, "figures/main")
GENESET_YAML <- file.path(ROOT, "config/gene_sets_rnaseq.yaml")

# Create directories
dir.create(INTERIM_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

message("\n========================================")
message("BULK RNA-SEQ PROCESSING PIPELINE")
message("WITH STATISTICAL VALIDATION")
message("========================================\n")

# Open validation report
validation_report <- file.path(PROCESSED_DIR, "rnaseq_validation_report.txt")
validation_log <- file(validation_report, "w")

write_validation <- function(msg) {
  message(msg)
  writeLines(msg, validation_log)
}

write_validation("==============================================")
write_validation("BULK RNA-SEQ VALIDATION REPORT")
write_validation(paste("Generated:", Sys.time()))
write_validation("==============================================\n")

#=============================================================================
# STEP 1: BUILD EXPRESSION MATRICES
#=============================================================================

message("[STEP 1/8] Building log2-transformed expression matrices")

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

message("  ✓ Expression matrices written")

write_validation("Expression matrices built:")
write_validation(sprintf("  Bl6: %d genes x %d samples", nrow(bl6_expr)-1, ncol(bl6_expr)-1))
write_validation(sprintf("  PAP_SCC: %d genes x %d samples", nrow(pap_expr)-1, ncol(pap_expr)-1))
write_validation(sprintf("  PDV: %d genes x %d samples", nrow(pdv_expr)-1, ncol(pdv_expr)-1))

#=============================================================================
# STEP 2: BUILD SAMPLE METADATA
#=============================================================================

message("\n[STEP 2/8] Building sample metadata")

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
message("  ✓ Metadata written")

write_validation("\nSample metadata:")
cond_table <- table(metadata$condition)
for (cond in names(cond_table)) {
  write_validation(sprintf("  %s: %d samples", cond, cond_table[cond]))
}

#=============================================================================
# STEP 3: COMPUTE MODULE SCORES
#=============================================================================

message("\n[STEP 3/8] Computing module scores")

# Load gene sets
gene_sets <- read_yaml(GENESET_YAML)
message("  Gene sets loaded: ", paste(names(gene_sets), collapse = ", "))

write_validation("\nModule gene sets:")
for (mod in names(gene_sets)) {
  write_validation(sprintf("  %s: %d genes", mod, length(gene_sets[[mod]])))
}

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
    
    # Report coverage
    message(sprintf("    %s: %d/%d genes present", mod, length(present), length(genes)))
    
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
# STEP 4: STATISTICAL VALIDATION
#=============================================================================

message("\n[STEP 4/8] Statistical validation of module scores")

scores_long <- all_scores %>%
  pivot_longer(
    cols = c(TGFb_bulk, mTOR_bulk, Angio_bulk, CSC_bulk),
    names_to = "module",
    values_to = "score"
  )

write_validation("\n==============================================")
write_validation("STATISTICAL TESTS: Module ~ Condition")
write_validation("==============================================\n")

# Test each module across conditions
module_stats <- list()

for (mod in c("TGFb_bulk", "mTOR_bulk", "Angio_bulk", "CSC_bulk")) {
  test_data <- scores_long %>%
    filter(module == mod) %>%
    mutate(condition = factor(condition, levels = c("Normal", "Papilloma", "SCC", "PDV_WT", "PDV_LeprKO")))
  
  # ANOVA
  anova_result <- aov(score ~ condition, data = test_data)
  anova_summary <- summary(anova_result)
  f_stat <- anova_summary[[1]]$`F value`[1]
  p_val <- anova_summary[[1]]$`Pr(>F)`[1]
  
  # Effect size (eta-squared)
  ss_condition <- anova_summary[[1]]$`Sum Sq`[1]
  ss_total <- sum(anova_summary[[1]]$`Sum Sq`)
  eta_squared <- ss_condition / ss_total
  
  write_validation(sprintf("%s:", mod))
  write_validation(sprintf("  ANOVA: F=%.2f, p=%.2e, η²=%.3f %s", 
                          f_stat, p_val, eta_squared,
                          ifelse(p_val < 0.05, "✓ SIGNIFICANT", "✗ NS")))
  
  # Post-hoc pairwise comparisons for progression
  pairwise_tests <- list()
  
  # Normal vs Papilloma
  norm_pap <- test_data %>% filter(condition %in% c("Normal", "Papilloma"))
  if (nrow(norm_pap) > 0) {
    t_test <- t.test(score ~ condition, data = norm_pap)
    cohens_d <- (mean(norm_pap$score[norm_pap$condition == "Papilloma"]) - 
                 mean(norm_pap$score[norm_pap$condition == "Normal"])) / 
                sd(norm_pap$score)
    pairwise_tests[["Normal→Papilloma"]] <- list(p = t_test$p.value, d = cohens_d)
    write_validation(sprintf("  Normal→Papilloma: p=%.2e, d=%.2f", t_test$p.value, cohens_d))
  }
  
  # Papilloma vs SCC
  pap_scc <- test_data %>% filter(condition %in% c("Papilloma", "SCC"))
  if (nrow(pap_scc) > 0) {
    t_test <- t.test(score ~ condition, data = pap_scc)
    cohens_d <- (mean(pap_scc$score[pap_scc$condition == "SCC"]) - 
                 mean(pap_scc$score[pap_scc$condition == "Papilloma"])) / 
                sd(pap_scc$score)
    pairwise_tests[["Papilloma→SCC"]] <- list(p = t_test$p.value, d = cohens_d)
    write_validation(sprintf("  Papilloma→SCC: p=%.2e, d=%.2f", t_test$p.value, cohens_d))
  }
  
  # PDV WT vs LeprKO
  pdv_comp <- test_data %>% filter(condition %in% c("PDV_WT", "PDV_LeprKO"))
  if (nrow(pdv_comp) > 0) {
    t_test <- t.test(score ~ condition, data = pdv_comp)
    cohens_d <- (mean(pdv_comp$score[pdv_comp$condition == "PDV_LeprKO"]) - 
                 mean(pdv_comp$score[pdv_comp$condition == "PDV_WT"])) / 
                sd(pdv_comp$score)
    pairwise_tests[["PDV_WT→LeprKO"]] <- list(p = t_test$p.value, d = cohens_d)
    write_validation(sprintf("  PDV_WT→LeprKO: p=%.2e, d=%.2f %s", 
                            t_test$p.value, cohens_d,
                            ifelse(abs(cohens_d) > 0.8, "✓ LARGE EFFECT", "")))
  }
  
  module_stats[[mod]] <- list(
    anova_p = p_val,
    anova_f = f_stat,
    eta_squared = eta_squared,
    pairwise = pairwise_tests
  )
  
  write_validation("")
}

#=============================================================================
# STEP 5: CROSS-DATASET VALIDATION
#=============================================================================

message("\n[STEP 5/8] Cross-dataset validation")

write_validation("==============================================")
write_validation("CROSS-DATASET CONSISTENCY")
write_validation("==============================================\n")

# Compare Bl6 SCC vs PAP_SCC SCC
bl6_scc <- scores_long %>% filter(dataset == "Bl6", condition == "SCC")
pap_scc <- scores_long %>% filter(dataset == "PAP_SCC", condition == "SCC")

if (nrow(bl6_scc) > 0 && nrow(pap_scc) > 0) {
  # Aggregate by module
  bl6_means <- bl6_scc %>% group_by(module) %>% summarise(score_bl6 = mean(score))
  pap_means <- pap_scc %>% group_by(module) %>% summarise(score_pap = mean(score))
  
  cross_comp <- inner_join(bl6_means, pap_means, by = "module")
  
  if (nrow(cross_comp) > 0) {
    cor_test <- cor.test(cross_comp$score_bl6, cross_comp$score_pap, method = "spearman")
    
    write_validation("SCC samples across datasets (Bl6 vs PAP_SCC):")
    write_validation(sprintf("  Spearman correlation: r=%.3f, p=%.3f %s",
                            cor_test$estimate, cor_test$p.value,
                            ifelse(cor_test$p.value < 0.05 & abs(cor_test$estimate) > 0.7, 
                                   "✓ CONSISTENT", "⚠ WEAK")))
    
    for (i in 1:nrow(cross_comp)) {
      write_validation(sprintf("    %s: Bl6=%.2f, PAP=%.2f, diff=%.2f",
                              cross_comp$module[i],
                              cross_comp$score_bl6[i],
                              cross_comp$score_pap[i],
                              abs(cross_comp$score_bl6[i] - cross_comp$score_pap[i])))
    }
  }
}

write_validation("")

#=============================================================================
# STEP 6: PCA FOR BATCH EFFECTS
#=============================================================================

message("\n[STEP 6/8] PCA and batch effect assessment")

# Prepare matrix for PCA (samples x modules)
pca_matrix <- all_scores %>%
  select(sample_id, dataset, condition, TGFb_bulk, mTOR_bulk, Angio_bulk, CSC_bulk) %>%
  column_to_rownames("sample_id") %>%
  select(TGFb_bulk, mTOR_bulk, Angio_bulk, CSC_bulk)

pca_meta <- all_scores %>%
  select(sample_id, dataset, condition)

# Run PCA
pca_result <- prcomp(pca_matrix, scale. = TRUE, center = TRUE)
pca_df <- as.data.frame(pca_result$x) %>%
  rownames_to_column("sample_id") %>%
  left_join(pca_meta, by = "sample_id")

# Variance explained
var_explained <- summary(pca_result)$importance[2, 1:2] * 100

write_validation("PCA Analysis:")
write_validation(sprintf("  PC1 variance: %.1f%%", var_explained[1]))
write_validation(sprintf("  PC2 variance: %.1f%%", var_explained[2]))

# Test if samples cluster by dataset (batch effect)
pc1_aov <- aov(PC1 ~ dataset, data = pca_df)
pc1_p <- summary(pc1_aov)[[1]]$`Pr(>F)`[1]

write_validation(sprintf("  PC1 ~ dataset: p=%.3f %s", pc1_p,
                        ifelse(pc1_p < 0.05, "⚠ BATCH EFFECT", "✓ NO BATCH EFFECT")))

# PCA plot
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, shape = dataset)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_brewer(palette = "Set1", name = "Condition") +
  scale_shape_manual(values = c(16, 17, 15), name = "Dataset") +
  labs(
    title = "PCA of bulk RNA-seq module scores",
    x = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
    y = sprintf("PC2 (%.1f%% variance)", var_explained[2])
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(FIG_DIR, "rnaseq_pca_module_scores.png"),
       p_pca, width = 7, height = 5, dpi = 300)

#=============================================================================
# STEP 7: PUBLICATION FIGURES
#=============================================================================

message("\n[STEP 7/8] Generating publication figures")

# Prepare data for plotting
scores_long_plot <- scores_long %>%
  mutate(
    condition = factor(condition, levels = c("Normal", "Papilloma", "SCC", "PDV_WT", "PDV_LeprKO")),
    module = factor(module, levels = c("TGFb_bulk", "Angio_bulk", "CSC_bulk", "mTOR_bulk"))
  )

# Figure 1: Violin plots with statistics
p_violins <- ggplot(scores_long_plot, aes(x = condition, y = score, fill = condition)) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6) +
  facet_wrap(~ module, scales = "free_y", ncol = 2) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Module scores across conditions",
    x = "Condition",
    y = "Module score (z-score within dataset)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  stat_compare_means(method = "anova", label.y.npc = 0.95, size = 3)

ggsave(file.path(FIG_DIR, "rnaseq_module_violins_with_stats.png"),
       p_violins, width = 9, height = 8, dpi = 400)

# Figure 2: Heatmap of condition means
condition_summary <- scores_long %>%
  group_by(condition, module) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

heatmap_matrix <- condition_summary %>%
  pivot_wider(names_from = module, values_from = mean_score) %>%
  column_to_rownames("condition") %>%
  as.matrix()

# Row order: progression
row_order <- c("Normal", "Papilloma", "SCC", "PDV_WT", "PDV_LeprKO")
heatmap_matrix <- heatmap_matrix[intersect(row_order, rownames(heatmap_matrix)), ]

png(file.path(FIG_DIR, "rnaseq_module_heatmap.png"), width = 2400, height = 1600, res = 300)
pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize = 12,
  fontsize_number = 10,
  main = "Module scores by condition (mean z-scores)",
  border_color = "grey60",
  cellwidth = 60,
  cellheight = 40
)
dev.off()

# Figure 3: Progression trajectory (Normal → Pap → SCC)
progression_data <- scores_long %>%
  filter(condition %in% c("Normal", "Papilloma", "SCC")) %>%
  mutate(
    condition_numeric = case_when(
      condition == "Normal" ~ 1,
      condition == "Papilloma" ~ 2,
      condition == "SCC" ~ 3
    ),
    module_clean = gsub("_bulk", "", module)
  )

p_trajectory <- ggplot(progression_data, aes(x = condition_numeric, y = score, color = module_clean)) +
  geom_point(alpha = 0.4, size = 2) +
  geom_smooth(method = "loess", se = TRUE, size = 1.2) +
  scale_x_continuous(breaks = c(1, 2, 3), labels = c("Normal", "Papilloma", "SCC")) +
  scale_color_brewer(palette = "Set1", name = "Module") +
  labs(
    title = "Module score progression (benign → malignant)",
    x = "Disease stage",
    y = "Module score (z-score)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(FIG_DIR, "rnaseq_progression_trajectory.png"),
       p_trajectory, width = 7, height = 5, dpi = 300)

#=============================================================================
# STEP 8: EXPORT CALIBRATION TARGETS WITH CONFIDENCE INTERVALS
#=============================================================================

message("\n[STEP 8/8] Exporting calibration targets with confidence intervals")

# Compute mean ± SE for each condition
targets_with_ci <- scores_long %>%
  group_by(dataset, condition, module) %>%
  summarise(
    n_samples = n(),
    mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    se_score = sd_score / sqrt(n_samples),
    ci_lower = mean_score - 1.96 * se_score,
    ci_upper = mean_score + 1.96 * se_score,
    .groups = "drop"
  )

# Pivot to wide format for model
targets <- targets_with_ci %>%
  select(dataset, condition, module, mean_score) %>%
  pivot_wider(names_from = module, values_from = mean_score) %>%
  mutate(
    C_target = CSC_bulk,
    A_target = Angio_bulk,
    T_target = TGFb_bulk,
    M_target = mTOR_bulk
  ) %>%
  select(dataset, condition, C_target, A_target, T_target, M_target, everything())

# Also save version with CIs for uncertainty quantification
targets_detailed <- targets_with_ci %>%
  select(dataset, condition, module, mean_score, se_score, ci_lower, ci_upper) %>%
  pivot_wider(
    names_from = module,
    values_from = c(mean_score, se_score, ci_lower, ci_upper)
  )

write_csv(targets, file.path(PROCESSED_DIR, "ras_csc_calibration_targets.csv"))
write_csv(targets_detailed, file.path(PROCESSED_DIR, "ras_csc_calibration_targets_with_ci.csv"))

write_validation("\n==============================================")
write_validation("CALIBRATION TARGETS EXPORTED")
write_validation("==============================================\n")

for (i in 1:nrow(targets)) {
  write_validation(sprintf("%s - %s:", targets$dataset[i], targets$condition[i]))
  write_validation(sprintf("  C=%.3f, A=%.3f, T=%.3f, M=%.3f",
                          targets$C_target[i], targets$A_target[i],
                          targets$T_target[i], targets$M_target[i]))
}

# Close validation report
close(validation_log)

message("\n========================================")
message("PIPELINE COMPLETE")
message("========================================")
message("\nCalibration targets written to:")
message("  ", file.path(PROCESSED_DIR, "ras_csc_calibration_targets.csv"))
message("  ", file.path(PROCESSED_DIR, "ras_csc_calibration_targets_with_ci.csv"))
message("\nValidation report:")
message("  ", validation_report)
message("\nFigures:")
message("  ", file.path(FIG_DIR, "rnaseq_*.png"))
message("\nConditions available for model fitting:")
print(table(targets$condition))
message("\n✓ READ VALIDATION REPORT FOR STATISTICAL TEST RESULTS\n")