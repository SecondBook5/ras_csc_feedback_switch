#!/usr/bin/env Rscript
# pipelines/data_processing/01_process_bulk_rnaseq.R
#
# COMPLETE BULK RNA-SEQ PIPELINE WITH RIGOROUS VALIDATION
# Runs all steps with statistical testing:
#   1. Build expression matrices
#   2. Build sample metadata (canonical, matching old 02 script)
#   3. Compute module scores with z-scoring
#   4. Statistical validation (ANOVA, effect sizes, correlations)
#   5. PCA and clustering for batch effects
#   6. Cross-dataset validation
#   7. Generate publication figures
#   7B. Advanced hypothesis-driven figures for Ras–CSC loop
#   8. Export calibration targets with confidence intervals
#
# Output:
#   data/processed/rnaseq/ras_csc_calibration_targets.csv
#   data/processed/rnaseq/ras_csc_calibration_targets_with_ci.csv
#   data/processed/rnaseq/rnaseq_validation_report.txt

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(yaml)
  library(ggplot2)
  library(ggpubr)
  library(RColorBrewer)
  library(ggrepel)
  # Advanced heatmap tools
  library(ComplexHeatmap)
  library(circlize)
})

# Yuan CSC 101 marker panel (used to enforce/augment CSC_bulk module)
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

ROOT <- "."
RAW_DIR <- file.path(ROOT, "data/raw/GSE190411")
INTERIM_DIR <- file.path(ROOT, "data/interim/rnaseq")
PROCESSED_DIR <- file.path(ROOT, "data/processed/rnaseq")
FIG_DIR <- file.path(ROOT, "figures/rnaseq")
GENESET_YAML <- file.path(ROOT, "config/gene_sets_rnaseq.yaml")

#-------------------------------------------------------------------------------
# Basic path checks and directory creation
#-------------------------------------------------------------------------------

if (!dir.exists(RAW_DIR)) {
  stop(
    "[ERROR] Raw RNA-seq directory not found: ", RAW_DIR,
    "\nRun this script from the project root."
  )
}
if (!file.exists(GENESET_YAML)) {
  stop("[ERROR] Gene set YAML not found at: ", GENESET_YAML)
}

dir.create(INTERIM_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

message("\n========================================")
message("BULK RNA-SEQ PROCESSING PIPELINE")
message("WITH STATISTICAL VALIDATION")
message("========================================\n")

#-------------------------------------------------------------------------------
# Validation report setup
#-------------------------------------------------------------------------------

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

# ============================================================================ #
# STEP 1: BUILD EXPRESSION MATRICES (MATCHING OLD 02 LOGIC)
# ============================================================================ #

message("[STEP 1/8] Building log2-transformed expression matrices")

log2p1 <- function(x) log2(x + 1)

f_bl6 <- file.path(RAW_DIR, "GSE190411_Yuan2021_Bl6_Norm_Pap_SCC_RNAseq_counts.csv.gz")
f_pap <- file.path(RAW_DIR, "GSE190411_Yuan2021_PAP_SCC_RNAseq_counts.txt.gz")
f_pdv <- file.path(RAW_DIR, "GSE190411_Yuan2021_PDV_WT_leprKO_RNAseq_counts.txt.gz")

if (!file.exists(f_bl6)) stop("[ERROR] Missing Bl6 file: ", f_bl6)
if (!file.exists(f_pap)) stop("[ERROR] Missing PAP/SCC file: ", f_pap)
if (!file.exists(f_pdv)) stop("[ERROR] Missing PDV file: ", f_pdv)

# 1A: Bl6 Normal/Papilloma/SCC
message("  Processing Bl6 dataset...")
bl6_raw <- read_csv(f_bl6, show_col_types = FALSE)

if (!"transcript_id" %in% colnames(bl6_raw)) {
  stop(
    "[ERROR] Expected column 'transcript_id' in Bl6 file, found: ",
    paste(colnames(bl6_raw)[1:min(5, ncol(bl6_raw))], collapse = ", ")
  )
}

bl6_expr <- bl6_raw %>%
  mutate(gene_symbol = str_replace(transcript_id, "^[^_]+_", "")) %>%
  select(gene_symbol, everything(), -transcript_id) %>%
  group_by(gene_symbol) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(across(where(is.numeric), log2p1))

write_tsv(bl6_expr, file.path(INTERIM_DIR, "bl6_expression.tsv"))

# 1B: PAP/SCC
message("  Processing PAP/SCC dataset...")
pap_raw <- read_tsv(f_pap, show_col_types = FALSE)

if (!"gene" %in% colnames(pap_raw)) {
  stop(
    "[ERROR] Expected column 'gene' in PAP/SCC file, found: ",
    paste(colnames(pap_raw)[1:min(5, ncol(pap_raw))], collapse = ", ")
  )
}

pap_expr <- pap_raw %>%
  rename(gene_symbol = gene) %>%
  group_by(gene_symbol) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(across(where(is.numeric), log2p1))

write_tsv(pap_expr, file.path(INTERIM_DIR, "pap_scc_expression.tsv"))

# 1C: PDV WT / LeprKO
message("  Processing PDV dataset...")
pdv_raw <- read_tsv(f_pdv, show_col_types = FALSE)

if (!"gene" %in% colnames(pdv_raw)) {
  stop(
    "[ERROR] Expected column 'gene' in PDV file, found: ",
    paste(colnames(pdv_raw)[1:min(5, ncol(pdv_raw))], collapse = ", ")
  )
}

pdv_expr <- pdv_raw %>%
  rename(gene_symbol = gene) %>%
  group_by(gene_symbol) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(across(where(is.numeric), log2p1))

write_tsv(pdv_expr, file.path(INTERIM_DIR, "pdv_expression.tsv"))

message("  ✓ Expression matrices written")

write_validation("Expression matrices built:")
write_validation(sprintf("  Bl6: %d genes x %d samples", nrow(bl6_expr) - 1, ncol(bl6_expr) - 1))
write_validation(sprintf("  PAP_SCC: %d genes x %d samples", nrow(pap_expr) - 1, ncol(pap_expr) - 1))
write_validation(sprintf("  PDV: %d genes x %d samples", nrow(pdv_expr) - 1, ncol(pdv_expr) - 1))

# ============================================================================ #
# STEP 2: BUILD SAMPLE METADATA (CANONICAL MAPPING FROM OLD 02)
# ============================================================================ #

message("\n[STEP 2/8] Building sample metadata")

ids_bl6 <- setdiff(colnames(bl6_expr), "gene_symbol")
ids_pap <- setdiff(colnames(pap_expr), "gene_symbol")
ids_pdv <- setdiff(colnames(pdv_expr), "gene_symbol")

metadata <- tibble(
  sample_id = c(ids_bl6, ids_pap, ids_pdv),
  dataset = c(
    rep("Bl6", length(ids_bl6)),
    rep("PAP_SCC", length(ids_pap)),
    rep("PDV", length(ids_pdv))
  )
) %>%
  mutate(
    sample_id_lower = str_to_lower(sample_id),
    condition = case_when(
      str_detect(sample_id_lower, "norm") ~ "Normal",
      str_detect(sample_id_lower, "pap") ~ "Papilloma",
      str_detect(sample_id_lower, "leprko") ~ "PDV_LeprKO",
      str_detect(sample_id_lower, "^ko_rep") ~ "PDV_LeprKO",
      str_detect(sample_id_lower, "pdv") ~ "PDV_WT",
      str_detect(sample_id_lower, "^wt_rep") ~ "PDV_WT",
      str_detect(sample_id_lower, "tumor") ~ "SCC",
      str_detect(sample_id_lower, "scc") ~ "SCC",
      str_detect(sample_id_lower, "tgfbr") ~ "SCC",
      str_detect(sample_id_lower, "tetohras") ~ "SCC",
      TRUE ~ "UNKNOWN"
    ),
    replicate = 1L
  ) %>%
  select(sample_id, dataset, condition, replicate)

unknowns <- metadata %>% filter(condition == "UNKNOWN")
if (nrow(unknowns) > 0) {
  stop(
    "[ERROR] Unknown conditions found. Update rules in case_when:\n  ",
    paste(unknowns$sample_id, collapse = ", ")
  )
}

write_csv(metadata, file.path(INTERIM_DIR, "sample_metadata_GSE190411.csv"))
message("  ✓ Metadata written")

write_validation("\nSample metadata (condition counts):")
cond_table <- table(metadata$condition)
for (cond in names(cond_table)) {
  write_validation(sprintf("  %s: %d samples", cond, cond_table[cond]))
}

write_validation("\nDataset x condition table:")
capture.output(print(table(metadata$dataset, metadata$condition))) |>
  writeLines(con = validation_log)

# ============================================================================ #
# STEP 3: COMPUTE MODULE SCORES (WITH DEFENSIVE CHECKS)
# ============================================================================ #

message("\n[STEP 3/8] Computing module scores")

gene_sets <- read_yaml(GENESET_YAML)
if (!is.list(gene_sets) || length(gene_sets) == 0L) {
  stop("[ERROR] Gene set YAML did not parse as a non-empty list: ", GENESET_YAML)
}

# Enforce/augment CSC module with Yuan CSC 101 markers
if ("CSC_bulk" %in% names(gene_sets)) {
  gene_sets$CSC_bulk <- union(gene_sets$CSC_bulk, YUAN_CSC_101_MARKERS)
} else {
  gene_sets$CSC_bulk <- YUAN_CSC_101_MARKERS
}

message(
  "  Gene sets loaded (CSC_bulk includes Yuan CSC 101): ",
  paste(names(gene_sets), collapse = ", ")
)

write_validation("\nModule gene sets:")
for (mod in names(gene_sets)) {
  write_validation(sprintf("  %s: %d genes", mod, length(gene_sets[[mod]])))
}

compute_modules <- function(expr_df, dataset_name, gene_sets, metadata) {
  message("  Computing modules for: ", dataset_name)

  expr_long <- expr_df %>%
    pivot_longer(-gene_symbol, names_to = "sample_id", values_to = "expr")

  meta_sub <- metadata %>%
    filter(dataset == dataset_name)

  if (nrow(meta_sub) == 0L) {
    stop("[ERROR] No metadata rows for dataset: ", dataset_name)
  }

  expr_samples <- unique(expr_long$sample_id)
  meta_samples <- unique(meta_sub$sample_id)

  missing_in_meta <- setdiff(expr_samples, meta_samples)
  if (length(missing_in_meta) > 0L) {
    stop(
      "[ERROR] Dataset ", dataset_name,
      " has expression samples with no metadata rows:\n  ",
      paste(missing_in_meta, collapse = ", ")
    )
  }

  missing_in_expr <- setdiff(meta_samples, expr_samples)
  if (length(missing_in_expr) > 0L) {
    warning(
      "[WARN] Dataset ", dataset_name,
      " has metadata samples with no expression columns:\n  ",
      paste(missing_in_expr, collapse = ", ")
    )
    meta_sub <- meta_sub %>%
      filter(sample_id %in% expr_samples)
  }

  expr_long <- expr_long %>%
    group_by(gene_symbol) %>%
    mutate(
      expr_z = (expr - mean(expr, na.rm = TRUE)) / sd(expr, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      expr_z = if_else(is.na(expr_z), 0, expr_z)
    )

  module_list <- list()

  for (mod in names(gene_sets)) {
    genes <- gene_sets[[mod]]
    present <- intersect(genes, unique(expr_long$gene_symbol))

    if (length(present) == 0L) {
      warning("  [WARN] Module ", mod, " has no genes in dataset ", dataset_name)
      next
    }

    message(sprintf("    %s: %d/%d genes present", mod, length(present), length(genes)))

    mod_scores <- expr_long %>%
      filter(gene_symbol %in% present) %>%
      group_by(sample_id) %>%
      summarise(score = mean(expr_z, na.rm = TRUE), .groups = "drop") %>%
      mutate(module = mod)

    module_list[[mod]] <- mod_scores
  }

  if (length(module_list) == 0L) {
    stop("[ERROR] No module scores computed for dataset: ", dataset_name)
  }

  scores <- bind_rows(module_list) %>%
    left_join(meta_sub, by = "sample_id")

  if (any(is.na(scores$dataset)) || any(is.na(scores$condition))) {
    stop(
      "[ERROR] After join, some rows have NA dataset/condition for dataset: ",
      dataset_name
    )
  }

  scores %>%
    select(sample_id, dataset, condition, module, score) %>%
    pivot_wider(names_from = module, values_from = score)
}

scores_bl6 <- compute_modules(bl6_expr, "Bl6", gene_sets, metadata)
scores_pap <- compute_modules(pap_expr, "PAP_SCC", gene_sets, metadata)
scores_pdv <- compute_modules(pdv_expr, "PDV", gene_sets, metadata)

all_scores <- bind_rows(scores_bl6, scores_pap, scores_pdv)

required_modules <- c("TGFb_bulk", "mTOR_bulk", "Angio_bulk", "CSC_bulk")
missing_mods <- setdiff(required_modules, colnames(all_scores))
if (length(missing_mods) > 0L) {
  stop(
    "[ERROR] The following required modules are missing from all_scores:\n  ",
    paste(missing_mods, collapse = ", ")
  )
}

write_csv(all_scores, file.path(PROCESSED_DIR, "module_scores_by_sample.csv"))
message("  ✓ Module scores written")

# ============================================================================ #
# STEP 4: STATISTICAL VALIDATION
# ============================================================================ #

message("\n[STEP 4/8] Statistical validation of module scores")

scores_long <- all_scores %>%
  pivot_longer(
    cols = all_of(required_modules),
    names_to = "module",
    values_to = "score"
  )

write_validation("\n==============================================")
write_validation("STATISTICAL TESTS: Module ~ Condition")
write_validation("==============================================\n")

module_stats <- list()

for (mod in required_modules) {
  test_data <- scores_long %>%
    filter(module == mod) %>%
    mutate(condition = factor(
      condition,
      levels = c("Normal", "Papilloma", "SCC", "PDV_WT", "PDV_LeprKO")
    ))

  anova_result <- aov(score ~ condition, data = test_data)
  anova_summary <- summary(anova_result)
  f_stat <- anova_summary[[1]]$`F value`[1]
  p_val <- anova_summary[[1]]$`Pr(>F)`[1]

  ss_condition <- anova_summary[[1]]$`Sum Sq`[1]
  ss_total <- sum(anova_summary[[1]]$`Sum Sq`)
  eta_squared <- ss_condition / ss_total

  write_validation(sprintf("%s:", mod))
  write_validation(sprintf(
    "  ANOVA: F=%.2f, p=%.2e, η²=%.3f %s",
    f_stat, p_val, eta_squared,
    ifelse(p_val < 0.05, "✓ SIGNIFICANT", "✗ NS")
  ))

  pairwise_tests <- list()

  norm_pap <- test_data %>% filter(condition %in% c("Normal", "Papilloma"))
  if (nrow(norm_pap) > 0) {
    t_test <- t.test(score ~ condition, data = norm_pap)
    cohens_d <- (mean(norm_pap$score[norm_pap$condition == "Papilloma"]) -
      mean(norm_pap$score[norm_pap$condition == "Normal"])) /
      sd(norm_pap$score)
    pairwise_tests[["Normal→Papilloma"]] <- list(p = t_test$p.value, d = cohens_d)
    write_validation(sprintf("  Normal→Papilloma: p=%.2e, d=%.2f", t_test$p.value, cohens_d))
  }

  pap_scc <- test_data %>% filter(condition %in% c("Papilloma", "SCC"))
  if (nrow(pap_scc) > 0) {
    t_test <- t.test(score ~ condition, data = pap_scc)
    cohens_d <- (mean(pap_scc$score[pap_scc$condition == "SCC"]) -
      mean(pap_scc$score[pap_scc$condition == "Papilloma"])) /
      sd(pap_scc$score)
    pairwise_tests[["Papilloma→SCC"]] <- list(p = t_test$p.value, d = cohens_d)
    write_validation(sprintf("  Papilloma→SCC: p=%.2e, d=%.2f", t_test$p.value, cohens_d))
  }

  pdv_comp <- test_data %>% filter(condition %in% c("PDV_WT", "PDV_LeprKO"))
  if (nrow(pdv_comp) > 0) {
    t_test <- t.test(score ~ condition, data = pdv_comp)
    cohens_d <- (mean(pdv_comp$score[pdv_comp$condition == "PDV_LeprKO"]) -
      mean(pdv_comp$score[pdv_comp$condition == "PDV_WT"])) /
      sd(pdv_comp$score)
    pairwise_tests[["PDV_WT→LeprKO"]] <- list(p = t_test$p.value, d = cohens_d)
    write_validation(sprintf(
      "  PDV_WT→LeprKO: p=%.2e, d=%.2f %s",
      t_test$p.value, cohens_d,
      ifelse(abs(cohens_d) > 0.8, "✓ LARGE EFFECT", "")
    ))
  }

  module_stats[[mod]] <- list(
    anova_p = p_val,
    anova_f = f_stat,
    eta_squared = eta_squared,
    pairwise = pairwise_tests
  )

  write_validation("")
}

# ============================================================================ #
# STEP 5: CROSS-DATASET VALIDATION
# ============================================================================ #

message("\n[STEP 5/8] Cross-dataset validation")

write_validation("==============================================")
write_validation("CROSS-DATASET CONSISTENCY")
write_validation("==============================================\n")

bl6_scc <- scores_long %>% filter(dataset == "Bl6", condition == "SCC")
pap_scc <- scores_long %>% filter(dataset == "PAP_SCC", condition == "SCC")

if (nrow(bl6_scc) > 0 && nrow(pap_scc) > 0) {
  bl6_means <- bl6_scc %>%
    group_by(module) %>%
    summarise(score_bl6 = mean(score), .groups = "drop")
  pap_means <- pap_scc %>%
    group_by(module) %>%
    summarise(score_pap = mean(score), .groups = "drop")

  cross_comp <- inner_join(bl6_means, pap_means, by = "module")

  if (nrow(cross_comp) > 0) {
    cor_test <- cor.test(cross_comp$score_bl6, cross_comp$score_pap, method = "spearman")

    write_validation("SCC samples across datasets (Bl6 vs PAP_SCC):")
    write_validation(sprintf(
      "  Spearman correlation: r=%.3f, p=%.3f %s",
      cor_test$estimate, cor_test$p.value,
      ifelse(cor_test$p.value < 0.05 & abs(cor_test$estimate) > 0.7,
        "✓ CONSISTENT", "⚠ WEAK"
      )
    ))

    for (i in seq_len(nrow(cross_comp))) {
      write_validation(sprintf(
        "    %s: Bl6=%.2f, PAP=%.2f, diff=%.2f",
        cross_comp$module[i],
        cross_comp$score_bl6[i],
        cross_comp$score_pap[i],
        abs(cross_comp$score_bl6[i] - cross_comp$score_pap[i])
      ))
    }
  }
}

write_validation("")

# ============================================================================ #
# STEP 6: PCA FOR BATCH EFFECTS
# ============================================================================ #

message("\n[STEP 6/8] PCA and batch effect assessment")

pca_matrix <- all_scores %>%
  select(sample_id, dataset, condition, all_of(required_modules)) %>%
  column_to_rownames("sample_id") %>%
  select(all_of(required_modules))

pca_meta <- all_scores %>%
  select(sample_id, dataset, condition)

pca_result <- prcomp(pca_matrix, scale. = TRUE, center = TRUE)
pca_df <- as.data.frame(pca_result$x) %>%
  rownames_to_column("sample_id") %>%
  left_join(pca_meta, by = "sample_id")

var_explained <- summary(pca_result)$importance[2, 1:2] * 100

write_validation("PCA Analysis:")
write_validation(sprintf("  PC1 variance: %.1f%%", var_explained[1]))
write_validation(sprintf("  PC2 variance: %.1f%%", var_explained[2]))

pc1_aov <- aov(PC1 ~ dataset, data = pca_df)
pc1_p <- summary(pc1_aov)[[1]]$`Pr(>F)`[1]

write_validation(sprintf(
  "  PC1 ~ dataset: p=%.3f %s",
  pc1_p,
  ifelse(pc1_p < 0.05, "⚠ BATCH EFFECT", "✓ NO BATCH EFFECT")
))

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
  p_pca,
  width = 7, height = 5, dpi = 300
)

# ============================================================================ #
# STEP 7: PUBLICATION FIGURES (BASIC)
# ============================================================================ #

message("\n[STEP 7/8] Generating publication figures")

scores_long_plot <- scores_long %>%
  mutate(
    condition = factor(condition,
      levels = c("Normal", "Papilloma", "SCC", "PDV_WT", "PDV_LeprKO")
    ),
    module = factor(module,
      levels = c("TGFb_bulk", "Angio_bulk", "CSC_bulk", "mTOR_bulk")
    )
  )

p_violins <- ggplot(scores_long_plot, aes(x = condition, y = score, fill = condition)) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6) +
  facet_wrap(~module, scales = "free_y", ncol = 2) +
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
  ggpubr::stat_compare_means(method = "anova", label.y.npc = 0.95, size = 3)

ggsave(file.path(FIG_DIR, "rnaseq_module_violins_with_stats.png"),
  p_violins,
  width = 9, height = 8, dpi = 400
)

# --------------------------------------------------------------------------- #
# 7B. Sample-level ComplexHeatmap with row annotations (fixed)
# --------------------------------------------------------------------------- #

if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
  requireNamespace("circlize", quietly = TRUE)) {
  message("  Generating ComplexHeatmap for sample-level module scores...")

  sample_mat <- all_scores %>%
    select(sample_id, dataset, condition, all_of(required_modules)) %>%
    distinct() %>%
    column_to_rownames("sample_id")

  sample_meta <- sample_mat %>%
    as.data.frame() %>%
    select(dataset, condition)

  sample_mat <- sample_mat[, required_modules, drop = TRUE]

  sample_mat_scaled <- scale(sample_mat)
  rownames(sample_mat_scaled) <- rownames(sample_mat)
  colnames(sample_mat_scaled) <- colnames(sample_mat)

  ht_col_fun <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("blue", "white", "red")
  )

  dataset_levels <- sort(unique(sample_meta$dataset))
  condition_levels <- c("Normal", "Papilloma", "SCC", "PDV_WT", "PDV_LeprKO")
  condition_levels <- condition_levels[condition_levels %in% unique(sample_meta$condition)]

  dataset_cols <- setNames(
    RColorBrewer::brewer.pal(max(3, length(dataset_levels)), "Set1")[seq_along(dataset_levels)],
    dataset_levels
  )

  condition_cols <- setNames(
    RColorBrewer::brewer.pal(max(5, length(condition_levels)), "Set2")[seq_along(condition_levels)],
    condition_levels
  )

  ann_df <- sample_meta %>%
    mutate(
      dataset = factor(dataset, levels = dataset_levels),
      condition = factor(condition, levels = condition_levels)
    )

  ha_row <- ComplexHeatmap::rowAnnotation(
    df = ann_df,
    col = list(
      Dataset = dataset_cols,
      Condition = condition_cols
    ),
    annotation_name_side = "top"
  )

  row_split_vec <- factor(
    ann_df$condition,
    levels = condition_levels
  )

  ht <- ComplexHeatmap::Heatmap(
    sample_mat_scaled,
    name = "Z-score",
    col = ht_col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    left_annotation = ha_row,
    row_split = row_split_vec,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_title = "Ras–CSC loop modules",
    row_title = "Samples",
    heatmap_legend_param = list(
      title = "Module\nz-score",
      at = c(-2, 0, 2),
      labels = c("Low", "Mid", "High")
    )
  )

  png(file.path(FIG_DIR, "rnaseq_module_complex_heatmap_samples.png"),
    width = 2600, height = 2200, res = 300
  )
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legend = TRUE
  )
  dev.off()
} else {
  message("  [WARN] ComplexHeatmap or circlize not available; skipping advanced heatmap.")
}

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
  scale_x_continuous(
    breaks = c(1, 2, 3),
    labels = c("Normal", "Papilloma", "SCC")
  ) +
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
  p_trajectory,
  width = 7, height = 5, dpi = 300
)

# ============================================================================ #
# STEP 7B: ADVANCED HYPOTHESIS-DRIVEN FIGURES (Ras–CSC LOOP)
# ============================================================================ #

message("\n[STEP 7B] Generating hypothesis-driven loop figures")

# 7B.1: Pairwise loop edges (A ↔ T ↔ M ↔ C) as scatter plots with correlations

loop_df <- all_scores %>%
  select(sample_id, dataset, condition, all_of(required_modules)) %>%
  mutate(
    condition = factor(
      condition,
      levels = c("Normal", "Papilloma", "SCC", "PDV_WT", "PDV_LeprKO")
    )
  )

loop_pairs <- tibble::tibble(
  x_mod = c("Angio_bulk", "Angio_bulk", "TGFb_bulk", "mTOR_bulk"),
  y_mod = c("TGFb_bulk", "CSC_bulk", "mTOR_bulk", "CSC_bulk"),
  edge_label = c("A → T", "C → A", "T → M", "M → C")
)

loop_plots <- list()

for (i in seq_len(nrow(loop_pairs))) {
  x_mod <- loop_pairs$x_mod[i]
  y_mod <- loop_pairs$y_mod[i]
  edge_label <- loop_pairs$edge_label[i]

  pair_df <- loop_df %>%
    filter(
      !is.na(.data[[x_mod]]),
      !is.na(.data[[y_mod]])
    )

  if (nrow(pair_df) < 4) {
    next
  }

  p_pair <- ggplot(
    pair_df,
    aes(
      x = .data[[x_mod]],
      y = .data[[y_mod]],
      color = condition,
      shape = dataset
    )
  ) +
    geom_point(size = 2.5, alpha = 0.85) +
    geom_smooth(method = "lm", se = FALSE, size = 0.8, linetype = "dashed") +
    ggpubr::stat_cor(
      method = "spearman",
      label.x.npc = "left",
      label.y.npc = "top",
      size = 3
    ) +
    scale_color_brewer(palette = "Set1", name = "Condition") +
    scale_shape_manual(values = c(16, 17, 15), name = "Dataset") +
    labs(
      title = paste0("Loop edge: ", edge_label),
      x = x_mod,
      y = y_mod
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  loop_plots[[i]] <- p_pair
}

if (length(loop_plots) > 0) {
  p_loop_grid <- ggpubr::ggarrange(
    plotlist = loop_plots,
    ncol = 2,
    nrow = ceiling(length(loop_plots) / 2),
    common.legend = TRUE,
    legend = "right"
  )

  ggsave(
    file.path(FIG_DIR, "rnaseq_loop_pairs_scatter.png"),
    p_loop_grid,
    width = 9,
    height = 7,
    dpi = 400
  )
}

# 7B.2: Condition-level Ras–CSC loop fingerprint (C, A, T, M) per dataset

loop_summary <- all_scores %>%
  select(sample_id, dataset, condition, all_of(required_modules)) %>%
  group_by(dataset, condition) %>%
  summarise(
    across(
      all_of(required_modules),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  mutate(
    condition = factor(
      condition,
      levels = c("Normal", "Papilloma", "SCC", "PDV_WT", "PDV_LeprKO")
    )
  ) %>%
  tidyr::pivot_longer(
    cols = all_of(required_modules),
    names_to = "module",
    values_to = "mean_score"
  ) %>%
  mutate(
    module = factor(
      module,
      levels = c("TGFb_bulk", "Angio_bulk", "CSC_bulk", "mTOR_bulk")
    )
  )

if (nrow(loop_summary) > 0) {
  p_parallel <- ggplot(
    loop_summary,
    aes(
      x = module,
      y = mean_score,
      group = condition,
      color = condition
    )
  ) +
    geom_line(size = 1.1, alpha = 0.9) +
    geom_point(size = 2) +
    facet_wrap(~dataset, nrow = 1) +
    scale_color_brewer(palette = "Set1", name = "Condition") +
    labs(
      title = "Ras–CSC loop fingerprint by condition and dataset",
      x = "Module (proxy for loop component)",
      y = "Mean module score (z-score)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(
    file.path(FIG_DIR, "rnaseq_loop_parallel_fingerprint.png"),
    p_parallel,
    width = 9,
    height = 4.5,
    dpi = 400
  )
}

# 7B.3: Module–module correlation structure (Spearman) as ComplexHeatmap

if (nrow(all_scores) > 3 &&
  requireNamespace("ComplexHeatmap", quietly = TRUE) &&
  requireNamespace("circlize", quietly = TRUE)) {
  corr_mat <- all_scores %>%
    select(all_of(required_modules)) %>%
    as.matrix()

  corr_spearman <- suppressWarnings(
    cor(
      corr_mat,
      use = "pairwise.complete.obs",
      method = "spearman"
    )
  )

  col_fun_corr <- circlize::colorRamp2(
    c(-1, 0, 1),
    c("blue", "white", "red")
  )

  ht_corr <- ComplexHeatmap::Heatmap(
    corr_spearman,
    name = "rho",
    col = col_fun_corr,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_title = "Spearman correlations between Ras–CSC loop modules"
  )

  png(
    file.path(FIG_DIR, "rnaseq_module_correlation_heatmap.png"),
    width = 2000,
    height = 2000,
    res = 300
  )
  ComplexHeatmap::draw(
    ht_corr,
    heatmap_legend_side = "right",
    merge_legend = TRUE
  )
  dev.off()
} else {
  message("  [WARN] Skipping correlation ComplexHeatmap; n<4 or ComplexHeatmap/circlize missing.")
}

# 7B.4: Cross-dataset SCC comparison (Bl6 vs PAP_SCC) for loop modules

scc_cross <- scores_long %>%
  filter(condition == "SCC") %>%
  group_by(dataset, module) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = dataset,
    values_from = mean_score
  )

if (all(c("Bl6", "PAP_SCC") %in% colnames(scc_cross))) {
  scc_plot_df <- scc_cross %>%
    filter(module %in% required_modules)

  if (nrow(scc_plot_df) > 0) {
    p_cross <- ggplot(
      scc_plot_df,
      aes(x = Bl6, y = PAP_SCC, label = module)
    ) +
      geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
      geom_point(size = 3) +
      ggrepel::geom_text_repel(size = 3.5) +
      coord_equal() +
      labs(
        title = "SCC module scores: Bl6 vs PAP_SCC",
        x = "Mean SCC module score (Bl6)",
        y = "Mean SCC module score (PAP_SCC)"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid.minor = element_blank()
      )

    ggsave(
      file.path(FIG_DIR, "rnaseq_scc_cross_dataset_scatter.png"),
      p_cross,
      width = 5,
      height = 5,
      dpi = 400
    )
  }
}

message("[STEP 7B] Advanced figures written to:")
message("  ", file.path(FIG_DIR, "rnaseq_loop_pairs_scatter.png"))
message("  ", file.path(FIG_DIR, "rnaseq_loop_parallel_fingerprint.png"))
message("  ", file.path(FIG_DIR, "rnaseq_module_correlation_heatmap.png"))
message("  ", file.path(FIG_DIR, "rnaseq_scc_cross_dataset_scatter.png"))

# ============================================================================ #
# STEP 8: EXPORT CALIBRATION TARGETS WITH CONFIDENCE INTERVALS
# ============================================================================ #

message("\n[STEP 8/8] Exporting calibration targets with confidence intervals")

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

targets_detailed <- targets_with_ci %>%
  select(dataset, condition, module, mean_score, se_score, ci_lower, ci_upper) %>%
  pivot_wider(
    names_from = module,
    values_from = c(mean_score, se_score, ci_lower, ci_upper)
  )

calib_plot_df <- targets_with_ci %>%
  filter(module %in% required_modules) %>%
  mutate(
    condition = factor(
      condition,
      levels = c("Normal", "Papilloma", "SCC", "PDV_WT", "PDV_LeprKO")
    ),
    module = factor(
      module,
      levels = c("TGFb_bulk", "Angio_bulk", "CSC_bulk", "mTOR_bulk")
    )
  )

if (nrow(calib_plot_df) > 0) {
  p_calib <- ggplot(
    calib_plot_df,
    aes(x = condition, y = mean_score, color = module, group = module)
  ) +
    geom_point(position = position_dodge(width = 0.4), size = 2.4) +
    geom_errorbar(
      aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.2,
      position = position_dodge(width = 0.4),
      alpha = 0.9
    ) +
    facet_wrap(~dataset, nrow = 1) +
    scale_color_brewer(palette = "Set1", name = "Module") +
    labs(
      title = "Calibration targets with 95% CI by dataset and condition",
      x = "Condition",
      y = "Module score (mean ± 95% CI)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )

  ggsave(
    file.path(FIG_DIR, "rnaseq_calibration_targets_ci_by_condition.png"),
    p_calib,
    width = 10,
    height = 4,
    dpi = 400
  )
}

write_csv(targets, file.path(PROCESSED_DIR, "ras_csc_calibration_targets.csv"))
write_csv(targets_detailed, file.path(PROCESSED_DIR, "ras_csc_calibration_targets_with_ci.csv"))

write_validation("\n==============================================")
write_validation("CALIBRATION TARGETS EXPORTED")
write_validation("==============================================\n")

for (i in seq_len(nrow(targets))) {
  write_validation(sprintf("%s - %s:", targets$dataset[i], targets$condition[i]))
  write_validation(sprintf(
    "  C=%.3f, A=%.3f, T=%.3f, M=%.3f",
    targets$C_target[i], targets$A_target[i],
    targets$T_target[i], targets$M_target[i]
  ))
}

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
