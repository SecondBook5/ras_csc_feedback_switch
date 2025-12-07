#!/usr/bin/env Rscript
# pipelines/data_processing/03_process_atac.R
#
# COMPREHENSIVE ATAC-SEQ PIPELINE WITH HYPOTHESIS TESTING
# Processes GSE190414 narrowPeak files for chromatin accessibility
#
# Pipeline steps:
#   1. Import narrowPeak files and annotate to genes (mm10)
#   2. Build gene-level accessibility matrix
#   3. Compute module scores with z-scoring
#   4. Statistical validation across conditions
#   5. Cross-omics comparison (ATAC vs RNA-seq)
#   6. Mechanistic hypothesis testing
#   7. Advanced publication figures
#   8. Export calibration targets
#   9. Summary statistics
#
# Outputs:
#   data/interim/atac/atac_gene_accessibility.tsv
#   data/processed/atac/atac_module_scores_by_sample.csv
#   data/processed/atac/atac_calibration_targets.csv
#   data/processed/atac/atac_validation_report.txt
#   figures/atac/atac_*.png

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(stringr)
    library(yaml)
    library(ggplot2)
    library(ggpubr)
    library(patchwork)
    library(RColorBrewer)
    library(GenomicRanges)
    library(rtracklayer)
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    library(org.Mm.eg.db)
    library(ChIPseeker)
    # Advanced heatmaps for adjacency / sample structure
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
})

# Hard-coded Yuan CSC 101 markers (used for CSC_Yuan101 module)
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

message("=== RUNNING pipelines/data_processing/03_process_atac.R ===")
condition <- NULL

message("\n========================================")
message("ATAC-SEQ PROCESSING PIPELINE")
message("WITH HYPOTHESIS TESTING")
message("========================================\n")

# Paths
ROOT <- "."
RAW_DIR <- file.path(ROOT, "data/raw/GSE190414")
INTERIM_DIR <- file.path(ROOT, "data/interim/atac")
PROCESSED_DIR <- file.path(ROOT, "data/processed/atac")
FIG_DIR <- file.path(ROOT, "figures/atac")
GENESET_YAML <- file.path(ROOT, "config/gene_sets_rnaseq.yaml")

# RNA-seq results for cross-omics comparison
RNA_SCORES <- file.path(ROOT, "data/processed/rnaseq/module_scores_by_sample.csv")

# Create directories
dir.create(INTERIM_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Validation report
validation_report <- file.path(PROCESSED_DIR, "atac_validation_report.txt")
validation_log <- file(validation_report, "w")

write_validation <- function(msg) {
    message(msg)
    writeLines(msg, validation_log)
}

write_validation("==============================================")
write_validation("ATAC-SEQ VALIDATION REPORT")
write_validation(paste("Generated:", Sys.time()))
write_validation("==============================================\n")

# =============================================================================
# STEP 1: DEFINE SAMPLES AND FILE MAPPING
# =============================================================================

message("[STEP 1/9] Defining ATAC-seq samples")

samples <- tibble::tribble(
    ~sample_id,       ~npeak_file,                                          ~condition,        ~stage,
    "HFSC_ATAC",      "GSE190414_HFSC_ATAC_Rep1.narrowPeak.gz",             "HFSC",            "Normal",
    "IFE_ATAC",       "GSE190414_IFE_ATAC_Rep1.narrowPeak.gz",              "IFE",             "Normal",
    "PAP_mneg_ATAC",  "GSE190414_PAP_mneg_ATAC_Reps_POOL.narrowPeak.gz",    "Papilloma_mneg",  "Papilloma",
    "PAP_mpos_ATAC",  "GSE190414_PAP_mpos_ATAC_Reps_POOL.narrowPeak.gz",    "Papilloma_mpos",  "Papilloma",
    "SCC_mneg_ATAC",  "GSE190414_SCC_mneg_ATAC_Reps_POOL.narrowPeak.gz",    "SCC_mneg",        "SCC",
    "SCC_mpos_ATAC",  "GSE190414_SCC_mpos_ATAC_Reps_POOL.narrowPeak.gz",    "SCC_mpos",        "SCC"
) %>%
    mutate(
        path = file.path(RAW_DIR, npeak_file),
        dataset = "ATAC",
        # m+ status (CSC enriched in Yuan's paper)
        csc_enriched = str_detect(condition, "mpos")
    )

# Verify files exist
for (i in seq_len(nrow(samples))) {
    if (!file.exists(samples$path[i])) {
        stop("[ERROR] Missing narrowPeak file: ", samples$path[i])
    }
}

write_validation("ATAC-seq samples:")
write_validation(sprintf("  Total samples: %d", nrow(samples)))
write_validation("\nConditions:")
for (cond in unique(samples$condition)) {
    n <- sum(samples$condition == cond)
    stage <- unique(samples$stage[samples$condition == cond])
    write_validation(sprintf("  %s (%s): %d sample(s)", cond, stage, n))
}

# Save metadata
write_csv(
    samples %>% dplyr::select(sample_id, dataset, condition, stage, csc_enriched),
    file.path(INTERIM_DIR, "sample_metadata_ATAC.csv")
)

# =============================================================================
# STEP 2: ANNOTATE PEAKS TO GENES
# =============================================================================

message("\n[STEP 2/9] Annotating ATAC peaks to genes (mm10)")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoter_upstream <- 2000
promoter_downstream <- 2000

annotate_peaks_to_genes <- function(sample_row, txdb) {
    sample_id <- sample_row$sample_id
    peak_file <- sample_row$path

    message("  Processing: ", sample_id)

    # Import narrowPeak
    gr <- rtracklayer::import(peak_file, format = "narrowPeak")

    if (length(gr) == 0L) {
        warning("[WARN] No peaks in: ", peak_file)
        return(NULL)
    }

    # Annotate to nearest genes
    peak_anno <- ChIPseeker::annotatePeak(
        gr,
        TxDb = txdb,
        tssRegion = c(-promoter_upstream, promoter_downstream),
        annoDb = "org.Mm.eg.db",
        verbose = FALSE
    )

    anno_df <- as.data.frame(peak_anno)

    # Check for required columns
    if (!"SYMBOL" %in% colnames(anno_df)) {
        stop("[ERROR] Missing SYMBOL column for: ", sample_id)
    }

    # Use signalValue as accessibility measure
    signal_col <- if ("signalValue" %in% colnames(anno_df)) {
        "signalValue"
    } else if ("score" %in% colnames(anno_df)) {
        "score"
    } else {
        stop("[ERROR] No signal column for: ", sample_id)
    }

    # Collapse to gene level (mean accessibility per gene)
    gene_access <- anno_df %>%
        dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
        dplyr::group_by(SYMBOL) %>%
        dplyr::summarise(
            n_peaks = n(),
            access_raw = mean(.data[[signal_col]], na.rm = TRUE),
            .groups = "drop"
        ) %>%
        dplyr::mutate(
            gene_symbol = SYMBOL,
            access_log2 = log2(access_raw + 1),
            sample_id = sample_id
        ) %>%
        dplyr::select(sample_id, gene_symbol, n_peaks, access_log2)

    gene_access
}

# Process all samples
access_list <- list()
peak_stats <- list()

for (i in seq_len(nrow(samples))) {
    res <- annotate_peaks_to_genes(samples[i, ], txdb)
    if (!is.null(res)) {
        access_list[[samples$sample_id[i]]] <- res

        # Track peak statistics
        peak_stats[[samples$sample_id[i]]] <- tibble(
            sample_id = samples$sample_id[i],
            n_genes = nrow(res),
            total_peaks = sum(res$n_peaks),
            mean_peaks_per_gene = mean(res$n_peaks)
        )
    }
}

peak_summary <- bind_rows(peak_stats)

write_validation("\nPeak annotation statistics:")
for (i in seq_len(nrow(peak_summary))) {
    write_validation(sprintf(
        "  %s: %d genes, %d peaks (%.1f peaks/gene)",
        peak_summary$sample_id[i],
        peak_summary$n_genes[i],
        peak_summary$total_peaks[i],
        peak_summary$mean_peaks_per_gene[i]
    ))
}

# Combine to long format
access_long <- bind_rows(access_list)

# Convert to wide format (genes × samples)
access_wide <- access_long %>%
    dplyr::select(gene_symbol, sample_id, access_log2) %>%
    tidyr::pivot_wider(
        id_cols = gene_symbol,
        names_from = sample_id,
        values_from = access_log2,
        values_fill = 0 # Genes without peaks = 0 accessibility
    )

write_tsv(access_wide, file.path(INTERIM_DIR, "atac_gene_accessibility.tsv"))
message("  ✓ Gene accessibility matrix written")

# =============================================================================
# STEP 3: COMPUTE MODULE SCORES
# =============================================================================

message("\n[STEP 3/9] Computing module scores")

gene_sets <- read_yaml(GENESET_YAML)
message("  Gene sets read from: ", GENESET_YAML)
message("  Gene sets available: ", paste(names(gene_sets), collapse = ", "))

write_validation("\nModule gene sets from YAML:")
for (mod in names(gene_sets)) {
    write_validation(sprintf("  %s: %d genes", mod, length(gene_sets[[mod]])))
}

# Attach Yuan CSC 101 markers as a separate module (new CSC definition for chromatin)
gene_sets$CSC_Yuan101 <- YUAN_CSC_101_MARKERS
write_validation(sprintf(
    "\nAdded module CSC_Yuan101 (Yuan CSC 101 markers): %d genes",
    length(gene_sets$CSC_Yuan101)
))

# Long format for z-scoring
expr_long <- access_wide %>%
    tidyr::pivot_longer(
        cols = -gene_symbol,
        names_to = "sample_id",
        values_to = "access"
    )

# Z-score per gene across all ATAC samples
expr_long <- expr_long %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::mutate(
        access_z = (access - mean(access, na.rm = TRUE)) /
            sd(access, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        access_z = ifelse(is.na(access_z), 0, access_z)
    )

# Compute module scores
module_list <- list()

for (mod in names(gene_sets)) {
    genes_mod <- gene_sets[[mod]]
    present <- intersect(genes_mod, unique(expr_long$gene_symbol))

    if (length(present) == 0L) {
        warning("  [WARN] Module ", mod, " has no genes in ATAC data")
        next
    }

    message(sprintf("  %s: %d/%d genes present", mod, length(present), length(genes_mod)))

    mod_scores <- expr_long %>%
        dplyr::filter(gene_symbol %in% present) %>%
        dplyr::group_by(sample_id) %>%
        dplyr::summarise(
            score = mean(access_z, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        dplyr::mutate(module = mod)

    module_list[[mod]] <- mod_scores

    write_validation(sprintf(
        "  %s: %d/%d genes present (%.1f%%)",
        mod, length(present), length(genes_mod),
        100 * length(present) / length(genes_mod)
    ))
}

scores_long <- bind_rows(module_list)

# Join metadata
scores_long <- scores_long %>%
    left_join(
        samples %>% dplyr::select(sample_id, dataset, condition, stage, csc_enriched),
        by = "sample_id"
    )

# Wide format
scores_wide <- scores_long %>%
    dplyr::select(sample_id, dataset, condition, stage, csc_enriched, module, score) %>%
    tidyr::pivot_wider(names_from = module, values_from = score)
write_csv(scores_wide, file.path(PROCESSED_DIR, "atac_module_scores_by_sample.csv"))
message("  ✓ Module scores written")

# =============================================================================
# STEP 4: STATISTICAL VALIDATION
# =============================================================================

message("\n[STEP 4/9] Statistical validation")

# Required modules for downstream comparisons / loop adjacency
required_modules <- c("TGFb_bulk", "mTOR_bulk", "Angio_bulk", "CSC_bulk", "CSC_Yuan101")

# Check that required modules exist in gene_sets
missing_mods <- setdiff(required_modules, names(gene_sets))
if (length(missing_mods) > 0) {
    stop(
        "[ERROR] Missing required modules in ", GENESET_YAML, " or local definitions: ",
        paste(missing_mods, collapse = ", ")
    )
}

# Sanity check that CSC_Yuan101 has 101 genes
if ("CSC_Yuan101" %in% names(gene_sets)) {
    if (length(gene_sets$CSC_Yuan101) != 101L) {
        warning(
            "[WARN] CSC_Yuan101 has length ",
            length(gene_sets$CSC_Yuan101),
            " (expected 101). Confirm this is intentional."
        )
    }
}

write_validation("\n==============================================")
write_validation("STATISTICAL TESTS: ATAC Module Scores")
write_validation("==============================================\n")

# Test 1: CSC and other modules enrichment in m+ vs m-
for (mod in required_modules) {
    test_data <- scores_long %>%
        filter(module == mod, stage %in% c("Papilloma", "SCC"))

    if (nrow(test_data) > 2) {
        t_test <- t.test(score ~ csc_enriched, data = test_data)

        mean_mpos <- mean(test_data$score[test_data$csc_enriched], na.rm = TRUE)
        mean_mneg <- mean(test_data$score[!test_data$csc_enriched], na.rm = TRUE)

        cohens_d <- (mean_mpos - mean_mneg) / sd(test_data$score)

        write_validation(sprintf("%s (m+ vs m-):", mod))
        write_validation(sprintf("  m+ mean: %.3f", mean_mpos))
        write_validation(sprintf("  m- mean: %.3f", mean_mneg))
        write_validation(sprintf(
            "  t-test: p=%.3f, d=%.2f %s",
            t_test$p.value, cohens_d,
            ifelse(t_test$p.value < 0.05 & cohens_d > 0, "✓ CSC ENRICHED",
                ifelse(t_test$p.value < 0.05 & cohens_d < 0, "⚠ DEPLETED", "✗ NS")
            )
        ))
        write_validation("")
    }
}

# Test 2: Progression trends (Normal -> Papilloma -> SCC)
progression_summary <- scores_long %>%
    filter(module %in% required_modules) %>%
    group_by(stage, module) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

write_validation("Progression trends (mean scores):")
for (mod in required_modules) {
    prog_data <- progression_summary %>%
        filter(module == mod) %>%
        arrange(stage)

    norm <- prog_data$mean_score[prog_data$stage == "Normal"]
    pap <- prog_data$mean_score[prog_data$stage == "Papilloma"]
    scc <- prog_data$mean_score[prog_data$stage == "SCC"]

    if (length(norm) > 0 && length(scc) > 0) {
        fold_change <- scc / (norm + 0.001)

        write_validation(sprintf(
            "  %s: Normal=%.2f -> Pap=%.2f -> SCC=%.2f (FC=%.2f)",
            mod,
            ifelse(length(norm) > 0, norm, NA),
            ifelse(length(pap) > 0, mean(pap), NA),
            ifelse(length(scc) > 0, mean(scc), NA),
            fold_change
        ))
    }
}

# =============================================================================
# STEP 5: CROSS-OMICS COMPARISON (ATAC vs RNA-seq)
# =============================================================================

message("\n[STEP 5/9] Cross-omics comparison (ATAC vs RNA-seq)")

if (file.exists(RNA_SCORES)) {
    # Load RNA module scores
    rna_scores <- readr::read_csv(RNA_SCORES, show_col_types = FALSE)

    message("  RNA score columns: ", paste(colnames(rna_scores), collapse = ", "))

    # Ensure we have a 'condition' column
    if (!"condition" %in% colnames(rna_scores)) {
        if ("Condition" %in% colnames(rna_scores)) {
            message("  [INFO] Found 'Condition' column; renaming to 'condition'")
            rna_scores <- dplyr::rename(rna_scores, condition = Condition)
        } else {
            stop(
                "[ERROR] RNA module score file (",
                RNA_SCORES,
                ") does not contain a 'condition' column.\n",
                "       Columns present: ",
                paste(colnames(rna_scores), collapse = ", "),
                "\n       Re-run 01_process_bulk_rnaseq.R to regenerate module_scores_by_sample.csv."
            )
        }
    }

    # New column used for mapping, preserve original condition
    rna_scores <- rna_scores %>%
        dplyr::mutate(cond_for_atac = condition)

    # Mapping between RNA conditions and ATAC stages
    condition_mapping <- tibble::tribble(
        ~rna_condition, ~atac_stage,
        "Normal",       "Normal",
        "Papilloma",    "Papilloma",
        "SCC",          "SCC"
    )

    # Only use modules that exist in RNA file
    modules_for_cross_omics <- intersect(required_modules, colnames(rna_scores))

    # RNA mean scores by condition
    rna_summary <- rna_scores %>%
        tidyr::pivot_longer(
            cols = dplyr::all_of(modules_for_cross_omics),
            names_to = "module",
            values_to = "score"
        ) %>%
        dplyr::group_by(cond_for_atac, module) %>%
        dplyr::summarise(
            rna_mean = mean(score, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        dplyr::rename(rna_condition = cond_for_atac)

    # ATAC mean scores by stage
    atac_summary <- scores_long %>%
        dplyr::filter(module %in% modules_for_cross_omics) %>%
        dplyr::group_by(stage, module) %>%
        dplyr::summarise(
            atac_mean = mean(score, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        dplyr::rename(atac_stage = stage)

    # Join and compute concordance
    cross_omics <- condition_mapping %>%
        dplyr::left_join(rna_summary, by = "rna_condition") %>%
        dplyr::left_join(atac_summary, by = c("atac_stage", "module")) %>%
        dplyr::filter(!is.na(rna_mean) & !is.na(atac_mean))

    write_validation("\n==============================================")
    write_validation("CROSS-OMICS CONCORDANCE (RNA vs ATAC)")
    write_validation("==============================================\n")

    for (mod in unique(cross_omics$module)) {
        mod_data <- cross_omics %>% dplyr::filter(module == mod)

        if (nrow(mod_data) >= 3) {
            cor_test <- cor.test(
                mod_data$rna_mean,
                mod_data$atac_mean,
                method = "spearman"
            )

            write_validation(sprintf("%s:", mod))
            write_validation(sprintf(
                "  Spearman r = %.3f, p = %.3f %s",
                cor_test$estimate, cor_test$p.value,
                ifelse(
                    cor_test$p.value < 0.05 & abs(cor_test$estimate) > 0.7,
                    "✓ CONCORDANT", "⚠ WEAK"
                )
            ))

            for (i in seq_len(nrow(mod_data))) {
                write_validation(sprintf(
                    "    %s: RNA=%.2f, ATAC=%.2f",
                    mod_data$rna_condition[i],
                    mod_data$rna_mean[i],
                    mod_data$atac_mean[i]
                ))
            }
            write_validation("")
        }
    }

    readr::write_csv(
        cross_omics,
        file.path(PROCESSED_DIR, "cross_omics_rna_atac_comparison.csv")
    )
} else {
    write_validation("\n⚠ RNA-seq scores not found; skipping cross-omics comparison")
}

# =============================================================================
# STEP 6: MECHANISTIC HYPOTHESIS TESTING
# =============================================================================

message("\n[STEP 6/9] Testing mechanistic hypotheses")

write_validation("\n==============================================")
write_validation("MECHANISTIC HYPOTHESIS TESTING")
write_validation("==============================================\n")

# Prepare wide table with both CSC definitions if available
test_df_wide <- scores_wide %>%
    filter(!is.na(TGFb_bulk))

# H1: TGFb -> CSC for both definitions
if (all(c("TGFb_bulk", "CSC_bulk") %in% colnames(test_df_wide))) {
    cor_T_C_legacy <- cor.test(
        test_df_wide$TGFb_bulk,
        test_df_wide$CSC_bulk,
        method = "spearman"
    )
    write_validation("H1a: TGFβ -> CSC (legacy CSC_bulk)")
    write_validation(sprintf(
        "  Correlation: r = %.3f, p = %.3f %s",
        cor_T_C_legacy$estimate, cor_T_C_legacy$p.value,
        ifelse(cor_T_C_legacy$p.value < 0.05 & cor_T_C_legacy$estimate > 0.3,
            "✓ SUPPORTED", "✗ NOT SUPPORTED"
        )
    ))
}

if (all(c("TGFb_bulk", "CSC_Yuan101") %in% colnames(test_df_wide))) {
    cor_T_C_yuan <- cor.test(
        test_df_wide$TGFb_bulk,
        test_df_wide$CSC_Yuan101,
        method = "spearman"
    )
    write_validation("\nH1b: TGFβ -> CSC (Yuan CSC 101 chromatin)")
    write_validation(sprintf(
        "  Correlation: r = %.3f, p = %.3f %s",
        cor_T_C_yuan$estimate, cor_T_C_yuan$p.value,
        ifelse(cor_T_C_yuan$p.value < 0.05 & cor_T_C_yuan$estimate > 0.3,
            "✓ SUPPORTED", "✗ NOT SUPPORTED"
        )
    ))
}

# H2: mTOR -> CSC for both definitions
if (all(c("mTOR_bulk", "CSC_bulk") %in% colnames(test_df_wide))) {
    cor_M_C_legacy <- cor.test(
        test_df_wide$mTOR_bulk,
        test_df_wide$CSC_bulk,
        method = "spearman"
    )
    write_validation("\nH2a: mTOR -> CSC (legacy CSC_bulk)")
    write_validation(sprintf(
        "  Correlation: r = %.3f, p = %.3f %s",
        cor_M_C_legacy$estimate, cor_M_C_legacy$p.value,
        ifelse(cor_M_C_legacy$p.value < 0.05 & cor_M_C_legacy$estimate > 0.3,
            "✓ SUPPORTED", "✗ NOT SUPPORTED"
        )
    ))
}

if (all(c("mTOR_bulk", "CSC_Yuan101") %in% colnames(test_df_wide))) {
    cor_M_C_yuan <- cor.test(
        test_df_wide$mTOR_bulk,
        test_df_wide$CSC_Yuan101,
        method = "spearman"
    )
    write_validation("\nH2b: mTOR -> CSC (Yuan CSC 101 chromatin)")
    write_validation(sprintf(
        "  Correlation: r = %.3f, p = %.3f %s",
        cor_M_C_yuan$estimate, cor_M_C_yuan$p.value,
        ifelse(cor_M_C_yuan$p.value < 0.05 & cor_M_C_yuan$estimate > 0.3,
            "✓ SUPPORTED", "✗ NOT SUPPORTED"
        )
    ))
}

# H3: Angio -> TGFb
if (all(c("Angio_bulk", "TGFb_bulk") %in% colnames(test_df_wide))) {
    cor_A_T <- cor.test(
        test_df_wide$Angio_bulk,
        test_df_wide$TGFb_bulk,
        method = "spearman"
    )

    write_validation("\nH3: Angio -> TGFβ pathway")
    write_validation(sprintf(
        "  Correlation: r = %.3f, p = %.3f %s",
        cor_A_T$estimate, cor_A_T$p.value,
        ifelse(cor_A_T$p.value < 0.05 & cor_A_T$estimate > 0.3,
            "✓ SUPPORTED", "✗ NOT SUPPORTED"
        )
    ))
}

# =============================================================================
# STEP 7: PUBLICATION FIGURES (UPGRADED, NO PHEATMAP)
# =============================================================================

message("\n[STEP 7/9] Generating publication figures")

# Theme for plots
theme_atac <- theme_minimal(base_size = 12) +
    theme(
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0),
        legend.position = "right"
    )

# Figure 1: Module scores by condition (barplot with error bars)
scores_summary <- scores_long %>%
    filter(module %in% required_modules) %>%
    group_by(condition, module) %>%
    summarise(
        mean_score = mean(score, na.rm = TRUE),
        se_score = sd(score, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
    ) %>%
    mutate(
        module = factor(module, levels = required_modules),
        condition = factor(condition, levels = c(
            "HFSC", "IFE",
            "Papilloma_mneg", "Papilloma_mpos", "SCC_mneg", "SCC_mpos"
        ))
    )

p1 <- ggplot(scores_summary, aes(x = condition, y = mean_score, fill = condition)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_score - se_score, ymax = mean_score + se_score),
        width = 0.3
    ) +
    facet_wrap(~module, scales = "free_y", ncol = 2) +
    scale_fill_brewer(palette = "Set3") +
    labs(
        title = "ATAC-seq module scores across conditions",
        x = "Condition",
        y = "Module score (z-score)"
    ) +
    theme_atac +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )

ggsave(file.path(FIG_DIR, "atac_module_scores_by_condition.png"),
    p1,
    width = 9, height = 8, dpi = 400
)

# Figure 2: CSC enrichment (m+ vs m- comparison)
scores_mcomp <- scores_long %>%
    filter(module %in% required_modules, stage %in% c("Papilloma", "SCC")) %>%
    mutate(
        module = factor(module, levels = required_modules),
        csc_status = ifelse(csc_enriched, "m+ (CSC-enriched)", "m- (CSC-depleted)")
    )

p2 <- ggplot(scores_mcomp, aes(x = csc_status, y = score, fill = csc_status)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
    facet_wrap(~module, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c(
        "m+ (CSC-enriched)" = "#E74C3C",
        "m- (CSC-depleted)" = "#3498DB"
    )) +
    labs(
        title = "Module scores: CSC-enriched (m+) vs CSC-depleted (m-)",
        x = "",
        y = "Module score (z-score)"
    ) +
    theme_atac +
    theme(legend.position = "top", legend.title = element_blank()) +
    stat_compare_means(method = "t.test", label = "p.format", size = 3)

ggsave(file.path(FIG_DIR, "atac_csc_enrichment_mpos_vs_mneg.png"),
    p2,
    width = 9, height = 8, dpi = 400
)

# Figure 3: Advanced sample-level ComplexHeatmap with row annotations
heatmap_matrix <- scores_wide %>%
    dplyr::select(sample_id, all_of(required_modules)) %>%
    column_to_rownames("sample_id") %>%
    as.matrix()

annotation_row <- samples %>%
    dplyr::select(sample_id, stage, csc_enriched) %>%
    dplyr::mutate(CSC = ifelse(csc_enriched, "m+", "m-")) %>%
    dplyr::select(sample_id, stage, CSC) %>%
    column_to_rownames("sample_id")

if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    requireNamespace("circlize", quietly = TRUE)) {
    message("  Generating ComplexHeatmap for ATAC sample-level modules...")

    # Column-wise scaling (per module)
    heatmap_scaled <- scale(heatmap_matrix)
    rownames(heatmap_scaled) <- rownames(heatmap_matrix)
    colnames(heatmap_scaled) <- colnames(heatmap_matrix)

    ann_df <- annotation_row %>% as.data.frame()

    ann_colors <- list(
        stage = c("Normal" = "#2ECC71", "Papilloma" = "#F39C12", "SCC" = "#E74C3C"),
        CSC   = c("m+" = "#8E44AD", "m-" = "#95A5A6")
    )

    ha_row <- ComplexHeatmap::rowAnnotation(
        df = ann_df,
        col = list(
            stage = ann_colors$stage,
            CSC   = ann_colors$CSC
        ),
        annotation_name_side = "top"
    )

    col_fun_atac <- circlize::colorRamp2(
        c(-2, 0, 2),
        c("blue", "white", "red")
    )

    ht_atac_samples <- ComplexHeatmap::Heatmap(
        heatmap_scaled,
        name = "Z-score",
        col = col_fun_atac,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_title = "ATAC Ras–CSC loop modules",
        row_title = "Samples",
        column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        row_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        left_annotation = ha_row,
        heatmap_legend_param = list(
            title = "Module\nz-score",
            at = c(-2, 0, 2),
            labels = c("Low", "Mid", "High")
        )
    )

    png(file.path(FIG_DIR, "atac_module_complex_heatmap_samples.png"),
        width = 2600, height = 2200, res = 300
    )
    ComplexHeatmap::draw(
        ht_atac_samples,
        heatmap_legend_side = "right",
        annotation_legend_side = "right",
        merge_legend = TRUE
    )
    dev.off()
} else {
    message("  [WARN] ComplexHeatmap or circlize not available; skipping advanced sample heatmap.")
}

# Figure 4: Progression trajectory
progression_plot_data <- scores_long %>%
    filter(module %in% required_modules) %>%
    mutate(
        stage_numeric = case_when(
            stage == "Normal" ~ 1,
            stage == "Papilloma" ~ 2,
            stage == "SCC" ~ 3
        ),
        module_clean = gsub("_bulk", "", module)
    )

p4 <- ggplot(progression_plot_data, aes(x = stage_numeric, y = score, color = module_clean)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "loess", se = TRUE, size = 1.2) +
    scale_x_continuous(breaks = c(1, 2, 3), labels = c("Normal", "Papilloma", "SCC")) +
    scale_color_brewer(palette = "Set1", name = "Module") +
    labs(
        title = "ATAC accessibility progression (chromatin landscape)",
        subtitle = "Module scores across disease stages",
        x = "Disease stage",
        y = "Module score (z-score)"
    ) +
    theme_atac

ggsave(file.path(FIG_DIR, "atac_progression_trajectory.png"),
    p4,
    width = 7, height = 5, dpi = 400
)

# Figure 5: Mechanistic correlations (if RNA-seq available)
if (exists("cross_omics") && nrow(cross_omics) > 0) {
    p5 <- ggplot(cross_omics, aes(x = rna_mean, y = atac_mean, color = rna_condition)) +
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
        facet_wrap(~module, scales = "free", ncol = 2) +
        scale_color_brewer(palette = "Dark2", name = "Condition") +
        labs(
            title = "Cross-omics concordance: RNA-seq vs ATAC-seq",
            x = "RNA-seq module score",
            y = "ATAC-seq module score"
        ) +
        theme_atac

    ggsave(file.path(FIG_DIR, "atac_rna_cross_omics_correlation.png"),
        p5,
        width = 9, height = 8, dpi = 400
    )
}

# Figure 6: ATAC module correlation matrix + RNA vs ATAC loop adjacency (ComplexHeatmap)
if (nrow(scores_wide) >= 3 &&
    requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    requireNamespace("circlize", quietly = TRUE)) {
    # Full ATAC module correlation (required modules)
    cor_matrix_atac <- scores_wide %>%
        dplyr::select(all_of(required_modules)) %>%
        cor(method = "spearman", use = "complete.obs")

    col_fun_corr <- circlize::colorRamp2(
        c(-1, 0, 1),
        c("blue", "white", "red")
    )

    ht_atac_corr <- ComplexHeatmap::Heatmap(
        cor_matrix_atac,
        name = "ATAC\nrho",
        col = col_fun_corr,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_title = "ATAC module correlation matrix (Spearman)",
        row_title = "",
        column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        heatmap_legend_param = list(
            title = "Spearman\nρ",
            at = c(-1, 0, 1),
            labels = c("-1", "0", "1")
        )
    )

    png(file.path(FIG_DIR, "atac_module_correlation_matrix_complex.png"),
        width = 2200, height = 2000, res = 300
    )
    ComplexHeatmap::draw(
        ht_atac_corr,
        heatmap_legend_side = "right",
        merge_legend = TRUE
    )
    dev.off()

    # Loop adjacency comparison (RNA vs ATAC) if RNA available
    if (exists("rna_scores")) {
        modules_for_corr <- intersect(
            c("TGFb_bulk", "Angio_bulk", "CSC_bulk", "mTOR_bulk"),
            colnames(rna_scores)
        )

        # Only proceed if we have the core loop modules in both omics
        if (length(modules_for_corr) >= 3 &&
            all(modules_for_corr %in% colnames(cor_matrix_atac))) {
            # RNA loop adjacency
            rna_corr <- rna_scores %>%
                dplyr::select(all_of(modules_for_corr)) %>%
                cor(method = "spearman", use = "complete.obs")

            # ATAC loop adjacency (subset to same modules)
            atac_corr_loop <- cor_matrix_atac[modules_for_corr, modules_for_corr]

            col_fun_loop <- circlize::colorRamp2(
                c(-1, 0, 1),
                c("blue", "white", "red")
            )

            ht_rna <- ComplexHeatmap::Heatmap(
                rna_corr,
                name = "RNA\nrho",
                col = col_fun_loop,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                column_title = "RNA-seq loop adjacency",
                row_title = "",
                column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
                heatmap_legend_param = list(
                    title = "Spearman\nρ",
                    at = c(-1, 0, 1),
                    labels = c("-1", "0", "1")
                )
            )

            ht_atac_loop <- ComplexHeatmap::Heatmap(
                atac_corr_loop,
                name = "ATAC\nrho",
                col = col_fun_loop,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                column_title = "ATAC loop adjacency",
                row_title = "",
                column_title_gp = grid::gpar(fontsize = 12, fontface = "bold")
            )

            png(file.path(FIG_DIR, "atac_rna_loop_adjacency_complex_heatmap.png"),
                width = 2600, height = 1400, res = 300
            )
            ComplexHeatmap::draw(
                ht_rna + ht_atac_loop,
                heatmap_legend_side = "right",
                merge_legend = TRUE
            )
            dev.off()
        } else {
            message("  [WARN] Not enough overlapping loop modules for RNA vs ATAC adjacency ComplexHeatmap.")
        }
    } else {
        message("  [WARN] rna_scores object missing; skipping RNA vs ATAC adjacency figure.")
    }
} else {
    message("  [WARN] Skipping correlation/adjacency ComplexHeatmap; either n < 3 or ComplexHeatmap/circlize missing.")
}

# =============================================================================
# STEP 8: EXPORT CALIBRATION TARGETS
# =============================================================================

message("\n[STEP 8/9] Exporting calibration targets")

calibration_targets <- scores_long %>%
    filter(module %in% required_modules) %>%
    group_by(condition, stage, module) %>%
    summarise(
        mean_score = mean(score, na.rm = TRUE),
        sd_score = sd(score, na.rm = TRUE),
        n_samples = n(),
        .groups = "drop"
    ) %>%
    pivot_wider(names_from = module, values_from = c(mean_score, sd_score)) %>%
    mutate(
        # New primary C target based on Yuan CSC 101 module
        C_target        = mean_score_CSC_Yuan101,
        # Legacy C target from original CSC_bulk
        C_target_legacy = if ("mean_score_CSC_bulk" %in% colnames(.)) mean_score_CSC_bulk else NA_real_,
        A_target        = mean_score_Angio_bulk,
        T_target        = mean_score_TGFb_bulk,
        M_target        = mean_score_mTOR_bulk
    )

write_csv(
    calibration_targets,
    file.path(PROCESSED_DIR, "atac_calibration_targets.csv")
)

write_validation("\n==============================================")
write_validation("CALIBRATION TARGETS (ATAC)")
write_validation("==============================================\n")

for (i in seq_len(nrow(calibration_targets))) {
    write_validation(sprintf(
        "%s (%s):",
        calibration_targets$condition[i],
        calibration_targets$stage[i]
    ))
    write_validation(sprintf(
        "  C_Yuan101=%.3f, C_legacy=%.3f, A=%.3f, T=%.3f, M=%.3f",
        calibration_targets$C_target[i],
        calibration_targets$C_target_legacy[i],
        calibration_targets$A_target[i],
        calibration_targets$T_target[i],
        calibration_targets$M_target[i]
    ))
}

# =============================================================================
# STEP 9: SUMMARY STATISTICS
# =============================================================================

message("\n[STEP 9/9] Computing summary statistics")

write_validation("\n==============================================")
write_validation("SUMMARY STATISTICS")
write_validation("==============================================\n")

write_validation(sprintf("Total samples processed: %d", nrow(samples)))
write_validation(sprintf("Total genes with accessibility: %d", nrow(access_wide)))
write_validation(sprintf("Modules computed (required set): %d", length(required_modules)))

write_validation("\nGene coverage per required module:")
for (mod in required_modules) {
    genes_in_set <- length(gene_sets[[mod]])
    genes_in_data <- sum(gene_sets[[mod]] %in% access_wide$gene_symbol)
    pct <- 100 * genes_in_data / genes_in_set
    write_validation(sprintf(
        "  %s: %d/%d genes (%.1f%%)",
        mod, genes_in_data, genes_in_set, pct
    ))
}

close(validation_log)

message("\n========================================")
message("ATAC-SEQ PIPELINE COMPLETE")
message("========================================")
message("\nOutputs:")
message("  ", file.path(INTERIM_DIR, "atac_gene_accessibility.tsv"))
message("  ", file.path(PROCESSED_DIR, "atac_module_scores_by_sample.csv"))
message("  ", file.path(PROCESSED_DIR, "atac_calibration_targets.csv"))
message("\nValidation report:")
message("  ", validation_report)
message("\nFigures:")
message("  ", file.path(FIG_DIR, "atac_*.png"))
message("\n✓ READ VALIDATION REPORT FOR STATISTICAL RESULTS\n")
