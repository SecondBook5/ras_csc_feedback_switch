#!/usr/bin/env Rscript

## R/install_packages.R
##
## Install all R + Bioconductor dependencies needed for:
##   - pipelines/data_processing/01_process_bulk_rnaseq.R
##   - pipelines/data_processing/02_process_scrna.R
##   - pipelines/data_processing/03_process_atacseq.R
##
## If you later hit "there is no package called 'X'",
## add X to one of the lists below and re-run this script.

message("=== Ras–CSC pipeline: installing R + Bioconductor dependencies ===")

## ----------------------------------------------------------------------
## 1. Make sure BiocManager exists
## ----------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

## ----------------------------------------------------------------------
## 2. CRAN packages used in the pipelines
## ----------------------------------------------------------------------
cran_pkgs <- c(
  # Core tidyverse-style data handling
  "readr",
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",

  # Config and utilities
  "yaml",
  "jsonlite",

  # Plotting
  "ggplot2",
  "ggpubr",
  "pheatmap",
  "RColorBrewer",
  "patchwork",
  "circlize",

  # Statistical helpers
  "statmod",

  # Matrix / linear algebra (needed for Seurat and friends)
  "Matrix",

  # Seurat is on CRAN as well as phateR
  "Seurat",
  "phateR"
)

## ----------------------------------------------------------------------
## 3. Bioconductor packages used in the pipelines
## ----------------------------------------------------------------------
bioc_pkgs <- c(
  # Bulk RNA-seq / DE
  "DESeq2",
  "genefilter",
  "ComplexHeatmap",

  # Single-cell + matrix infrastructure
  "DelayedMatrixStats",
  "SingleCellExperiment",
  "slingshot",

  # Annotation and Ensembl access
  "biomaRt",

  # Genomic ranges and ATAC handling
  "GenomicRanges",
  "rtracklayer",
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "org.Mm.eg.db",
  "ChIPseeker"
)

## ----------------------------------------------------------------------
## 4. Helpers: install-if-missing
## ----------------------------------------------------------------------
install_cran_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("[CRAN] Installing: %s", pkg))
      tryCatch(
        {
          install.packages(pkg, repos = "https://cloud.r-project.org")
        },
        error = function(e) {
          message(sprintf("[CRAN] ERROR installing %s: %s", pkg, e$message))
        }
      )
    } else {
      message(sprintf("[CRAN] Already installed: %s", pkg))
    }
  }
}

install_bioc_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("[Bioc] Installing: %s", pkg))
      tryCatch(
        {
          BiocManager::install(pkg, ask = FALSE, update = FALSE)
        },
        error = function(e) {
          message(sprintf("[Bioc] ERROR installing %s: %s", pkg, e$message))
        }
      )
    } else {
      message(sprintf("[Bioc] Already installed: %s", pkg))
    }
  }
}

## ----------------------------------------------------------------------
## 5. Run installers
## ----------------------------------------------------------------------
install_cran_if_missing(cran_pkgs)
install_bioc_if_missing(bioc_pkgs)

message("=== Finished installing R dependencies ===")
message("If a script still fails with 'there is no package called ...',")
message("add that package name to cran_pkgs or bioc_pkgs above and re-run.")
