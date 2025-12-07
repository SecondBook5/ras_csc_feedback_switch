#!/usr/bin/env Rscript

# ===========================================================
# update_gene_sets_from_moesm8.R
#
# Purpose:
#   1) Read Yuan et al. MOESM8 (Fig. 1d GO table), sheet "Fig. 1d".
#   2) Extract gene lists for selected GO terms:
#        - GO:0001525~angiogenesis
#        - GO:0016477~cell migration
#        - GO:0030036~actin cytoskeleton organization
#        - GO:0098609~cell-cell adhesion
#        - GO:0006468~protein phosphorylation
#        - GO:0035556~intracellular signal transduction
#   3) Map these to three bulk modules:
#        - Angio_bulk
#        - CSC_bulk   (migration + actin + adhesion)
#        - mTOR_bulk  (phosphorylation + signal transduction)
#   4) Convert symbols to mouse-style case (PDGFRA -> Pdgfra).
#   5) Expand TGFb_bulk using:
#        - GO terms whose Term contains "TGF"
#        - MSigDB GOMF_TRANSFORMING_GROWTH_FACTOR_BETA_BINDING
#        - a small canonical TGFβ core
#        - any existing TGFb_bulk genes in YAML
#   6) Update config/gene_sets_rnaseq.yaml in place and back up
#      the original YAML to config/gene_sets_rnaseq.yaml.bak.
# ===========================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(yaml)
})

#---------------- paths ------------------------------------

# Path to MOESM8 from Yuan 2022
moesm_path <- "data/raw/yuan2022_moesm/41586_2022_5475_MOESM8_ESM.xlsx"

# Sheet containing Fig. 1d GO table
moesm_sheet <- "Fig. 1d"

# YAML with current gene sets
yaml_path <- "config/gene_sets_rnaseq.yaml"
yaml_backup <- paste0(yaml_path, ".bak")

#---------------- sanity checks ----------------------------

if (!file.exists(moesm_path)) {
  stop("[ERROR] MOESM8 file not found at: ", moesm_path)
}

if (!file.exists(yaml_path)) {
  stop("[ERROR] gene_sets_rnaseq.yaml not found at: ", yaml_path)
}

#---------------- read MOESM8 sheet ------------------------

message("[INFO] Reading MOESM8 sheet '", moesm_sheet, "'")
df_go <- readxl::read_excel(moesm_path, sheet = moesm_sheet)

required_cols <- c("Category", "Term", "Genes")
if (!all(required_cols %in% colnames(df_go))) {
  stop(
    "[ERROR] MOESM8 sheet '", moesm_sheet,
    "' is missing required columns: ",
    paste(setdiff(required_cols, colnames(df_go)), collapse = ", ")
  )
}

# Optionally restrict to BP terms if Category exists
if ("Category" %in% colnames(df_go)) {
  df_go <- df_go %>% filter(str_detect(Category, "GOTERM_BP"))
}

#---------------- module <-> GO mapping --------------------

# We use exact Term strings from the paper
module_go_terms <- list(
  Angio_bulk = c("GO:0001525~angiogenesis"),
  CSC_bulk = c(
    "GO:0016477~cell migration",
    "GO:0030036~actin cytoskeleton organization",
    "GO:0098609~cell-cell adhesion"
  ),
  mTOR_bulk = c(
    "GO:0006468~protein phosphorylation",
    "GO:0035556~intracellular signal transduction"
  )
)

#---------------- helpers ----------------------------------

parse_genes_column <- function(genes_vec) {
  genes_vec %>%
    discard(is.na) %>%
    str_split(",\\s*") %>%
    unlist() %>%
    str_trim() %>%
    discard(~ .x == "") %>%
    unique()
}

to_mouse_case <- function(genes) {
  genes %>%
    str_to_lower() %>%
    str_to_sentence() %>%
    sort() %>%
    unique()
}

# Extract genes for Angio / CSC / mTOR modules
extract_genes_for_module <- function(go_df, go_terms, module_name) {
  sub <- go_df %>%
    filter(Term %in% go_terms)

  if (nrow(sub) == 0L) {
    warning(
      "[WARN] No rows found in MOESM8 for module ", module_name,
      " with GO terms: ", paste(go_terms, collapse = ", ")
    )
    return(character(0))
  }

  genes <- parse_genes_column(sub$Genes)

  if (length(genes) == 0L) {
    warning("[WARN] Genes column empty after parsing for module ", module_name)
    return(character(0))
  }

  genes_mouse <- to_mouse_case(genes)

  message(
    "[INFO] Module ", module_name, ": extracted ",
    length(genes_mouse), " genes from GO terms."
  )

  genes_mouse
}

# Extract TGFb-related genes from MOESM8 based on Term containing "TGF"
extract_tgfb_genes <- function(go_df) {
  sub <- go_df %>%
    filter(str_detect(Term, regex("TGF", ignore_case = TRUE)))

  if (nrow(sub) == 0L) {
    warning("[WARN] No GO terms containing 'TGF' found in MOESM8.")
    return(character(0))
  }

  genes <- parse_genes_column(sub$Genes)

  if (length(genes) == 0L) {
    warning("[WARN] Genes column empty after parsing TGF-related GO terms.")
    return(character(0))
  }

  genes_mouse <- to_mouse_case(genes)

  message(
    "[INFO] TGFb_bulk: extracted ",
    length(genes_mouse), " genes from TGF-related GO terms."
  )

  genes_mouse
}

#---------------- build new module gene sets ---------------

new_sets <- list()

for (mod in names(module_go_terms)) {
  go_terms <- module_go_terms[[mod]]
  new_sets[[mod]] <- extract_genes_for_module(df_go, go_terms, mod)
}

tgfb_from_go <- extract_tgfb_genes(df_go)

#---------------- load existing YAML -----------------------

message("[INFO] Reading existing YAML: ", yaml_path)
gene_sets <- yaml::read_yaml(yaml_path)

if (!is.list(gene_sets)) {
  stop("[ERROR] gene_sets_rnaseq.yaml did not parse as a list.")
}

if (is.null(gene_sets$TGFb_bulk)) {
  warning("[WARN] TGFb_bulk not present in YAML; it will be created.")
  gene_sets$TGFb_bulk <- character(0)
}

#---------------- backup YAML -----------------------------

if (!file.exists(yaml_backup)) {
  ok <- file.copy(yaml_path, yaml_backup, overwrite = FALSE)
  if (ok) {
    message("[INFO] Backed up original YAML to: ", yaml_backup)
  } else {
    warning("[WARN] Failed to create YAML backup at: ", yaml_backup)
  }
} else {
  message("[INFO] Backup YAML already exists at: ", yaml_backup)
}

#---------------- update Angio / CSC / mTOR  ---------------

for (mod in names(new_sets)) {
  genes_mod <- new_sets[[mod]]

  if (length(genes_mod) == 0L) {
    warning(
      "[WARN] Skipping update for module ", mod,
      " because no genes were extracted."
    )
    next
  }

  gene_sets[[mod]] <- genes_mod
  message(
    "[INFO] Updated module ", mod,
    " with ", length(genes_mod), " genes from MOESM8."
  )
}

#---------------- expand TGFb_bulk -------------------------

# MSigDB GOMF_TRANSFORMING_GROWTH_FACTOR_BETA_BINDING gene set (mouse)
tgfb_msigdb_core <- c(
  "Acvr1", "Acvr2b", "Acvrl1", "Agrn", "Cd109", "Cd36", "Chrdl1", "Eng",
  "Hyal2", "Itgav", "Lrrc32", "Ltbp1", "Ltbp3", "Ltbp4", "Nrros", "Tgfb3",
  "Tgfbr1", "Tgfbr2", "Tgfbr3", "Tgfbr3l", "Thbs1", "Tsku", "Twsg1",
  "Vasn", "Wfikkn1", "Wfikkn2"
)

# Manual canonical TGFβ core
tgfb_manual_core <- c(
  "Tgfb1", "Tgfb2", "Tgfb3",
  "Smad2", "Smad3", "Smad4", "Smad7",
  "Tgfbr1", "Tgfbr2",
  "Acvr1",
  "Serpine1", # PAI-1
  "Ctgf",
  "Id1"
)

# Union:
#  - existing YAML TGFb_bulk
#  - TGF-related GO genes from MOESM8
#  - MSigDB TGF-beta binding genes
#  - manual canonical core
tgfb_new <- unique(c(
  gene_sets$TGFb_bulk,
  tgfb_from_go,
  tgfb_msigdb_core,
  tgfb_manual_core
))

tgfb_new <- sort(unique(tgfb_new))
gene_sets$TGFb_bulk <- tgfb_new

message(
  "[INFO] Updated TGFb_bulk with ",
  length(tgfb_new),
  " genes (YAML existing + TGF GO + MSigDB MF + core)."
)

#---------------- write updated YAML -----------------------

yaml::write_yaml(gene_sets, yaml_path)
message("[DONE] Updated gene_sets_rnaseq.yaml written to: ", yaml_path)
