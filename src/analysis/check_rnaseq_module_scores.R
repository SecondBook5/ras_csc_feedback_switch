#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

scores <- read_csv("data/processed/rnaseq/module_scores_by_sample.csv",
                   show_col_types = FALSE)

# Look at the structure
print(head(scores, 10))

# Condition-level means per dataset
cond_means <- scores %>%
  group_by(dataset, condition) %>%
  summarise(
    across(where(is.numeric), mean, .names = "mean_{.col}"),
    .groups = "drop"
  )

print(cond_means)

# Optional: write out so you can use it directly in model calibration
write_csv(cond_means,
          "data/processed/rnaseq/module_scores_condition_means.csv")
