library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

df <- read_csv("data/processed/model_fits/model_vs_data.csv")

make_panel <- function(var) {
     ggplot(df, aes_string(
          x = paste0(var, "_target"),
          y = paste0(var, "_model"),
          color = "condition",
          shape = "dataset"
     )) +
          geom_point(size = 3) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
          theme_bw() +
          labs(
               title = var,
               x = paste0(var, " (target)"),
               y = paste0(var, " (model)")
          ) +
          coord_equal()
}

pC <- make_panel("C")
pA <- make_panel("A")
pT <- make_panel("T")
pM <- make_panel("M")

out <- (pC | pA) / (pT | pM)

ggsave("figures/main/fig_model_vs_data.png", out, width = 9, height = 9, dpi = 300)
ggsave("figures/main/fig_model_vs_data.pdf", out, width = 9, height = 9)
