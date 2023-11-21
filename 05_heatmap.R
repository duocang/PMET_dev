suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(hrbrthemes)

  source("scripts/R_utils/process_pmet_result.R")
  source("scripts/R_utils/motif_pair_gene_diagonal.R")
  source("scripts/R_utils/motif_pair_diagonal.R")
  source("scripts/R_utils/motif_pair_plot_hetero.R")
  source("scripts/R_utils/motif_pair_plot_homog.R")
  source("scripts/R_utils/heatmap.R")
})

args <- commandArgs(trailingOnly = TRUE)
# 1. method
# 2. filename
# 3. pmet.out
if (length(args) != 3) {
  stop("You need to provide exactly 3 arguments: method, filename, and pmet.out.")
}
method   <- args[1]
filename <- args[2]
pmet.out <- args[3]

# # method       <- "All"
# method       <- "Overlap"
# filename     <- "results/05_plot/heatmap.png"
# pmet.out     <- "data/motif_output.txt"



heatmap.func(filename = filename,
             method = method,
             p.adj.threshold =0.05,
             p.adj.method = "Adjusted p-value (Bonf)",
             pmet.out = pmet.out,
             draw.histgram = FALSE,
             unique.cmbination = TRUE)
