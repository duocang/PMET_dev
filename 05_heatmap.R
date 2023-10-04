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

# args <- commandArgs(trailingOnly = TRUE)
# # 1. method
# # 2. filename
# # 3. pmet.out
# if (length(args) != 3) {
#   stop("You need to provide exactly 3 arguments: method, filename, and pmet.out.")
# }
# method   <- args[1]
# filename <- args[2]
# pmet.out <- args[3]


method       <- "Overlap"
# method       <- "All"
filename     <- "results/05_plot/heatmap_old_fimo_plus_pmetindex.png"
pmet.out <- "results/04_new_fimo_vs_old_fimo_plus_pmetindex/old_fimo_plus_pmetindex.txt"



heatmap.func(filename, method, pmet.out, TRUE)
