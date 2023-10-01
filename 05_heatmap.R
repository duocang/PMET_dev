suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  
  source("scripts/R_utils/process_pmet_result.R")
  source("scripts/R_utils/motif_pair_gene_diagonal.R")
  source("scripts/R_utils/motif_pair_diagonal.R")
  source("scripts/R_utils/motif_pair_plot_hetero.R")
  source("scripts/R_utils/motif_pair_plot_homog.R")
})

args <- commandArgs(trailingOnly = TRUE)
# 1. method
# 2. filename
# 3. pmet.out.dir
if (length(args) != 3) {
  stop("You need to provide exactly 3 arguments: method, filename, and pmet.out.dir.")
}

method       <- args[1]
filename     <- args[2]
pmet.out.dir <- args[3]


# method       <- "Overlap"
# # method       <- "All"
# filename     <- "results/05_plot/heatmap_old_fimo_plus_pmetindex.png"
# pmet.out.dir <- "results/04_new_fimo_vs_old_fimo_plus_pmetindex/old_fimo_plus_pmetindex.txt"

pmet.result <- data.table::fread(pmet.out.dir,
  select = c(
    "Cluster", "Motif 1", "Motif 2",
    "Number of genes in cluster with both motifs",
    "Adjusted p-value (BH)", "Genes"
  ), verbose = FALSE) %>%
  setNames(c("cluster", "motif1", "motif2", "gene_num", "p_adj", "genes")) %>%
  # dplyr::filter(gene_num > 0) %>%
  arrange(desc(p_adj)) %>%
  mutate(`motif_pair` = paste0(motif1, "^^", motif2))



pmet.result.processed <- ProcessPmetResult( pmet_result       = pmet.result,
                                            p_adj_limt        = 0.05,
                                            gene_portion      = 0.05,
                                            topn              = 20,
                                            unique_cmbination = TRUE)



results   <- pmet.result.processed
clusters <- names(results$pmet_result) %>% sort()

if (method == "Overlap") {
  
  motifs <- TopMotifsGenerator(pmet.result.processed$motifs, by.cluster = FALSE, exclusive.motifs = TRUE)
  num.motifs <- length(motifs)
  
  # expend ggplot with genes for hover information
  dat_list   <- MotifPairGeneDiagonal(pmet.result.processed$pmet_result, motifs, counts = "p_adj")
  clusters   <- names(dat_list) %>% sort()
  
  # merge data into DF[[1]]
  dat <- dat_list[[1]]
  # move all non-NA values from other DFs to DF[[1]]
  for (i in 2:length(dat_list)) {
    indx                 <- which(!is.na(dat_list[[i]][, "cluster"]))
    dat[indx,          ] <- dat_list[[i]][indx, ]
    dat[indx, "cluster"] <- names(dat_list)[i]
  }
  
  p <- MotifPairPlotHetero(dat,  "p_adj", motifs, clusters)
} else {
  if (method == "All") {
    respective.plot <- FALSE
  } else if (method %in% clusters) {
    respective.plot <- TRUE
  }
  
  p <- MotifPairPlotHomog(results$pmet_result,
                          results$motifs,
                          counts            = "p_adj",
                          exclusive.motifs  = TRUE,
                          by.cluster        = FALSE,
                          show.cluster      = FALSE,
                          legend.title      = "-log10(p.adj)",
                          nrow_             = ceiling(length(clusters)/2),
                          ncol_             = 2,
                          axis.lables       = "",
                          show.axis.text    = TRUE,
                          diff.colors       = TRUE,
                          respective.plot   = respective.plot
  )
  
  if (method %in% clusters) {
    p <- p[[method]]
  }
}

# set size of saved plot
if (method == "Overlap") {
  wid <- 20
  hei <- 20
} else {
  if (method == "All") {
    wid <- 20
    hei <- 10 * ceiling(length(clusters)/2)
  } else if (method %in% clusters) {
    wid <- 20
    hei <- 20
  }
}
ggsave(filename, p, width = wid, height = hei, dpi = 320, units = "in")
