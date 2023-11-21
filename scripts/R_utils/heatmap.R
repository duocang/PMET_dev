heatmap.func <- function(filename = NULL,
                         method = NULL,
                         pmet.out = NULL,
                         p.adj.threshold = 0.05,
                         p.adj.method = "Adjusted p-value (Bonf)",
                         draw.histgram = TRUE,
                         unique.cmbination = TRUE)
{
  if (draw.histgram) {
    histgram.path <- filename %>%
    tools::file_path_sans_ext() %>%
    stringr::str_replace("heatmap_", "histogram_")
  } else {
    histgram.path <- NULL
  }

  pmet.result <- data.table::fread(pmet.out,
    select = c(
      "Cluster",
      "Motif 1",
      "Motif 2",
      "Number of genes in cluster with both motifs",
      p.adj.method,
      "Genes"
    ), verbose = FALSE) %>%
    setNames(c("cluster", "motif1", "motif2", "gene_num", "p_adj", "genes")) %>%
    arrange(desc(p_adj)) %>%
    mutate(`motif_pair` = paste0(motif1, "^^", motif2))

  pmet.result.processed <- ProcessPmetResult( pmet_result       = pmet.result,
                                              p_adj_limt        = p.adj.threshold,
                                              gene_portion      = 0.05,
                                              topn              = 20,
                                              histgram_dir      = histgram.path,
                                              unique_cmbination = TRUE)

  if (is.null(pmet.result.processed)) {
    cat("No meaningfull data left after filtering!\n")
    return(NULL)
  }

  results   <- pmet.result.processed
  clusters <- names(results$pmet_result) %>% sort()

  if (method == "Overlap") {
    # after filtering, some clusters may be gone and only one cluster is left
    if (length(clusters) == 1) {

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
                        respective.plot   = TRUE)
      p <- p[[1]]

    } else {
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
    }
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
}
