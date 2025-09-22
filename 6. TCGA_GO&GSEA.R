# ------------------- Libraries -------------------
library(TCGAplot)
library(openxlsx)
library(ggplot2)
library(gridExtra)

# ------------------- Gene set -------------------
geneset <- c("SIM2", "OTX1", "MYBL2", "HJURP", "IQGAP3", "KIF18B",
             "TPX2", "RRM2", "SKA3", "RHPN1", "SLC7A11", "TNS4", 
             "FOLR1", "MYL9", "ABCG2", "SORBS1", "NRG1", "CRYAB", 
             "PYGM", "SYNM","CAV1", "MASP1", "FHL1", "RXRG", "TCEAL2", "DES")

# ------------------- Cancer types -------------------
cancer_types <- c("BRCA", "PRAD", "COAD", "LUAD")

# ------------------- Function to collect plots + results -------------------
save_cancer_plots <- function(cancer, genes, method = c("GO", "KEGG"), ncol = 5) {
  method <- match.arg(method)
  plot_list <- list()
  results_list <- list()
  
  for (gene in genes) {
    res <- tryCatch({
      if (method == "GO") {
        gene_gsea_go(cancer, gene)
      } else {
        gene_gsea_kegg(cancer, gene)
      }
    }, error = function(e) NULL)
    
    # store plot
    if (!is.null(res)) {
      if ("gg" %in% class(res) | inherits(res, "ggplot")) {
        # case 1: only a plot returned
        plot_list[[paste0(cancer, "_", gene)]] <- res
      } else if (is.list(res)) {
        # case 2: result is list with plot + table
        if (!is.null(res$plot)) {
          plot_list[[paste0(cancer, "_", gene)]] <- res$plot
        }
        if (!is.null(res$table)) {
          df <- res$table
          df$cancer <- cancer
          df$gene <- gene
          results_list[[paste0(cancer, "_", gene)]] <- df
        }
      }
    }
  }
  
  # save plots in PDF
  if (length(plot_list) > 0) {
    nrow <- ceiling(length(plot_list) / ncol)
    pdf(paste0(cancer, "_", method, "_GSEA.pdf"), width = 16, height = 9)
    do.call(grid.arrange, c(plot_list, ncol = ncol, nrow = nrow))
    dev.off()
  }
  
  # save results in CSV
  if (length(results_list) > 0) {
    all_res <- do.call(rbind, results_list)
    write.csv(all_res, paste0(cancer, "_", method, "_GSEA.csv"), row.names = FALSE)
  }
}

# ------------------- Run for all cancers -------------------
for (cancer in cancer_types) {
  # Save GO plots & results
  save_cancer_plots(cancer, geneset, method = "GO", ncol = 5)
  
  # Save KEGG plots & results
  save_cancer_plots(cancer, geneset, method = "KEGG", ncol = 5)
}
