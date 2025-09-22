#  LINK::: https://github.com/tjhwangxiong/TCGAplot?tab=readme-ov-file
###  LINK-2 : ONLINE ::: cSurvival (https://tau.cmmt.ubc.ca/cSurvival)

# NOTE :: # Geneset based analysis--> GO to line 84

library(TCGAplot)
#-----------------------------------------------------------------------------------------
## YOUR QUERY GENE/s
gene <- "DNAJC15"

#------------------------------------------PART-1------------------------------------------------
# Create a pan-cancer box plot for a single gene with symbols indicating statistical significance.

pan_boxplot(gene,palette="lancet",legend="right",method="wilcox.test") # or direct ---> pan_boxplot("KLF7") :::: the color palette to be used for filling by groups. e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".

# Create a pan-cancer paired box plot for a single gene with symbols indicating statistical significance.
pan_paired_boxplot(gene,legend = "none")
pan_paired_boxplot(gene,palette="lancet",legend="right",method="wilcox.test")


#Create a pan-cancer box plot for a single gene in tumor samples.
pan_tumor_boxplot(gene)

# Pan-cancer correlation analysis:: Create a pan-cancer radar chart for gene expression and TMB correlation.
gene_TMB_radar(gene,method = "pearson")
#Pan-cancer gene expression and MSI correlation radar chart
gene_MSI_radar(gene,method = "pearson")

# Create a pan-cancer heatmap with symbols indicating statistical significance to reveal the correlation between the expression of 
# a single gene and ICGs (immune checkpoint genes).
gene_checkpoint_heatmap(gene,method="pearson",lowcol="blue",highcol="red",cluster_row=T,cluster_col=T,legend=T)

##### IMPORTANT
## Pan-cancer gene expression and chemokine correlation
gene_chemokine_heatmap(gene,method="pearson",lowcol="blue",highcol="red",cluster_row=T,cluster_col=T,legend=T)


## Pan-cancer gene expression and chemokine receptor correlation
gene_receptor_heatmap(gene,method="pearson",lowcol="blue",highcol="red",cluster_row=T,cluster_col=T,legend=T)

## Pan-cancer gene expression and immune stimulator correlation
gene_immustimulator_heatmap(gene,method="pearson",lowcol="blue",highcol="red",cluster_row=T,cluster_col=T,legend=T)

## Pan-cancer gene expression and immune inhibitor correlation
gene_immuinhibitor_heatmap(gene,method="pearson",lowcol="blue",highcol="red",cluster_row=T,cluster_col=T,legend=T)

## Pan-cancer gene expression and immune infiltration correlation
gene_immucell_heatmap(gene,method="pearson",lowcol="blue",highcol="red",cluster_row=T,cluster_col=T,legend=T)

## Pan-cancer gene expression and immune score correlation
gene_immunescore_heatmap(gene,method="pearson",lowcol="blue",highcol="red",cluster_row=T,cluster_col=T,legend=T)

## Create a pan-cancer triangle reveals the correlation between the expression of a single gene and immune scores, including Stromal score, immune score, and ESTIMATE score.
gene_immunescore_triangle(gene,method="pearson")

## Create a pan-cancer Cox regression forest plot for a specific gene.
pan_forest(gene,adjust=F)
#--------------------------------------------------------PART-1.1: Above same Analysis with automation ------------
# Create directory if it doesn't exist
if (!dir.exists("0.Expression_Analysis")) {
  dir.create("0.Expression_Analysis")
}

# Function to save plot with consistent dimensions
save_plot <- function(plot_func, gene, filename, ...) {
  pdf(file = file.path("0.Expression_Analysis", filename), 
      width = 15, height = 6)
  print(plot_func(gene, ...))
  dev.off()
}

# Save all plots
gene <- "DNAJC15"  # Replace with your gene of interest

# Boxplots
save_plot(pan_boxplot, gene, "pan_boxplot.pdf", palette="lancet", legend="right", method="wilcox.test")
save_plot(pan_paired_boxplot, gene, "pan_paired_boxplot.pdf", palette="lancet", legend="right", method="wilcox.test")
save_plot(pan_tumor_boxplot, gene, "pan_tumor_boxplot.pdf")

# Radar charts
save_plot(gene_TMB_radar, gene, "gene_TMB_radar.pdf", method = "pearson")
save_plot(gene_MSI_radar, gene, "gene_MSI_radar.pdf", method = "pearson")

# Heatmaps
save_plot(gene_checkpoint_heatmap, gene, "gene_checkpoint_heatmap.pdf", 
          method="pearson", lowcol="blue", highcol="red", cluster_row=T, cluster_col=T, legend=T)
save_plot(gene_chemokine_heatmap, gene, "gene_chemokine_heatmap.pdf", 
          method="pearson", lowcol="blue", highcol="red", cluster_row=T, cluster_col=T, legend=T)
save_plot(gene_receptor_heatmap, gene, "gene_receptor_heatmap.pdf", 
          method="pearson", lowcol="blue", highcol="red", cluster_row=T, cluster_col=T, legend=T)
save_plot(gene_immustimulator_heatmap, gene, "gene_immustimulator_heatmap.pdf", 
          method="pearson", lowcol="blue", highcol="red", cluster_row=T, cluster_col=T, legend=T)
save_plot(gene_immuinhibitor_heatmap, gene, "gene_immuinhibitor_heatmap.pdf", 
          method="pearson", lowcol="blue", highcol="red", cluster_row=T, cluster_col=T, legend=T)
save_plot(gene_immucell_heatmap, gene, "gene_immucell_heatmap.pdf", 
          method="pearson", lowcol="blue", highcol="red", cluster_row=T, cluster_col=T, legend=T)
save_plot(gene_immunescore_heatmap, gene, "gene_immunescore_heatmap.pdf", 
          method="pearson", lowcol="blue", highcol="red", cluster_row=T, cluster_col=T, legend=T)

# Triangle plot
save_plot(gene_immunescore_triangle, gene, "gene_immunescore_triangle.pdf", method="pearson")

# Forest plot
save_plot(pan_forest, gene, "pan_forest.pdf", adjust=F)


#-------------------------------------------------------PART-2------------------------------------------------------------------
####### Cancer type specific analysis #################

# Create a tumor-normal box plot for a single gene with symbols indicating statistical significance in a specific type of cancer.
# # cancer name likes "BRCA
# cancer <-"DNAJC15"
# tcga_boxplot(cancer,gene, add = "jitter",palette="jco",legend="none",label="p.signif",method="wilcox.test")
# 
#----------------------------------------------------------------------------
# ## for Multiple CAncers
#_---------------------------------------------------------------------------

# List of TCGA cancer types you want to analyze ( It is better to consider only signifiacnt one from boxplot)
## ALL CANCER TYPES
# cancer_types <- c(
#   "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC",
#   "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV",
#   "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
#   "UCEC", "UCS", "UVM"
# )

## Significant Cancer types as per Box plot result
cancer_types <- c("UCEC", "STAD", "READ", "PRAD", 
                  "LUAD", "LUSC", "KIRC","KICH", 
                  "HNSC", "COAD", "CHOL")

# Your gene of interest
gene <- "DNAJC15"  # Replace with your gene of interest

#---------------------------------------
## SINGLE GENE BOXPLOTS
#---------------------------------------
# Create a directory for the plots (if it doesn't exist)
output_dir <- "1.Boxplots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Your gene of interest
gene <- "DNAJC15"  # Replace with your gene of interest

# Create PDFs for each cancer type
for (cancer in cancer_types) {
  # Create PDF file path with directory
  pdf_file <- file.path(output_dir, paste0(cancer, "_boxplot.pdf"))
  
  # Create PDF
  pdf(file = pdf_file, width = 4, height = 4)
  
  # Generate the boxplot
  print(tcga_boxplot(cancer, gene, 
                     add = "jitter", 
                     palette = "lancet", 
                     legend = "none", 
                     label = "p.signif", 
                     method = "wilcox.test"))
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  message("Boxplot for ", cancer, " saved to ", pdf_file)
}

#---------------------------------------
## PAIRED BOXPLOTS
#---------------------------------------
output_dir <- "2.Paired_boxplot"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Create PDFs for each cancer type
for (cancer in cancer_types) {
  # Create PDF file path with directory
  pdf_file <- file.path(output_dir, paste0(cancer, "_Paired_boxplot.pdf"))
  
  # Create PDF
  pdf(file = pdf_file, width = 6, height = 6)
  
  # Generate the boxplot
  print(paired_boxplot(cancer, gene,
                       palette = "jco",
                       legend = "none",
                       label = "p.signif",
                       method = "wilcox.test"))
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  message("Paired boxplot for ", cancer, " saved to ", pdf_file)
}

#---------------------------------------
## STAGEWISE BOXPLOTS
#---------------------------------------
output_dir <- "3.Stagewise_boxplot"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Create PDFs for each cancer type
for (cancer in cancer_types) {
  # Create PDF file path with directory
  pdf_file <- file.path(output_dir, paste0(cancer, "_Stagewise_boxplot.pdf"))
  
  # Create PDF
  pdf(file = pdf_file, width = 6, height = 6)
  
  # Generate the boxplot
  print(gene_stage(cancer, gene,
                   add = "jitter",
                   palette = "jco",
                   legend = "none",
                   label = "p.signif",
                   method = "wilcox.test"))
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  message("Stagewise boxplot for ", cancer, " saved to ", pdf_file)
}
#---------------------------------------
## HEATMAPS
#---------------------------------------
## Expression analysis grouped by the expression of a spcecific gene
## Create a heatmap for differentially expressed genes grouped by the expression of a single gene in a specific type of cancer.
## checking how other genes behave in samples with high vs low expression of a single gene (like TP53) within a certain cancer type (like breast cancer).
# also produce full list of DEGS in csv

# Create a directory for the plots (if it doesn't exist)
output_dir <- "4.Heatmaps"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

gene <- "DNAJC15"  # Replace with your gene of interest

# Create PDFs for each cancer type
for (cancer in cancer_types) {
  # Create PDF file path with directory
  pdf_file <- file.path(output_dir, paste0(cancer, "_Heatmap.pdf"))
  
  # Create PDF
  pdf(file = pdf_file, width = 10, height = 6)
  
  # Generate the heatmap with error handling
  tryCatch({
    print(gene_deg_heatmap(cancer, gene, top_n = 20))
  }, error = function(e) {
    message("Error generating heatmap for ", cancer, ": ", e$message)
  })
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  message("Heatmap for ", cancer, " saved to ", pdf_file)
}
#------------------------------------------PART-3---------------------------------------------------------------
                              # GSEA-GO, GSEA_KEGG
##----------------------------------------------------------------------------
# Tumor samples are divided into DNAJC15-high and DNAJC15-low groups based on a predefined
# finally, we analyze what biological processes those genes are involved in using GSEA-GO.”

gene <- "DNAJC15"  # Replace with your gene of interest  # gene_gsea_go(cancer,gene) ## single cancer/gene -> gene_gsea_go(cancer,gene)

# Create output directory (if it doesn't exist)
output_dir <- "5.GSEA_mapGO"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)  # Added recursive=TRUE for nested directories
}

# Create PDFs for each cancer type
for (cancer in cancer_types) {
  # Create proper file path using file.path()
  pdf_file <- file.path(output_dir, paste0(cancer, "_GSEA_GO.pdf"))
  
  # Create PDF with appropriate dimensions
  pdf(file = pdf_file, width = 10, height = 8)  # Increased dimensions for better visualization
  
  # Generate the plot with error handling
  tryCatch({
    print(gene_gsea_go(cancer, gene))
    
    # Print success message only if no error
    message("Successfully created GSEA GO plot for ", cancer, " saved to ", pdf_file)
  }, error = function(e) {
    message("ERROR generating GSEA GO plot for ", cancer, ": ", e$message)
    dev.off()  # Ensure device is closed if error occurs
    file.remove(pdf_file)  # Remove empty/incomplete PDF
  })
  
  # Close the PDF device (if not already closed by error)
  if (dev.cur() > 1) dev.off()
}

# Final summary message
message("\nGSEA GO analysis completed for ", length(cancer_types), " cancer types")

#---------------------------------------------------------------------
## GSEA-KEGG grouped by the expression of a spcecific gene
#---------------------------------------------------------------------
# Gene of interest
gene <- "DNAJC15"  # Replace with your gene of interest

# Create output directory (if it doesn't exist)
output_dir <- "6.GSEA_mapKegg"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)  # Added recursive=TRUE for nested directories
}

# Create PDFs for each cancer type
for (cancer in cancer_types) {
  # Create proper file path using file.path()
  pdf_file <- file.path(output_dir, paste0(cancer, "_GSEA_KEGG.pdf"))
  
  # Create PDF with appropriate dimensions
  pdf(file = pdf_file, width = 10, height = 8)  # Increased dimensions for better visualization
  
  # Generate the plot with error handling
  tryCatch({
    print(gene_gsea_kegg(cancer, gene))
    
    # Print success message only if no error
    message("Successfully created GSEA KEGG plot for ", cancer, " saved to ", pdf_file)
  }, error = function(e) {
    message("ERROR generating GSEA KEGG plot for ", cancer, ": ", e$message)
    dev.off()  # Ensure device is closed if error occurs
    file.remove(pdf_file)  # Remove empty/incomplete PDF
  })
  
  # Close the PDF device (if not already closed by error)
  if (dev.cur() > 1) dev.off()
}

# Final summary message
message("\nGSEA KEGG analysis completed for ", length(cancer_types), " cancer types")


#-----------------------------------------PART-4: Diagnostic ROC Curve------------------------------------------

# A function to draw the ROC curve and calculate the AUC of a diagnostic model using the expression of a single gene in a specific type of cancer.

gene <- "DNAJC15"  # Replace with your gene of interest

# Create a directory for the plots (if it doesn't exist)
output_dir <- "7.ROC_Curves"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Create PDFs for each cancer type
for (cancer in cancer_types) {
  # Create PDF file path with directory
  pdf_file <- file.path(output_dir, paste0(cancer, "_ROC_Curve.pdf"))
  
  # Create PDF
  pdf(file = pdf_file, width = 5, height = 5)
  
  # Generate the ROC curve
  tryCatch({
    print(tcga_roc(cancer, gene))
  }, error = function(e) {
    message("Error generating ROC curve for ", cancer, ": ", e$message)
  })
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  message("ROC curve for ", cancer, " saved to ", pdf_file)
}

##--------------------------------------PART-5 :: Cancer type specific correlation analysis -------------------------------------

# # Gene-gene correlation scatter
# # Scatter plot of gene and gene correlation in a specific type cancer.
# gene_gene_scatter(cancer,gene1,gene2,density="F") ##  gene_gene_scatter("BLCA","CBX2","CBX3") OR gene_gene_scatter("BLCA","CBX2","CBX3",density="T")

gene <- "DNAJC15"  # Replace with your gene of interest

# Create a directory for the plots (if it doesn't exist)
output_dir <- "8.Gene-Gene_CORR"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Create PDFs for each cancer type
for (cancer in cancer_types) {
  # Create PDF file path with directory
  pdf_file <- file.path(output_dir, paste0(cancer, "_CORR.pdf"))
  
  # Create PDF
  pdf(file = pdf_file, width = 5, height = 5)
  
  # Generate the ROC curve
  tryCatch({
    print(gene_gene_scatter(cancer, gene, "TP53", density = "T"))
  }, error = function(e) {
    message("Error generating Corr.plot for ", cancer, ": ", e$message)
  })
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  message("Corr.plot for ", cancer, " saved to ", pdf_file)
}


##--------------------------------------PART-6 :: Gene-promoter methylation correlation scatter -------------------------------------


gene <- "DNAJC15"  # Replace with your gene of interest

# Create a directory for the plots (if it doesn't exist)
output_dir <- "9.Gene_Methylation_Corr"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Create PDFs for each cancer type
for (cancer in cancer_types) {
  # Create PDF file path with directory
  pdf_file <- file.path(output_dir, paste0(cancer, "_gmCORR.pdf"))
  
  # Create PDF
  pdf(file = pdf_file, width = 5, height = 5)
  
  # Generate the ROC curve
  tryCatch({
    print(gene_methylation_scatter(cancer,gene))
  }, error = function(e) {
    message("Error generating Gene_Methylation_Corr for ", cancer, ": ", e$message)
  })
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  message("Gene_Methylation_Corr plot for ", cancer, " saved to ", pdf_file)
}

## It can give message if gene is not exist --> "This gene is not included in the methylation matrix!"

#-----------------------------------PART-7 :Expression heatmap of significantly correlated genes and GO analysis----------------------------

# Heatmap and Go enrichment of the positive and negative co-expressed genes of a single gene in a specific type of cancer.
# Example = gene_coexp_heatmap("STAD","KLF7")

gene <- "DNAJC15"  # Replace with your gene of interest

# Create a directory for the plots (if it doesn't exist)
output_dir <- "10.Gene_Coexp_Heatmap_BP"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Create PDFs for each cancer type
for (cancer in cancer_types) {
  # Create PDF file path with directory
  pdf_file <- file.path(output_dir, paste0(cancer, "_genCORRHM.pdf"))
  
  # Create PDF
  pdf(file = pdf_file, width = 15, height = 15)
  
  # Generate the ROC curve
  tryCatch({
    print(gene_coexp_heatmap(cancer,gene,top_n=20, method="pearson"))
  }, error = function(e) {
    message("Error generating Gene_Coexp_Heatmap_BP for ", cancer, ": ", e$message)
  })
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  message("GGene_Coexp_Heatmap_BP plot for ", cancer, " saved to ", pdf_file)
}

#--------------------------PART-8: Survavial analysis based on the expression of a single gene----------------------
# K_M survival plot for a single gene in a specific type of cancer.
# Example :: tcga_kmplot("COAD","KLF7")


gene <- "DNAJC15"  # Replace with your gene of interest

# Create a directory for the plots (if it doesn't exist)
output_dir <- "11.K_M_survival_plot"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Create PDFs for each cancer type
for (cancer in cancer_types) {
  # Create PDF file path with directory
  pdf_file <- file.path(output_dir, paste0(cancer, "_km_surv.pdf"))
  
  # Create PDF
  pdf(file = pdf_file, width = 6, height = 6)
  
  # Generate the ROC curve
  tryCatch({
    print(tcga_kmplot(cancer,gene,palette='lancet'))
  }, error = function(e) {
    message("Error generating K_M_survival_plot for ", cancer, ": ", e$message)
  })
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  message("K_M_survival_plot for ", cancer, " saved to ", pdf_file)
}



##----------------------------------------------------------------------------------------------------------------
                                           # Geneset based analysis
# ##-------------------------------------------------------------------------------------------------------------
# Both geneset listed in MSigDB and user defined geneset in the form of character vector were supported to perform geneset based pan-cancer and cancer type specific analysis. get_geneset() function could extract the whole built in geneset list from MSigDB.
## Create a pan-cancer box plot for a geneset with symbols indicating statistical significance.

# Define gene sets
geneset_U <- c("SIM2", "OTX1", "MYBL2", "HJURP", "IQGAP3", "KIF18B",
               "TPX2", "RRM2", "SKA3", "RHPN1", "SLC7A11")

geneset_D <- c("MYL9", "ABCG2", "SORBS1", "NRG1", "CRYAB", "PYGM", "SYNM",
               "CAV1", "MASP1", "FHL1", "RXRG", "TCEAL2", "DES")

geneset_alias <- "Key Signature Genes"
#-----------------------------------------------------------------------------------------
# 1. Chemokine heatmaps - Separate PDFs for Up and Down
pdf("1.Chemokine_Heatmap_UPREGULATED1.pdf", width = 10, height = 10)
print("Generating chemokine heatmap for UPREGULATED genes...")
gs_chemokine_heatmap(geneset_U, method = "pearson", lowcol = "blue", highcol = "red", 
                     cluster_row = F, cluster_col = F, legend = T)
dev.off()

pdf("1.Chemokine_Heatmap_DOWNREGULATED1.pdf", width = 10, height = 10)
print("Generating chemokine heatmap for DOWNREGULATED genes...")
gs_chemokine_heatmap(geneset_D, method = "pearson", lowcol = "blue", highcol = "red", 
                     cluster_row = F, cluster_col = F, legend = T)
dev.off()
print("Chemokine heatmaps saved to separate PDF files")
#-----------------------------------------------------------------------------------------
# 2. Receptor heatmaps - Separate PDFs for Up and Down
pdf("2.Receptor_Heatmap_UPREGULATED.pdf", width = 10, height = 10)
print("Generating receptor heatmap for UPREGULATED genes...")
gs_receptor_heatmap(geneset_U, method = "pearson", lowcol = "blue", highcol = "red", 
                    cluster_row = F, cluster_col = F, legend = T)
dev.off()

pdf("2.Receptor_Heatmap_DOWNREGULATED.pdf", width = 10, height = 10)
print("Generating receptor heatmap for DOWNREGULATED genes...")
gs_receptor_heatmap(geneset_D, method = "pearson", lowcol = "blue", highcol = "red", 
                    cluster_row = F, cluster_col = F, legend = T)
dev.off()
print("Receptor heatmaps saved to separate PDF files")
#-----------------------------------------------------------------------------------------
# 3. Immune stimulator heatmaps - Separate PDFs for Up and Down
pdf("3.ImmuneStimulator_Heatmap_UPREGULATED.pdf", width = 10, height = 10)
print("Generating immune stimulator heatmap for UPREGULATED genes...")
gs_immustimulator_heatmap(geneset_U, method = "pearson", lowcol = "blue", highcol = "red", 
                          cluster_row = F, cluster_col = F, legend = T)
dev.off()

pdf("3.ImmuneStimulator_Heatmap_DOWNREGULATED.pdf", width = 10, height = 10)
print("Generating immune stimulator heatmap for DOWNREGULATED genes...")
gs_immustimulator_heatmap(geneset_D, method = "pearson", lowcol = "blue", highcol = "red", 
                          cluster_row = F, cluster_col = F, legend = T)
dev.off()
print("Immune stimulator heatmaps saved to separate PDF files")
#-----------------------------------------------------------------------------------------
# 4. Immune inhibitor heatmaps - Separate PDFs for Up and Down
pdf("4.ImmuneInhibitor_Heatmap_UPREGULATED.pdf", width = 10, height = 10)
print("Generating immune inhibitor heatmap for UPREGULATED genes...")
gs_immuinhibitor_heatmap(geneset_U, method = "pearson", lowcol = "blue", highcol = "red", 
                         cluster_row = F, cluster_col = F, legend = T)
dev.off()

pdf("4.ImmuneInhibitor_Heatmap_DOWNREGULATED.pdf", width = 10, height = 10)
print("Generating immune inhibitor heatmap for DOWNREGULATED genes...")
gs_immuinhibitor_heatmap(geneset_D, method = "pearson", lowcol = "blue", highcol = "red", 
                         cluster_row = F, cluster_col = F, legend = T)
dev.off()
print("Immune inhibitor heatmaps saved to separate PDF files")
#-----------------------------------------------------------------------------------------
# 5. Immune score heatmaps - Separate PDFs for Up and Down
pdf("5.ImmuneScore_Heatmap_UPREGULATED.pdf", width = 10, height = 10)
print("Generating immune score heatmap for UPREGULATED genes...")
gs_immunescore_heatmap(geneset_U, method = "pearson", lowcol = "blue", highcol = "red", 
                       cluster_row = F, cluster_col = F, legend = T)
dev.off()

pdf("5.ImmuneScore_Heatmap_DOWNREGULATED.pdf", width = 10, height = 10)
print("Generating immune score heatmap for DOWNREGULATED genes...")
gs_immunescore_heatmap(geneset_D, method = "pearson", lowcol = "blue", highcol = "red", 
                       cluster_row = F, cluster_col = F, legend = T)
dev.off()
print("Immune score heatmaps saved to separate PDF files")
#-----------------------------------------------------------------------------------------
# 6. Immune infiltration heatmaps - Separate PDFs for Up and Down
pdf("6.ImmuneInfiltration_Heatmap_UPREGULATED.pdf", width = 10, height = 10)
print("Generating immune infiltration heatmap for UPREGULATED genes...")
gs_immucell_heatmap(geneset_U, method = "pearson", lowcol = "blue", highcol = "red", 
                    cluster_row = F, cluster_col = F, legend = T)
dev.off()

pdf("6.ImmuneInfiltration_Heatmap_DOWNREGULATED.pdf", width = 10, height = 10)
print("Generating immune infiltration heatmap for DOWNREGULATED genes...")
gs_immucell_heatmap(geneset_D, method = "pearson", lowcol = "blue", highcol = "red", 
                    cluster_row = F, cluster_col = F, legend = T)
dev.off()
print("Immune infiltration heatmaps saved to separate PDF files")
#-----------------------------------------------------------------------------------------
pdf("7.immune_checkpoint_genes_UPREGULATED.pdf", width = 10, height = 10)
print("Generating immune_checkpoint_gene heatmap for DOWNREGULATED genes...")
gs_checkpoint_heatmap(geneset_U, method = "pearson", lowcol = "blue", highcol = "red", 
                    cluster_row = F, cluster_col = F, legend = T)
dev.off()
print("immune_checkpoint_gene heatmaps saved to separate PDF files")


pdf("7.immune_checkpoint_gene_Heatmap_DOWNREGULATED.pdf", width = 10, height = 10)
print("Generating immune_checkpoint_gene heatmap for DOWNREGULATED genes...")
gs_checkpoint_heatmap(geneset_D, method = "pearson", lowcol = "blue", highcol = "red", 
                    cluster_row = F, cluster_col = F, legend = T)
dev.off()
print("immune_checkpoint_gene heatmaps saved to separate PDF files")
#-----------------------------------------------------------------------------------------
print("All analyses completed! 12 separate PDF files have been created:")
print("UPREGULATED GENE SET:")
print("1. Chemokine_Heatmap_UPREGULATED.pdf")
print("2. Receptor_Heatmap_UPREGULATED.pdf")
print("3. ImmuneStimulator_Heatmap_UPREGULATED.pdf")
print("4. ImmuneInhibitor_Heatmap_UPREGULATED.pdf")
print("5. ImmuneScore_Heatmap_UPREGULATED.pdf")
print("6. ImmuneInfiltration_Heatmap_UPREGULATED.pdf")
print("7. immune_checkpoint_genes_UPREGULATED.pdf")
print("")
print("DOWNREGULATED GENE SET:")
print("8. Chemokine_Heatmap_DOWNREGULATED.pdf")
print("9. Receptor_Heatmap_DOWNREGULATED.pdf")
print("10. ImmuneStimulator_Heatmap_DOWNREGULATED.pdf")
print("11. ImmuneInhibitor_Heatmap_DOWNREGULATED.pdf")
print("12. ImmuneScore_Heatmap_DOWNREGULATED.pdf")
print("13. ImmuneInfiltration_Heatmap_DOWNREGULATED.pdf")
print("14. immune_checkpoint_genes_DOWNREGULATED.pdf")


# ### Genset based cancer type specific analysis
# # Genset based expression analysis grouped by clinical information (Genset based tumor-normal boxplot)
# gs_boxplot(cancer,geneset,geneset_alias,add = "jitter",
#            palette="jco",legend="none",label="p.signif",method="wilcox.test")
# 
# #----------------------------------------------------------------
# # Survavial analysis
# # Survavial analysis based on the expression of a single gene
# #----------------------------------------------------------------
# # Create a directory to save the plots if it doesn't exist
# output_dir <- "km_plots"
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir)
# }
# 
# # Loop through each gene and cancer type, and save the Kaplan-Meier curve as a PDF
# for (gene in genes) {
#   for (cancer in cancers) {
#     # Define the file name and path
#     file_name <- paste0(output_dir, "/", cancer, "_", gene, "_kmplot.pdf")
# 
#     # Open a PDF device
#     pdf(file_name)
# 
#     # Check if the function is generating a plot
#     print(paste("Plotting and saving", gene, "for", cancer))
# 
#     # Capture the plot object
#     km_plot <- tcga_kmplot(cancer, gene, palette = 'jco')
# 
#     # Check if the plot was generated
#     if (!is.null(km_plot)) {
#       print(km_plot)  # Explicitly print the plot
#     } else {
#       message(paste("No plot generated for", gene, "in", cancer))
#     }
# 
#     # Close the PDF device
#     dev.off()
#   }
# }
# ###################################### Stagewsie boxplot ######################################################
library(TCGAplot)
library(openxlsx)
library(ggplot2)
library(gridExtra)

# Define gene sets
geneset_U <- c("SIM2", "OTX1", "MYBL2", "HJURP", "IQGAP3", "KIF18B",
               "TPX2", "RRM2", "SKA3", "RHPN1", "SLC7A11")

geneset_D <- c("MYL9", "ABCG2", "SORBS1", "NRG1", "CRYAB", "PYGM", "SYNM",
               "CAV1", "MASP1", "FHL1", "RXRG", "TCEAL2", "DES")

geneset_alias <- "Key Signature Genes"



library(gridExtra)

# Save plots into variables-Up Regulated
p1 <- gs_stage("BRCA", geneset_U, "geneset_alias")
p2 <- gs_stage("LUAD", geneset_U, "geneset_alias")
p3 <- gs_stage("PRAD", geneset_U, "geneset_alias")
p4 <- gs_stage("COAD", geneset_U, "geneset_alias")

# Open PDF
pdf("geneset_U.pdf", width = 5, height = 5)
# Arrange plots in 2x2 grid
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()

## Geneset based Stagewise _Down Regulated
p5 <-gs_stage("BRCA", geneset_D, "geneset_alias")
p6 <-gs_stage("COAD", geneset_D, "geneset_alias")
p7 <-gs_stage("LUAD", geneset_D, "geneset_alias")
p8 <-gs_stage("PRAD", geneset_D, "geneset_alias")

# Open PDF
pdf("geneset_D.pdf", width = 5, height = 5)
# Arrange plots in 2x2 grid
grid.arrange(p5, p6, p7, p8, ncol = 2)
dev.off()

# Geneset based diagnostic ROC Curve_Up Regulated
p9 <-gs_roc("BRCA",geneset_U,"geneset_alias")
p10 <-gs_roc("COAD",geneset_U,"geneset_alias")
p11 <-gs_roc("LUAD",geneset_U,"geneset_alias")
p12 <-gs_roc("PRAD",geneset_U,"geneset_alias")
# Open PDF
pdf("geneset_ROC_UP.pdf", width = 5, height = 5)
# Arrange plots in 2x2 grid
grid.arrange(p9, p10, p11, p12, ncol = 2)
dev.off()

# Geneset based diagnostic ROC Curve_Down Regulated
p13 <-gs_roc("BRCA",geneset_D,"geneset_alias")
p14 <-gs_roc("COAD",geneset_D,"geneset_alias")
p15 <-gs_roc("LUAD",geneset_D,"geneset_alias")
p16 <-gs_roc("PRAD",geneset_D,"geneset_alias")

# Open PDF
pdf("geneset_ROC_DOWN.pdf", width = 5, height = 5)
# Arrange plots in 2x2 grid
grid.arrange(p13, p14, p15, p16, ncol = 2)
dev.off()

#--------------------------------------------------------------------------------------------------------
## NORMAL EXPRESSION of GENESET BETWEEN TWO GROUPS (NORMAL vs TUMOR)

p17<-gs_boxplot("BRCA",geneset_U,"geneset_alias",add = "jitter",
                palette="lancet",legend="none",label="p.signif",method="wilcox.test")

p18<-gs_boxplot("BRCA",geneset_D,"geneset_alias",add = "jitter",
                palette="lancet",legend="none",label="p.signif",method="wilcox.test")

p19<-gs_boxplot("COAD",geneset_U,"geneset_alias",add = "jitter",
                palette="lancet",legend="none",label="p.signif",method="wilcox.test")

p20<-gs_boxplot("COAD",geneset_D,"geneset_alias",add = "jitter",
                palette="lancet",legend="none",label="p.signif",method="wilcox.test")

p21<-gs_boxplot("LUAD",geneset_U,"geneset_alias",add = "jitter",
                palette="lancet",legend="none",label="p.signif",method="wilcox.test")

p22<-gs_boxplot("LUAD",geneset_D,"geneset_alias",add = "jitter",
                palette="lancet",legend="none",label="p.signif",method="wilcox.test")

p23<-gs_boxplot("PRAD",geneset_U,"geneset_alias",add = "jitter",
                palette="lancet",legend="none",label="p.signif",method="wilcox.test")

p24<-gs_boxplot("PRAD",geneset_D,"geneset_alias",add = "jitter",
                palette="lancet",legend="none",label="p.signif",method="wilcox.test")

# Open PDF
pdf("geneset_N.vs.T_expression.pdf", width = 4, height = 4)
# Arrange plots in 2x2 grid
grid.arrange(p17, p18, p19, p20, p21, p22, p23, p24, ncol = 2)
dev.off()





#--------------------------------------------------------------------------------------------------------
# 2. Genset based pan-cancer Cox regression analysis (Genset based pan-cancer Cox regression forest plot)

library(TCGAplot)

# Define gene sets
geneset <- c("SIM2", "OTX1", "MYBL2", "HJURP", "IQGAP3", "KIF18B",
             "TPX2", "RRM2", "SKA3", "RHPN1", "SLC7A11", "TNS4", 
             "FOLR1", "MYL9", "ABCG2", "SORBS1", "NRG1", "CRYAB", 
             "PYGM", "SYNM","CAV1", "MASP1", "FHL1", "RXRG", "TCEAL2", "DES")

geneset_alias <- "Key Signature Genes"

# Genset based pan-cancer Cox regression forest plot
pdf("HAZARD_ratio.pdf", width = 10, height = 10)
print("HAZARD_ratio Calculation.......")
gs_pan_forest(geneset,geneset_alias,adjust=F)
dev.off()

print("HAZARD_ratio heatmaps saved to PDF files")


















