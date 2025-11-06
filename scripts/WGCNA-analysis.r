# Load necessary libraries
library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(dplyr)
library(tibble)
library(biomaRt)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
library(foreach)
library(doParallel)
library(parallel)

# Set strings not as factors
options(stringsAsFactors = FALSE)
# Enable multithreading
enableWGCNAThreads()

# Read the raw data (rows are the sample and columns are the genes)
expressiondata = read.csv("extracted_DEGs_L.csv")
### step-1 duplicate entries
expressiondata<- expressiondata[!duplicated(expressiondata$gene_names), ]
#remove any NA values
expressiondata <- na.omit(expressiondata)

# Create a new format for expression data - remove gene name column
expression = as.data.frame(expressiondata[, -c(1)]) 
expression = t(expression)

# Set gene names as column names
colnames(expression) = expressiondata$gene_names
rownames(expression) = names(expressiondata)[-c(1)]

# Group data in a dendrogram to check outliers
sampleTree = hclust(dist(expression), method = "average")

# Read phenotypic data
traitData <- read.csv("LC_manifest.csv", sep = ",")
# Assuming 'expression' is a matrix or data frame with row names
Samples <- rownames(expression)
# Match samples in expression data with traitData
traitRows <- match(Samples, traitData$ID)
# Subset traitData to match the samples in expression data, excluding the ID column
datTraits <- traitData[traitRows, -1]
# Convert datTraits to a data frame if it's not already (just to ensure)
datTraits <- as.data.frame(datTraits)

rownames(datTraits) = traitData[traitRows, 1]
# Regroup samples
sampleTree2 = hclust(dist(expression), method = "average")
# Check unique values in datTraits
unique(datTraits$datTraits)
# Assign disease_state_bin based on datTraits values
datTraits$disease_state_bin <- ifelse(datTraits$datTraits == 'Primary_Solid_Tumour', 1, 0)

# Function to convert values to colors (adjust as per your numbers2colors function)
numbers2colors <- function(data, signed = FALSE) {
  # Example implementation, replace with your actual function
  ifelse(data == "Primary_Solid_Tumour", "red", "blue")
}

# Add traitColors column based on datTraits values
traitColors <- numbers2colors(datTraits$datTraits)

# Plot a sample dendrogram with the colors below
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")

# Creating the network
par(mfrow = c(1,2))
cex1 = 0.9



# Set gene_names as row names
rownames(expressiondata) <- expressiondata$gene_names
expressiondata <- expressiondata[ , -1]

data <- as.data.frame(expressiondata)
#---Needed when using rNA-seq count data

# data <- data %>%
#   filter(rowSums(.) >= 0)
#--------------------------------------
# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = round(data),
                              colData = traitData,
                              design = ~ dex)
dds$dex <- factor(dds$dex, levels = c("Solid_Tumour_Normal, Primary_Solid_Tumour"))

# Perform variance stabilization
dds_norm <- vst(dds)

# Get normalized counts
expression1 <- assay(dds_norm) %>% t() 

##if you data already normalized
expression1 <- assay(dds) %>% t() 


# -------------------------------
# 1. Choose soft-thresholding power
# -------------------------------
powers <- c(1:10, seq(12, 30, by = 2))  # Slightly cleaner sequence
sft <- pickSoftThreshold(expression,
                         powerVector = powers,
                         networkType = "signed",
                         verbose = 5)
#==================================================================================================
# # Automated selection
# softPower <- sft$powerEstimate
# selection_method <- "WGCNA estimate"
# 
# if (is.na(softPower)) {
#   softPower <- which(sft$fitIndices[,2] > 0.8)[1]
#   selection_method <- "First power with R² > 0.8"
# }
# 
# if (is.na(softPower)) {
#   max_r2_index <- which.max(sft$fitIndices[,2])
#   softPower <- sft$fitIndices[max_r2_index, 1]
#   selection_method <- "Power with maximum R²"
# }
# 
# # FIX: Find the correct row index for the selected power
# selected_row <- which(sft$fitIndices[,1] == softPower)
# selected_r2 <- sft$fitIndices[selected_row, 2]
# selected_connectivity <- sft$fitIndices[selected_row, 5]
# 
# # Plot results with selection highlighted
# par(mfrow = c(1,2))
# cex1 = 0.9
# 
# # Plot 1: Scale independence
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R²",
#      type="n", main = paste("Scale independence\n", selection_method))
# 
# # Plot all points
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=power, cex=cex1, col="red")
# 
# # Highlight selected power in BLUE
# text(sft$fitIndices[selected_row, 1], 
#      -sign(sft$fitIndices[selected_row, 3])*sft$fitIndices[selected_row, 2],
#      labels=softPower, cex=cex1, col="blue", font=2)
# 
# abline(h=0.85, col="red")
# 
# # Plot 2: Mean connectivity
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)", ylab="Mean Connectivity", 
#      type="n", main = paste("Mean connectivity"))
# 
# # Plot all points
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=power, cex=cex1, col="red")
# 
# # Highlight selected power in BLUE
# text(sft$fitIndices[selected_row, 1], sft$fitIndices[selected_row, 5],
#      labels=softPower, cex=cex1, col="blue", font=2)
# 
# abline(h=40, col="red")
# 
# # Print summary
# cat("=== SOFT THRESHOLDING POWER SELECTION ===\n")
# cat("Selection method:", selection_method, "\n")
# cat("Selected power:", softPower, "\n")
# cat("R² at selected power:", round(selected_r2, 4), "\n")
# cat("Mean connectivity:", round(selected_connectivity, 2), "\n")
#===================================================================================
# -------------------------------
# 2. Plot scale-free topology and mean connectivity
# -------------------------------
par(mfrow = c(1, 2))
cex1 <- 0.9

# Scale-free topology plot
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.80, col = "red", lty = 2)

# Mean connectivity plot
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")
abline(h = 40, col = "red", lty = 2)

# -------------------------------
# 3. Construct network
# -------------------------------
softPower <- 16  # Chosen based on above plots
adjacency <- adjacency(expression, power = softPower, type = "signed")
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM  # distance matrix for clustering

# -------------------------------
# 4. Hierarchical clustering of genes
# -------------------------------
geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree, main = "Clustering of genes", xlab = "", sub = "")

# -------------------------------
# 5. Dynamic tree cut for modules
# -------------------------------
dynamicColors <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                               deepSplit = 2, pamRespectsDendro = FALSE,
                               minClusterSize = 30)
table(dynamicColors)

# -------------------------------
# 6. Calculate module eigengenes
# -------------------------------
MEList <- moduleEigengenes(expression, colors = dynamicColors)
MEs <- MEList$eigengenes

# Optional: remove first column if necessary (usually not needed)
# MEs <- MEs[, -1]

# -------------------------------
# 7. Cluster module eigengenes
# -------------------------------
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
MEDissThres <- 0.55
abline(h = MEDissThres, col = "red", lty = 2)

# -------------------------------
# 8. Merge similar modules
# -------------------------------
merge <- mergeCloseModules(expression, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

# -------------------------------
# 9. Plot dendrogram with merged modules
# -------------------------------
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)




# 2. Calculate the number of genes for each module
module_gene_counts <- as.data.frame(table(merge$colors))
colnames(module_gene_counts) <- c("Module", "GeneCount")
module_gene_counts

# 3. Function to create a heatmap for the module eigengenes
library(ComplexHeatmap)
library(circlize)

# Assuming traitData$dex is a factor with levels "normal" and "tumor"
# Define a color map for the annotations
annotation_colors <- list(dex = c("Solid_Tumour_Normal" = "blue", "Primary_Solid_Tumour" = "red"))

# 3. Function to create a heatmap for the module eigengenes
plot_eigengene_heatmap <- function() {
  # Transform module eigengenes data for heatmap
  module_eigengenes_matrix <- as.matrix(MEs)
  
  # Ensure that the row names of the matrix are set to sample IDs or something similar
  rownames(module_eigengenes_matrix) <- rownames(traitData)
  
  # Create a row annotation based on traitData$dex
  row_annotation <- rowAnnotation(
    dex = traitData$dex,
    col = annotation_colors,
    show_legend = TRUE,
    annotation_legend_param = list(dex = list(title = "Group"))
  )
  
  # Define the color function for the heatmap
  color_func <- colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
  
  # Create the heatmap
  heatmap <- Heatmap(module_eigengenes_matrix,
                     name = "Eigengenes",
                     col = color_func,
                     cluster_rows = FALSE,
                     cluster_columns = TRUE,
                     show_row_names = FALSE,
                     show_column_names = TRUE,
                     column_names_gp = gpar(fontsize = 10),
                     heatmap_legend_param = list(title = "Eigengenes"),
                     left_annotation = row_annotation)  # Use left_annotation for row annotations
  draw(heatmap)
}

# Call the function to create the heatmap
plot_eigengene_heatmap()


# Additional analysis can be performed here, such as relating module eigengenes to external traits
# Calculate correlations between module eigengenes and traits
moduleTraitCor <- cor(mergedMEs, datTraits$disease_state_bin, use = "p")   ## if want to remove additioan ncolum that not match with matrix just used mergedMEs[, -colno.]
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(mergedMEs))

# Prepare text matrix for heatmap annotation
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# Define xLabels based on the number of columns in moduleTraitCor
ncol_moduleTraitCor <- ncol(moduleTraitCor)
xLabels <- names(datTraits)[1:ncol_moduleTraitCor]

# Adjusted margins to fit within the plotting device
par(mar = c(3, 3, 1, 1))  # Example: bottom, left, top, right

# Plotting the heatmap
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = xLabels,
  yLabels = names(mergedMEs),  # Assuming names(mergedMEs) is correctly defined
  ySymbols = names(mergedMEs),  # Assuming names(mergedMEs) is correctly defined
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = "Module-trait relationships"
)


# Display the correlation values and p-values
print(moduleTraitCor)
print(moduleTraitPvalue)

# Identify significant modules (e.g., p-value < 0.05)
significantModules = names(mergedMEs)[apply(moduleTraitPvalue, 1, function(p) any(p < 0.05))]
print(significantModules)

#---------------------------------------------------------------------------------------------
##Defining the variable state containing the column state of datTrait
state = as.data.frame(datTraits$disease_state_bin)
names(state) = "state"
# Define numbers of genes and samples
nSamples <- nrow(expression1)
nGenes <- ncol(expression1)
Genes <- row.names(expressiondata)
# Calculate module membership (kME)
moduleMembership <- as.data.frame(cor(expression1, mergedMEs, use = "p"))
names(moduleMembership) <- paste("p.MM", names(mergedMEs), sep = "")
# Calculate p-values for module membership measures
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(moduleMembership), nSamples))
# Calculate gene significance (GS)
geneTraitSignificance <- as.data.frame(cor(expression1, datTraits$disease_state_bin, use = "p"))
names(geneTraitSignificance) <- paste("GS", names(datTraits$disease_state_bin), sep = "")
# Calculate p-values for gene significance
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
write.csv(GSPvalue, file ="GSPvalue_C.csv", row.names = T)

#make a dataframe of genes and associated module colour
module_df <- data.frame(
  gene_id = colnames(expression),
  colors = labels2colors(merge$colors),
  stringsAsFactors = FALSE  # Ensure colors are treated as characters
)

# Combine the results
geneInfo = data.frame(Genes,
                      moduleMembership,module_df,
                      geneTraitSignificance,GSPvalue)
write.csv(geneInfo, file ="geneInfo.csv", row.names = T)


# Ensure moduleGenes correctly subsets rows based on moduleColors and module
module <- 6
moduleColors <- merge$colors
moduleGenes <- moduleColors == module
# Ensure column correctly identifies the column index corresponding to the module
column <- match(paste("p.MMME", module, sep=""), names(moduleMembership))
# Plotting the scatterplot
par(mfrow = c(1, 1))  # Reset the plot layout to 1x1
Scplot <- verboseScatterplot(
  abs(moduleMembership[moduleGenes, column]),
  abs(geneTraitSignificance[moduleGenes, 1]),
  xlab = paste("Module Membership in", module, "module"),
  ylab = "Gene significance for Tumour",
  main = "Module membership vs. gene significance",
  cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,
  pch = 19, col = module
)
# Overlay with unfilled dots to create border effect
points(
  abs(moduleMembership[moduleGenes, column]),
  abs(geneTraitSignificance[moduleGenes, 1]),
  col = "black",  # Set border color
  pch = 1  # Set pch = 1 for unfilled dots (circles)
)

#-------------------------------------# Extract gene names for modules --------------------------
#-----------------------------------GENE ID CONVERTER----------------------------------------------
library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(dplyr)
library(tibble)
library(biomaRt)
library(org.Hs.eg.db)
library(tidyverse)
library(dplyr)

geneInfo <- geneInfo[, -1]
Input <- as.data.frame(geneInfo)
# Convert row names back to a column
Input <- Input %>%
  rownames_to_column("Gene")

# input list of Ensembl IDs
ensembl_ids <- Input$Gene

# Use biomaRt to connect to Ensembl
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)

# Use the human dataset ('hsapiens_gene_ensembl')
ensembl.con <- useMart("ensembl", dataset = 'mmusculus_gene_ensembl')

# Retrieve the gene names based on Ensembl IDs
gene_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                   filters = "ensembl_gene_id",
                   values = ensembl_ids,
                   mart = ensembl.con)

# Merge with the original data frame
Input <- Input %>%
  left_join(gene_info, by = c("Gene" = "ensembl_gene_id"))

# View the result
print(Input)
geneInfo <- Input
rm(Input)
write.csv(geneInfo, file ="geneInfo__FINAL2.csv", row.names = F)
   #if you have ensemblegene
#--------------------------------------------EXPORT NETWORK--------------------------------------------------------
#Select the gene modules
module_colors <- merge$colors
genes = colnames(expression1)

#if you want export specific colors, substitute the second modulecolors by above modules
inModule = is.finite(match(module_colors, module_colors))
modGenes = genes[inModule]

#Select the corresponding topologic overlap 
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modGenes, modGenes)
modTOMSignificantes = which(modTOM>0.0)

#Create the dataframe since the beginning
geneInfo0 = data.frame(ESTs = genes,
                       moduleColor = merge$colors,
                       geneTraitSignificance,
                       GSPvalue)

#Export the network in list files os n edges that cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = "CytoscapeEdgeFile7.txt",
                               nodeFile = "CytoscapeNodeFile7.txt",
                               weighted = TRUE,
                               threshold = 0.1,
                               nodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])

