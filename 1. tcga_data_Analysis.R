###D:\3.Integerated_Results for Manuscripts\Aftab\scripts

# Load packages
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")

#################################################################################################################
## TCGA DATA AQUISATION                   PART-1
#################################################################################################################
# We can check all the available projects at TCGA with the command bellow.
GDCprojects = getGDCprojects()
head(GDCprojects[c("project_id", "name")])

# Creat A folder (Directory)
dir.create("BRCA")

# We can use the following function to get details on all data deposited for TCGA-BRCA.
TCGAbiolinks:::getProjectSummary("TCGA-BRCA")           ## TCGA-BRCA, TCGA-COAD, TCGA-LUAD

# Note that performing this query will take a few of minutes
query_TCGA = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts")

# To visualize the query results in a more readable way :make results as table
lihc_res = getResults(query_TCGA)
head(lihc_res) # data of the first 6 patients.
colnames(lihc_res) # columns present in the table

## One interesting question is the tissue type measured at an experiment (normal, solid tissue, cell line). This information is present at column “sample_type”.
head(lihc_res$sample_type) # first 6 types of tissue

# We can visualize it better with the summary function.
summary(factor(lihc_res$sample_type))

# build a query to retrieve gene expression data ------------
query_TCGA <- GDCquery(project = 'TCGA-BRCA',
                       data.category = 'Transcriptome Profiling',    ## Proteome Profiling
                       experimental.strategy = 'RNA-Seq',
                       data.type = 'Gene Expression Quantification',
                       workflow.type = 'STAR - Counts',
                       access = 'open',
                       data.format = 'TSV',
                       sample.type = c('Primary Tumor', 'Solid Tissue Normal'))     ## "Primary Tumor", "Solid Tissue Normal"

getResults(query_TCGA)


# Let us now download the files specified in the query.
GDCdownload(query_TCGA)


## Finally, lets load the actual RNASeq data into R. Remember that the output directory set must be the same to where you downloaded the data.
tcga_data = GDCprepare(query_TCGA)

#####
# There are 3 functions that allow us to access to most important data present in this object, these are: colData(), rowData(), assays(). 
# Use the command ?SummarizedExperiment to find more details. colData() allows us to access the clinical data associated with our samples. 
# The functions colnames() and rownames() can be used to extract the column and rows names from a given table respectively.

#----------------------------------------------------------------------------------------------------------------------
# In R (and other programming languages) you can often
# chain functions to save time and space
colnames(colData(tcga_data))

## This will provides a basic explanation about tcga_data
table(tcga_data@colData$vital_status)
table(tcga_data@colData$definition)
table(tcga_data@colData$tissue_or_organ_of_origin)
table(tcga_data@colData$gender)
table(tcga_data@colData$race)
dim(assay(tcga_data)) 

# gene expression matrices
dim(assay(tcga_data))

# expression of first 6 genes and first 10 samples
head(assay(tcga_data)[,1:10])
# ensembl id and gene id of the first 6 genes.
head(rowData(tcga_data))    


# Save the data as a file, if you need it later, you can just load this file
# instead of having to run the whole pipeline again
saveRDS(object = tcga_data,
        file = "tcga_data.RDS",
        compress = FALSE)

# The data can be loaded with the following command
# tcga_data = readRDS(file = "tcga_data.RDS")


#################################################################################################################
## RNASeq Normalization                    PART-2
#################################################################################################################
###  USE THIS FUNCTION

limma_pipeline = function(
    tcga_data,
    condition_variable,
    reference_group=NULL){
  
  design_factor = colData(tcga_data)[, condition_variable, drop=T]
  
  group = factor(design_factor)
  if(!is.null(reference_group)){group = relevel(group, ref=reference_group)}
  
  design = model.matrix(~ group)
  
  dge = DGEList(counts=assay(tcga_data),
                samples=colData(tcga_data),
                genes=as.data.frame(rowData(tcga_data)))
  
  # filtering
  keep = filterByExpr(dge,design)
  dge = dge[keep,,keep.lib.sizes=FALSE]
  rm(keep)
  
  # Normalization (TMM followed by voom)
  dge = calcNormFactors(dge)
  v = voom(dge, design, plot=TRUE)
  
  # Fit model to data given design
  fit = lmFit(v, design)
  fit = eBayes(fit)
  
  # Show top genes
  topGenes = topTable(fit, coef=ncol(design), number=3700, sort.by="p")    ## otherwisw 100 or 200 top DEGs is OK
  
  return(
    list(
      voomObj=v, # normalized data
      fit=fit, # fitted model and statistics
      topGenes=topGenes # the 100 most differentially expressed genes
    )
  )
}
########################################################################################
## With the following command, we can obtain the DE analysis comparing samples
#-------------------------------------------------------------------------------------------------------------------
## option-1 : when you have two group of samples (Normal +controle) then defination will be other samples

limma_res = limma_pipeline(
  tcga_data=tcga_data,
  condition_variable="definition",
  reference_group="Solid Tissue Normal"   ## Control group
)
#SAVE OBJECT
saveRDS(object = limma_res,
        file = "limma_res.RDS",
        compress = FALSE)
#-------------------------------------------------------------------------------------------------------------------
## option-2 : within samples comparison

# As an example, we also show here how we can use limma_pipeline to perform DE analysis by grouping patients by gender instead of by tissue type.
limma_res = limma_pipeline(
  tcga_data=tcga_data,
  condition_variable="gender",
  reference_group="female"
)

#SAVE OBJECT
saveRDS(object = limma_res,
        file = "limma_res.RDS",
        compress = FALSE)
#-------------------------------------------------------------------------------------------------------------------
### Visualization
plot_PCA = function(voomObj, condition_variable){
  group = factor(voomObj$targets[, condition_variable])
  pca = prcomp(t(voomObj$E))
  # Take PC1 and PC2 for the plot
  plot(pca$x[,1:2],col=group, pch=19)
  # include a legend for points
  legend("bottomleft", inset=.01, levels(group), pch=19, col=1:length(levels(group)))
  return(pca)
}
res_pca = plot_PCA(limma_res$voomObj, "definition")

# Produce a PCA plot using gender as condition (hint: remember that we already run the pipeline for gender). 
# emember that you can see all available clinical features/phenotypes with the following command:
#   colnames(colData(tcga_data))
########################################################################################
## PART-2 (OPTIONAL)
# We will select a particular clinical feature of the data to use as class for grouping the samples as either tumor vs normal tissue. 
# so we convert it as such
#----------------------------------------------------------------------------------------------------------------------

clinical_data = colData(tcga_data)
group = factor(clinical_data$definition)
 
#define Solid Tissue Normal as being our base or reference level.
group = relevel(group, ref="Solid Tissue Normal")   ## (our control samples)

## we need to create a design matrix, which will indicate the conditions to be compared by the DE analysis
design = model.matrix(~group)
head(design)

## Before performing DE analysis, we remove genes, which have low amount of counts. We transform our tcga_data object as DGEList, which provides functions for filtering. 
# Create DGEList object and filter genes
dge <- DGEList(counts=assay(tcga_data), samples=colData(tcga_data), genes=as.data.frame(rowData(tcga_data)))
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
rm(keep) #  use rm() to remove objects from memory if you don't need them anymore

# Normalize data and apply VOOM transformation
dge <- calcNormFactors(dge, method="TMM")
v <- voom(dge, design, plot=TRUE)

# Fit linear model and apply empirical Bayes smoothing
fit <- lmFit(v, design)
fit <- eBayes(fit)
#---------------------------------------------------------------------------------------------------------------------
# Calculate average expression for each group
group_labels <- clinical_data$definition
tumor_samples <- which(group_labels == "Primary solid Tumor")
normal_samples <- which(group_labels == "Solid Tissue Normal")
average_tumor <- rowMeans(assay(tcga_data)[, tumor_samples])
average_normal <- rowMeans(assay(tcga_data)[, normal_samples])

# Extract differential expression results
results <- topTable(fit, coef=2, sort.by = "P", number=Inf)
results$Gene <- rownames(results)  # Ensure 'Gene' column is present

# Create a data frame with average expressions
average_expression <- data.frame(
  Gene = rownames(assay(tcga_data)),
  Average_Tumor = average_tumor,
  Average_Normal = average_normal
)

# Merge the average expression data with the DE results
results_with_averages <- merge(results, average_expression, by = "Gene")
write.csv(results_with_averages, file = "full_genes_list.csv", row.names = FALSE)

results1 <- results_with_averages[which(abs(results_with_averages$logFC) > 1.5 & results_with_averages$adj.P.Val < 0.05),]
write.csv(results1, file = "l2ogFC1.5.csv", row.names = FALSE)

#---------------------------------------------------------------------------------------------------------------------
#### Extract full sample information and merge
# sample_info <- as.data.frame(colData(tcga_data))
# results_with_sample_info <- cbind(results_with_averages, sample_info)
# 
# # View the final results
# head(results_with_sample_info)
#---------------------------------------------------------------------------------------------------------------------
# #option-3#  Using the function topTable we can check the top10 genes classified as being differentially expressed. 
# topGenes = topTable(fit, coef=1, sort.by="p")
# print(topGenes)
# write.csv(topGenes, file = "topGenes2_list.csv", row.names = FALSE)
#---------------------------------------------------------------------------------------------------------------------


#################################################################################################################
##  Classification(Train and test paradigm)                 PART-3    https://www.costalab.org/wp-content/uploads/2020/11/R_class_D3.html
#################################################################################################################
## Now, we will explore the use of a machine learning method to classify an unseen sample as being a tumor or not. 
#To achieve this goal we are going to build first a simple linear model (with feature selection), and then an Elastic Net model. 
#----------------------------------------------------------------------------------------------------------------------
load("BRCA/BRCA1.RData")
# Transpose and make it into a matrix object
d_mat = as.matrix(t(limma_res$voomObj$E))

# As before, we want this to be a factor
d_resp = as.factor(limma_res$voomObj$targets$definition)


# Divide data into training and testing set

# Set (random-number-generator) seed so that results are consistent between runs
set.seed(42)
train_ids = createDataPartition(d_resp, p=0.75, list=FALSE)

x_train = d_mat[train_ids, ]
x_test  = d_mat[-train_ids, ]

y_train = d_resp[train_ids]
y_test  = d_resp[-train_ids]

### Train Elastic Net model
###########################

# Train model on training dataset using cross-validation
res = cv.glmnet(
  x = x_train,
  y = y_train,
  alpha = 0.5,
  family = "binomial"
)

###  Model Evaluation
###########################

# Test/Make prediction on test dataset
y_pred = predict(res, newx=x_test, type="class", s="lambda.min")

confusion_matrix = table(y_pred, y_test)  ## A confusion matrix is a simple table that compares the predictions from our model against their real values. Therefore, it shows us the true positives, true negatives, false positives and false negatives. We can use it to compute a number of accuracy metrics that we can then use to understand how good our model actually is.

# Evaluation statistics
print(confusion_matrix)

print(paste0("Sensitivity: ",sensitivity(confusion_matrix)))   #  ## As you see = 1, for this data and under the current parameters, Elastic Net has great accuracy.
print(paste0("Specificity: ",specificity(confusion_matrix)))
print(paste0("Precision: ",precision(confusion_matrix)))

# We can now take a look at the genes (coefficients), that Elastic Net selected to build its model.
# Getting genes that contribute for the prediction
res_coef = coef(res, s="lambda.min") # the "coef" function returns a sparse matrix
dim(res_coef)
head(res_coef) # in a sparse matrix the "." represents the value of zero

## Of course, the number of coefficients is large (there are many genes!). 
#We only want to consider coefficients with non-zero values, as these represent variables (genes) selected by the Elastic Net.

# get coefficients with non-zero values
res_coef = res_coef[res_coef[,1] != 0,]
# note how performing this operation changed the type of the variable
head(res_coef)

# remove first coefficient as this is the intercept, a variable of the model itself
res_coef = res_coef[-1]

relevant_genes = names(res_coef) # get names of the (non-zero) variables.
length(relevant_genes) # number of selected genes
head(relevant_genes) # few select genes

##ensembl to gene name correspondence.
head(limma_res$voomObj$genes)

relevant_gene_names = limma_res$voomObj$genes[relevant_genes,"gene_name"]

head(relevant_gene_names) # few select genes (with readable names now)

# Create a data frame with gene names
gene_df <- data.frame(GeneName = relevant_gene_names)

# Save the data frame to a CSV file
write.csv(gene_df, "relevant_genes.csv", row.names = FALSE)

# Did limma and Elastic Net select some of the same genes? 
# We can check the common genes between the two results by using the intersect function.

print(intersect(limma_res$topGenes$gene_id, relevant_genes)) ## NOTE: we do not expect a high overlap between genes selected by limma and Elastic net. The reason for this is the fact Elastic Net criteria bias the selection of genes, which are not highly correlated against each other, while not such bias is present in limma.

#################################################################################################################
## Hierarchical clustering                 PART-4
#################################################################################################################
## Finally we can take a look at how our samples cluster together by running an hierarchical clustering algorithm. 
#We will only be looking at the genes Elastic Net found and use these to cluster the samples. 
#The genes highlighted in green are the ones that limma had also selected as we’ve seen before. 
#The samples highlighted in red are Solid Tissue Normal, the samples highlighted in black are Primary solid Tumor
#----------------------------------------------------------------------------------------------------------------------
# define the color palette for the plot
hmcol = colorRampPalette(rev(brewer.pal(9, "RdBu")))(256)

# perform complete linkage clustering
clust = function(x) hclust(x, method="complete")
# use the inverse of correlation as distance.
dist = function(x) as.dist((1-cor(t(x)))/2)

# Show green color for genes that also show up in DE analysis
colorLimmaGenes = ifelse(
  # Given a vector of boolean values
  (relevant_gene_names %in% limma_res$topGenes$gene_name),
  "green", # if true, return green for that value
  "white" # if false, return white for that value
)
###################### OPTIONAL ##############################   ## If you want to comapre with your DE gene list
# Show green color for genes that also show up in DE analysis
colorLimmaGenes = ifelse(
  # Given a vector of boolean values
  (relevant_gene_names %in% results1$gene_name),
  "green", # if true, return green for that value
  "white" # if false, return white for that value
)
############################################################


# Define the PDF output file
pdf("gene_heatmap.pdf", width = 8, height = 12)

# As you've seen a good looking heatmap involves a lot of parameters
gene_heatmap = heatmap.2(
  t(d_mat[,relevant_genes]),
  scale="row",          # Scale the values for each gene (row)
  density.info="none",  # Turn off density plot inside color legend
  trace="none",         # Turn off trace lines inside the heat map
  col=hmcol,            # Use the generated color map
  labRow=relevant_gene_names, # Use gene names instead of Ensembl annotation
  RowSideColors=colorLimmaGenes,
  labCol=FALSE,         # Not showing column labels
  ColSideColors=as.character(as.numeric(d_resp)), # Show colors for each response class
  dendrogram="both",    # Show dendrograms for both axes
  hclust = clust,       # Define hierarchical clustering method
  distfun = dist,       # Using correlation coefficient for distance function
  cexRow=.8,            # Resize row labels
  margins=c(1,5)        # Define margin spaces
)
dev.off()

#--------------------------------------------------------------------------
# Using the same method as in Day-2, get the dendrogram from the heatmap
# and cut it to get the 2 classes of genes

# Extract the hierarchical cluster from heatmap to class "hclust"
hc = as.hclust(gene_heatmap$rowDendrogram)

# Cut the tree into 2 groups, up-regulated in tumor and up-regulated in control
clusters = cutree(hc, k=2)
table(clusters)

# Extract genes in each cluster
cluster_genes <- split(relevant_gene_names, clusters)
# Get the maximum length of genes among clusters
max_len <- max(length(cluster_genes[[1]]), length(cluster_genes[[2]]))
# Pad shorter clusters with NA to make them equal length
cluster1_padded <- c(cluster_genes[[1]], rep(NA, max_len - length(cluster_genes[[1]])))
cluster2_padded <- c(cluster_genes[[2]], rep(NA, max_len - length(cluster_genes[[2]])))
# Create a data frame with two columns for cluster genes
cluster_df <- data.frame(Cluster1 = cluster1_padded, Cluster2 = cluster2_padded)
# Write the data frame to a tab-delimited text file
write.table(cluster_df, "cluster_genes.txt", sep = "\t", row.names = FALSE)

###--------------------------------------------------------------------------------------------------------------

# selecting just a few columns so that its easier to visualize the table
gprofiler_cols = c("significant","p.value","overlap.size","term.id","term.name")

# make sure the URL uses https
set_base_url("https://biit.cs.ut.ee/gprofiler")

# Group 1, up in tumor
gprofiler(names(clusters[clusters %in% 1]))[, gprofiler_cols]


# Group 2, up in control
gprofiler(names(clusters[clusters %in% 2]))[, gprofiler_cols]

#------------------ If gprofiler does not work for you, you can also use the web interface: https://biit.cs.ut.ee/gprofiler/gost

## ----------------------------------------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------Excersize (Optional)--------------------------------------------------------------------------
## Another way to reduce dimensionality is through the use of variance filtering.

# retain only a small subset of the genes (see documentation for ?varFilter)
d_mat = varFilter(limma_res$voomObj$E, var.func=IQR, var.cutoff=0.95)
# transpose the matrix, so that it has the same shape as the d_mat we used at the beginning
d_mat = t(d_mat)
#
print(dim(d_mat))


## This function takes the normalized RNASeq data stored in the limma_res object, and returns a matrix of 
#the same form What varFilter does is calculate the interquartile range (IQR) for all genes. 
#It then removes 95% of the genes, starting with those with lower IQR. Therefore, only 5% of the genes will be retained - those with the highest dispersion, 
#and therefore possibly the most useful information. This can easily be understood that remembering that a gene that has the exact same expression across all samples, has a variance (and IQR) of zero. We definitely want to remove those.

# size before
print(dim(x_train))
print(dim(x_test))


x_train = d_mat[train_ids, ]
x_test  = d_mat[-train_ids, ]

# size after
print(dim(x_train))
print(dim(x_test))


