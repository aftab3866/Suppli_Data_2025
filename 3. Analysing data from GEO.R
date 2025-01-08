###D:\3.Integerated_Results for Manuscripts\Aftab\scripts

## Load  packages before starting !
library("forcats")
library("stringr")
library("ggplot2")
library("ggrepel")
library("tidyr")
library("survminer")
library("GEOquery")
library("survminer")
library("limma")
library("org.Hs.eg.db")
library("pheatmap")
library("dplyr")
library("impute")


#------------------------------------------------------------------------------------------------
#                                  1. DATA PreProcessing
#------------------------------------------------------------------------------------------------

## change my_id to be the dataset that you want.
my_id <- "GSE6919"
gse <- getGEO(my_id)

## check how many platforms used
length(gse)

gse <- gse[[1]]      ## ## if more than one dataset is present, you can analyse the other dataset by changing the number inside the [[...]]
                     ## e.g. gse <- gse[[2]]
gse

pData(gse) ## print the sample information
fData(gse) ## print the gene annotation
exprs(gse) ## print the expression data


## Check the normalisation and scales used (The methods we will use assume the data are on a log2 scale; typically in the range of 0 to 16.)
## exprs get the expression levels as a data frame and get the distribution
summary(exprs(gse))

#------------------------------Check if values are higher------------------------------------------------------------------
### From this output if we that the values go beyond 16, we should  perform a log2 transformation otherwise directly make boxplot. 
exprs(gse) <- log2(exprs(gse))
summary(exprs(gse))
#------------------------------------------------------------------------------------------------
### Normalize using quantile normalization
gse_normalized <- normalizeQuantiles(exprs(gse))

# Before normalization
boxplot(exprs(gse), outline = FALSE, main = "Before Normalization")

# After normalization
boxplot(gse_normalized, outline = FALSE, main = "After Normalization")

#------------------------------------------------------------------------------------------------
###                               2. Inspect the clinical variables
#------------------------------------------------------------------------------------------------
sampleInfo <- pData(gse)
colnames(sampleInfo)

## source_name_ch1 and characteristics_ch1 (Normal & TUmour) seem to contain factors we might need for the analysis. Let's pick just those columns
sampleInfo <- select(sampleInfo, characteristics_ch1, geo_accession)
## Optionally, rename to more convenient column names
sampleInfo <- rename(sampleInfo, group = characteristics_ch1, samples=geo_accession)

#--------------------------------------OptiONAL-1--------------------------------------------------
# #### if you want to include only "Normal" and "Tumour" samples:
sampleInfo <- filter(sampleInfo, characteristics_ch1 %in% c("Tissue: primary prostate tumor", "Tissue:normal prostate tissue adjacent to tumor"))
# # rename columns to more convenient names
sampleInfo <- rename(sampleInfo, group = characteristics_ch1, samples = geo_accession)
#------------------------------------------------------------------------------------------------

## Sample clustering and Principal Components Analysis
corMatrix <- cor(exprs(gse),use="c")   ## argument use="c" stops an error if there are any missing data points

pheatmap(corMatrix) 

####################################### OPTINAL-2 ###############################################
## Print the rownames of the sample information and check it matches the correlation matrix
# rownames(sampleInfo)
# colnames(corMatrix)
### If not, force the rownames to match the columns
# rownames(sampleInfo) <- colnames(corMatrix)
# pheatmap(corMatrix,
#          annotation_col=sampleInfo)
###############################################################################################

## Principal Components Analysis (PCA)

## MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX
pca <- prcomp(t(exprs(gse)))
#---------------------------------------------------------------
# if yoy found error
# Load expression data
exprs_data <- exprs(gse)
# Replace infinite values with NA
exprs_data[is.infinite(exprs_data)] <- NA
# Impute missing values using k-nearest neighbors imputation
imputed_data <- impute.knn(exprs_data)$data
pca <- prcomp(t(imputed_data), scale. = TRUE)
#-----------------------------------------------------------------

## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=group,label=paste(samples))) + geom_point() + geom_text_repel()

#-----------------------------Optional: If you dont want label only DOTs-----------------------
## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y = PC2, col = group)) + 
  geom_point()
#----------------------------------------------------------------------------------------------
## Exporting the data
library(readr)
full_output <- cbind(fData(gse),exprs(gse))
write_csv(full_output, path="GSE6919//gse_full_output.csv")

#--------------------------------------------OPTIONAL-3--------------------------------------------------
# Or Select only few important column
### Look at the features data frame and decide the names of the columns you want to keep
features <- fData(gse)
colnames(features) ## Check the total numver of Columns

features <- features[, c("Gene Title", "ENTREZ_GENE_ID", "Gene Symbol")]
features <- features %>%
  tibble::rownames_to_column(var = "ID")

full_output <- cbind(features,exprs(gse))
write_csv(full_output, path="GSE6919/gse_full_output_ANNO.csv") # ## FULL DATA with expression (save it if needed)
#--------------------------------------------------------------------------------------------------------------------------------------
# Extract expression matrix from gse
expr_matrix <- exprs(gse)
#If yiu want to convert row name to colums
# expr_matrix <- expr_matrix %>%
#   tibble::rownames_to_column(var = "ID")

#---------------------------------OPTIONAL-4----------------------------------------------
# Get the sample IDs of the selected samples (see OPTIONAL-1)
selected_samples <- sampleInfo$samples
# Filter the expression matrix based on the selected samples
expr_matrix <- expr_matrix[, colnames(expr_matrix) %in% selected_samples]
#---------------------------------------------------------------------------------------

# Create a data frame from the expression matrix with row names as "ID"
expr_df <- data.frame(ID = rownames(expr_matrix), expr_matrix, stringsAsFactors = FALSE)

write.csv(expr_df, file = "GSE6919/gse_full_output_expression.csv", row.names = T)

#----------------------------------------------------------------------------------------------
###                                         Differential Expression:
#----------------------------------------------------------------------------------------------
library(limma)
design <- model.matrix(~0 + sampleInfo$group)
design

## the column names are a bit ugly, so we will rename
colnames(design) <- c("Tumour", "Normal")
design

summary(exprs(gse))                    ## replace "exprs(gse)" with "expr_matrix" if you vae selected samples from data

## calculate median expression level
cutoff <- median(exprs(gse))          ## replace "exprs(gse)" with "expr_matrix" if you vae selected samples from data

## TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed <- exprs(gse) > cutoff    ## replace "exprs(gse)" with "expr_matrix" if you vae selected samples from data

## Identify genes expressed in more than 5samples

keep <- rowSums(is_expressed) > 5

## check how many genes are removed / retained.
table(keep)

## subset to just those expressed genes
gse <- gse[keep,]

## The lmFit function is used to fit the model to the data. 
##The result of which is to estimate the expression level in each of the groups that we specified.
fit <- lmFit(exprs(gse), design)            ## replace "exprs(gse)" with "expr_matrix" if you vae selected samples from data
head(fit$coefficients)                      

# Define contrast matrix
contrast.matrix <- makeContrasts(Tumour - Normal, levels=design)

## can define multiple contrasts
## e.g. makeContrasts(Group1 - Group2, Group2 - Group3,....levels=design)

# Fit contrasts and apply empirical Bayes moderation
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

##check top list
topTable(fit2)

## If you want to see results for other contrasts
topTable(fit2, coef=1)

### to see the results of the second contrast (if it exists)
# topTable(fit2, coef=2)
#----------------------------------------------------------------------------------------------
## If we want to know how many genes are differentially-expressed overall we can use the decideTests function.
decideTests(fit2)
table(decideTests(fit2))

## Coping with outliers  {Ignore this step}
#----------------------------------------------------------------------------------------------
## It is tempting to discard any arrays which seem to be outliers prior to differential expressions.
## calculate relative array weights
# aw <- arrayWeights(exprs(gse),design)
# aw
# ## The lmFit function can accept weights, and the rest of the code proceeds as above.
# fit <- lmFit(exprs(gse), design,
#              weights = aw)
# contrasts <- makeContrasts(Normal - Tumour, levels=design)
# fit2 <- contrasts.fit(fit, contrasts)
# fit2 <- eBayes(fit2)
#----------------------------------------------------------------------------------------------
### Further processing and visualization of DE results
#----------------------------------------------------------------------------------------------
anno <- fData(gse)
anno

# Select specific columns from the annotation data
anno <- anno[, c("Gene Title", "ENTREZ_GENE_ID", "Gene Symbol")]

# Add annotation data to 'fit2' object
fit2$genes <- anno

# Retrieve top table from 'fit2' object
topTable(fit2)
# Retrieve full results
full_results <- topTable(fit2, number = Inf)
##----------------------------------------------------------------------------------------------------------------------------------
##                                                 AVERAGE EXPRESSION
##----------------------------------------------------------------------------------------------------------------------------------
# First Check the sample name
sampleInfo$group

normal_samples <- filter(sampleInfo, group == "Tissue:normal prostate tissue adjacent to tumor")
cancer_samples <- filter(sampleInfo, group == "Tissue: primary prostate tumor")
# 
# # Get the indices of normal and cancer samples
normal_indices <- match(normal_samples$samples, colnames(exprs(gse)))    ## replace "exprs(gse)" with "expr_matrix" if you vae selected samples from data
cancer_indices <- match(cancer_samples$samples, colnames(exprs(gse)))    ## replace "exprs(gse)" with "expr_matrix" if you vae selected samples from data
# 
# # Calculate the average expression for normal and cancer samples
normal_avg_expression <- rowMeans(exprs(gse)[, normal_indices])          ## replace "exprs(gse)" with "expr_matrix" if you vae selected samples from data
cancer_avg_expression <- rowMeans(exprs(gse)[, cancer_indices])          ## replace "exprs(gse)" with "expr_matrix" if you vae selected samples from data
# 
# # Convert normal_avg_expression and cancer_avg_expression to data frames
normal_avg_expression_df <- data.frame(avg_Normal_exp = normal_avg_expression)
normal_avg_expression_df <- normal_avg_expression_df %>%
  tibble::rownames_to_column(var = "ID")
cancer_avg_expression_df <- data.frame(avg_Cancer_exp = cancer_avg_expression)
cancer_avg_expression_df <- cancer_avg_expression_df %>%
  tibble::rownames_to_column(var = "ID")
# 
# # Merge normal_avg_expression_df, cancer_avg_expression_df, and features by ID
Avg_expression_mtrx <- merge(normal_avg_expression_df, cancer_avg_expression_df, by = "ID", all = F)

##----------------------------------------------------------------------------------------------------------------------------------
# Retrieve full results from 'fit2' object :
#----------------------------------------------------------------------------------------------
## make rownames_to_column as "ID"
full_results <- full_results %>%
  tibble::rownames_to_column(var = "ID")

   ## Also add the avg expression colums of samples
full_DEGS_results <- merge(full_results, Avg_expression_mtrx, by = "ID", all = T) ### If you want to add all samples colums, then replace the "Avg_expression_mtrx" with "expr_matrix"

# Write the full results to a CSV file
write.csv(full_DEGS_results, "full_DEGS_results.csv", row.names = FALSE)
#----------------------------------------------------------------------------------------------
# Filter results based on adjusted p-values and absolute log fold change

filtered_results <- subset(full_DEGS_results, adj.P.Val < 0.05 & abs(logFC) > 1)
write.csv(filtered_results, "filtered_results_file.csv", row.names = FALSE)


save.image("D:/3. Data-Analysis-TCGA/3. Validation with other independent dataset/PRAD/GSE6919/GSE6919.RData")
#----------------------------------------------------------------------------------------------
# ## Get the results for particular gene of interest
# filter(full_results, Gene.Symbol == "TP53")

# ## Get results for genes with TP53 in the name
# filter(full_results, grepl("TP53", Gene.Symbol))

# ## Get results for one chromosome
# filter(full_results, Chromosome==20)
#----------------------------------------------------------------------------------------------

#                                  Visualisation
#                                 HEATMAPS, VOLCANO
#----------------------------------------------------------------------------------------------
## VOLCANO
library(ggrepel)        
p_cutoff <- 0.05
fc_cutoff <- 1.5
topN <- 200

full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>%     ## rename col names if gave error ->:: colnames(full_results)[colnames(full_results) == "Gene Symbol"] <- "Gene_Symbol"
  mutate(Rank = 1:n(), Label = ifelse(Rank < topN, Gene_Symbol,"")) %>% 
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black")


## Heatmap
# Load necessary libraries
library(dplyr)
library(pheatmap)

# Set the number of top genes to visualize
topN <- 200             ### top 20, 100, 200 or Inf

# # Convert rownames to a column named "ID" in the full_results dataframe
#full_results <- full_results %>%
#tibble::rownames_to_column(var = "ID")

# Extract the top N gene IDs
ids_of_interest <- mutate(full_results, Rank = 1:n()) %>% 
  filter(Rank < topN) %>% 
  pull(ID)

# Extract the gene symbols for the top N genes
gene_names <- mutate(full_results, Rank = 1:n()) %>% 
  filter(Rank < topN) %>% 
  pull(Gene_Symbol)

## Get the rows corresponding to ids_of_interest and all columns
gene_matrix <- exprs(gse)[ids_of_interest,]             # expr_matrix

# Generate and plot the heatmap
pheatmap(gene_matrix,
         labels_row = gene_names,
         scale = "row")
#----------------------------------------------------------------------------------------------




