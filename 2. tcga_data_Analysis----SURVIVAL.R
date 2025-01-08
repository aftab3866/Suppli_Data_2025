###D:\3.Integerated_Results for Manuscripts\Aftab\scripts

#### Survival Analysis

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

# NB: make sure you set the working directory of RStudio correctly

tcga_data = readRDS(file = "tcga_data.RDS")
limma_res = readRDS(file = "limma_res.RDS")

# extract clinical data
clinical = tcga_data@colData
dim(clinical)

# we are only interested in the "Primary solid Tumor" cases for survival
clin_df = clinical[clinical$definition == "Primary solid Tumor",
                   c("patient",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up",
                     "gender",
                     "tumor_stage")]

# Now we have a new dataframe, clin_df, containing only the information that is relevant to survival analysis. 
#In addition to gender, we have added vital_status (whether patient is alive or dead), 
#tumor_stage (from stage 1 to 4) and two important variables: days_to_death, 
#that is the number of days passed from the initial diagnosis to the patient’s death (clearly, this is only relevant for dead patients),
#and days_to_last_follow_up that is the number of days passed from the initial diagnosis to the last visit.

# create a new boolean variable that has TRUE for dead patients
# and FALSE for live patients
clin_df$deceased = clin_df$vital_status == "Dead"

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)

# show first 10 samples
head(clin_df)

## Kaplan-Meier plots
########################
#As a first step, we need to define a survival formula with the help of the Surv object.

#The survival package provides us with an object, Surv, to form a dependent variable out of the overall_survival and deceased information:
Surv(clin_df$overall_survival, clin_df$deceased)

# Now that the survival time has been tagged with the censoring, we can add the categorical independent variable gender, and effectively create a formula
Surv(clin_df$overall_survival, clin_df$deceased) ~ clin_df$gender

##We now have a survival formula that can be passed to the survfit function to fit the survival model, 
# and then to another function to produce the Kaplan-Meier plots.
# fit a survival model
fit = survfit(Surv(overall_survival, deceased) ~ gender, data=clin_df)
print(fit)

# we produce a Kaplan Meier plot
ggsurvplot(fit, data=clin_df)

# This Kaplan-Meier plot shows two very similar trends until almost the 2000-day mark, where females seem to have a worse survival probability. 
#  But is there a significant difference?

# The difference between two such “event curves” is best tested via P value
ggsurvplot(fit, data=clin_df, pval=T)

# Can we see the number of patients dying (or being “censored”) as Time increases? Indeed we can, with what is called the “at risk table

ggsurvplot(fit, data=clin_df, pval=T, risk.table=T, risk.table.col="strata")

#### Another question could be: how does tumor stage affect survival?
# remove any of the letters "a", "b" or "c", but only if they are at the end
# of the name, eg "stage iiia" would become simply "stage iii"
clin_df$tumor_stage = gsub("[abc]$", "", clin_df$tumor_stage)

# we remove those with stage "not reported", since they are unknown
clin_df[which(clin_df$tumor_stage == "not reported"), "tumor_stage"] = NA
# finally, we also remove those with tumor stage 4, since they are too few
clin_df[which(clin_df$tumor_stage == "stage iv"), "tumor_stage"] = NA
table(clin_df$tumor_stage)

# We can now fit a new survival model with the tumor stage groups (one to four, plus the “not reported”):
  
fit = survfit(Surv(overall_survival, deceased) ~ tumor_stage, data=clin_df)
# we can extract the survival p-value and print it
pval = surv_pvalue(fit, data=clin_df)$pval
print(pval)

# we produce a Kaplan-Meier plot from the fitted model
ggsurvplot(fit, data=clin_df, pval=T, risk.table=T)


#########################################################################################
# Gene expression and survival
#########################################################################################
# We already have the top differentially expressed genes, ordered by significance, in the limma_res$topGenes dataframe, so we just have to take the first one.

# let's extract the table of differential expression we got earlier
expr_df = limma_res$topGenes

# print the first row, to see the gene name, the logFC value and the p-value
print(expr_df[1, ])

# get the ensembl gene id of the first row
gene_id = expr_df[1, "ensembl_gene_id"]

# also get the common gene name of the first row
gene_name = expr_df[1, "external_gene_name"]

#  We now have selected a gene. Let’s visualize how much differentially expressed it is:
  
# visualize the gene expression distribution on the diseased samples (in black)
# versus the healthy samples (in red)
expr_diseased = d_mat[rownames(clin_df), gene_id]
expr_healthy = d_mat[setdiff(rownames(d_mat), rownames(clin_df)), gene_id]

boxplot(expr_diseased, expr_healthy,
        names=c("Diseased", "Healthy"), main="Distribution of gene expression")
# get the expression values for the selected gene
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]

# find the median value of the gene and print it
median_value = median(clin_df$gene_value)
print(median_value)

# divide patients in two groups, up and down regulated.
# if the patient expression is greater or equal to them median we put it
# among the "up-regulated", otherwise among the "down-regulated"
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")

# we can fit a survival model, like we did in the previous section
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)

# we can extract the survival p-value and print it
pval = surv_pvalue(fit, data=clin_df)$pval
print(pval)

# and finally, we produce a Kaplan-Meier plot
ggsurvplot(fit, data=clin_df, pval=T, risk.table=T, title=paste(gene_name))

