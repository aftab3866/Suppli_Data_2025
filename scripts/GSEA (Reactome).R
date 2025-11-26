library ("ReactomePA")
library ("clusterProfiler")
library ("org.Hs.eg.db")
library ("plyr")
library ("splitstackshape")
library("readxl")

##################################
## Selecting DE proteins
##################################
data.human<- read.csv("PRAD-199.csv", header = TRUE, sep = ",")
#OR
data.human<- read_excel("199gebes_GSEA.xlsx", sheet = 4, col_names = TRUE)

input <- data.human
#-------------------------------------------------------------------------
uniprot2entrezID <- read.csv("D:/3. Data-Analysis-TCGA/Scripts/file/uniprotkb_fulllist_Hm_ETZ.csv")

attach.entregene <- function (x,uniprotkb_col) {
  output <- merge (x, uniprot2entrezID, by.x=uniprotkb_col,by.y="Gene_Names",all.x=TRUE,sort=F)
  output <- output [!is.na (output$ENTREZ_GENE),]
  return (output)
}
input <- attach.entregene(input,"Genes_P")

#------------------------ if you select yoy DE data from global list ------
# data.human.up   <- input [input$Pvalue <= 0.25 & input$Fold_change > 1,]
# data.human.down <- input [input$Pvalue <= 0.25 & input$Fold_change < 1,]
# input<- NULL
#-------------------------------------------------------------------------
input.sorted <- input [input$adj.P.Val <= 0.05,]

genelist <- input.sorted$logFC
names (genelist) <- input.sorted$ENTREZ_GENE
genelist <- genelist [!duplicated (genelist)]
genelist <- genelist [!duplicated (names (genelist))]
genelist <- sort (genelist, decreasing = TRUE)

#library(UniProt.ws)
#up <- UniProt.ws(taxId=9606)
set.seed(12345)
y <- gsePathway(genelist, nPerm=2000,
                pvalueCutoff=1, minGSSize = 2, maxGSSize = 700,
                eps = 1e-10,exponent = 1,
                pAdjustMethod="BH", by = "DOSE", verbose=FALSE)     ##by one of 'fgsea' or 'DOSE' 

res <- as.data.frame(y)

y.readable <- setReadable(y, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

res.readable <- as.data.frame(y.readable)
write.csv (res.readable,"P_gsea-199.csv", row.names = f)



#Mapping completed|||

########################################################## PLOTS ############################################################
#source("https://bioconductor.org/biocLite.R")
#biocLite("enrichplot")
library (enrichplot)
library (ggplot2)

heatplot (y.readable, foldChange = genelist)
pdf ("test3-all.pdf", width = 120, height = 30,useDingbats=FALSE)
pp <- heatplot (y.readable, foldChange = genelist, showCategory = 5000,)
print (pp)
dev.off()






