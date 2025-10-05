# Explanation of Analysis Scripts

## 1. **tcga_data_Analysis.R**
This script performs data analysis using TCGA (The Cancer Genome Atlas) datasets. Key functionalities include:  
- Loading and pre-processing TCGA data.  
- Identifying and visualizing Differentially Expressed Genes (DEGs).  
- Generating summary statistics and exploratory plots.  

## 2. **tcga_data_Analysis----SURVIVAL.R**
This script extends the analysis of TCGA datasets by incorporating survival analysis. It includes:  
- Kaplan-Meier survival plots to assess the impact of specific gene expressions on patient survival.  
- Cox proportional hazards regression models for advanced survival analysis.  

## 3. **Analysing data from GEO.R**
This script is tailored to analyze data from the Gene Expression Omnibus (GEO). Its key features include:  
- Importing and normalizing GEO datasets.  
- Identifying DEGs across experimental conditions.  
- Integrating GEO data analysis with findings from TCGA datasets for broader insights.  

## 4. **WGCNA-final script__Tested for DEGs profile__Perfectly working.R**
This is a comprehensive script for conducting Weighted Gene Co-expression Network Analysis (WGCNA). Highlights include:  
- Constructing co-expression networks using DEGs profiles.  
- Identifying gene modules correlated with clinical traits or phenotypes.  
- Exporting WGCNA results for further pathway and topological analysis.  
- Ensured to be tested and optimized for accuracy and reproducibility.

  ## 5. PAN-CANCER ANALYSIS TCGAPLOT.R**
 This script allows analysis of any single gene or multiple genes of interest for expression patterns across pan-cancer datasets. It can also perform many other important tasks. For details, visit: https://github.com/tjhwangxiong/TCGAplot
 [cited: PMID: 38105215]

  ## 6. TCGA_GO&GSEA.R**
 This script allows GSEA-GO to be grouped by the expression of a specific gene. For details, visit: https://github.com/tjhwangxiong/TCGAplot
 [cited: PMID: 38105215]

  ## 7-8. PART-1.1__TCGA-OS-COX_HR---## Perfreable (gene by gene).R** &&&& PART-1.2__TCGA-OS-COX_HR---## Perfreable ( use when satge and gender not defined).R
 This script evaluates the prognostic value of candidate genes while adjusting for clinical factors and creates publication-ready visualizations and reports
 
# Supplementary Data and Figures

### Supplementary Data
- **1. Suppli-1 data ( total_DEGs and Common):**  
  This section presents information on differentially expressed genes (DEGs) identified in selected cancer types. (A) It includes a list of common DEGs and a comparison with an independent dataset (GEO), visualized through a heatmap of log2 fold changes and a scatter plot showing the correlation between the TCGA and GEO datasets. (B) cross validated with Promics data. (C) Total no. of DEGs across SCTs. (D) List of Shared DEGs among the SCTs. (E) Geneset Enrichment Analysis (GSEA). (F) Comparison expression in P. Tumoro vs Metastases across SCTs.
  
- **2. Suppli-2 data (WGCNA result combined):**  
  Comprehensive results of the weighted gene co-expression network analysis.

- **3. Suppli-3 data (TOPO+PATHWAYS):**  
  (Sheet-1).PPI (Protein-Protein Interaction) network along with its topological properties. (Sheet 3-4)Inferred genes list and their Network topological properties.(Sheet-4) Detailed information about 11 complexities and their pathway enrichment analysis (ORA)

- **4. Suppli-4 data (Cluster info):**  
   Full details of the hub genes, their interactions, and functional roles within each of the 11 clusters.

- **5. Suppli-4 data (Immunotherapy response):**  
   Full details on the hub genes and their corresponding immunotherapy response, including predictive values and expression changes.

### Cytoscape files
- **PPI network and clusters.cys:**  
  Cytoscape file with PPI Network and cluster information.

### Supplementary Figures

  - **Supplementary Figure 1:**  
  A plot showing the Principal Component (PC) analysis for sample variation in selected cancer types (SCT) and a volcano plot displaying the differentially expressed genes in SCT.

- **Supplementary Figure 2:**  
  Commonly differentially expressed genes across SCT.

- **Supplementary Figure 3:**  
  Machine learning model using a linear model with feature selection, followed by an Elastic Net model. Limma also selected genes highlighted in green. Samples in red represent Solid Tissue Normal, while samples in black represent Primary Solid Tumor.

- **Supplementary Figure 4:**
  Gene set enrichment analysis (GSEA) mapped onto Reactome pathways, detailing the functional annotation of all differentially expressed genes from four cancer types. The colour of the nodes represents the normalized enrichment score (NES); positive NES values indicate increased pathway activity, while negative NES values indicate decreased activity. Each node represents a pathway. The significance of each path is highlighted by the node border, with darker colours indicating greater significance. Grey nodes indicate irrelevant pathways.

- **Supplementary Figure 5:**
  The expression patterns of all key genes (26 genes × 4 cancer types) across Tumor and Normal tissues were analyzed using TNMplot.

- **Supplementary Figure 6:**  
  This figure illustrates the Pan-Cancer Expression Landscape of 26 signature genes across 33 cancer types. The box plot highlights the expression pattern for significantly Up and down-regulated signature genes in 33 Types of cancer.


