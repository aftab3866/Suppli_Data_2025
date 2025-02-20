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
This script is tailored for analyzing data from the Gene Expression Omnibus (GEO). Its key features include:  
- Importing and normalizing GEO datasets.  
- Identifying DEGs across experimental conditions.  
- Integrating GEO data analysis with findings from TCGA datasets for broader insights.  

## 4. **WGCNA-final script__Tested for DEGs profile__Perfectly working.R**
This is a comprehensive script for conducting Weighted Gene Co-expression Network Analysis (WGCNA). Highlights include:  
- Constructing co-expression networks using DEGs profiles.  
- Identifying gene modules correlated with clinical traits or phenotypes.  
- Exporting WGCNA results for further pathway and topological analysis.  
- Ensured to be tested and optimized for accuracy and reproducibility.  

# Supplementary Data and Figures

### Supplementary Data
- **Suppli-1 data (DEGs_TCGA):**  
  Information about DEGs (Differentially Expressed Genes) from selected cancer types.
  
- **Suppli-2 data (WGCNA result combined):**  
  Comprehensive results of the weighted gene co-expression network analysis.

- **Suppli-3 data (TOPO+PATHWAYS):**  
  PPI (Protein-Protein Interaction) network along with its topological properties.

- **Suppli-4 data (Integrated Drug-Target):**  
  Details about drug-target interactions from various databases.

- **Suppli-5 data (Detailed Clusters info):**
- A comprehensive discussion of all 11 functional clusters, along with relevant references, is presented in


### Cytoscape files
- **Drug-Target Interaction Network.cys:**  
  Cytoscape file containing the Drug-Target Interaction Network.

- **PPI network and clusters.cys:**  
  Cytoscape file with PPI Network and cluster information.

---

### Supplementary Figures
- **Supplementary Figure 1:**  
  A plot showing the Principal Component (PC) analysis for sample variation in selected cancer types (SCT) and a volcano plot displaying the differentially expressed genes in SCT.

- **Supplementary Figure 2:**  
  Common differentially expressed genes across SCT.

- **Supplementary Figure 3:**  
  Machine learning model using a linear model with feature selection followed by an Elastic Net model. Limma also selected genes highlighted in green. Samples in red represent Solid Tissue Normal, while samples in black represent Primary Solid Tumor.

- **Supplementary Figure 4:**
  The expression patterns of all key genes (26 genes × 4 cancer types) across Tumor and Normal tissues were analyzed using TNMplot. The heatmap shows that hub genes exhibit higher expression in metastatic tumors compared to primary tumors in colon, breast, and prostate cancers.

- **Supplementary Figure 5:**  
  A figure showing drug-target interactions of non-approved drugs (right side) alongside approved drugs (left side).

- **Supplementary Figure 6:**  
  This figure illustrates the expression patterns of nine genes across 33 cancer types. The heatmap highlights similar expression trends among different cancers, while the bar plot depicts consistent expression patterns of each gene across multiple cancer types.


