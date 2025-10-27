# Cancer Biomarker Discovery Pipeline

A comprehensive bioinformatics workflow for identifying and validating hub genes across multiple cancer types using TCGA data.

## Workflow Overview

### 1. Data Acquisition & Pre-processing
- **Data Sources**: TCGA datasets (BRCA, LUAD, COAD, PRAD)
- **Sample Types**: 'Primary Tumor' and 'Solid Tissue Normal'
- **Tools**: `TCGAbiolinks` R package
- **Pre-processing**: Read Count normalization, `calcNormFactors()`, `model.matrix()`

### 2. Differential Expression Analysis
- **Methods**: `lmFit()` and `eBayes()` from limma package
- **Thresholds**: 
  - `|logFC| > 1.5`
  - `padj < 0.05`
- **Output**: Total DEGs identification

### 3. Multi-Method Gene Discovery
- **DEG Analysis**: Limma-based differential expression
- **WGCNA Analysis**: TOM, MEs, MM for co-expression networks
- **Machine Learning**: ElasticNet pipeline using `glmnet` and `caret`
- **Integration**: Common genes identification across all methods

### 4. Network Analysis & Functional Enrichment
- **PPI Construction**: BioGrid and HINT databases
- **Clustering**: MCL-clustering (degree cut-off ≥ 3, inflation parameter ≥ 2)
- **Enrichment**: 
  - `gsePathway()` analysis
  - ClusterProfiler package
  - MSigDB, KEGG, EPC, Reactome databases

### 5. Hub Gene Identification & Validation
- **Hub Gene Selection**: Central nodes in PPI networks
- **Expression Validation**: TNMPlot tool
- **ROC Analysis**: Diagnostic power assessment using TCGAplot
- **Protein-Level Validation**: UALCAN database
- **Pan-Cancer Analysis**: 33 cancer types consistency check

### 6. Clinical & Immunological Analysis
- **Survival Analysis**: Cox Proportional Hazards Model
- **Clinical Features**: age, gender, tumor stage integration
- **Packages**: `survival`, `forestmodel`, `survminer`, `rms`
- **Immune Analysis**: 
  - Immune-correlation analysis
  - Immune-checkpoint inhibitors (ICIs) response
  - CEPIA3.0 for OS & PFI analysis

## Key R Packages
```r
TCGAbiolinks    # Data download
limma           # DEG analysis
ClusterProfiler # Enrichment analysis
WGCNA           # Co-expression networks
glmnet, caret   # Machine Learning
survival        # Survival analysis
TCGAplot        # Visualization
