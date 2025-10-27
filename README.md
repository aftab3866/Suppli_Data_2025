# 🧬 Cancer Biomarker Discovery Pipeline

A comprehensive bioinformatics workflow for identifying and validating hub genes across multiple cancer types using TCGA data.

![Workflow Diagram](workflow.png)
*Figure 1: Comprehensive analysis workflow from data acquisition to clinical validation*

## 📋 Overview

This pipeline integrates multi-omics data analysis, machine learning, and network biology to identify robust cancer biomarkers with clinical relevance across 33 cancer types.

## 🔄 Workflow Steps

### 1. 📊 Data Acquisition & Pre-processing
- **Data Sources**: TCGA datasets (BRCA, LUAD, COAD, PRAD)
- **Sample Types**: 'Primary Tumor' and 'Solid Tissue Normal'
- **Tools**: `TCGAbiolinks` R package
- **Pre-processing**: Read Count normalization using `calcNormFactors()`, `model.matrix()`

### 2. 🔍 Differential Expression Analysis
- **Methods**: `lmFit()` and `eBayes()` from limma package
- **Thresholds**: 
  - `|logFC| > 1.5`
  - `padj < 0.05`
- **Output**: Identification of significantly differentially expressed genes

### 3. 🎯 Multi-Method Gene Discovery
| Method | Approach | Key Features |
|--------|----------|--------------|
| **DEG Analysis** | Limma-based | Traditional differential expression |
| **WGCNA** | Co-expression networks | TOM, MEs, module membership |
| **Machine Learning** | ElasticNet pipeline | Feature selection with `glmnet` and `caret` |

### 4. 🌐 Network Analysis & Functional Enrichment
- **PPI Construction**: BioGrid and HINT databases
- **Clustering**: MCL-clustering (degree cut-off ≥ 3, inflation parameter ≥ 2)
- **Pathway Analysis**: 
  - `gsePathway()` enrichment
  - ClusterProfiler package
  - MSigDB, KEGG, EPC, Reactome databases

### 5. ✅ Hub Gene Validation
- **Expression Level**: TNMPlot validation
- **Protein Level**: UALCAN database
- **Diagnostic Power**: ROC analysis with TCGAplot
- **Pan-Cancer**: Consistency across 33 TCGA cancer types

### 6. 🏥 Clinical & Immunological Analysis
- **Survival Analysis**: Cox Proportional Hazards Model
- **Clinical Integration**: age, gender, tumor stage
- **Immune Profiling**: 
  - Immune-correlation analysis
  - Immune-checkpoint inhibitors response
  - CEPIA3.0 for overall and progression-free survival

## 🛠️ Installation & Usage

```bash
# Clone repository
git clone https://github.com/yourusername/cancer-biomarker-workflow.git
cd cancer-biomarker-workflow

# Install required R packages
Rscript install_dependencies.R
