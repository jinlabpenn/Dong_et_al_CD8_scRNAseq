# Tumor-associated Microbiota Suppresses Anti-tumor Immunity Analysis

## Project Description
Analysis of single-cell RNA sequencing data from CD8 T cells isolated from germ-free (GF) and specific-pathogen-free (SPF) tumor-bearing mice to investigate how tumor-associated microbiota affects anti-tumor immunity through cDC1 dysfunction in lung adenocarcinoma.

## Author
QIANG DONG

## Dependencies
- **Core Analysis**: Seurat (v5.2.1), harmony (v1.2.3), ProjecTILs (v3.5.0)
- **Visualization**: ggplot2 (v3.5.1), ComplexHeatmap (v2.22.0), EnhancedVolcano (v1.13.2)
- **Data Manipulation**: tidyverse (v2.0.0), dplyr (v1.1.4), reshape2 (v1.4.4)
- **Other**: patchwork (v1.3.0), ggpubr (v0.6.0), Matrix (v1.6-4), cowplot (v1.1.3), limma (v3.62.1), rstatix (v0.7.2), RColorBrewer (v1.1-3), writexl (v1.5.1), conflicted

## Required Data Files (Direct Download Links)

| File | Description | Download |
| ---- | ----------- | -------- |
| KP_sc_data.RData | Single-cell RNA sequencing data | [📥 Download](https://www.dropbox.com/scl/fi/n89uzggv42zl8tpbhoujx/KP_sc_data.RData?rlkey=qq2m1ec2enqtza0z2asg9eni4&st=5hdx88r2&dl=1) |
| ref_TILAtlas_mouse_v1.rds | ProjecTILs reference dataset | [📥 Download](https://www.dropbox.com/scl/fi/3la672wfyzkj2fd3en31k/ref_TILAtlas_mouse_v1.rds?rlkey=p3da4v6he9bigfxb17hguvit7&st=pzob2xbg&dl=1) |

**Note**: The download links are set to direct download (dl=1). After downloading, place KP_sc_data.RData in the RData folder and ref_TILAtlas_mouse_v1.rds in the ProjecTILs folder.

## Project Setup
The analysis script automatically sets up the following:

1. Creates necessary directories (Figures, Tables, rdsData, RData, ProjecTILs)
2. Downloads required data files from Dropbox
3. Sets up a user-writeable R library path (R_libs folder)
4. Installs and loads all required packages
5. Resolves common function conflicts between packages

## Directory Structure
```
project-root/
├── main.R                 # Main analysis script
├── Figures/               # Output figures (PDFs)
├── Tables/                # Output data tables (CSVs)
├── rdsData/               # Processed Seurat objects (RDS files)
├── RData/                 # Raw data files
│   └── KP_sc_data.RData   # Initial scRNA-seq data
├── ProjecTILs/            # ProjecTILs reference data
│   └── ref_TILAtlas_mouse_v1.rds # Reference for ProjecTILs analysis
└── R_libs/                # Local R package library
```

## Analysis Workflow

### 1. Data Integration
- Merge data from multiple sequencing lanes
- Add hash tag oligo (HTO) information
- Annotate experimental groups (GF, SPF, Healthy)

### 2. Data Processing
- Quality control filtering (nFeature, nCount, percent.mt)
- Normalization and scaling
- Dimensionality reduction with PCA
- Batch correction using Harmony
- UMAP visualization and clustering

### 3. T and NK Cell Analysis
- Subset and classify T/NK cells
- Create density and enrichment plots
- Perform ProjecTILs analysis
- Calculate module scores

### 4. CD8 T Cell Analysis
- Focus on CD8 T cell subpopulations
- Analyze cluster distributions between GF and SPF
- Identify marker genes with balloon plots and heatmaps
- Feature visualization of exhaustion markers

## Expected Output Files

### Figures
- **Figure_1A_Densityplot.pdf** - T/NK cell density by experimental group
- **Figure_1A_enrichment_plot.pdf** - Cell population enrichment analysis
- **Figure_1C.pdf** - CD8 T cell UMAP visualization
- **Figure_1D.pdf** - CD8 T cell cluster distribution comparison
- **Figure_1E.pdf** - Marker gene balloon plot
- **Figure_S1A.pdf** - T and NK cell populations UMAP
- **Figure_S1D.pdf** - CD8 T cell marker gene heatmap
- **Figure_S1E.pdf** - Exhaustion marker feature plots

### Tables
- **Figure_1A_enrichment_analysis.csv** - Enrichment statistics
- **CD8_markers_complete.csv** - CD8 T cell cluster markers
- **CD8_heatmap_gene_list.csv** - Genes for heatmap visualization
- **CD8_selected_genes.csv** - Key genes highlighted in analyses

### Processed Data
- **KP_merged.rds** - Merged Seurat object
- **T_NK_cells_final.rds** - Processed T and NK cell data
- **CD8_T_cells_final.rds** - Processed CD8 T cell data

## Running the Analysis
1. Clone this repository
2. Download the required data files using the direct links above
3. Place the files in their respective directories
4. Open `main.R` in RStudio 
5. Run the script

## Troubleshooting
- If automatic downloads fail, use the direct download links above
- If package installation fails, check write permissions and consider using the local R_libs directory
- For ProjecTILs installation issues, refer to https://github.com/carmonalab/ProjecTILs

## Contact
For questions or issues: jinlabpenn@gmail.com
