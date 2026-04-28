# Code for: Chemerin promotes obesity-accelerated orthodontic tooth movement by activating IRE1α-dependent ER stress in PDL fibroblasts

## Description
This repository contains R scripts used for re-analysis of publicly available transcriptome and single-cell RNA-seq datasets in the study: "Chemerin promotes obesity-accelerated orthodontic tooth movement by activating IRE1α-dependent ER stress in PDL fibroblasts".

The scripts perform:
- Differential expression, KEGG enrichment, and GSEA of transcriptome data from hPDL fibroblasts under intermittent compressive force (GSE112122)
- Single-cell RNA-seq analysis and visualization for mouse OTM (GSE287729) and human periodontitis (GSE171213)
- Identification of cell types expressing *Rarres2*, differential expression between *Rarres2*+ and *Rarres2*- cells, and AUCell scoring of ER stress pathway activity

## File structure
├── README.md
├── 1_Transcriptome_GSE112122_with_annotation.R # Transcriptome re-analysis (Fig. 4A–C)
├── 2_scRNA_OTM_mouse_GSE287729.R # Mouse OTM scRNA-seq analysis (Fig. 1G, 4K–M)
└── 3_scRNA_PD_human_GSE171213.R # Human periodontitis scRNA-seq analysis (Fig. 1G)


## Data availability
Raw data are publicly available at Gene Expression Omnibus (GEO):
- **GSE112122**: Intermittent compressive force on hPDL fibroblasts (transcriptome)
- **GSE287729**: Mouse orthodontic tooth movement model (scRNA-seq)
- **GSE171213**: Human periodontitis (scRNA-seq)

This repository only contains analysis scripts; users must download the data from GEO before running the code.

## System requirements & Dependencies
- **R** version 4.2.1
- **R packages**:
  - Bioconductor: DESeq2, clusterProfiler, org.Hs.eg.db, org.Mm.eg.db, scran, SingleCellExperiment, multtest, KEGGREST, DOSE
  - CRAN: Seurat (v4.3.0), dplyr, ggplot2, tidyr, patchwork, harmony, msigdbr, AUCell (v1.25.2), ggsignif, pheatmap, cowplot, ggrepel, ggpubr, ggraph, tidyverse, Rcpp, forcats, Matrix, readr, scales, stringr, reshape2, ggsci, GseaVis (GitHub)
- The scripts contain automatic installation commands (`if(!require(package))...`), but manual installation may be necessary if dependencies fail.

## Usage
1. **Download the datasets** from GEO and place the data files in the appropriate directory structure as indicated in each script.
2. **Set the working directory (setwd) paths** to match your local file locations. The current paths are relative (`./`); you must adjust them.
3. **Run the scripts** in the following order or independently:
   - `1_Transcriptome_GSE112122_with_annotation.R` – performs DESeq2 analysis, KEGG enrichment, GSEA, and violin plots for ER stress markers.
   - `2_scRNA_OTM_mouse_GSE287729.R` – processes mouse single-cell data, annotates cell types, generates UMAP plots, differential expression, and AUCell scoring.
   - `3_scRNA_PD_human_GSE171213.R` – processes human periodontitis single-cell data, generates UMAP and identifies RARRES2-expressing cells.

## Important notes
- The scripts have been translated from their original Chinese comments into English for public sharing.
- The ER stress gene set used for AUCell scoring is derived from the KEGG pathway "Protein processing in endoplasmic reticulum" (mmu04141).
- The violin plots in `1_Transcriptome_GSE112122_with_annotation.R` produce the UPR marker gene expression figures (Fig. 4C).
- **GseaVis** must be installed from GitHub: `devtools::install_github("junjunlab/GseaVis")`.

## License
This code is provided under the MIT License. You are free to use, modify, and distribute it with appropriate attribution.

## Contact
For questions or issues with the code, please contact:
- Yuanyuan Yin,  The Affiliated Stomatological Hospital of Chongqing Medical University, 
- Email: yinyuanyuan@hospital.cqmu.edu.cn
