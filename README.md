{\rtf1\ansi\ansicpg936\cocoartf2580
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fnil\fcharset0 HelveticaNeue;\f1\fnil\fcharset134 PingFangSC-Regular;\f2\froman\fcharset0 TimesNewRomanPSMT;
\f3\fnil\fcharset134 STSongti-SC-Regular;}
{\colortbl;\red255\green255\blue255;\red13\green14\blue17;\red255\green255\blue255;\red0\green0\blue0;
}
{\*\expandedcolortbl;;\cssrgb\c5882\c6667\c8235;\cssrgb\c100000\c100000\c100000;\cssrgb\c0\c0\c0;
}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # Code for: Chemerin promotes obesity-accelerated orthodontic tooth movement by activating IRE1\uc0\u945 -dependent ER stress in PDL fibroblasts\
\
## Description\
This repository contains R scripts used for re-analysis of publicly available transcriptome and single-cell RNA-seq datasets in the study: "Chemerin promotes obesity-accelerated orthodontic tooth movement by activating IRE1\uc0\u945 -dependent ER stress in PDL fibroblasts".\
\
The scripts perform:\
- Differential expression, KEGG enrichment, and GSEA of transcriptome data from hPDL fibroblasts under intermittent compressive force (GSE112122)\
- Single-cell RNA-seq analysis and visualization for mouse OTM (GSE287729) and human periodontitis (GSE171213)\
- Identification of cell types expressing *Rarres2*, differential expression between *Rarres2*+ and *Rarres2*- cells, and AUCell scoring of ER stress pathway activity\
\
## File structure\
\pard\pardeftab720\partightenfactor0

\f1 \cf2 \cb3 \expnd0\expndtw0\kerning0
\'a9\'c0
\f0 \uc0\u9472 \u9472  README.md\

\f1 \'a9\'c0
\f0 \uc0\u9472 \u9472  1_Transcriptome_GSE112122_with_annotation.R # Transcriptome re-analysis (Fig. 4A\'96C)\

\f1 \'a9\'c0
\f0 \uc0\u9472 \u9472  2_scRNA_OTM_mouse_GSE287729.R # Mouse OTM scRNA-seq analysis (Fig. 1G, 4K\'96M)\
\uc0\u9492 \u9472 \u9472  3_scRNA_PD_human_GSE171213.R # Human periodontitis scRNA-seq analysis (Fig. 1G)\
\
\
## Data availability\
Raw data are publicly available at Gene Expression Omnibus (GEO):\
- **GSE112122**: Intermittent compressive force on hPDL fibroblasts (transcriptome)\
- **GSE287729**: Mouse orthodontic tooth movement model (scRNA-seq)\
- **GSE171213**: Human periodontitis (scRNA-seq)\
\
This repository only contains analysis scripts; users must download the data from GEO before running the code.\
\
## System requirements & Dependencies\
- **R** version 4.2.1\
- **R packages**:\
  - Bioconductor: DESeq2, clusterProfiler, org.Hs.eg.db, org.Mm.eg.db, scran, SingleCellExperiment, multtest, KEGGREST, DOSE\
  - CRAN: Seurat (v4.3.0), dplyr, ggplot2, tidyr, patchwork, harmony, msigdbr, AUCell (v1.25.2), ggsignif, pheatmap, cowplot, ggrepel, ggpubr, ggraph, tidyverse, Rcpp, forcats, Matrix, readr, scales, stringr, reshape2, ggsci, GseaVis (GitHub)\
- The scripts contain automatic installation commands (`if(!require(package))...`), but manual installation may be necessary if dependencies fail.\
\
## Usage\
1. **Download the datasets** from GEO and place the data files in the appropriate directory structure as indicated in each script.\
2. **Set the working directory (setwd) paths** to match your local file locations. The current paths are relative (`./`); you must adjust them.\
3. **Run the scripts** in the following order or independently:\
   - `1_Transcriptome_GSE112122_with_annotation.R` \'96 performs DESeq2 analysis, KEGG enrichment, GSEA, and violin plots for ER stress markers.\
   - `2_scRNA_OTM_mouse_GSE287729.R` \'96 processes mouse single-cell data, annotates cell types, generates UMAP plots, differential expression, and AUCell scoring.\
   - `3_scRNA_PD_human_GSE171213.R` \'96 processes human periodontitis single-cell data, generates UMAP and identifies RARRES2-expressing cells.\
\
## Important notes\
- The scripts have been translated from their original Chinese comments into English for public sharing.\
- The ER stress gene set used for AUCell scoring is derived from the KEGG pathway "Protein processing in endoplasmic reticulum" (mmu04141).\
- The violin plots in `1_Transcriptome_GSE112122_with_annotation.R` produce the UPR marker gene expression figures (Fig. 4C).\
- **GseaVis** must be installed from GitHub: `devtools::install_github("junjunlab/GseaVis")`.\
\
## License\
This code is provided under the MIT License. You are free to use, modify, and distribute it with appropriate attribution.\
\
## Contact\
For questions or issues with the code, please contact:\
- Yuanyuan Yin,  
\f2\fs26\fsmilli13333 \cf4 \cb1 The Affiliated Stomatological Hospital of Chongqing Medical University,
\f3\fs24  
\f0 \cf2 \cb3 \
- Email: yinyuanyuan@hospital.cqmu.edu.cn\
}