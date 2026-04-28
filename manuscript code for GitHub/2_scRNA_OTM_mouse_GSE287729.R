
#########Single-cell RNA-seq analysis: Mouse OTM dataset (GSE287729) ###########
# NOTE: Modify file paths according to your local directory structure.
# Load required packages
if(!require(BiocManager))install.packages("BiocManager")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))install.packages("mindr")
if(!require(tidyverse))install.packages("tidyverse")
if(!require(patchwork))install.packages("patchwork")
if(!require(R.utils))install.packages("R.utils")
if(!require(cowplot))install.packages("cowplot")
if(!require(colorspace))install.packages("colorspace")
if(!require(devtools))install.packages("devtools")
if(!require(ggthemes))install.packages("ggthemes")
if(!require(ggraph))install.packages("ggraph")
if(!require(readr))install.packages("readr")
if(!require(scales))install.packages("scales")
if(!require(msigdbr))install.packages("msigdbr")
if(!require(harmony))install.packages("harmony")
if(!require(enrichplot))install.packages("enrichplot")
if(!require(org.Mm.eg.db))BiocManager::install("org.Mm.eg.db")
if(!require(org.Hs.eg.db))BiocManager::install("org.Hs.eg.db")
if(!require(scran))BiocManager::install("scran", force = T)
if(!require(clusterProfiler))BiocManager::install("clusterProfiler")
if(!require(DOSE))BiocManager::install("DOSE")
if(!require(SingleCellExperiment))BiocManager::install("SingleCellExperiment")
if(!require(multtest))BiocManager::install("multtest")
if(!require(ggplot2))install.packages("ggplot2")
if(!require(colorspace))install.packages("colorspace")
if(!require(ggthemes))install.packages("ggthemes")
if(!require(ggraph))install.packages("ggraph")
if(!require(readr))install.packages("readr")
if(!require(scales))install.packages("scales")
if(!require(forcats))install.packages("forcats")
if(!require(Matrix))install.packages("Matrix")
if(!require(Rcpp))install.packages("Rcpp")
if(!require(rvcheck))install.packages("rvcheck")
if(!require(patchwork))install.packages("patchwork")
if(!require(rvcheck))install.packages("rvcheck")

# Set working directory
setwd('./Singlecell/GSE287729_OTM_mouse')
getwd( )

# Import data
OTM.data <- Read10X(data.dir = "./GSE287729_OTM_mouse/GSE287729_RAW/GSM8750486_OTM") 
CON.data <- Read10X(data.dir = "./GSE287729_OTM_mouse/GSE287729_RAW/GSM8750487_CON") 

# Create Seurat objects
CON <- CreateSeuratObject(counts = CON.data, project = "CON", min.cells = 50, min.features = 500)
CON@meta.data$Group <- "CON"
OTM <- CreateSeuratObject(counts = OTM.data, project = "OTM", min.cells = 50, min.features = 500)
OTM@meta.data$Group <- "OTM"

# Merge datasets
mouse_OTM.combined <- merge(CON, y = c(OTM), add.cell.ids = c("CON", "OTM"), project = "mouse_OTM.combined")
mouse_OTM.combined

# Calculate mitochondrial gene percentage
mouse_OTM.combined[["percent.mt"]] <- PercentageFeatureSet(mouse_OTM.combined, pattern = "^mt-") 

# Quality control filtering
mouse_OTM.combined <- subset(mouse_OTM.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &percent.mt < 10)   
ncol(as.data.frame(mouse_OTM.combined[["RNA"]]@counts))

# Normalize data
mouse_OTM.combined <- NormalizeData(mouse_OTM.combined, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable genes
mouse_OTM.combined <- FindVariableFeatures(mouse_OTM.combined, selection.method = "vst", nfeatures = 5000)


# Scale data
mouse_OTM.combined <- ScaleData(mouse_OTM.combined, features = rownames(mouse_OTM.combined))

# Perform PCA
mouse_OTM.combined <- RunPCA(mouse_OTM.combined, features = VariableFeatures(object = mouse_OTM.combined))
print(mouse_OTM.combined[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(mouse_OTM.combined, dims = 1:2, reduction = "pca")
DimPlot(mouse_OTM.combined, reduction = "pca")
DimHeatmap(mouse_OTM.combined, dims = 1:15, cells = 500, balanced = TRUE)

# UMAP and clustering (before batch correction)
mouse_OTM.combined <- FindNeighbors(mouse_OTM.combined, dims = 1:10) 
mouse_OTM.combined <- FindClusters(mouse_OTM.combined, resolution = 0.2)
mouse_OTM.combined <- RunUMAP(mouse_OTM.combined, dims = 1:10)#与上面的PC要一致
DimPlot(mouse_OTM.combined,reduction = "umap") + plot_annotation(title = "mouse_OTM.combined, before integration")
before <- DimPlot(mouse_OTM.combined,reduction = "umap") + plot_annotation(title = "mouse_OTM.combined, before integration")

# Batch correction using Harmony
mouse_OTM.harmony <- mouse_OTM.combined %>%
  RunHarmony(group.by.vars = 'orig.ident', plot_convergence = T)

# UMAP and clustering using Harmony embeddings
mouse_OTM.harmony <- mouse_OTM.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.1)
after <- DimPlot(mouse_OTM.harmony,reduction = "umap") + plot_annotation(title = "mouse_OTM.harmony, after integration")
after1 <- DimPlot(mouse_OTM.harmony,reduction = "umap", group.by = 'orig.ident') + plot_annotation(title = "mouse_OTM.harmony, after integration")
mouse_OTM.harmony <- SetIdent(mouse_OTM.harmony,value = "seurat_clusters")
DimPlot(mouse_OTM.harmony,label = T) + NoLegend()
after <- DimPlot(mouse_OTM.harmony,label = T) + NoLegend()
after


# Cell type annotation
names(mouse_OTM.harmony@meta.data)
unique(mouse_OTM.harmony$Group)
# Annotate clusters based on canonical markers
new.cluster.ids <- c("Neutrophils", "Fibroblasts", "B cells", "HSC", "Fibroblasts", 
                     "Macrophages", "Endothelial", "Erythroblasts", "Megakaryocytes", "T cells", 
                     "HSC", "Plasmacytoid dendritic cells", "Fibroblasts", "Schwann cells", "Fibroblasts")

names(new.cluster.ids) <- levels(mouse_OTM.harmony)
mouse_OTM.harmony <- RenameIdents(mouse_OTM.harmony, new.cluster.ids)
DimPlot(mouse_OTM.harmony, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
mouse_OTM.harmony$celltype <- Idents(mouse_OTM.harmony)
mouse_OTM.harmony$celltype.Group <- paste(mouse_OTM.harmony$celltype, mouse_OTM.harmony$Group, sep = "_")
Idents(mouse_OTM.harmony) <- "celltype"
DimPlot(mouse_OTM.harmony,label = T) + NoLegend()
DimPlot(mouse_OTM.harmony,label = T,   split.by = "Group") + NoLegend()


######Fig.1G, UMAP visualization with Rarres2 expression######
df_complete <- mouse_OTM.harmony@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(
    cluster = mouse_OTM.harmony@meta.data$seurat_clusters,
    celltype = mouse_OTM.harmony@meta.data$celltype
  )

rarres2_data <- FetchData(mouse_OTM.harmony, vars = "rna_Rarres2")
df_complete$Rarres2 <- rarres2_data$Rarres2

#Gradient coloring for Rarres2 expression
final_plot <- ggplot(df_complete, aes(x = UMAP_1, y = UMAP_2)) +
  stat_ellipse(aes(fill = mouse_OTM.harmony@meta.data$celltype), 
               geom = "polygon", linetype = 2, alpha = 0.15,
               linewidth = 0.5, show.legend = FALSE, level = 0.95) +
  geom_point(color = "grey90", size = 0.2, alpha = 0.1) +
  geom_point(aes(color = Rarres2), size = 0.6, alpha = 1) +
  scale_fill_manual(values = colors2) +
  scale_color_gradientn(
    colors = c("lightgrey", "blue", "purple", "red"),
    name = "RARRES2 Expression"
  ) +
  theme_classic() +
  labs(title = "Cell Types with RARRES2 Expression Overlay") +
  guides(color = guide_colorbar(barwidth = 0.3,   # 更窄的宽度
                                barheight = 3,     # 高度
                                title.position = "top",
                                title.hjust = 0.5)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 8))

print(final_plot)


label_positions <- df_complete %>%
  group_by(celltype = mouse_OTM.harmony@meta.data$celltype) %>%
  summarize(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = 'drop'
  )

final_plot_with_labels <- final_plot +
  geom_label(data = label_positions,
             aes(x = UMAP_1, y = UMAP_2, label = celltype),
             size = 3, alpha = 0.8, label.size = 0,
             fill = NA,color = "black") +
  labs(title = "Cell Types with RARRES2 Expression and Labels")

print(final_plot_with_labels)

pdf(file = 'celltype_UMAP_RARRES2_in_mouse_OTM.harmony2.pdf',width = 5.5,height = 3.3)
print(final_plot_with_labels)
dev.off()

####Fig.4K，UMAP without annotation and with Rarres2 expression####
#Fig.4K (left)
Idents(mouse_OTM.harmony) <- "celltype"
pdf(file = 'mouse_OTM.harmony_cluster2.pdf',width = 4,height = 3)
DimPlot(mouse_OTM.harmony,label = F) + NoLegend()
dev.off()

#Fig.4K (right)
pdf(file = 'mouse_OTM.harmony_Rarres2_3.pdf',width = 3,height = 3)
FeaturePlot(mouse_OTM.harmony, features = c("rna_Rarres2"))
dev.off()

#######Fig.4L，KEGG #######
# Differential expression between Rarres2+ and Rarres2- cells (all cells)
mouse_OTM_Rarres2 <- subset(x = mouse_OTM.harmony, subset = rna_Rarres2 > 0)
mouse_OTM_Rarres2_negative <- subset(x = mouse_OTM.harmony, subset = rna_Rarres2 > 0, invert = TRUE)
mouse_OTM.harmony$Rarres2_group  <- "YES"
mouse_OTM.harmony$Rarres2_group[Cells(mouse_OTM_Rarres2)] <- paste("Rarres2_YES")
mouse_OTM.harmony$Rarres2_group[Cells(mouse_OTM_Rarres2_negative)] <- paste("Rarres2_NO")
Idents(mouse_OTM.harmony) <- "Rarres2_group"
Rarres2.markers <- FindMarkers(mouse_OTM.harmony, ident.1 = "Rarres2_YES", ident.2 = "Rarres2_NO", min.pct = 0.1)
head(Rarres2.markers, n = 10) 
write.csv(Rarres2.markers,"GSE287729_mouse_OTM_all_cells_Rarres2.markers.csv", row.names = TRUE)
Rarres2.markers <- read.table("GSE287729_mouse_OTM_all_cells_Rarres2.markers.csv", sep=",", header = T, row.names = 1) #2750个基因
#Filter upregulated genes: avg_log2FC > 0 & pvalue < 0.05
mydeg <- read.table("GSE287729_mouse_OTM_all_cells_Rarres2.markers.csv", sep=",", header = T, row.names = 1)
mygene1 <- subset(mydeg, avg_log2FC > 0 & p_val_adj < 0.05)
mygene1 <- mygene1 %>% rownames()
mygene1[1:10]
gene.df <- bitr(mygene1,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                OrgDb = org.Mm.eg.db)
# Helper function for enrichment dot plot
erich2plot <- function(data4plot){
  library(ggplot2)
  data4plot <- data4plot[order(data4plot$qvalue,decreasing = F)[1:10],]
  data4plot$BgRatio<-
    apply(data4plot,1,function(x){
      as.numeric(strsplit(x[3],'/')[[1]][1])
    })/apply(data4plot,1,function(x){
      as.numeric(strsplit(x[4],'/')[[1]][1])
    })
  
  p <- ggplot(data4plot,aes(BgRatio,Description))
  p<-p + geom_point()
  
  pbubble <- p + geom_point(aes(size=Count,color=-1*log10(qvalue)))
  
  pr <- pbubble + scale_colour_gradient(low="blue",high="red") + 
    labs(color=expression(-log[10](qvalue)),size="observed.gene.count", 
         x="Richfactor", y="term.description",title="Enrichment Process")
  
  pr <- pr + theme_bw()
  pr
}
# KEGG enrichment analysis
ekegg <- enrichKEGG(unique(gene.df$ENTREZID), organism='mmu',
                    pvalueCutoff=0.05,pAdjustMethod='BH',qvalueCutoff=0.05,
                    minGSSize=10,maxGSSize=500,use_internal_data=F)
ekegg <- setReadable(ekegg,'org.Mm.eg.db','ENTREZID')
ekegg_result <- ekegg@result
erich2plot(ekegg@result)
pdf(file = 'GSE287729_mouse_OTM_all_cells_Rarres2_padj0.05.pdf',width = 7.8,height = 3.2)
erich2plot(ekegg@result)
dev.off()


####Fig.4M，AUCell analysis: ER stress signature####
setwd("./GSE287729_OTM_mouse")
getwd(   )
library(Seurat)
library(stringr)
library(dplyr)
library(future)
library(future.apply)
library(msigdbr)
library(clusterProfiler)
library(devtools)
library(harmony)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(org.Hs.eg.db)
library(tidyverse)
library(dplyr)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(DESeq2)
library(AUCell)
library(reshape2)
library(ggsci)
cors <- pal_igv()(15)

# Build complete ER stress gene set from KEGG mmu04141
if (!requireNamespace("KEGGREST", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("KEGGREST")
}
library(KEGGREST)

pathway_id <- "mmu04141"
gene_info <- KEGGREST::keggGet(pathway_id)
gene_list <- gene_info[[1]]$GENE
gene_ids <- gene_list[seq(1, length(gene_list), 2)]
gene_descriptions <- gene_list[seq(2, length(gene_list), 2)]

extract_gene_symbol <- function(description) {
  symbol <- gsub(";.*", "", description)
  symbol <- trimws(symbol)
  return(symbol)
}

gene_symbols <- sapply(gene_descriptions, extract_gene_symbol)

gene_df <- data.frame(
  gene_symbol = gene_symbols,
  function. = "ER stress",
  stringsAsFactors = FALSE
)
rownames(gene_df) <- gene_symbols

head(gene_df)
print(paste("Table dimensions:", nrow(gene_df), "rows", ncol(gene_df), "columns"))

ER_sig <- gene_df# Note: gene_df is created later from KEGG mmu04141, replace as needed

# Build gene rankings for AUCell
cells_rankings <- AUCell_buildRankings(mouse_OTM.harmony@assays$RNA@data, splitByBlocks = TRUE)
# Calculate AUC scores (top 10% expressed genes)
cells_AUC_ER <- AUCell_calcAUC(ER_sig, #geneSets or signatures
                               cells_rankings, 
                               aucMaxRank = nrow(cells_rankings)*0.1)

# Explore thresholds
cells_assignment <- AUCell_exploreThresholds(cells_AUC_ER, plotHist = TRUE, assign = TRUE)
dev.off()

# Add AUC scores to Seurat metadata
cells_AUC_ER <- as.numeric(getAUC(cells_AUC_ER))
mouse_OTM.harmony$cells_AUC_ER <- cells_AUC_ER #添加至metadata中
head(mouse_OTM.harmony@meta.data)

# Violin plot
VlnPlot(mouse_OTM.harmony,features = 'cells_AUC_ER', #features也可改为AUCell
        pt.size = 0, group.by = "Rarres2",col = cors) #按细胞类型分组
dev.off()

# Box plot
my_comparisons = list(c("Rarres2_NO", "Rarres2_YES"))
p1 <- ggboxplot(mouse_OTM.harmony@meta.data, x="Rarres2", y="cells_AUC_ER", width = 0.6, #按group分组
                color = "black",#轮廓颜色
                fill="Rarres2",#填充
                palette = cors,
                xlab = F, #不显示x轴的标签
                bxp.errorbar=T,#显示误差条
                bxp.errorbar.width=0.5, #误差条大小
                size=0.5, #箱型图边线的粗细
                outlier.shape=NA, #不显示outlier
                legend = "right") + 
  #ylim(0,0.15) + 
  #facet_wrap(~celltype,ncol=6,scales="free_y") + # 按照细胞类型分面
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test") # 添加t检验
dev.off()

pdf(file = 'AUCell_RARRES2_in_mouse_OTM.harmony3.pdf',width = 3.5,height = 3.3)
p1
dev.off()




