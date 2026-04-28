
##########Single-cell RNA-seq analysis: Human periodontitis dataset (GSE171213) ##############
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
setwd('./Singlecell/Theranostics2022')
getwd( )

# Load pre-processed Seurat objects
hcS.big <- readRDS("hcS.big.rds")
pdS.big <- readRDS("pdS.big.rds")
hcS.big@meta.data$Group <- "HC"
pdS.big@meta.data$Group <- "PD"

# Merge datasets
combined_HCPD <- merge(hcS.big, y = c(pdS.big), add.cell.ids = c("hc", "pd"), project = "combined_HC_PD")
combined_HCPD
saveRDS(combined_HCPD,"combined_HCPD.rds")

head(colnames(combined_HCPD))
tail(colnames(combined_HCPD))

# Calculate mitochondrial gene percentage
combined_HCPD[["percent.mt"]] <- PercentageFeatureSet(combined_HCPD, pattern = "^MT-") 
head(combined_HCPD@meta.data,5)
VlnPlot(combined_HCPD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

plot1 <- FeatureScatter(combined_HCPD, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(combined_HCPD, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2)) 

# Quality control filtering
combined_HCPD <- subset(combined_HCPD, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 40)   
ncol(as.data.frame(combined_HCPD[["RNA"]]@counts))

# Normalize data
combined_HCPD <- NormalizeData(combined_HCPD, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable genes
combined_HCPD <- FindVariableFeatures(combined_HCPD, selection.method = "vst", nfeatures = 3000)
#取出前10个高变基因# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(combined_HCPD), 10)
top10

# Scale data
combined_HCPD <- ScaleData(combined_HCPD, features = rownames(combined_HCPD))

# Perform PCA
combined_HCPD <- RunPCA(combined_HCPD, features = VariableFeatures(object = combined_HCPD))
print(combined_HCPD[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(combined_HCPD, dims = 1:2, reduction = "pca")
DimPlot(combined_HCPD, reduction = "pca")
DimHeatmap(combined_HCPD, dims = 1:15, cells = 500, balanced = TRUE)

# Advanced PCA (JackStraw)
combined_HCPD <- JackStraw(combined_HCPD, num.replicate = 100)
combined_HCPD <- ScoreJackStraw(combined_HCPD, dims = 1:20)
JackStrawPlot(combined_HCPD, dims = 1:20)#可视化前20个PC
ElbowPlot(combined_HCPD) #肘型图

# UMAP and clustering (before batch correction)
combined_HCPD <- FindNeighbors(combined_HCPD, dims = 1:10) 
combined_HCPD <- FindClusters(combined_HCPD, resolution = 0.2)
combined_HCPD <- RunUMAP(combined_HCPD, dims = 1:10)#与上面的PC要一致
DimPlot(combined_HCPD,reduction = "umap") + plot_annotation(title = "combined_HCPD, before integration")
 
# Batch correction using Harmony
combined_HCPD <- combined_HCPD %>%
  RunHarmony(group.by.vars = 'orig.ident', plot_convergence = T)

# UMAP and clustering using Harmony embeddings
combined_HCPD <- combined_HCPD %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.2)
after <- DimPlot(combined_HCPD,reduction = "umap") + plot_annotation(title = "combined_HCPD, after integration")
after1 <- DimPlot(combined_HCPD,reduction = "umap", group.by = 'orig.ident') + plot_annotation(title = "combined_HCPD.harmony, after integration")

combined_HCPD <- SetIdent(combined_HCPD, value = "seurat_clusters")
DimPlot(combined_HCPD, label = T) + NoLegend()
after <- DimPlot(combined_HCPD,label = T) + NoLegend()
after


# Cell type annotation
names(combined_HCPD@meta.data)
unique(combined_HCPD$Group)
new.cluster.ids <- c("T cells", "T cells", "Endothelial", "Neutrophil", "Plasma cells", 
                     "B cells", "Monocytes", "Fibroblasts", "Mast cells", "Epithelial", 
                     "Neural progenitors", "Pericytes", "Others")

names(new.cluster.ids) <- levels(combined_HCPD)
combined_HCPD <- RenameIdents(combined_HCPD, new.cluster.ids)
DimPlot(combined_HCPD, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
combined_HCPD$celltype <- Idents(combined_HCPD)
combined_HCPD$celltype.Group <- paste(combined_HCPD$celltype, combined_HCPD$Group, sep = "_")
Idents(combined_HCPD) <- "celltype"
DimPlot(combined_HCPD,label = T) + NoLegend()
DimPlot(combined_HCPD,label = T,   split.by = "Group") + NoLegend()

Idents(combined_HCPD) <- "celltype"
pdf(file = 'VlnPlot_RARRES2_in_combined_HCPD.pdf',width = 5.2,height = 3.3)
VlnPlot(combined_HCPD, features = c("rna_RARRES2"), pt.size = 0.1)
dev.off()


#####Fig.1G,UMAP visualization with RARRES2 expression#####
df_complete <- combined_HCPD@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(
    cluster = combined_HCPD@meta.data$seurat_clusters,
    celltype = combined_HCPD@meta.data$celltype
  )

rarres2_data <- FetchData(combined_HCPD, vars = "RARRES2")
df_complete$RARRES2 <- rarres2_data$RARRES2

# Gradient coloring
final_plot <- ggplot(df_complete, aes(x = UMAP_1, y = UMAP_2)) +
  stat_ellipse(aes(fill = combined_HCPD@meta.data$celltype), 
               geom = "polygon", linetype = 2, alpha = 0.15,
               linewidth = 0.5, show.legend = FALSE, level = 0.95) +
  geom_point(color = "grey90", size = 0.2, alpha = 0.1) +
  geom_point(aes(color = RARRES2), size = 0.6, alpha = 1) +
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
  group_by(celltype = combined_HCPD@meta.data$celltype) %>%
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

pdf(file = 'celltype_UMAP_RARRES2_in_combined_HCPD 3.pdf',width = 6,height = 3.3)
print(final_plot_with_labels)
dev.off()
