

#################################Transcriptome data re-analysis: GSE112122 #################################
setwd('./Transcriptome/GSE112122_PDL_OTM/2_Intermittent')
getwd( )
library(dplyr)

#Read probe-to-gene mapping file
probe2gene <- read.table("probe2gene.csv", sep=",", header = T)

#Read raw count matrix
probes_expr <- read.csv("GSE112122_raw_counts_Intermittent.csv", header = TRUE, row.names = 1)
probes_expr <- as.matrix(probes_expr)
dim(probes_expr)
max(probes_expr)

#Aggregate duplicate gene symbols by mean expression
filterEM2 <- function(probes_expr, probe2gene,method="mean"){
  colnames(probe2gene) <- c("probeid","symbol")
  probe2gene$probeid=as.character(probe2gene$probeid)
  probe2gene$symbol=trimws(probe2gene$symbol)
  # head(probe2gene)
  
  message(paste0('input expression matrix is ',nrow(probes_expr),' rows(genes or probes) and ',ncol(probes_expr),' columns(samples).\n'))
  message(paste0('input probe2gene is ',nrow(probe2gene),' rows(genes or probes)\n'))
  
  probe2gene=na.omit(probe2gene)
  # if one probe mapped to many genes, we will only keep one randomly.
  probe2gene=probe2gene[!duplicated(probe2gene$probeid),]
  # 这个地方是有问题的，随机挑选一个注释进行后续分析。
  probe2gene = probe2gene[probe2gene$probeid %in% rownames(probes_expr),]
  
  message(paste0('after remove NA or useless probes for probe2gene, ',nrow(probe2gene),' rows(genes or probes) left\n'))
  
  #probes_expr <- exprs(eSet);dim(probes_expr)
  probes_expr <- as.data.frame(probes_expr)
  genes_expr <- tibble::rownames_to_column(probes_expr,var = "probeid") %>%
    merge(probe2gene,.,by= "probeid")
  genes_expr <- genes_expr[-1]
  ##remove duplicates symbol:method = mean, also median, max ,min
  #https://www.jingege.wang/2021/08/28/geo/
  message("remove duplicate symbols, it will take a while, be patient!")
  genes_expr<-aggregate(x=genes_expr[,2:ncol(genes_expr)],by=list(genes_expr$symbol),FUN=method,na.rm=T) 
  genes_expr<- tibble::column_to_rownames(genes_expr,var = "Group.1")
  
  message(paste0('output expression matrix is ',nrow(genes_expr),' rows(genes or probes) and ',ncol(genes_expr),' columns(samples).'))
  return(genes_expr)
}
genes_expr <- filterEM2(probes_expr, probe2gene )

#Remove genes with mean count <= 100 across all samples
filtered_genes_expr = genes_expr[rowMeans(genes_expr)>100,] 
#Round counts to integers (required by DESeq2)
filtered_genes_expr <- round(filtered_genes_expr)
write.csv(filtered_genes_expr, file = "genes_expr_GSE112122_raw_counts_Intermittent.csv", row.names = T)

#Read group information
group_list <- read.csv("group_list.csv", stringsAsFactors=T)
#Check sample order alignment
colnames(filtered_genes_expr) ==group_list$X

#Differential expression analysis using DESeq2
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = filtered_genes_expr, 
                              colData = group_list, 
                              design= ~ group)
dds$group <- relevel(dds$group, ref = "Control_Intermittent")

dds <- DESeq(dds)
res <- results(dds)
head(res)
class(res)
res_1 <- data.frame(res)
res_1 %>% 
  mutate(group = case_when(
    log2FoldChange >= 0.5 & padj <= 0.05 ~ "UP",
    log2FoldChange <= -0.5 & padj <= 0.05 ~ "DOWN",
    TRUE ~ "NOT_CHANGE"
  )) -> res_2

table(res_2$group)
write.csv(res_2, file = "diff_result_GSE112122_Intermittent.csv", row.names = T, quote = F)#保存文件

#####Fig4A，KEGG#####
# Filter upregulated genes: log2FoldChange > 0.5 and padj < 0.05
mydeg <- read.table("diff_result_GSE112122_Intermittent.csv", sep=",", header = T, row.names = 1)
mygene1 <- subset(mydeg, log2FoldChange > 0.5 & padj < 0.05)
mygene1 <- mygene1 %>% rownames()

mygene1[1:10]
gene.df <- bitr(mygene1,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                OrgDb = org.Hs.eg.db)

# Helper function for enrichment dot plot
erich2plot <- function(data4plot){
  library(ggplot2)
  data4plot <- data4plot[order(data4plot$qvalue,decreasing = F)[1:15],]
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
ekegg <- enrichKEGG(unique(gene.df$ENTREZID), organism='hsa',
                    pvalueCutoff=0.05,pAdjustMethod='BH',qvalueCutoff=0.05,
                    minGSSize=10,maxGSSize=500,use_internal_data=F)
ekegg <- setReadable(ekegg,'org.Hs.eg.db','ENTREZID')
ekegg_result <- ekegg@result
erich2plot(ekegg@result)
pdf(file = 'KEGG_OTM_vs_Con_log2FoldChange0.5且padj0.05小.pdf',width = 5.8,height = 3)
erich2plot(ekegg@result)
dev.off()

#####Fig4B，GSEA of unfolded protein response#####
# GSEA analysis
library(GSEABase) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)
library(enrichplot)
library(msigdbr)

diff_result_GSE112122_Intermittent <- read.csv("diff_result_GSE112122_Intermittent.csv", row.names=1)
deg <- diff_result_GSE112122_Intermittent
ge = deg$log2FoldChange
names(ge) = rownames(deg)
ge = sort(ge,decreasing = T)
head(ge)

# Obtain human Hallmark gene sets from MSigDB
geneset <- msigdbr(species = "Homo sapiens",category = "H") %>%
  dplyr::select(gs_name,gene_symbol)
geneset[1:4,]

# Format gene set names
geneset$gs_name = geneset$gs_name %>%
  str_split("_",simplify = T,n = 2)%>%
  .[,2]%>%
  str_replace_all("_"," ") %>% 
  str_to_sentence()

# Run GSEA
em <- GSEA(ge, TERM2GENE = geneset)

# Plot GSEA results
gseaplot2(em, geneSetID = 9, pvalue_table = T, title = em$Description[9])

# GseaVis for visualization
devtools::install_github("junjunlab/GseaVis")
library(GseaVis)
gseaNb(object = em,
       geneSetID = 'Unfolded protein response')

mygene <- c("XBP1","ERN1","HSPA5")
gseaNb(object = em,
       geneSetID = 'Unfolded protein response',
       addGene = mygene,
       subPlot = 2,
       geneSize = 4,
       addPval = T)
pdf(file = 'GSEA_OTM_vs_Con_Unfolded protein response.pdf',width = 5.2,height = 3.3)
dev.off()




####Fig.4C  Violin plots of UPR marker genes (Ire1a, Xbp1, and Bip)####
library(ggplot2)
library(ggsignif)
library(DESeq2)

# Read differential expression results
diff_result <- read.csv("diff_result_GSE112122_Intermittent.csv", row.names = 1)

# Read expression matrix and group information
filtered_genes_expr <- read.table("genes_expr_GSE112122_raw_counts_Intermittent.csv", 
                                  sep = ",", header = TRUE, row.names = 1)
group_list <- read.csv("group_list.csv", stringsAsFactors = TRUE)

# Normalize using DESeq2
dds <- DESeqDataSetFromMatrix(countData = filtered_genes_expr,
                              colData = group_list,
                              design = ~ group)
dds$group <- relevel(dds$group, ref = "Control_Intermittent")
dds <- DESeq(dds)

# Extract normalized count data
vsd <- vst(dds, blind = FALSE)
expr_norm <- assay(vsd)

# Generate violin plot for XBP1 and save as PDF
if("XBP1" %in% rownames(expr_norm)) {
  plot_data_xbp1 <- data.frame(
    expression = expr_norm["XBP1", ],
    group = group_list$group
  )
  
  gene_diff_xbp1 <- diff_result["XBP1", ]
  padj_xbp1 <- gene_diff_xbp1$padj
  
  significance_xbp1 <- ifelse(padj_xbp1 < 0.001, "***",
                              ifelse(padj_xbp1 < 0.01, "**",
                                     ifelse(padj_xbp1 < 0.05, "*", "ns")))
  
  p_xbp1 <- ggplot(plot_data_xbp1, aes(x = group, y = expression, fill = group)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    labs(title = "Expression of XBP1",
         y = "Normalized Expression",
         x = "") +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.5)
    ) +
    scale_fill_manual(values = c("lightblue", "lightcoral"))
  
  if(significance_xbp1 != "ns") {
    p_xbp1 <- p_xbp1 + geom_signif(
      comparisons = list(c("Control_Intermittent", "Intermittent")),
      annotations = significance_xbp1,
      y_position = max(plot_data_xbp1$expression) * 1.05,
      tip_length = 0.01,
      textsize = 4
    )
  }
  
  ggsave("violin_plot_XBP1.pdf", plot = p_xbp1, width = 2.5, height = 3)
  print("XBP1 violin plot saved as PDF")
} else {
  message("Gene XBP1 not found in expression matrix")
}

# Generate violin plot for ERN1 and save as PDF
if("ERN1" %in% rownames(expr_norm)) {
  plot_data_ern1 <- data.frame(
    expression = expr_norm["ERN1", ],
    group = group_list$group
  )
  
  gene_diff_ern1 <- diff_result["ERN1", ]
  padj_ern1 <- gene_diff_ern1$padj
  
  significance_ern1 <- ifelse(padj_ern1 < 0.001, "***",
                              ifelse(padj_ern1 < 0.01, "**",
                                     ifelse(padj_ern1 < 0.05, "*", "ns")))
  
  p_ern1 <- ggplot(plot_data_ern1, aes(x = group, y = expression, fill = group)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    labs(title = "Expression of ERN1",
         y = "Normalized Expression",
         x = "") +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.5)
    ) +
    scale_fill_manual(values = c("lightblue", "lightcoral"))
  
  if(significance_ern1 != "ns") {
    p_ern1 <- p_ern1 + geom_signif(
      comparisons = list(c("Control_Intermittent", "Intermittent")),
      annotations = significance_ern1,
      y_position = max(plot_data_ern1$expression) * 1.05,
      tip_length = 0.01,
      textsize = 4
    )
  }
  
  ggsave("violin_plot_ERN1.pdf", plot = p_ern1, width = 2.5, height = 3)
  print("ERN1 violin plot saved as PDF")
} else {
  message("Gene ERN1 not found in expression matrix")
}

# Generate violin plot for HSPA5 and save as PDF
if("HSPA5" %in% rownames(expr_norm)) {
  plot_data_hspa5 <- data.frame(
    expression = expr_norm["HSPA5", ],
    group = group_list$group
  )
  
  gene_diff_hspa5 <- diff_result["HSPA5", ]
  padj_hspa5 <- gene_diff_hspa5$padj
  
  significance_hspa5 <- ifelse(padj_hspa5 < 0.001, "***",
                               ifelse(padj_hspa5 < 0.01, "**",
                                      ifelse(padj_hspa5 < 0.05, "*", "ns")))
  
  p_hspa5 <- ggplot(plot_data_hspa5, aes(x = group, y = expression, fill = group)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    labs(title = "Expression of HSPA5",
         y = "Normalized Expression",
         x = "") +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.5)
    ) +
    scale_fill_manual(values = c("lightblue", "lightcoral"))
  
  if(significance_hspa5 != "ns") {
    p_hspa5 <- p_hspa5 + geom_signif(
      comparisons = list(c("Control_Intermittent", "Intermittent")),
      annotations = significance_hspa5,
      y_position = max(plot_data_hspa5$expression) * 1.05,
      tip_length = 0.01,
      textsize = 4
    )
  }
  
  ggsave("violin_plot_HSPA5.pdf", plot = p_hspa5, width = 2.5, height = 3)
  print("HSPA5 violin plot saved as PDF")
} else 
  message("Gene HSPA5 not found in expression matrix")
  