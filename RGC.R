################################################################################
#...RGC...
################################################################################
#loading R packages
library(Matrix)
library(Seurat)
library(tidyverse)
library(rliger)
library(SeuratDisk)
library(SeuratWrappers)
library(patchwork)
library(cowplot)
library(ggplot2)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(MySeuratWrappers)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(DoubletFinder)
library(magrittr)
library(stringr)
library(multtest)
library(mindr)
library(scRNAtoolVis)
library(Scillus)

# extracting RGC cells
table(scRNA@meta.data[["celltype"]])
scRNA = scRNA[,scRNA@meta.data[["celltype"]] %in% c("RGC")]
# data integration
count <- scRNA@assays$RNA@counts
meta <- subset(scRNA@meta.data, select= c("orig.ident","percent.mt","age","gender"))
scRNA <- CreateSeuratObject(counts = count, meta.data = meta , project = "TS", min.cells = 3, min.features = 200)
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, split.by="orig.ident", do.center=FALSE)
nFactors=20 
scRNA <- RunOptimizeALS(scRNA, k=nFactors, split.by="orig.ident")
scRNA <- RunQuantileNorm(scRNA, split.by="orig.ident")
scRNA$clusters <- factor(scRNA$clusters, levels=1:length(levels(scRNA$clusters)))
save(scRNA,file="scRNA-1.RData")
scRNA <- FindNeighbors(scRNA, reduction="iNMF", dims=1:20)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims=1:nFactors, reduction="iNMF")
scRNA <- RunTSNE(scRNA, dims = 1:nFactors, reduction = "iNMF")
#
table(scRNA@meta.data$orig.ident)
table(scRNA@meta.data$seurat_clusters)
table(scRNA@meta.data[["seurat_clusters"]], scRNA@meta.data[["age"]])
table(scRNA@meta.data[["seurat_clusters"]], scRNA@meta.data[["age"]])
table(scRNA@meta.data[["seurat_clusters"]], scRNA@meta.data[["orig.ident"]])
save(scRNA,file="scRNA.RData")

# visualization
p7 <- DimPlot(scRNA, reduction = "umap", cols = c('#b6e14f','#28bd6a','#090ba2'), group.by = "age", pt.size=0.01)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(scRNA, reduction = "umap", group.by = "ident", label.size = 8, pt.size=0.01, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap_Cluster1.pdf", plot = umap2, device = 'pdf', width = 41, height = 16, units = 'cm')

markers.to.plot <- c("OPN1LW","OPN1SW","ARR3","GRK7","PDE6H","PDE6C","GNAT2","RCVRN","CRX","PDE6B","PDE6A","ROM1","SAG","CA10","GRIK1li1","ISL1","TRPM1","NETO1","ONECUT1","ONECUT2","NDRG1","RLBP1","WIF1","CHRDL1","TF","GLUL","C1QA","C1QB","CSF1R","CX3CR1","PTPRC","AQP4","SLC1A3","SLC14A1","NRN1","SLC17A6","RBPMS","THY1","PAX6","SYNPR","NRXN1","MEIS2","TFAP2A","TFAP2B","C1QL2","SLC6A1","TCF4","SLC6A9","SLC17A8","NFIA","NFIB","SLC18A3","MEGF11","ROBO1","SLIT2","SPP1","RUNX1","FOXP2")
dp2 <- DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2, dot.scale = 6) + RotatedAxis() #+ggplot2:::coord_flip()
ggsave(filename = "DotPlot1.pdf", plot = dp2, device = 'pdf', width = 37, height = 30, units = 'cm')

# Poor quality and double cells were excluded
scRNA <- scRNA[,scRNA@meta.data[["seurat_clusters"]] %in% c("1",	"2",	"3",	"4",	"5",	"6",	"7",	"8", "9",	"10",	"11",	"12",	"13",	"14",	"15",	"16",	"17",	"18",	"19",	"20", "21",	"22",	"23",	"24",	"25",	"26")]

# Integrate again
count <- scRNA@assays$RNA@counts
meta <- subset(scRNA@meta.data, select= c("orig.ident","percent.mt","age","gender"))
scRNA <- CreateSeuratObject(counts = count, meta.data = meta , project = "TS", min.cells = 3, min.features = 200)
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, split.by="orig.ident", do.center=FALSE)
nFactors=15
scRNA <- RunOptimizeALS(scRNA, k=nFactors, split.by="orig.ident")
scRNA <- RunQuantileNorm(scRNA, split.by="orig.ident")
scRNA$clusters <- factor(scRNA$clusters, levels=1:length(levels(scRNA$clusters)))
scRNA <- FindNeighbors(scRNA, reduction="iNMF", dims=1:15)
scRNA <- FindClusters(scRNA, resolution = 0.4)
scRNA <- RunUMAP(scRNA, dims=1:nFactors, reduction="iNMF")
scRNA <- RunTSNE(scRNA, dims = 1:nFactors, reduction = "iNMF")
#
table(scRNA@meta.data$orig.ident)
table(scRNA@meta.data$seurat_clusters)
table(scRNA@meta.data[["seurat_clusters"]], scRNA@meta.data[["age"]])
table(scRNA@meta.data[["seurat_clusters"]], scRNA@meta.data[["gender"]])
save(scRNA,file="scRNA.RData")

# visualization
p7 <- DimPlot(scRNA, reduction = "umap", cols = c('#b6e14f','#28bd6a','#090ba2'), group.by = "age", pt.size=0.01)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(scRNA, reduction = "umap", group.by = "ident", label.size = 8, pt.size=0.01, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap_Cluster2.pdf", plot = umap2, device = 'pdf', width = 32, height = 13, units = 'cm')

markers.to.plot <- c("OPN1LW","OPN1SW","ARR3","GRK7","PDE6H","PDE6C","GNAT2","RCVRN","CRX","PDE6B","PDE6A","ROM1","SAG","CA10","GRIK1li1","ISL1","TRPM1","NETO1","ONECUT1","ONECUT2","NDRG1","RLBP1","WIF1","CHRDL1","TF","GLUL","C1QA","C1QB","CSF1R","CX3CR1","PTPRC","AQP4","SLC1A3","SLC14A1","NRN1","SLC17A6","RBPMS","THY1","PAX6","SYNPR","NRXN1","MEIS2","TFAP2A","TFAP2B","C1QL2","SLC6A1","TCF4","SLC6A9","SLC17A8","NFIA","NFIB","SLC18A3","MEGF11","ROBO1","SLIT2","SPP1","RUNX1","FOXP2")
dp2 <- DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2, dot.scale = 7) + RotatedAxis() #+ggplot2:::coord_flip()
ggsave(filename = "DotPlot2.pdf", plot = dp2, device = 'pdf', width = 37, height = 20, units = 'cm')

#FindAllMarkers
markers <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.25, only.pos = F, test.use = "wilcox")  
write.table(markers,file="markers_cluster.txt",quote=F,sep="\t",row.names=F,col.names=T)

# Add cell labels
current.cluster.ids <- c("0",	"1",	"2",	"3",	"4",	"5",	"6",	"7",	"8",	"9",	"10")
new.cluster.ids <- c("C0", "C1",	"C2",	"C3",	"C4",	"C5",	"C6",	"C7",	"C8",	"C9",	"C10") 
scRNA@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA@meta.data$celltype)
Idents(scRNA) <- "celltype"
Idents(scRNA) <- factor(Idents(scRNA), levels = c("C0","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10"))
DimPlot(scRNA, reduction = "umap", label = TRUE)
table(scRNA@meta.data[["celltype"]], scRNA@meta.data[["gender"]])
table(scRNA@meta.data[["celltype"]], scRNA@meta.data[["age"]])
table(scRNA@meta.data[["celltype"]], scRNA@meta.data[["orig.ident"]])
save(scRNA,file="RGC.RData")

#DimPlot
p7 <- DimPlot(scRNA, reduction = "umap", cols = c('#b6e14f','#28bd6a','#090ba2'), group.by = "age", pt.size=0.01)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(scRNA, reduction = "umap", cols = c('#6E90D0', '#b4b47b', '#93B727', '#6d9150', '#6c584c', '#e4c2ce', '#595757','#E999C6', '#EDDD8E', '#BCDD9D', '#C25A7F', "#FBBD23","#F0B371",'#6f4fa3', '#ff9a94','#C2A3DA','#CD6F44','#FFE27E',"#f0dbb4",'#b6e14f','#28bd6a','#090ba2',"#168c90","#77c99b"), group.by = "celltype", label.size = 7, pt.size=0.01, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap_Cluster分群后.pdf", plot = umap2, device = 'pdf', width = 31, height = 13, units = 'cm')

#DotPlot
markers.to.plot <- c("SLC17A6","RBPMS","ROBO1","SLIT2","EBF2","HTR2C","FOXP2","CDH8","RUNX1","WLS","TBR1","SPP1","IL1RAPL2","CALCRL")
dp2 <- DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "celltype", col.min = 0,col.max = 2, dot.scale = 8) + RotatedAxis() #+ggplot2:::coord_flip()
ggsave(filename = "DotPlot分群后.pdf", plot = dp2, device = 'pdf', width = 15, height = 4.5, units = 'cm')

#jjplot
jjplot <- jjDotPlot(object = scRNA,
                    gene = c("SLC17A6","RBPMS","ROBO1","SLIT2","EBF2","HTR2C","FOXP2","CDH8","RUNX1","WLS","TBR1","SPP1","IL1RAPL2","CALCRL"),
                    gene.order = c("SLC17A6","RBPMS","ROBO1","SLIT2","EBF2","HTR2C","FOXP2","CDH8","RUNX1","WLS","TBR1","SPP1","IL1RAPL2","CALCRL"),
                    id = 'celltype',
                    cluster.order = c("C0","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10"),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 8)#+ggplot2:::coord_flip()
jjplot
ggsave(filename = "jjDotPlot-celltype.pdf", plot = jjplot, device = 'pdf', width = 19, height = 16, units = 'cm')

#FindAllMarkers
markers <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.25, only.pos = F, test.use = "wilcox") 
write.table(markers,file="markers_celltype.txt",quote=F,sep="\t",row.names=F,col.names=T)

# Violin diagram of differential genes
library(ggpubr)
library(ggplot2)
library(ggsignif)
library(patchwork)
#
my9color <- c('#fbd3b9','#a1c9e5','#417bb9')
my_comparisons <- list(c('Child','Adult'),c('Adult','Old'),c('Child','Old'))
#
VlnPlot(scRNA, idents = c("C1"), features = c("MDGA2"), group.by = "age",  pt.size=0) +
  NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") +
  ylim(0, 8.5)+
  geom_boxplot(fill="#a6a7ac",color="#a6a7ac",
               width = 0.1,outlier.shape = NA)+
  scale_color_manual(values = my9color) + 
  scale_fill_manual(values = my9color)
ggsave('MDGA2_vlo.pdf',width = 5,height = 4)

c("SLC35F1","MDGA2","DYSF")
#
VlnPlot(scRNA, idents = c("C0"), features = c("NELL1",	"LRRN3",	"CNTNAP2-AS-12",	"CTNNA3",	"PCDH10",	"TMEM163",	"ATP2B1",	"PDE3A-AS-3"), group.by = "age",  stacked = T,pt.size=0,combine = T) +
  NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") +
  ylim(0, 8.5)+
  geom_boxplot(fill="#a6a7ac",color="#a6a7ac",
               width = 0.1,outlier.shape = NA)+
  scale_color_manual(values = my9color) + 
  scale_fill_manual(values = my9color)
ggsave('RGC3_vlo.pdf',width = 5,height = 4)

# correlation analysis
table(scRNA$age)  
av<-AverageExpression(scRNA,group.by = "celltype", assays = "RNA")
av=av[[1]]
head(av)
cg=names(tail(sort(apply(av,1,sd)),1000))
View(av[cg,])
View(cor(av[cg,],method = "spearman"))
write.csv(cor(av[cg,],method = "spearman"),"cor_celltype.csv")
pdf("cor_celltype.pdf", width = 6,height = 5.5)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),treeheight_row = 15,treeheight_col = 15,fontsize = 15)
dev.off()


################################################################################
#...SCENIC...
################################################################################
library(Seurat)
library(AUCell)
library(SCENIC)
library(RcisTarget)
library(GENIE3)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)

# load data
setwd("/home/test-user/TS/SCENIC/RGC")
getwd()
load('scRNA.RData')
AB <- subset(scRNA, downsample=100)
exprMat<-GetAssayData(object = AB)
dim(exprMat)
exprMat[1:4,1:10]
cellInfo <- AB@meta.data[,c("nCount_RNA","nFeature_RNA","age","celltype")]
head(cellInfo)
cellTypeColumn <- "celltype"
colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "CellType"
cbind(table(cellInfo$CellType))
head(cellInfo)
# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(CellType=c("C0"="forestgreen","C1"="darkorange","C2"="magenta4","C3"="hotpink","C4"="red3","C5"="skyblue","C6"="darkblue","C7"="peru","C8"="salmon","C9"="violet","C10"="Pink4"))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
org="hgnc" # or mgi-mouse, or dmel-fly
dbDir="/home/test-user/TS/SCENIC/RGC/cisTarget_databases" # RcisTarget databases location
myDatasetTitle="TS" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir,  nCores=16)
#Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
#Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
saveRDS(cellInfo, file="int/cellInfo.Rds")
save(exprMat,file = "int/exprMat.Rds")
saveRDS(colVars, file="int/colVars.Rds")
exprMat<-as.matrix(exprMat)
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
exprMat_filtered_log <- na.omit(exprMat_filtered_log)
runGenie3(exprMat_filtered_log, scenicOptions)

scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 16
scenicOptions@settings$seed <- 123

scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] 
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget"))
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log,skipTsne = T)
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)

# celltype
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[-c(4,9,10,11,12,13,15,16,17),]
##选择
pdf("RGC-TFs.pdf", width = 6, height = 8)
pheatmap::pheatmap(regulonActivity_byCellType_Scaled,
                   color = colorRampPalette(c("#330066","#66CC66","#FFCC33"))(100),
                   breaks = seq(-1, 2, length.out = 100),
                   border_color = "NA")
dev.off()

# filename="regulonActivity_byCellType.pdf", width=10, height=20)
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)
write.csv(topRegulators,"topRegulators_cell.csv")


minPerc <- .5
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$CellType),
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
pheatmap::pheatmap(binaryActPerc_subset, # fontsize_row=5,
                   color = colorRampPalette(c("white","pink","red"))(100), breaks=seq(0, 1, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
# group
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$age),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
pheatmap::pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)

# filename="regulonActivity_byCellType.pdf", width=10, height=20)
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "age", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)
write.csv(topRegulators,"topRegulators_age.csv")

minPerc <- .1
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$age),
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
pheatmap::pheatmap(binaryActPerc_subset, # fontsize_row=5,
                   color = colorRampPalette(c("white","pink","red"))(100), breaks=seq(0, 1, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)


################################################################################
# velocity
################################################################################
# analysis in python
import os
os.chdir("/home/xll-zxcv/TS/loom/Seurat_dataset/")
os.getcwd()

import scanpy as sc
import anndata 
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

type="RGC"

#(1) load sparse matrix:
X = io.mmread(f"/home/xll-zxcv/TS/loom/Seurat_dataset/RGC/counts.mtx")

# create anndata object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)

#(2) load cell metadata:
cell_meta = pd.read_csv(f"/home/xll-zxcv/TS/loom/Seurat_dataset/RGC/metadata.csv")

# load gene names:
with open(f"/home/xll-zxcv/TS/loom/Seurat_dataset/RGC/gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']

adata.var.index = gene_names

#(3) load dimensional reduction:
pca = pd.read_csv(f"/home/xll-zxcv/TS/loom/Seurat_dataset/RGC/pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:
sc.pl.umap(adata, color='seurat_clusters', frameon=False, save=True)

#(4) load loom data
import scvelo as scv
import scanpy as sc
#import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

import matplotlib.pyplot as plt

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, 
                               fontsize=6, color_map='viridis',
                               frameon=False)

# load loom files for spliced/unspliced matrices for each sample:
sample_loom_1="/home/xll-zxcv/TS/loom/Child1_count.loom"
ldata1 = scv.read( sample_loom_1, cache=True)
#ldata1.obs['batch'] = '1'
#ldata1.obs['BATCH'] = '1'
ldata1

sample_loom_2="/home/xll-zxcv/TS/loom/Child2_count.loom"
ldata2 = scv.read( sample_loom_2, cache=True)
#ldata2.obs['batch'] = '2'
#ldata2.obs['BATCH'] = '2'
ldata2

sample_loom_3="/home/xll-zxcv/TS/loom/Child3_count.loom"
ldata3 = scv.read( sample_loom_3, cache=True)
#ldata3.obs['batch'] = '3'
#ldata3.obs['BATCH'] = '3'
ldata3

sample_loom_4="/home/xll-zxcv/TS/loom/Child4_count.loom"
ldata4 = scv.read( sample_loom_4, cache=True)
#ldata4.obs['batch'] = '4'
#ldata4.obs['BATCH'] = '4'
ldata4

sample_loom_5="/home/xll-zxcv/TS/loom/Child5_count.loom"
ldata5 = scv.read( sample_loom_5, cache=True)
#ldata5.obs['batch'] = '5'
#ldata5.obs['BATCH'] = '5'
ldata5

sample_loom_6="/home/xll-zxcv/TS/loom/Adult1_count.loom"
ldata6 = scv.read( sample_loom_6, cache=True)
#ldata6.obs['batch'] = '6'
#ldata6.obs['BATCH'] = '6'
ldata6

sample_loom_7="/home/xll-zxcv/TS/loom/Adult2_count.loom"
ldata7 = scv.read( sample_loom_7, cache=True)
#ldata7.obs['batch'] = '7'
#ldata7.obs['BATCH'] = '7'
ldata7

sample_loom_8="/home/xll-zxcv/TS/loom/Adult3_count.loom"
ldata8 = scv.read( sample_loom_8, cache=True)
#ldata8.obs['batch'] = '8'
#ldata8.obs['BATCH'] = '8'
ldata8

sample_loom_9="/home/xll-zxcv/TS/loom/Adult4_count.loom"
ldata9 = scv.read( sample_loom_9, cache=True)
#ldata9.obs['batch'] = '9'
#ldata9.obs['BATCH'] = '9'
ldata9

sample_loom_10="/home/xll-zxcv/TS/loom/Adult5_count.loom"
ldata10 = scv.read( sample_loom_10, cache=True)
#ldata10.obs['batch'] = '10'
#ldata10.obs['BATCH'] = '10'
ldata10

sample_loom_11="/home/xll-zxcv/TS/loom/Old1_count.loom"
ldata11 = scv.read( sample_loom_11, cache=True)
#ldata11.obs['batch'] = '11'
#ldata11.obs['BATCH'] = '11'
ldata11

sample_loom_12="/home/xll-zxcv/TS/loom/Old2_count.loom"
ldata12 = scv.read( sample_loom_12, cache=True)
#ldata12.obs['batch'] = '12'
#ldata12.obs['BATCH'] = '12'
ldata12

sample_loom_13="/home/xll-zxcv/TS/loom/Old3_count.loom"
ldata13 = scv.read( sample_loom_13, cache=True)
#ldata13.obs['batch'] = '13'
#ldata13.obs['BATCH'] = '13'
ldata13

sample_loom_14="/home/xll-zxcv/TS/loom/Old4_count.loom"
ldata14 = scv.read( sample_loom_14, cache=True)
#ldata14.obs['batch'] = '14'
#ldata14.obs['BATCH'] = '14'
ldata14

sample_loom_15="/home/xll-zxcv/TS/loom/Old5_count.loom"
ldata15 = scv.read( sample_loom_15, cache=True)
#ldata15.obs['batch'] = '15'
#ldata15.obs['BATCH'] = '15'
ldata15

# (5) merge loom data
ldata1.obs.index[0:2]
ldata2.obs.index[0:2]
ldata3.obs.index[0:2]
ldata4.obs.index[0:2]
ldata5.obs.index[0:2]
ldata6.obs.index[0:2]
ldata7.obs.index[0:2]
ldata8.obs.index[0:2]
ldata9.obs.index[0:2]
ldata10.obs.index[0:2]
ldata11.obs.index[0:2]
ldata12.obs.index[0:2]
ldata13.obs.index[0:2]
ldata14.obs.index[0:2]
ldata15.obs.index[0:2]

# ldata1
barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata1.obs.index = barcodes
ldata1.obs.index[0:5]

# ldata2
barcodes = [bc.split(':')[1] for bc in ldata2.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata2.obs.index = barcodes
ldata2.obs.index[0:5]

# ldata3
barcodes = [bc.split(':')[1] for bc in ldata3.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata3.obs.index = barcodes
ldata3.obs.index[0:5]

# ldata4
barcodes = [bc.split(':')[1] for bc in ldata4.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata4.obs.index = barcodes
ldata4.obs.index[0:5]

# ldata5
barcodes = [bc.split(':')[1] for bc in ldata5.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata5.obs.index = barcodes
ldata5.obs.index[0:5]

# ldata6
barcodes = [bc.split(':')[1] for bc in ldata6.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata6.obs.index = barcodes
ldata6.obs.index[0:5]

# ldata7
barcodes = [bc.split(':')[1] for bc in ldata7.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata7.obs.index = barcodes
ldata7.obs.index[0:5]

# ldata8
barcodes = [bc.split(':')[1] for bc in ldata8.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata8.obs.index = barcodes
ldata8.obs.index[0:5]

# ldata9
barcodes = [bc.split(':')[1] for bc in ldata9.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata9.obs.index = barcodes
ldata9.obs.index[0:5]

# ldata10
barcodes = [bc.split(':')[1] for bc in ldata10.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata10.obs.index = barcodes
ldata10.obs.index[0:5]

# ldata11
barcodes = [bc.split(':')[1] for bc in ldata11.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata11.obs.index = barcodes
ldata11.obs.index[0:5]

# ldata12
barcodes = [bc.split(':')[1] for bc in ldata12.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata12.obs.index = barcodes
ldata12.obs.index[0:5]

# ldata13
barcodes = [bc.split(':')[1] for bc in ldata13.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata13.obs.index = barcodes
ldata13.obs.index[0:5]

# ldata14
barcodes = [bc.split(':')[1] for bc in ldata14.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata14.obs.index = barcodes
ldata14.obs.index[0:5]

# ldata15
barcodes = [bc.split(':')[1] for bc in ldata15.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata15.obs.index = barcodes
ldata15.obs.index[0:5]

ldata1.var.head()
ldata2.var.head()
ldata3.var.head()
ldata4.var.head()
ldata5.var.head()
ldata6.var.head()
ldata7.var.head()
ldata8.var.head()
ldata9.var.head()
ldata10.var.head()
ldata11.var.head()
ldata12.var.head()
ldata13.var.head()
ldata14.var.head()
ldata15.var.head()

ldata1.var_names_make_unique()
ldata2.var_names_make_unique()
ldata3.var_names_make_unique()
ldata4.var_names_make_unique()
ldata5.var_names_make_unique()
ldata6.var_names_make_unique()
ldata7.var_names_make_unique()
ldata8.var_names_make_unique()
ldata9.var_names_make_unique()
ldata10.var_names_make_unique()
ldata11.var_names_make_unique()
ldata12.var_names_make_unique()
ldata13.var_names_make_unique()
ldata14.var_names_make_unique()
ldata15.var_names_make_unique()

# merge
ldata = ldata1.concatenate([ldata2,ldata3,ldata4,ldata5,ldata6,ldata7,ldata8,ldata9,ldata10,ldata11,ldata12,ldata13,ldata14,ldata15])
ldata = sc.AnnData.concatenate(ldata1,ldata2,ldata3,ldata4,ldata5,ldata6,ldata7,ldata8,ldata9,ldata10,ldata11,ldata12,ldata13,ldata14,ldata15,batch_key = 'BATCH')
ldata.obs

# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)
adata.obs

adata.write_h5ad(f'/home/xll-zxcv/TS/loom/Seurat_dataset/RGC/RGC.adata_ldata.h5ad')

# (6) scVelo
adata = sc.read(f'/home/xll-zxcv/TS/loom/Seurat_dataset/RGC/RGC.adata_ldata.h5ad')
scv.pp.moments(adata, n_neighbors=30, n_pcs=30)

import gc
gc.collect()
#
temp_pre= f"RGC_nue.in_process2"
if False==os.path.exists(f"{temp_pre}.velo.gz.h5ad"):
    scv.tl.recover_dynamics(adata, var_names='all', n_jobs=18)
    scv.tl.velocity(adata, mode='dynamical')
    adata.write(f"{temp_pre}.velo.gz.h5ad", compression='gzip')
    print(">>Write to file")
else:
    adata = sc.read(f"{temp_pre}.velo.gz.h5ad", compression='gzip', ext="h5ad")
    print(">>read from file")

scv.tl.velocity_graph(adata, n_jobs=18)
scv.pl.velocity_embedding_stream(adata, basis="umap")
scv.pl.velocity_embedding(adata, basis="umap", save='embedding_RGC.pdf')

scv.pl.velocity_embedding_grid(adata, basis='umap', color='celltype', figsize=
(8,8), arrow_size=1.8,  arrow_length=2.5, size=100, alpha=0.1,
save='embedding_grid_RGC.pdf', title='', scale=0.1)

scv.pl.velocity_embedding_stream(adata, basis='umap', color=
['seurat_clusters'], figsize=(6,6),palette = ("#f26f66","#4962a1","#8572ad",'#3082b4',"#539075",'#63b9bd',"#ccaf74",
"#2F4858","#D5EC66","#86cc9d","#94755c","#92B674"), alpha=0.1, density=2.5, size=140,
arrow_size=1.2, linewidth=0.5, legend_fontsize=17, dpi=900,
save='embedding_stream_seurat_clusters_RGC.svg', title='')

scv.pl.velocity_embedding_stream(adata, basis='umap', color=
['celltype'], figsize=(6,6), palette = ("#f26f66","#4962a1","#8572ad",'#3082b4',"#539075",'#63b9bd',"#ccaf74",
"#2F4858","#D5EC66","#86cc9d","#94755c","#92B674"), alpha=0.1, density=2.5, size=140,
arrow_size=1.2, linewidth=0.5, legend_fontsize=17, dpi=900,
save='embedding_stream_celltype_RGC.svg', title='')

scv.pl.velocity_embedding(adata, arrow_length=2, arrow_size=1.5, dpi=900,
figsize=(8,8), save='embedding_stream2_RGC.pdf')

#visualization method 1
scv.pl.velocity(adata,
["NRN1","SLC17A6","RBPMS","THY1"],
figsize=(5,5), dpi=600, save='Interprete velocity_RGC1.png', ncols=2)

scv.pl.velocity(adata,
["ROBO1","SLIT2","EBF2","HTR2C"],
figsize=(5,5), dpi=600, save='Interprete velocity_RGC2.png', ncols=2)

scv.pl.velocity(adata,
["FOXP2","CDH8","RUNX1","WLS"], figsize=(5,5), dpi=600,
save='Interprete velocity_RGC3.png', ncols=2)

scv.pl.velocity(adata,
["TBR1","SPP1","IL1RAPL2","CALCRL"],
figsize=(5,5), dpi=600, save='Interprete velocity_RGC4.png', ncols=2)

#velocity and consistency
scv.tl.velocity_confidence(adata)
keys = 'velocity_length','velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], dpi=900, figsize=(6,5), save='Speed and coherence_RGC.pdf')

df = adata.obs.groupby('celltype')[keys].mean().T
df.style.background_gradient(cmap='coolwarm', axis=1)

scv.pl.velocity_graph(adata, threshold=.1, color='celltype', figsize=(8,8), dpi=900, save='velocity_graph_RGC.pdf')

#visualization
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', figsize=(5,5), dpi=300 ,save='velocity_pseudotime_RGC.pdf')

###PAGA
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='celltype')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.pl.paga(adata, basis='umap', size=80, alpha=.15,
            min_edge_width=2, node_size_scale=1.5, figsize=(6,5),arrowsize=20, node_size_power=0.5,dpi=900, save='paga_RGC.svg')

##Latent time
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, dpi=900, save='latent_time_RGC.pdf')
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]

scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='celltype', n_convolve=100, save='heatmap_RGC.pdf')

