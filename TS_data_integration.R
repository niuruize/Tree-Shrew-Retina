#loading R packages
library(Matrix)
library(Seurat)
library(tidyverse)
library(rliger)
library(SeuratDisk)
library(SeuratWrappers)
#remotes::install_github("mojaveazure/seurat-disk")
#remotes::install_github('satijalab/seurat-wrappers')
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
library(tidyverse)

################################################################################
#loading data
#Pre-processing and quality control
################################################################################
Child1 <- Read10X(data.dir = "Data/Child1")
Child2 <- Read10X(data.dir = "Data/Child2")
Child3 <- Read10X(data.dir = "Data/Child3")
Child4 <- Read10X(data.dir = "Data/Child4")
Child5 <- Read10X(data.dir = "Data/Child5")
Adult1 <- Read10X(data.dir = "Data/Adult1")
Adult2 <- Read10X(data.dir = "Data/Adult2")
Adult3 <- Read10X(data.dir = "Data/Adult3")
Adult4 <- Read10X(data.dir = "Data/Adult4")
Adult5 <- Read10X(data.dir = "Data/Adult5")
Old1 <- Read10X(data.dir = "Data/Old1")
Old2 <- Read10X(data.dir = "Data/Old2")
Old3 <- Read10X(data.dir = "Data/Old3")
Old4 <- Read10X(data.dir = "Data/Old4")
Old5 <- Read10X(data.dir = "Data/Old5")
Child1 <- CreateSeuratObject(counts = Child1, project = "Child1", min.cells = 3, min.features = 200)
Child2 <- CreateSeuratObject(counts = Child2, project = "Child2", min.cells = 3, min.features = 200)
Child3 <- CreateSeuratObject(counts = Child3, project = "Child3", min.cells = 3, min.features = 200)
Child4 <- CreateSeuratObject(counts = Child4, project = "Child4", min.cells = 3, min.features = 200)
Child5 <- CreateSeuratObject(counts = Child5, project = "Child5", min.cells = 3, min.features = 200)
Adult1 <- CreateSeuratObject(counts = Adult1, project = "Adult1", min.cells = 3, min.features = 200)
Adult2 <- CreateSeuratObject(counts = Adult2, project = "Adult2", min.cells = 3, min.features = 200)
Adult3 <- CreateSeuratObject(counts = Adult3, project = "Adult3", min.cells = 3, min.features = 200)
Adult4 <- CreateSeuratObject(counts = Adult4, project = "Adult4", min.cells = 3, min.features = 200)
Adult5 <- CreateSeuratObject(counts = Adult5, project = "Adult5", min.cells = 3, min.features = 200)
Old1 <- CreateSeuratObject(counts = Old1, project = "Old1", min.cells = 3, min.features = 200)
Old2 <- CreateSeuratObject(counts = Old2, project = "Old2", min.cells = 3, min.features = 200)
Old3 <- CreateSeuratObject(counts = Old3, project = "Old3", min.cells = 3, min.features = 200)
Old4 <- CreateSeuratObject(counts = Old4, project = "Old4", min.cells = 3, min.features = 200)
Old5 <- CreateSeuratObject(counts = Old5, project = "Old5", min.cells = 3, min.features = 200)
Child1[["percent.mt"]] <- PercentageFeatureSet(Child1, pattern = "^MT-")
Child2[["percent.mt"]] <- PercentageFeatureSet(Child2, pattern = "^MT-")
Child3[["percent.mt"]] <- PercentageFeatureSet(Child3, pattern = "^MT-")
Child4[["percent.mt"]] <- PercentageFeatureSet(Child4, pattern = "^MT-")
Child5[["percent.mt"]] <- PercentageFeatureSet(Child5, pattern = "^MT-")
Adult1[["percent.mt"]] <- PercentageFeatureSet(Adult1, pattern = "^MT-")
Adult2[["percent.mt"]] <- PercentageFeatureSet(Adult2, pattern = "^MT-")
Adult3[["percent.mt"]] <- PercentageFeatureSet(Adult3, pattern = "^MT-")
Adult4[["percent.mt"]] <- PercentageFeatureSet(Adult4, pattern = "^MT-")
Adult5[["percent.mt"]] <- PercentageFeatureSet(Adult5, pattern = "^MT-")
Old1[["percent.mt"]] <- PercentageFeatureSet(Old1, pattern = "^MT-")
Old2[["percent.mt"]] <- PercentageFeatureSet(Old2, pattern = "^MT-")
Old3[["percent.mt"]] <- PercentageFeatureSet(Old3, pattern = "^MT-")
Old4[["percent.mt"]] <- PercentageFeatureSet(Old4, pattern = "^MT-")
Old5[["percent.mt"]] <- PercentageFeatureSet(Old5, pattern = "^MT-")
VlnPlot(Child1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.001, ncol = 3)
ggsave("VlnPlot/Child1.pdf", width = 35, height = 25, units = "cm")
VlnPlot(Child2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.001, ncol = 3)
ggsave("VlnPlot/Child2.pdf", width = 35, height = 25, units = "cm")
VlnPlot(Child3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.001, ncol = 3)
ggsave("VlnPlot/Child3.pdf", width = 35, height = 25, units = "cm")
VlnPlot(Child4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.001, ncol = 3)
ggsave("VlnPlot/Child4.pdf", width = 35, height = 25, units = "cm")
VlnPlot(Child5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.001, ncol = 3)
ggsave("VlnPlot/Child5.pdf", width = 35, height = 25, units = "cm")

VlnPlot(Adult1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.001, ncol = 3)
ggsave("VlnPlot/Adult1.pdf", width = 35, height = 25, units = "cm")
VlnPlot(Adult2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.001, ncol = 3)
ggsave("VlnPlot/Adult2.pdf", width = 35, height = 25, units = "cm")
VlnPlot(Adult3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.001, ncol = 3)
ggsave("VlnPlot/Adult3.pdf", width = 35, height = 25, units = "cm")
VlnPlot(Adult4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.001, ncol = 3)
ggsave("VlnPlot/Adult4.pdf", width = 35, height = 25, units = "cm")
VlnPlot(Adult5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.001, ncol = 3)
ggsave("VlnPlot/Adult5.pdf", width = 35, height = 25, units = "cm")

VlnPlot(Old1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.001, ncol = 3)
ggsave("VlnPlot/Old1.pdf", width = 35, height = 25, units = "cm")
VlnPlot(Old2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.001, ncol = 3)
ggsave("VlnPlot/Old2.pdf", width = 35, height = 25, units = "cm")
VlnPlot(Old3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.001, ncol = 3)
ggsave("VlnPlot/Old3.pdf", width = 35, height = 25, units = "cm")
VlnPlot(Old4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.001, ncol = 3)
ggsave("VlnPlot/Old4.pdf", width = 35, height = 25, units = "cm")
VlnPlot(Old5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.001, ncol = 3)
ggsave("VlnPlot/Old5.pdf", width = 35, height = 25, units = "cm")
VlnPlot(Old6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.001, ncol = 3)
ggsave("VlnPlot/Old6.pdf", width = 35, height = 25, units = "cm")
VlnPlot(AD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.001, ncol = 3)
ggsave("VlnPlot/AD.pdf", width = 35, height = 25, units = "cm")
Child1 <- subset(Child1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
Child2 <- subset(Child2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
Child3 <- subset(Child3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
Child4 <- subset(Child4, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
Child5 <- subset(Child5, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
Adult1 <- subset(Adult1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
Adult2 <- subset(Adult2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
Adult3 <- subset(Adult3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
Adult4 <- subset(Adult4, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
Adult5 <- subset(Adult5, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
Old1 <- subset(Old1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
Old2 <- subset(Old2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
Old3 <- subset(Old3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
Old4 <- subset(Old4, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
Old5 <- subset(Old5, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)

################################################################################
#Removal of double cells
################################################################################
scRNA <- merge(Child1,y=c(Child2,Child3,Child4,Child5,Adult1,Adult2,Adult3,Adult4,Adult5,Old1,Old2,Old3,Old4,Old5)) -> scRNA.orig
samples <- c("Child1", "Child2", "Child3", "Child4", "Child5", "Adult1", "Adult2", "Adult3", "Adult4", "Adult5", "Old1", "Old2", "Old3", "Old4", "Old5")

for (sample in samples) {
  test.seu <- scRNA[, scRNA@meta.data[["orig.ident"]] %in% c(sample)]
  test.seu <- RunPCA(test.seu, features = VariableFeatures(test.seu), npcs = 50)
  sweep.res.list <- paramSweep_v3(test.seu, PCs = 1:10, sct = FALSE)
  
  for(i in 1:length(sweep.res.list)){
    if(length(sweep.res.list[[i]]$pANN[is.nan(sweep.res.list[[i]]$pANN)]) != 0){
      if(i != 1){
        sweep.res.list[[i]] <- sweep.res.list[[i - 1]]
      }else{
        sweep.res.list[[i]] <- sweep.res.list[[i + 1]]
      }
    }
  }
  
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pk_v <- as.numeric(as.character(bcmvn$pK))
  pk_good <- pk_v[bcmvn$BCmetric == max(bcmvn$BCmetric)]
  
  DoubletRate <- ncol(test.seu) * 8 * 1e-6
  nExp_poi <- round(DoubletRate * length(colnames(test.seu)))
  test.seu <- doubletFinder_v3(test.seu, PCs = 1:10, pN = 0.25, pK = pk_good, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  colnames(test.seu@meta.data)[ncol(test.seu@meta.data)] <- "DoubletFinder"
  DF_df <- test.seu@meta.data[, c("orig.ident", "DoubletFinder")]
  file_name <- paste("Singlet/", sample, "_DoubletFinder_result.txt", sep = "")
  write.table(DF_df, file = file_name, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  DimPlot(test.seu, reduction = "umap", pt.size = 0.3, group.by = "DoubletFinder")
  
  table(test.seu@meta.data[["DoubletFinder"]])
  assign(sample, test.seu[, test.seu@meta.data$DoubletFinder %in% c("Singlet")])
  
  save(test.seu, file = paste("Singlet/", sample, "-Singlet.RData", sep = ""))
}

################################################################################
#data integration---liger
################################################################################

scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, split.by="orig.ident", do.center=FALSE)
nFactors=20
scRNA <- RunOptimizeALS(scRNA, k=nFactors, split.by="orig.ident")
scRNA <- RunQuantileNorm(scRNA, split.by="orig.ident")
scRNA$clusters <- factor(scRNA$clusters, levels=1:length(levels(scRNA$clusters)))
save(scRNA,file="scRNA.RData")

################################################################################
#data dimension reduction; find clusters
################################################################################

scRNA <- FindNeighbors(scRNA, reduction="iNMF", dims=1:20)
scRNA <- FindClusters(scRNA, resolution = 0.4)
scRNA <- RunUMAP(scRNA, dims=1:nFactors, reduction="iNMF")
table(scRNA@meta.data$orig.ident)
table(scRNA@meta.data$seurat_clusters)
scRNA$gender <- str_replace_all(scRNA$orig.ident, c("Child1"="Female","Child2"="Female","Child3"="Female","Adult1"="Female","Adult2"="Female","Adult3"="Female","Old1"="Female","Old2"="Female","Old3"="Female","Child4"="Male","Child5"="Male","Adult4"="Male","Adult5"="Male","Old4"="Male","Old5"="Male"))
scRNA@meta.data[["age"]]<-factor(scRNA@meta.data[["age"]], levels=c("Child","Adult","Old"))
table(scRNA@meta.data[["seurat_clusters"]], scRNA@meta.data[["orig.ident"]])
#
save(scRNA,file="scRNA.RData")

################################################################################
#visualization
################################################################################

p7 <- DimPlot(scRNA, reduction = "umap", cols = c('#b6e14f','#28bd6a','#090ba2'), group.by = "age", pt.size=0.001)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(scRNA, reduction = "umap", group.by = "ident", label.size = 8, pt.size=0.001, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap_Cluster1.pdf", plot = umap2, device = 'pdf', width = 62, height = 28, units = 'cm')

markers.to.plot <- c("OPN1LW","OPN1SW","ARR3","GRK7","PDE6H","PDE6C","GNAT2","RCVRN","CRX","PDE6B","PDE6A","ROM1","SAG","CA10",
  "GRIK1","GRIK1li1","ISL1","TRPM1","NETO1","ONECUT1","ONECUT2","NDRG1","RLBP1","WIF1","CHRDL1","TF","GLUL","C1QA","C1QB","CSF1R",
  "CX3CR1","PTPRC","AQP4","SLC1A3","SLC14A1","NRN1","SLC17A6","RBPMS","THY1","PAX6","SYNPR","NRXN1","MEIS2","TFAP2A","TFAP2B",
  "C1QL2","SLC6A1","TCF4","SLC6A9","SLC17A8","NFIA","NFIB","SLC18A3","MEGF11","RPE65","IGFBP7","PLP1","OTX2","MPZ","VSX1","CHAT",
  "BEST2","CRHBP","COL18A1","ZIC1","OPTC","COL9A1","CPAMD8","ZIC2","HCN1","SYT1")
dp2 <- DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2, dot.scale = 8) + RotatedAxis() #+ggplot2:::coord_flip()
ggsave(filename = "DotPlot1.pdf", plot = dp2, device = 'pdf', width = 45, height = 30, units = 'cm')

#FindAllMarkers
markers <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.25, only.pos = F, test.use = "wilcox") 
write.table(markers,file="markers_cluster.txt",quote=F,sep="\t",row.names=F,col.names=T)

################################################################################
#"NPE","Lens" and the cluster where the double cells are located are removed
################################################################################
scRNA <- scRNA[,scRNA@meta.data[["seurat_clusters"]] %in% c("0",	"1",	"2",	"3",	"4",	"5",	"6",	"7",	"8",	"9",	"10",	"11",	"12",	"13",	"14",	"15",	"16",	"17",	"18",	"19",	"21",	"22",	"24",	"25")]
count <- scRNA@assays$RNA@counts
meta <- subset(scRNA@meta.data, select= c("orig.ident","percent.mt","age","gender"))
scRNA <- CreateSeuratObject(counts = count, meta.data = meta , project = "TS", min.cells = 3, min.features = 200)

#Add cell labels
current.cluster.ids <- c("0",	"1",	"2",	"3",	"4",	"5",	"6",	"7",	"8",	"9",	"10",	"11",	"12",	"13",	"14",	"15",	"16",	"17",	"18",	"19",	"20",	"21",	"22",	"23",	"24",	"25",	"26",	"27",	"28")
new.cluster.ids <- c("Cone",	"Amacrine",	"Muller",	"Bipolar",	"Cone",	"RGC",	"Amacrine",	"Bipolar",	"Amacrine",	"Bipolar",	"Bipolar",	"Bipolar",	"Bipolar",	"Amacrine",	"Amacrine",	"Amacrine",	"Amacrine",	"Amacrine",	"Bipolar",	"Amacrine",	"Amacrine",	"Amacrine",	"Bipolar",	"Horizontal",	"Astrocytes",	"Rod",	"Amacrine",	"Microglia",	"RGC") 
scRNA@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA@meta.data$celltype)

Idents(scRNA) <- "celltype"
Idents(scRNA) <- factor(Idents(scRNA), levels = c("Cone","Rod","Horizontal","Bipolar","Amacrine","RGC","Muller","Astrocytes","Microglia"))
DimPlot(scRNA, reduction = "umap", label = TRUE)

table(scRNA@meta.data[["celltype"]], scRNA@meta.data[["gender"]])
table(scRNA@meta.data[["celltype"]], scRNA@meta.data[["age"]])
table(scRNA@meta.data[["celltype"]], scRNA@meta.data[["orig.ident"]])
###
save(scRNA,file="scRNA.RData")

#FindAllMarkers
markers <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.25, only.pos = F, test.use = "wilcox") 
write.table(markers,file="markers_celltype.txt",quote=F,sep="\t",row.names=F,col.names=T)


