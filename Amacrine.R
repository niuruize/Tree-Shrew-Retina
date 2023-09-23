################################################################################
#...Amacrine...
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

# extracting amacrine cells
table(scRNA@meta.data[["celltype"]])
scRNA = scRNA[,scRNA@meta.data[["celltype"]] %in% c("Amacrine")]
# data integration
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, split.by="orig.ident", do.center=FALSE)
nFactors=30 
scRNA <- RunOptimizeALS(scRNA, k=nFactors, split.by="orig.ident")
scRNA <- RunQuantileNorm(scRNA, split.by="orig.ident")
scRNA$clusters <- factor(scRNA$clusters, levels=1:length(levels(scRNA$clusters)))
save(scRNA,file="scRNA.RData")
# data dimension reduction; find clusters
scRNA <- FindNeighbors(scRNA, reduction="iNMF", dims=1:30)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims=1:nFactors, reduction="iNMF")
scRNA <- RunTSNE(scRNA, dims = 1:nFactors, reduction = "iNMF")
# cell number
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

markers.to.plot <- c("OPN1LW","OPN1SW","ARR3","GRK7","PDE6H","PDE6C","GNAT2","RCVRN","CRX","PDE6B","PDE6A","ROM1","SAG","CA10","GRIK1li1","ISL1","TRPM1","NETO1","ONECUT1","ONECUT2","NDRG1","RLBP1","WIF1","CHRDL1","TF","GLUL","C1QA","C1QB","CSF1R","CX3CR1","PTPRC","AQP4","SLC1A3","SLC14A1","NRN1","SLC17A6","RBPMS","THY1","PAX6","SYNPR","NRXN1","MEIS2","TFAP2A","TFAP2B","C1QL2","SLC6A1","TCF4","SLC6A9","SLC17A8","NFIA","NFIB","SLC18A3","MEGF11")
dp2 <- DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2, dot.scale = 6) + RotatedAxis() #+ggplot2:::coord_flip()
ggsave(filename = "DotPlot1.pdf", plot = dp2, device = 'pdf', width = 37, height = 30, units = 'cm')

# Removal of double cells
scRNA <- scRNA[,scRNA@meta.data[["seurat_clusters"]] %in% c("0",	"1",	"2",	"3",	"4",	"5",	"6",	"7",	"9",	"10",	"11",	"12",	"13",	"14",	"15",	"16",	"17",	"18",	"19",	"21",	"22",	"23",	"24",	"25",	"26",	"27",	"28",	"29",	"30",	"31",	"33",	"34",	"36",	"37",	"38",	"39",	"40",	"41",	"42",	"44",	"45",	"47",	"48",	"50",	"51")]
# Integrate again
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
save(scRNA,file="scRNA-3.RData")
scRNA <- FindNeighbors(scRNA, reduction="iNMF", dims=1:20)
scRNA <- FindClusters(scRNA, resolution = 0.4)
scRNA <- RunUMAP(scRNA, dims=1:nFactors, reduction="iNMF")
scRNA <- RunTSNE(scRNA, dims = 1:nFactors, reduction = "iNMF")
p7 <- DimPlot(scRNA, reduction = "umap", cols = c('#b6e14f','#28bd6a','#090ba2'), group.by = "age", pt.size=0.01)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(scRNA, reduction = "umap", group.by = "ident", label.size = 8, pt.size=0.01, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap_Cluster2.pdf", plot = umap2, device = 'pdf', width = 32, height = 13, units = 'cm')

markers.to.plot <- c("CA10","GRIK1li1","PAX6","SYNPR","NRXN1","MEIS2","TFAP2A","TFAP2B","C1QL2","SLC6A1","TCF4","SLC6A9","SLC17A8","NFIA","NFIB","SLC18A3","ISL1","MEGF11","CHAT","PDGFRA","CBLN1","TACR1","EBF3","EBF2","NELL1",	"MAF","RELN","NPY","VIP","MEF2C","CALB2","KCNT2","KCND2","GRM1","GPC5","IDO1")
dp2 <- DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2, dot.scale = 7) + RotatedAxis() #+ggplot2:::coord_flip()
ggsave(filename = "DotPlot2.pdf", plot = dp2, device = 'pdf', width = 37, height = 20, units = 'cm')

# Removal of unknown cells ‘7’，‘21’
scRNA <- scRNA[,scRNA@meta.data[["seurat_clusters"]] %in% c("0",	"1",	"2",	"3",	"4",	"5",	"6",	"8",  "9",	"10",	"11",	"12",	"13",	"14",	"15",	"16",	"17",	"18",	"19",	"20", "22",	"23",	"24",	"25",	"26",	"27")]
# Integrate again
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
save(scRNA,file="scRNA-5.RData")
scRNA <- FindNeighbors(scRNA, reduction="iNMF", dims=1:20)
scRNA <- FindClusters(scRNA, resolution = 0.4)
scRNA <- RunUMAP(scRNA, dims=1:nFactors, reduction="iNMF")
scRNA <- RunTSNE(scRNA, dims = 1:nFactors, reduction = "iNMF")
save(scRNA,file="scRNA.RData")

# FindAllMarkers
markers <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.25, only.pos = F, test.use = "wilcox")
write.table(markers,file="markers分群前2.txt",quote=F,sep="\t",row.names=F,col.names=T)

# Add cell labels
current.cluster.ids <- c("0",	"1",	"2",	"3",	"4",	"5",	"6",	"7",	"8",	"9",	"10",	"11",	"12",	"13",	"14",	"15",	"16",	"17",	"18",	"19",	"20",	"21",	"22",	"23",	"24",	"25",	"26",	"27")
new.cluster.ids <- c("NFIB",	"SAC",	"GlyT1",	"GABA",	"TAC1",	"VGluT3",	"GlyT1",	"FAT4",	"GABA",	"GABA",	"GABA",	"GABA",	"GlyT1",	"DRD2",	"RELN",	"GABA",	"GlyT1",	"PDGFRα",	"GABA",	"GlyT1",	"GABA",	"GABA",	"GlyT1",	"GABA",	"GABA",	"GABA",	"GABA",	"VGluT3") 
scRNA@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA@meta.data$celltype)

Idents(scRNA) <- "celltype"
Idents(scRNA) <- factor(Idents(scRNA), levels = c("NFIB",	"SAC",	"GlyT1",	"GABA",	"TAC1",	"VGluT3",	"FAT4",	"DRD2",	"RELN",	"PDGFRα"))
DimPlot(scRNA, reduction = "umap", label = TRUE)

table(scRNA@meta.data[["celltype"]], scRNA@meta.data[["gender"]])
table(scRNA@meta.data[["celltype"]], scRNA@meta.data[["age"]])
table(scRNA@meta.data[["celltype"]], scRNA@meta.data[["orig.ident"]])
save(scRNA,file="AC.RData")
# FindAllMarkers
markers <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.25, only.pos = F, test.use = "wilcox")
write.table(markers,file="markers_celltype.txt",quote=F,sep="\t",row.names=F,col.names=T)

# visualization
p7 <- DimPlot(scRNA, reduction = "umap", cols = c('#b6e14f','#28bd6a','#090ba2'), group.by = "age", pt.size=0.01)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(scRNA, reduction = "umap", cols = c('#C75764','#4995C6','#6E90D0','#FF8C97',"#e4c2ce",'#1663A9','#A6C6FF','#8481BA',"#B4B4D5","#936674"), group.by = "celltype", label.size = 4, pt.size=0.01, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap_Cluster中群.pdf", plot = umap2, device = 'pdf', width = 38, height = 13, units = 'cm')

markers.to.plot <- c("PAX6","SYNPR","NRXN1","MEIS2","TFAP2A","TFAP2B","C1QL2","SLC6A1","TCF4","SLC6A9","SLC17A8","NFIA","NFIB","SLC18A3","ISL1","MEGF11","CHAT","PDGFRA","CBLN1","TACR1","EBF3","EBF2","NELL1",	"MAF","RELN","NPY","VIP","EPHA3","MEF2C","CALB2","TRPC4","FAT4","KCNQ5","EFEMP1","ITGA6","DRD2","CBLN2","EGFEM1","CDH9","AKAP2","FBN2","NTF3","IL1RAPL2","TRHDE","COL24A1","VEGFC","POSTN","NABP1","BCL11B","TAC1")
dp2 <- DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "celltype", col.min = 0,col.max = 2, dot.scale = 7) + RotatedAxis() #+ggplot2:::coord_flip()
ggsave(filename = "DotPlot中群.pdf", plot = dp2, device = 'pdf', width = 37, height = 17, units = 'cm')

jjplot <- jjDotPlot(object = scRNA,
                    gene = c("PAX6","NRXN1","TFAP2A","SLC6A1","GAD1","MEIS2","TCF4","SLC6A9","EBF3","EBF2","SLC18A3","ISL1","MEGF11","CHAT","SLC17A8","TAC1","VIP","EPHA3","NFIA","CBLN1","PDGFRA","TACR1","NFIB","FAT4","DRD2","RELN"),
                    gene.order = c("PAX6","NRXN1","TFAP2A","SLC6A1","GAD1","MEIS2","TCF4","SLC6A9","EBF3","EBF2","SLC18A3","ISL1","MEGF11","CHAT","SLC17A8","TAC1","VIP","EPHA3","NFIA","CBLN1","PDGFRA","TACR1","NFIB","FAT4","DRD2","RELN"),
                    id = 'celltype',
                    cluster.order = c("GABA","GlyT1","SAC","VGluT3","TAC1","PDGFRα","NFIB","FAT4","DRD2","RELN"),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 8)#+ggplot2:::coord_flip()
jjplot
ggsave(filename = "jjDotPlot1.pdf", plot = jjplot, device = 'pdf', width = 28, height = 16, units = 'cm')


jjplot <- jjDotPlot(object = scRNA,
                    gene = c("PAX6","NRXN1","AKAP2","BCL11B","CALB2","CDH9","COL24A1","DRD2",	"EBF3","EBF2","VIP","EPHA3","FAT4","FBN2","IL1RAPL2","ITGA6","MAF","NABP1","NELL1","NFIB","NPY","PDGFRA","POSTN","RELN","MEGF11","CHAT","TAC1","TRPC4","VEGFC","SLC17A8"),
                    gene.order = c("PAX6","NRXN1","AKAP2","BCL11B","CALB2","CDH9","COL24A1","DRD2",	"EBF3","EBF2","VIP","EPHA3","FAT4","FBN2","IL1RAPL2","ITGA6","MAF","NABP1","NELL1","NFIB","NPY","PDGFRA","POSTN","RELN","MEGF11","CHAT","TAC1","TRPC4","VEGFC","SLC17A8"),
                    id = 'celltype',
                    cluster.order = c("AKAP2",	"BCL11B",	"CALB2",	"CDH9",	"COL24A1",	"DRD2",	"EBF3",	"EPHA3",	"FAT4",	"FBN2",	"IL1RAPL2",	"ITGA6",	"MAF",	"NABP1",	"NELL1",	"NFIB",	"NPY",	"PDGFRα",	"POSTN",	"RELN",	"SAC",	"TRPC4",	"VEGFC",	"VGluT3"),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 7)#+ggplot2:::coord_flip()
jjplot
ggsave(filename = "jjDotPlot2.pdf", plot = jjplot, device = 'pdf', width = 32, height = 19, units = 'cm')

# SAC
markers.to.plot <- c("MEGF11","CHAT")
DotPlot(scRNA, idents = "SAC", features = markers.to.plot, cols = c("lightgrey", "red"), group.by = "age", col.min = -1,col.max = 1,dot.scale = 10) #+ggplot2:::coord_flip()
DotPlot(scRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "celltype", dot.scale = 8) + RotatedAxis()
DotPlot(scRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "seurat_clusters", dot.scale = 8) +  RotatedAxis()
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "celltype", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()

ggsave(filename = "DotPlot.pdf", device = 'pdf', width = 11, height = 7, units = 'cm')

################################################################################
#...cellchat...
################################################################################
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(Seurat)
library(ggpubr)
library(ggalluvial)

####
scRNA@meta.data[["age"]]<-factor(scRNA@meta.data[["age"]], levels=c("Child","Adult","Old"))
table(scRNA@meta.data[["celltype"]], scRNA@meta.data[["age"]])
age.list<-SplitObject(scRNA, split.by = "age")

# Create a CellChat object for the old, adult and child groups separately
data.input  <- age.list[["Old"]]@assays[["RNA"]]@data
celltype  <- age.list[["Old"]]@meta.data[["celltype"]]
data.input[1:2,1:2]
identity = data.frame(group = age.list[["Old"]]$celltype, row.names = names(age.list[["Old"]]$celltype)) # create a dataframe consisting of the cell labels
head(identity)
unique(identity$group) # check the cell labels
table(identity$group)
meta <- data.frame(labels = age.list[["Old"]]$celltype, row.names = names(identify))
cellchat <- createCellChat(object = data.input,meta = meta, group.by = "labels")
cellchat
summary(cellchat)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
head(cellchat@meta)
cellchat <- setIdent(cellchat, ident.use = "labels") 
CellChatDB <- CellChatDB.human 
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB) #"Secreted Signaling" "ECM-Receptor" "Cell-Cell Contact"
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human) 
cellchat <- computeCommunProb(cellchat,type = "triMean", population.size = F)  
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

Child <- cellchat
Adult <- cellchat
Old <- cellchat

saveRDS(Child,"Child.rds")
saveRDS(Adult,"Adult.rds")
saveRDS(Old,"Old.rds")

load("Child.rds")
load("Adult.rds")
load("Old.rds")

# merge CellChat data
con.list <- list(Child=Child,Adult=Adult,Old=Old)
cellchat <- mergeCellChat(con.list,add.names = names(con.list),cell.prefix = T)

# visualization
gg1 <- compareInteractions(cellchat, show.legend = F, color.use = c("#F3E0DF","#AFD8AD","#89ACC7") ,group = c(1,2,3),measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, color.use = c("#F3E0DF","#AFD8AD","#89ACC7") ,group = c(1,2,3),measure = "weight")
p <- gg1+gg2
ggsave("overview_number_strength.pdf", p, width = 6,height = 4)
dev.off()

par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat ,comparison = c(1, 2))
h2 <- netVisual_heatmap(cellchat, comparison = c(1, 2), measure = "weight")
pdf("Diff_number_strength_heatmap1_2.pdf", width = 10,height = 5.5)
h1+h2
dev.off()

par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat ,comparison = c(2, 3))
h2 <- netVisual_heatmap(cellchat, comparison = c(2, 3), measure = "weight")
pdf("Diff_number_strength_heatmap2_3.pdf", width = 10,height = 5.5)
h1+h2
dev.off()

gg1 <- rankNet(cellchat, mode = "comparison", color.use = c("#F3E0DF","#AFD8AD","#89ACC7"), comparison = c(1,2,3),stacked = T, do.stat = T)
gg2 <- rankNet(cellchat, mode = "comparison", color.use = c("#F3E0DF","#AFD8AD","#89ACC7"), comparison = c(1,2,3),stacked = F, do.stat = T)
p <- gg1 + gg2
ggsave("compare_pathway_strength.pdf", p,width = 10,height = 5.5)

library(ComplexHeatmap)
cellchat <- sham###sham h1 h4 h12 d3 d7
groupSize <- as.numeric(table(cellchat@idents))

pdf("sham_Number_weight.pdf", width = 20,height = 11)
par(mfrow = c(1,2), xpd=T)
s1 <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, edge.width.max = 13,title.name = "sham_Number of interactions")
s2 <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, edge.width.max = 13,title.name = "sham_Interaction weights")
dev.off()

library(ComplexHeatmap)
pathway.union <- union(con.list[[1]]@netP$pathways, c(con.list[[2]]@netP$pathways, con.list[[3]]@netP$pathways))

ht1 = netAnalysis_signalingRole_heatmap(con.list[[1]], pattern = "all", signaling = pathway.union,title = names(con.list)[1], width = 8, height = 10)

ht2 = netAnalysis_signalingRole_heatmap(con.list[[2]], pattern = "all", signaling = pathway.union,title = names(con.list)[2], width = 8, height = 10)

ht3 = netAnalysis_signalingRole_heatmap(con.list[[3]], pattern = "all", signaling = pathway.union,title = names(con.list)[3], width = 8, height = 10)

pdf("compare_signal_pattern_all.pdf", width = 15,height = 5.5)
draw(ht1 + ht2 + ht3 ,ht_gap = unit(0.5, "cm"))
garbage <- dev.off()

ht1 = netAnalysis_signalingRole_heatmap(con.list[[1]], pattern = "outgoing", signaling = pathway.union,title = names(con.list)[1], width = 8, height = 10)

ht2 = netAnalysis_signalingRole_heatmap(con.list[[2]], pattern = "outgoing", signaling = pathway.union,title = names(con.list)[2], width = 8, height = 10)

ht3 = netAnalysis_signalingRole_heatmap(con.list[[3]], pattern = "outgoing", signaling = pathway.union,title = names(con.list)[3], width = 8, height = 10)

pdf("compare_signal_pattern_outgoing.pdf", width = 15,height = 5.5)
draw(ht1 + ht2 + ht3 ,ht_gap = unit(0.5, "cm"))
garbage <- dev.off()

ht1 = netAnalysis_signalingRole_heatmap(con.list[[1]], pattern = "incoming", signaling = pathway.union,title = names(con.list)[1], width = 8, height = 10)

ht2 = netAnalysis_signalingRole_heatmap(con.list[[2]], pattern = "incoming", signaling = pathway.union,title = names(con.list)[2], width = 8, height = 10)

ht3 = netAnalysis_signalingRole_heatmap(con.list[[3]], pattern = "incoming", signaling = pathway.union,title = names(con.list)[3], width = 8, height = 10)

pdf("compare_signal_pattern_incoming.pdf", width = 15,height = 5.5)
draw(ht1 + ht2 + ht3 ,ht_gap = unit(0.5, "cm"))
garbage <- dev.off()

pathway.union
a <- cellchat@netP[["Child"]][["pathways"]]
b <- cellchat@netP[["Adult"]][["pathways"]]
c <- cellchat@netP[["Old"]][["pathways"]]
a_list <- list(a,b,c)
x <- do.call(cbind,lapply(lapply(a_list, unlist), `length<-`, max(lengths(a_list))))
colnames(x) <- c('Child','Adult','Old')
write.csv(x,file = "pathways.csv",row.names = F)

#
pathways.show <-c ("MAG")
weight.max <- getMaxWeight(con.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(6,1), xpd=TRUE)
pdf("MAG_circle.pdf", width = 5.5,height = 5.5)
for (i in 1:length(con.list)) {
  netVisual_aggregate(con.list[[i]],signaling = pathways.show, layout = "circle",
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      vertex.weight = as.numeric(table(con.list[[i]]@idents)),
                      signaling.name = paste(pathways.show, names(con.list)[i]))  
}
dev.off()

pdf("Compare_signal_net.pdf", width = 10,height = 10)
par(mfrow = c(2,2))
h1 <- netVisual_diffInteraction(cellchat,edge.weight.max = 10,comparison = c(1, 2), weight.scale = T)
h2 <- netVisual_diffInteraction(cellchat, measure = "weight",comparison = c(1, 2),vertex.size.max = 15, weight.scale = T)
h3 <- netVisual_diffInteraction(cellchat,edge.weight.max = 10,comparison = c(2, 3), weight.scale = T)
h4 <- netVisual_diffInteraction(cellchat, measure = "weight",comparison = c(2, 3),vertex.size.max = 15, weight.scale = T)
dev.off()

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

type="Amacrine"

# (1) load Seurat data for Amacrine
X = io.mmread(f"/home/xll-zxcv/TS/loom/Seurat_dataset/Amacrine/counts.mtx")
# create anndata object
adata = anndata.AnnData(X=X.transpose().tocsr())

# (2) load cell metadata:
cell_meta = pd.read_csv(f"/home/xll-zxcv/TS/loom/Seurat_dataset/Amacrine/metadata.csv")

# load gene names:
with open(f"/home/xll-zxcv/TS/loom/Seurat_dataset/Amacrine/gene_names.csv", 'r') as f:
  gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# (3) load dimensional reduction:
pca = pd.read_csv(f"/home/xll-zxcv/TS/loom/Seurat_dataset/Amacrine/pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:
sc.pl.umap(adata, color='seurat_clusters', frameon=False, save=True)

# (4) load loom data
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

adata.write_h5ad(f'/home/xll-zxcv/TS/loom/Seurat_dataset/Amacrine/Amacrine.adata_ldata.h5ad')

# (6) scVelo
adata = sc.read(f'/home/xll-zxcv/TS/loom/Seurat_dataset/Amacrine/Amacrine.adata_ldata.h5ad')

#scv.pp.filter_and_normalize(adata, min_shared_counts=5, min_shared_cells=3, log=True)
scv.pp.filter_and_normalize(adata)

scv.pp.moments(adata, n_neighbors=30, n_pcs=30)

import gc
gc.collect()
temp_pre= f"Amacrine_nue.in_process2"
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
scv.pl.velocity_embedding(adata, basis="umap", save='embedding_Amacrine.pdf')

scv.pl.velocity_embedding_stream(adata, basis='umap', color=
                                   ['seurat_clusters'], figsize=(8,8), palette =
                                   ('#85A0D2','#68A9E3',"#716EB5","#63AF70","#A2C260"),
                                 arrow_size=2, linewidth=1.5, legend_fontsize=25, dpi=900,
                                 save='embedding_stream_seurat_clusters_Amacrine.svg', title='')

scv.pl.velocity_embedding_stream(adata, basis='umap', color=
                                   ['celltype'], figsize=(8,8), palette = ('#85A0D2','#68A9E3',"#716EB5","#63AF70","#A2C260"), alpha=0.1, 
                                 arrow_size=1.2, linewidth=0.5, density=2.5, size=80, legend_fontsize=17, dpi=900,
                                 save='embedding_stream_celltype_Amacrine.svg', title='')

scv.pl.velocity_embedding_grid(adata, basis='umap', color='celltype', figsize=
                                 (8,8), arrow_size=1.8,  arrow_length=2.5, size=100, alpha=0.1,
                               save='embedding_grid_Amacrine.pdf', title='', scale=0.1)

scv.pl.velocity_embedding(adata, arrow_length=2, arrow_size=1.5, dpi=900,
                          figsize=(8,8), save='embedding_stream2_Amacrine.pdf')

scv.tl.velocity_confidence(adata)
keys = 'velocity_length','velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], dpi=900, figsize=(6,5), save='Speed and coherence_Amacrine.pdf')

df = adata.obs.groupby('celltype')[keys].mean().T
df.style.background_gradient(cmap='coolwarm', axis=1)

scv.pl.velocity_graph(adata, threshold=.1, color='celltype', figsize=(8,8), dpi=900, save='velocity_graph_Amacrine.pdf')

scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', figsize=(5,5), dpi=300 ,save='velocity_pseudotime_Amacrine.pdf')

###PAGA
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='celltype')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.pl.paga(adata, basis='umap', size=80, alpha=.15,
            min_edge_width=2, node_size_scale=1.5, figsize=(6,5),arrowsize=20, node_size_power=0.5,dpi=900, save='paga_Amacrine.svg')

##Latent time
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, dpi=900, save='latent_time_Amacrine.pdf')
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]

scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='celltype', n_convolve=100, save='heatmap_Amacrine.pdf')

