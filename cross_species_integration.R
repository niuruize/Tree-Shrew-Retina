################################################################################
#cross_species analysis
################################################################################
#loading R packages
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggplot2)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(MySeuratWrappers)

setwd("/home/fingerstyle/TS/cor_species")
################################################################################
# load data
################################################################################
###**********************TS*********************###
load("01TS/scRNA.RData")
count <- scRNA@assays$RNA@counts
meta <- scRNA@meta.data[,c("orig.ident","age","gender","celltype")]
TS = CreateSeuratObject(counts = count, assay = 'RNA', meta.data = meta, project = "TS", min.cells = 3, min.features = 200)
TS$Species="TS"
save(TS,file = "01TS/TS.Rdata")

###**********************Human*********************###
counts <-read.csv('02Human/GSE148077_count_mat_donor_H1.csv.gz', row.names=1,header = TRUE)
H1= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "H1", min.cells = 3, min.features = 600)
H1 <- subset(H1, subset = nFeature_RNA > 600 & nFeature_RNA < 7000)

counts <-read.csv('02Human/GSE148077_count_mat_donor_H2.csv.gz', row.names=1,header = TRUE)
H2= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "H2", min.cells = 3, min.features = 600)
H2 <- subset(H2, subset = nFeature_RNA > 600 & nFeature_RNA < 7000)

counts <-read.csv('02Human/GSE148077_count_mat_donor_H3.csv.gz', row.names=1,header = TRUE)
H3= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "H3", min.cells = 3, min.features = 600)
H3 <- subset(H3, subset = nFeature_RNA > 600 & nFeature_RNA < 7000)

counts <-read.csv('02Human/GSE148077_count_mat_donor_H4.csv.gz', row.names=1,header = TRUE)
H4= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "H4", min.cells = 3, min.features = 600)
H4 <- subset(H4, subset = nFeature_RNA > 600 & nFeature_RNA < 7000)

counts <-read.csv('02Human/GSE148077_count_mat_donor_H5.csv.gz', row.names=1,header = TRUE)
H5= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "H5", min.cells = 3, min.features = 600)
H5 <- subset(H5, subset = nFeature_RNA > 600 & nFeature_RNA < 7000)

counts <-read.csv('02Human/GSE148077_count_mat_donor_H9.csv.gz', row.names=1,header = TRUE)
H9= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "H9", min.cells = 3, min.features = 600)
H9 <- subset(H9, subset = nFeature_RNA > 600 & nFeature_RNA < 7000)

counts <-read.csv('02Human/GSE148077_count_mat_donor_H11.csv.gz', row.names=1,header = TRUE)
H11= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "H11", min.cells = 3, min.features = 600)
H11 <- subset(H11, subset = nFeature_RNA > 600 & nFeature_RNA < 7000)

Human <- merge(H1,list(H2,H3,H4,H5,H9,H11))
Human$Species="Human"
save(Human,file = "02Human/Human.Rdata")

###**********************Macaque*********************###
counts <-read.csv('03Macaque/GSE118852_CountMatrix_M4perCD73.csv', row.names=1,header = TRUE)
Per4= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "Per4", min.cells = 3, min.features = 500)
Per4 <- subset(Per4, subset = nFeature_RNA > 500 & nFeature_RNA < 6000)

counts <-read.csv('03Macaque/GSE118852_CountMatrix_M5perCD73.csv', row.names=1,header = TRUE)
Per5a= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "Per5a", min.cells = 3, min.features = 500)
Per5a <- subset(Per5a, subset = nFeature_RNA > 500 & nFeature_RNA < 6000)

counts <-read.csv('03Macaque/GSE118852_CountMatrix_M5perCD90.csv', row.names=1,header = TRUE)
Per5b= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "Per5b", min.cells = 3, min.features = 500)
Per5b <- subset(Per5b, subset = nFeature_RNA > 500 & nFeature_RNA < 6000)

counts <-read.csv('03Macaque/GSE118852_CountMatrix_M5perPNA.csv', row.names=1,header = TRUE)
Per5c= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "Per5c", min.cells = 3, min.features = 500)
Per5c <- subset(Per5c, subset = nFeature_RNA > 500 & nFeature_RNA < 6000)

counts <-read.csv('03Macaque/GSE118852_CountMatrix_M6perCD73.csv', row.names=1,header = TRUE)
Per6a= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "Per6a", min.cells = 3, min.features = 500)
Per6a <- subset(Per6a, subset = nFeature_RNA > 500 & nFeature_RNA < 6000)

counts <-read.csv('03Macaque/GSE118852_CountMatrix_M6perCD90.csv', row.names=1,header = TRUE)
Per6b= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "Per6b", min.cells = 3, min.features = 500)
Per6b <- subset(Per6b, subset = nFeature_RNA > 500 & nFeature_RNA < 6000)

counts <-read.csv('03Macaque/GSE118852_CountMatrix_M6perMixed.csv', row.names=1,header = TRUE)
Per6c= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "Per6c", min.cells = 3, min.features = 500)
Per6c <- subset(Per6c, subset = nFeature_RNA > 500 & nFeature_RNA < 6000)

counts <-read.csv('03Macaque/GSE118852_CountMatrix_M7perCD90.csv', row.names=1,header = TRUE)
Per7= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "Per7", min.cells = 3, min.features = 500)
Per7 <- subset(Per7, subset = nFeature_RNA > 500 & nFeature_RNA < 6000)

counts <-read.csv('03Macaque/GSE118546_CountMatrix_M1fovea.csv.gz', row.names=1,header = TRUE)
Fovea1= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "Fovea1", min.cells = 3, min.features = 500)
Fovea1 <- subset(Fovea1, subset = nFeature_RNA > 500 & nFeature_RNA < 6000)

counts <-read.csv('03Macaque/GSE118546_CountMatrix_M2fovea.csv.gz', row.names=1,header = TRUE)
Fovea2= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "Fovea2", min.cells = 3, min.features = 500)
Fovea2 <- subset(Fovea2, subset = nFeature_RNA > 500 & nFeature_RNA < 6000)

counts <-read.csv('03Macaque/GSE118546_CountMatrix_M3fovea.csv.gz', row.names=1,header = TRUE)
Fovea3= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "Fovea3", min.cells = 3, min.features = 500)
Fovea3 <- subset(Fovea3, subset = nFeature_RNA > 500 & nFeature_RNA < 6000)

counts <-read.csv('03Macaque/GSE118546_CountMatrix_M4fovea.csv.gz', row.names=1,header = TRUE)
Fovea4= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "Fovea4", min.cells = 3, min.features = 500)
Fovea4 <- subset(Fovea4, subset = nFeature_RNA > 500 & nFeature_RNA < 6000)

Macaque <- merge(Per4,list(Per5a,Per5b,Per5c,Per6a,Per6b,Per6c,Per7,Fovea1,Fovea2,Fovea3,Fovea4))
Macaque$Species="Macaque"
save(Macaque,file = "03Macaque/Macaque.Rdata")

###**********************Mouse*********************###
Mouse <- read.table("04Mouse/GSE63472_P14Retina_merged_digital_expression.txt", sep="", header = T,row.names = 1)
Mouse <- CreateSeuratObject(counts = Mouse, project = "Mouse", min.cells = 3, min.features = 200)
Mouse <- subset(Mouse, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)
Mouse$Species="Mouse"
save(Mouse,file = "04Mouse/Mouse.Rdata")

###**********************Chick*********************###
counts <-read.csv('05Chick/GSE159107_E16chick_count.matrix.csv.gz', row.names=1,header = TRUE)
E16= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "E16", min.cells = 3, min.features = 600)
counts <-read.csv('05Chick/GSE159107_E18chick_count.matrix.csv.gz', row.names=1,header = TRUE)
E18= CreateSeuratObject(counts = counts, assay = 'RNA',  project = "E18", min.cells = 3, min.features = 600)
Chick <- merge(E16,list(E18))
Chick$Species="Chick"
save(Chick,file = "05Chick/Chick.Rdata")

################################################################################
# Removal of double cells
################################################################################
#Human---
table(Human@meta.data$orig.ident)
samples <- c("H11FoveaS1","H1CD73dpS1","H1CD90S1","H2Fovea1S1","H2Fovea2S1","H3CD73dpS1","H3CD73dpS2","H3CD90S1","H3CD90S2","H3FoveaS1","H3FoveaS2","H4FoveaS1","H5FoveaS1","H5FoveaS2","H5FoveaS3","H5FoveaS4","H5FoveaS5","H9FoveaS1")
scRNA <- Human
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, split.by="orig.ident", do.center=FALSE)
nFactors=20
scRNA <- RunOptimizeALS(scRNA, k=nFactors, split.by="orig.ident")
scRNA <- RunQuantileNorm(scRNA, split.by="orig.ident")
scRNA$clusters <- factor(scRNA$clusters,levels=1:length(levels(scRNA$clusters)))
save(scRNA,file="Singlet/Human/scRNA-1.RData")
scRNA <- FindNeighbors(scRNA, reduction="iNMF", dims=1:20)
scRNA <- FindClusters(scRNA, resolution = 0.4)
scRNA <- RunUMAP(scRNA, dims=1:nFactors, reduction="iNMF")
save(scRNA,file="Singlet/Human/scRNA-2.RData")
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
  file_name <- paste("Singlet/Human", sample, "_DoubletFinder_result.txt", sep = "")
  write.table(DF_df, file = file_name, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  DimPlot(test.seu, reduction = "umap", pt.size = 0.3, group.by = "DoubletFinder")
  
  table(test.seu@meta.data[["DoubletFinder"]])
  assign(sample, test.seu[, test.seu@meta.data$DoubletFinder %in% c("Singlet")])
  
  save(test.seu, file = paste("Singlet/Human/", sample, "-Singlet.RData", sep = ""))
}

Human_single <- merge(H11FoveaS1, y=c(H1CD73dpS1, H1CD90S1, H2Fovea1S1, H2Fovea2S1, H3CD73dpS1, H3CD73dpS2, H3CD90S1, H3CD90S2, H3FoveaS1,  H3FoveaS2,  H3FoveaS3,  H4FoveaS1,  H5FoveaS1,  H5FoveaS2,  H5FoveaS3,  H5FoveaS4,  H5FoveaS5,  H9FoveaS1))
save(Human_single,file="Singlet/Human/Human-single.RData")

#Macaque----
table(Macaque@meta.data$orig.ident)
samples <- c("Fovea4S1","Fovea4S2","Fovea4S3","M1Fovea1","M1Fovea2","M1Fovea3","M1Fovea4","M1Fovea5","M1Fovea6","M1Fovea7","M1Fovea8","M2Fovea1","M2Fovea2","M2Fovea3","M2Fovea4","M2Fovea5","M2Fovea6","M2Fovea7","M2Fovea8","M3Fovea1","M3Fovea2","M3Fovea3","MacaqueCD73DP2S1","MacaqueCD73DP2S2","PerCd73S1","PerCd73S2","PerCd73S3","PerCd73S4","PerCd90PNAS1","PerCd90S1","PerCd90S2","PerCd90S3","PerCd90S4","PerCd90S5","PerCd90S6","PerCd90S7","PerCd90S8","PerCd90S9","PerMixedS1")
scRNA <- Macaque
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, split.by="orig.ident", do.center=FALSE)
nFactors=20  
scRNA <- RunOptimizeALS(scRNA, k=nFactors, split.by="orig.ident")
scRNA <- RunQuantileNorm(scRNA, split.by="orig.ident")
scRNA$clusters <- factor(scRNA$clusters,levels=1:length(levels(scRNA$clusters)))
save(scRNA,file="Singlet/Macaque/scRNA-1.RData")
scRNA <- FindNeighbors(scRNA, reduction="iNMF", dims=1:20)
scRNA <- FindClusters(scRNA, resolution = 0.4)
scRNA <- RunUMAP(scRNA, dims=1:nFactors, reduction="iNMF")
save(scRNA,file="Singlet/Macaque/scRNA-2.RData")
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
  file_name <- paste("Singlet/Macaque/", sample, "_DoubletFinder_result.txt", sep = "")
  write.table(DF_df, file = file_name, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  DimPlot(test.seu, reduction = "umap", pt.size = 0.3, group.by = "DoubletFinder")
  
  table(test.seu@meta.data[["DoubletFinder"]])
  assign(sample, test.seu[, test.seu@meta.data$DoubletFinder %in% c("Singlet")])
  
  save(test.seu, file = paste("Singlet/Macaque/", sample, "-Singlet.RData", sep = ""))
}

Macaque_single <- merge(Fovea4S1, y=c(Fovea4S2, Fovea4S3, M1Fovea1, M1Fovea2, M1Fovea3, M1Fovea4, M1Fovea5, M1Fovea6, M1Fovea7, M1Fovea8, M2Fovea1, M2Fovea2, M2Fovea3, M2Fovea4, M2Fovea5, M2Fovea6, M2Fovea7, M2Fovea8, M3Fovea1, M3Fovea2, M3Fovea3, MacaqueCD73DP2S1, MacaqueCD73DP2S2, PerCd73S1,  PerCd73S2,  PerCd73S3,  PerCd73S4,  PerCd90PNAS1, PerCd90S1,  PerCd90S2,  PerCd90S3,  PerCd90S4,  PerCd90S5,  PerCd90S6,  PerCd90S7,  PerCd90S8,  PerCd90S9,  PerMixedS1))

save(Macaque_single,file="Singlet/Macaque/Macaque-single.RData")

#Mouse----
table(Mouse@meta.data$orig.ident)
samples <- c("p1","r1","r2","r3","r4","r5","r6")
scRNA <- Mouse
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, split.by="orig.ident", do.center=FALSE)
nFactors=20 
scRNA <- RunOptimizeALS(scRNA, k=nFactors, split.by="orig.ident")
scRNA <- RunQuantileNorm(scRNA, split.by="orig.ident")
scRNA$clusters <- factor(scRNA$clusters,levels=1:length(levels(scRNA$clusters)))
save(scRNA,file="Singlet/Mouse/scRNA-1.RData")
scRNA <- FindNeighbors(scRNA, reduction="iNMF", dims=1:20)
scRNA <- FindClusters(scRNA, resolution = 0.4)
scRNA <- RunUMAP(scRNA, dims=1:nFactors, reduction="iNMF")
save(scRNA,file="Singlet/Mouse/scRNA-2.RData")
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
  file_name <- paste("Singlet/Mouse/", sample, "_DoubletFinder_result.txt", sep = "")
  write.table(DF_df, file = file_name, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  DimPlot(test.seu, reduction = "umap", pt.size = 0.3, group.by = "DoubletFinder")
  
  table(test.seu@meta.data[["DoubletFinder"]])
  assign(sample, test.seu[, test.seu@meta.data$DoubletFinder %in% c("Singlet")])
  
  save(test.seu, file = paste("Singlet/Mouse/", sample, "-Singlet.RData", sep = ""))
}

Mouse_single <- merge(p1, y=c(r1, r2, r3, r4, r5, r6))
save(Mouse_single,file="Singlet/Mouse/Mouse-single.RData")

#Chick----
table(Chick@meta.data$orig.ident)
samples <- c("Chicken1A","Chicken1B","Chicken1C","Chicken1D","ChickendRGC1","ChickendRGC2","ChickenvRGC1","ChickenvRGC2")
scRNA <- Chick
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, split.by="orig.ident", do.center=FALSE)
nFactors=20
scRNA <- RunOptimizeALS(scRNA, k=nFactors, split.by="orig.ident")
scRNA <- RunQuantileNorm(scRNA, split.by="orig.ident")
scRNA$clusters <- factor(scRNA$clusters,levels=1:length(levels(scRNA$clusters)))
save(scRNA,file="Singlet/Chick/scRNA-1.RData")
scRNA <- FindNeighbors(scRNA, reduction="iNMF", dims=1:20)
scRNA <- FindClusters(scRNA, resolution = 0.4)
scRNA <- RunUMAP(scRNA, dims=1:nFactors, reduction="iNMF")
save(scRNA,file="Singlet/Chick/scRNA-2.RData")
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
  file_name <- paste("Singlet/Chick/", sample, "_DoubletFinder_result.txt", sep = "")
  write.table(DF_df, file = file_name, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  DimPlot(test.seu, reduction = "umap", pt.size = 0.3, group.by = "DoubletFinder")
  
  table(test.seu@meta.data[["DoubletFinder"]])
  assign(sample, test.seu[, test.seu@meta.data$DoubletFinder %in% c("Singlet")])
  
  save(test.seu, file = paste("Singlet/Chick/", sample, "-Singlet.RData", sep = ""))
}

Chick_single <- merge(Chicken1A, y=c(Chicken1B, Chicken1C,  Chicken1D,  ChickendRGC1, ChickendRGC2, ChickenvRGC1, ChickenvRGC2))
save(Chick_single,file="Singlet/Chick/Chick-single.RData")

#
load("Singlet/TS/TS-single.RData")
TS <- TS
load("Singlet/Human/Human-single.RData")
Human <- Human_single
load("Singlet/Macaque/Macaque-single.RData")
Macaque <- Macaque_single
load("Singlet/Mouse/Mouse-single.RData")
Mouse <- Mouse_single
load("Singlet/Chick/Chick-single.RData")
Chick <- Chick_single

################################################################################
# Homologous gene conversion
################################################################################
gene_set1 <- data.frame(gene=row.names(TS))
gene_set2 <- data.frame(gene=row.names(Human))
gene_set3 <- data.frame(gene=row.names(Macaque))
gene_set4 <- data.frame(gene=row.names(Mouse))
gene_set4$gene <- str_to_title(gene_set4$gene)
gene_set5 <- data.frame(gene=row.names(Chick))

TS_human = homologene(gene_set1$gene, inTax = 9606, outTax = 9606)
human_human = homologene(gene_set2$gene, inTax = 9606, outTax = 9606)
Macaque_human = homologene(gene_set3$gene, inTax = 9544, outTax = 9606)
Mouse_human = homologene(gene_set4$gene, inTax = 10090, outTax = 9606)
Chick_human = homologene(gene_set5$gene, inTax = 9031, outTax = 9606)

#Remove duplicate symbols
duplicated(TS_human[,1])
TS_human[duplicated(TS_human[,1]) == F,]
TS_human_single <- TS_human[duplicated(TS_human[,1]) == F,]

duplicated(human_human[,1])
human_human[duplicated(human_human[,1]) == F,]
human_human_single <- human_human[duplicated(human_human[,1]) == F,]

duplicated(Macaque_human[,1])
Macaque_human[duplicated(Macaque_human[,1]) == F,]
Macaque_human_single <- Macaque_human[duplicated(Macaque_human[,1]) == F,]

duplicated(Mouse_human[,1])
Mouse_human[duplicated(Mouse_human[,1]) == F,]
Mouse_human_single <- Mouse_human[duplicated(Mouse_human[,1]) == F,]
Mouse_human_single[,1] <- toupper(Mouse_human_single[,1])

duplicated(Chick_human[,1])
Chick_human[duplicated(Chick_human[,1]) == F,]
Chick_human_single <- Chick_human[duplicated(Chick_human[,1]) == F,]

RenameGenesSeurat <- function(obj,newnames,gene.use=NULL,de.assay="RNA") { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}

TS <- RenameGenesSeurat(obj=TS,newnames=TS_human_single[,2],gene.use=TS_human_single[,1])
Human <- RenameGenesSeurat(obj=Human,newnames=human_human_single[,2],gene.use=human_human_single[,1])
Macaque <- RenameGenesSeurat(obj=Macaque,newnames=Macaque_human_single[,2],gene.use=Macaque_human_single[,1])
Mouse <- RenameGenesSeurat(obj=Mouse,newnames=Mouse_human_single[,2],gene.use=Mouse_human_single[,1])
Chick <- RenameGenesSeurat(obj=Chick,newnames=Chick_human_single[,2],gene.use=Chick_human_single[,1])

# Extract the list of intersection genes
gene_intersect <- Reduce(intersect, list(TS_human_single[,2], human_human_single[,2], Macaque_human_single[,2], Mouse_human_single[,2], Chick_human_single[,2]))

#Recreate the seurat object
#TS
TS_data <- TS@assays$RNA@data[TS@assays$RNA@data@Dimnames[[1]] %in% gene_intersect,]
meta <- TS@meta.data[,c("orig.ident","age","gender","celltype","Species")]
TS = CreateSeuratObject(counts = TS_data, assay = 'RNA', meta.data = meta, project = "TS", min.cells = 3, min.features = 200)
#Human
Human_data <- Human@assays$RNA@data[Human@assays$RNA@data@Dimnames[[1]] %in% gene_intersect,]
meta <- Human@meta.data[,c("orig.ident","Species")]
Human = CreateSeuratObject(counts = Human_data, assay = 'RNA', meta.data = meta, project = "Human", min.cells = 3, min.features = 200)
#Macaque
Macaque_data <- Macaque@assays$RNA@data[Macaque@assays$RNA@data@Dimnames[[1]] %in% gene_intersect,]
meta <- Macaque@meta.data[,c("orig.ident","Species")]
Macaque = CreateSeuratObject(counts = Macaque_data, assay = 'RNA', meta.data = meta, project = "Macaque", min.cells = 3, min.features = 200)
#Mouse
Mouse_data <- Mouse@assays$RNA@data[Mouse@assays$RNA@data@Dimnames[[1]] %in% gene_intersect,]
meta <- Mouse@meta.data[,c("orig.ident","Species")]
Mouse = CreateSeuratObject(counts = Mouse_data, assay = 'RNA', meta.data = meta, project = "Mouse", min.cells = 3, min.features = 200)
#Chick
Chick_data <- Chick@assays$RNA@data[Chick@assays$RNA@data@Dimnames[[1]] %in% gene_intersect,]
meta <- Chick@meta.data[,c("orig.ident","Species")]
Chick = CreateSeuratObject(counts = Chick_data, assay = 'RNA', meta.data = meta, project = "Chick", min.cells = 3, min.features = 200)

#save data
save(TS,file = "06homologene/09CCA/TS.Rdata")
save(Human,file = "06homologene/09CCA/Human.Rdata")
save(Macaque,file = "06homologene/09CCA/Macaque.Rdata")
save(Mouse,file = "06homologene/09CCA/Mouse.Rdata")
save(Chick,file = "06homologene/09CCA/Chick.Rdata")

#merge data
scRNA <- merge(TS,y=c(Human,Macaque,Mouse,Chick))
table(scRNA@meta.data$orig.ident)
table(scRNA@meta.data$Species)
save(scRNA,file = "06homologene/09CCA/Merge.Rdata")

################################################################################
# cross_species integration
################################################################################
scRNAlist <- SplitObject(scRNA, split.by = "orig.ident")
scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
})
scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
  x <- ScaleData(x, features = scRNA.features, verbose = FALSE)
  x <- RunPCA(x, features = scRNA.features, verbose = FALSE)
})

scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, reference = c(1, 6, 11, 16, 17,18,37,45,41,53,56,58,62,63,72,73,74,83,81,86,87),  normalization.method = "LogNormalize",reduction = "rpca", anchor.features = scRNA.features, dims = 1:50)  ##耗时久
scRNA <- IntegrateData(scRNA.anchors, normalization.method = "LogNormalize", dims = 1:50) 
DefaultAssay(scRNA) <- "integrated"
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, npcs = 100, verbose = FALSE)
save(scRNA,file = "09CCA/scRNA-2.RData")
ElbowPlot(scRNA, ndims = 50)
scRNA <- FindNeighbors(scRNA, dims = 1:30)
scRNA <- FindClusters(scRNA, resolution = 0.4)
scRNA <- RunUMAP(scRNA, dims = 1:30)
scRNA@meta.data[["Species"]]<-factor(scRNA@meta.data[["Species"]], levels=c("TS","Human","Macaque","Mouse","Chick"))
DefaultAssay(scRNA) <- "RNA"
save(scRNA,file = "09CCA/scRNA-3.RData")

#
table(scRNA@meta.data$orig.ident)
table(scRNA@meta.data$Species)
table(scRNA@meta.data$seurat_clusters)
table(scRNA@meta.data[["seurat_clusters"]], scRNA@meta.data[["Species"]])

#visualization
p10 <- DimPlot(scRNA, reduction = "umap", group.by = "ident", split.by = "Species", pt.size=0.01, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "09CCA/umap_Cluster.pdf", plot = p10, device = 'pdf', width = 65, height = 15, units = 'cm')

markers.to.plot <- c("OPN1LW","OPN1SW","ARR3","GRK7","PDE6H","PDE6C","GNAT2","RCVRN","CRX","PDE6B","PDE6A","ROM1","SAG","CA10","GRIK1","ISL1","TRPM1","NETO1","ONECUT1","ONECUT2","NDRG1","RLBP1","WIF1","CHRDL1","TF","GLUL","C1QA","LGMN","CSF1R","CX3CR1","PTPRC","AQP4","GFAP","SLC1A3","SLC14A1","NRN1","SLC17A6","RBPMS","THY1","PAX6","SYNPR","NRXN1","MEIS2","TFAP2A","TFAP2B","C1QL2","SLC6A1","TCF4","SLC6A9","SLC17A8","NFIA","NFIB","SLC18A3","MEGF11","RPE65","IGFBP7","VWF","PDGFRB","PLP1","CSPG4","PDGFRA","OTX2","MPZ","VSX1","CHAT","BEST2","CRHBP","COL18A1","ZIC1","OPTC","COL9A1","CPAMD8","ZIC2","HCN1","SYT1")
dp2 <- DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2, dot.scale = 8) + RotatedAxis() #+ggplot2:::coord_flip()
ggsave(filename = "09CCA/DotPlot.pdf", plot = dp2, device = 'pdf', width = 45, height = 30, units = 'cm')

#hclust
library(dendextend)
library(circlize)
library(Seurat)
scRNA$Species.cls <- paste(scRNA$Species, Idents(scRNA), sep = "-")
cluster.averages <- AverageExpression(scRNA,group.by = "Species.cls")
B.exp <- cluster.averages[["integrated"]]
hc <- as.dendrogram(hclust(dist(t(B.exp)), "ave"))
hc <- hc %>% color_branches(k = 60) %>% color_labels(k = 60)
pdf("09CCA/H-Cluster-Species.pdf", width = 35,height = 8)
plot(hc)
dev.off()
#
pdf("09CCA/C-Cluster-Species.pdf", width = 14,height = 14)
circlize_dendrogram(hc,labels_track_height = NA,dend_track_height = 0.4)
dev.off()

## Add cell labels
current.cluster.ids <- c("0", "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",  "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27")
new.cluster.ids <- c("Muller",  "RGC",  "Amacrine", "Cone", "Bipolar",  "Rod",  "Bipolar",  "Bipolar",  "Amacrine", "RGC",  "Bipolar",  "Amacrine", "Amacrine", "Amacrine", "Horizontal", "Bipolar",  "Amacrine", "Amacrine", "Bipolar",  "Pericytes",  "Amacrine", "Astrocytes", "Amacrine", "Amacrine", "Amacrine", "Microglia",  "Endothelial",  "Rod") 
scRNA@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA@meta.data$celltype)
Idents(scRNA) <- "celltype"
Idents(scRNA) <- factor(Idents(scRNA), levels = c("Cone", "Rod",  "Bipolar",  "Horizontal", "Amacrine", "RGC",  "Muller", "Astrocytes","Microglia", "Endothelial","Pericytes"))
DimPlot(scRNA, reduction = "umap", label = TRUE)
table(scRNA@meta.data[["celltype"]], scRNA@meta.data[["Species"]])
save(scRNA,file="09CCA/scRNA.RData")

#hclust
scRNA$Species.celltype <- paste(scRNA$Species, Idents(scRNA), sep = "-")
#
cluster.averages <- AverageExpression(scRNA,group.by = "Species.celltype")
B.exp <- cluster.averages[["integrated"]]
hc <- as.dendrogram(hclust(dist(t(B.exp)), "ave"))
hc <- hc %>% color_branches(k = 20) %>% color_labels(k = 20)
pdf("09CCA/H-celltype-Species.pdf", width = 18,height = 14)
plot(hc)
dev.off()
#
pdf("09CCA/C-celltype-Species.pdf", width = 8,height = 8)
circlize_dendrogram(hc,labels_track_height = NA,dend_track_height = 0.4)
dev.off()

#correlation analysis
table(scRNA$Species.celltype)  
av<-AverageExpression(scRNA,group.by = "Species.celltype", assays = "integrated")
av=av[[1]]
head(av)
cg=names(tail(sort(apply(av,1,sd)),1000))
View(av[cg,])
View(cor(av[cg,],method = "spearman"))
#
pdf("09CCA/cor_5Species-celltype.pdf", width = 12,height = 11.5)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),treeheight_row = 12,treeheight_col = 12,fontsize = 10)
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"09CCA/cor_Species-celltype.csv")

#DimPlot
p10 <- DimPlot(scRNA, cols = c('#8cb883','#643d2c','#FF8066','#c94052','#b5b736','#caa795','#80d1c8','#F7ab08',"#ddc07c",'#2C73D2','#845EC2'), reduction = "umap", group.by = "celltype", split.by = "Species", pt.size=0.01, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "09CCA/umap_celltype.pdf", plot = p10, device = 'pdf', width = 65, height = 15, units = 'cm')

p9 <- DimPlot(scRNA, cols = c('#8cb883','#643d2c','#FF8066','#c94052','#b5b736','#caa795','#80d1c8','#F7ab08',"#ddc07c",'#2C73D2','#845EC2'), reduction = "umap", group.by = "celltype", pt.size=0.01, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "09CCA/umap_celltype-all.pdf", plot = p9, device = 'pdf', width = 18.5, height = 15, units = 'cm')

#Remove the pericytes, endothelial,astrocyte and microglia
scRNA <- scRNA[,scRNA@meta.data[["celltype"]] %in% c("Cone",  "Rod",  "Bipolar",  "Horizontal", "Amacrine", "RGC",  "Muller")]
DefaultAssay(scRNA) <- "integrated"
ElbowPlot(scRNA, ndims = 50)
scRNA <- FindNeighbors(scRNA, dims = 1:25)
scRNA <- FindClusters(scRNA, resolution = 0.4)
scRNA <- RunUMAP(scRNA, dims = 1:25)
#Dimplot
p12 <- DimPlot(scRNA, reduction = "umap", group.by = "ident", split.by = "Species", pt.size=0.01, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "09CCA/umap_Cluster(no-glia).pdf", plot = p12, device = 'pdf', width = 65, height = 15, units = 'cm')
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("OPN1LW","OPN1SW","ARR3","GRK7","PDE6H","PDE6C","GNAT2","RCVRN","CRX","PDE6B","PDE6A","ROM1","SAG","CA10","GRIK1","ISL1","TRPM1","NETO1","ONECUT1","ONECUT2","NDRG1","RLBP1","WIF1","CHRDL1","TF","GLUL","C1QA","LGMN","CSF1R","CX3CR1","PTPRC","AQP4","GFAP","SLC1A3","SLC14A1","NRN1","SLC17A6","RBPMS","THY1","PAX6","SYNPR","NRXN1","MEIS2","TFAP2A","TFAP2B","C1QL2","SLC6A1","TCF4","SLC6A9","SLC17A8","NFIA","NFIB","SLC18A3","MEGF11","RPE65","IGFBP7","VWF","PDGFRB","PLP1","CSPG4","PDGFRA","OTX2","MPZ","VSX1","CHAT","BEST2","CRHBP","COL18A1","ZIC1","OPTC","COL9A1","CPAMD8","ZIC2","HCN1","SYT1")
dp2 <- DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2, dot.scale = 8) + RotatedAxis() #+ggplot2:::coord_flip()

#Re-label the cell type
current.cluster.ids <- c("0", "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",  "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24")
new.cluster.ids <- c("Muller",  "RGC",  "Amacrine", "Cone", "Bipolar",  "Rod",  "Bipolar",  "Bipolar",  "Amacrine", "Amacrine", "RGC",  "Bipolar",  "Bipolar",  "Amacrine", "Amacrine", "Horizontal", "Amacrine", "Bipolar",  "Amacrine", "Amacrine", "Amacrine", "Amacrine", "Amacrine", "Amacrine", "Rod") 

scRNA@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA@meta.data$celltype)
Idents(scRNA) <- "celltype"
Idents(scRNA) <- factor(Idents(scRNA), levels = c("Cone", "Rod",  "Bipolar",  "Horizontal", "Amacrine", "RGC",  "Muller"))
DimPlot(scRNA, reduction = "umap", label = TRUE)
table(scRNA@meta.data[["celltype"]], scRNA@meta.data[["Species"]])
save(scRNA,file="09CCA/scRNA(no-glia).RData")

#hclust
scRNA$Species.celltype <- paste(scRNA$Species, Idents(scRNA), sep = "-")
cluster.averages <- AverageExpression(scRNA,group.by = "Species.celltype")
B.exp <- cluster.averages[["integrated"]]
hc <- as.dendrogram(hclust(dist(t(B.exp)), "ave"))
hc <- hc %>% color_branches(k = 20) %>% color_labels(k = 20)
pdf("09CCA/H-celltype-Species(no-glia).pdf", width = 18,height = 14)
plot(hc)
dev.off()
pdf("09CCA/C-celltype-Species(no-glia).pdf", width = 8,height = 8)
circlize_dendrogram(hc,labels_track_height = NA,dend_track_height = 0.4)
dev.off()

#
cluster.averages <- AverageExpression(scRNA,group.by = "Species")
B.exp <- cluster.averages[["integrated"]]
hc <- as.dendrogram(hclust(dist(t(B.exp)), "ave"))

pdf("09CCA/H-Species(no-glia).pdf", height = 4, width = 6)
par(mar=c(9,3,6,2))
hc %>%
  set("labels_col", value = c('#eca8a9','#74aed4','#d3e2b7','#cfafd4','#f7c97e','#78a040','#E39A35',"#6778AE","#67adb7","#afacb7","#0d6e56","#97bfae",'#58A4C3', '#E4C755', '#B53E2B','#AA9A59'), k=4) %>%
  set("branches_k_color", value = c('#eca8a9','#74aed4','#d3e2b7','#cfafd4','#f7c97e','#78a040','#E39A35',"#6778AE","#67adb7","#afacb7","#0d6e56","#97bfae",'#58A4C3', '#E4C755', '#B53E2B','#AA9A59'), k = 4) %>%
  plot(axes=FALSE)
rect.dendrogram(hc, k=4, lty = 20,lwd = 3,x=c(2),prop_k_height = 0.1,col=rgb(0.1, 0.2, 0.4, 0.1)) 
dev.off()

#correlation analysis
scRNA$Species.celltype <- paste(scRNA$Species, Idents(scRNA), sep = "-")
table(scRNA$Species.celltype)  
av<-AverageExpression(scRNA,group.by = "Species.celltype", assays = "integrated")
av=av[[1]]
head(av)
cg=names(tail(sort(apply(av,1,sd)),1000))
View(av[cg,])
View(cor(av[cg,],method = "spearman"))
pdf("09CCA/cor_5Species-celltype(no-glia).pdf", width = 12,height = 11.5)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),treeheight_row = 12,treeheight_col = 12,fontsize = 10)
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"09CCA/cor_Species-celltype(no-glia).csv")

################################################################################
# Analysis of subtype across species
################################################################################
load("09CCA/scRNA(no-glia).RData")
scRNA_new2 <- scRNA
scRNA <- scRNA_new2[,scRNA_new2@meta.data[["celltype"]] %in% c("Cone","Rod")]
save(scRNA,file="10Species/PR.RData")
#Bipolar
scRNA <- scRNA_new2[,scRNA_new2@meta.data[["celltype"]] %in% c("Bipolar")]
save(scRNA,file="10Species/Bipolar.RData")
#Amacrine
scRNA <- scRNA_new2[,scRNA_new2@meta.data[["celltype"]] %in% c("Amacrine")]
save(scRNA,file="10Species/Amacrine.RData")
#RGC
scRNA <- scRNA_new2[,scRNA_new2@meta.data[["celltype"]] %in% c("RGC")]
save(scRNA,file="10Species/RGC.RData")
#Horizontal
scRNA <- scRNA_new2[,scRNA_new2@meta.data[["celltype"]] %in% c("Horizontal")]
save(scRNA,file="10Species/Horizontal.RData")
#Macroglia
scRNA <- scRNA_new2[,scRNA_new2@meta.data[["celltype"]] %in% c("Muller")]
save(scRNA,file="10Species/Muller.RData")

#The expression matrix was extracted for too-many-cells analysis
Idents(scRNA) <- "Species.celltype"
#ALL
scRNA_F=subset(scRNA, downsample=500)
count = as.data.frame(as.matrix(scRNA_F@assays$RNA@counts))
write.csv(count,"07too-many cells/allcells/count.csv")
label <- data.frame(item=colnames(count), label=scRNA_F$Species)
write.csv(label,"07too-many cells/allcells/labels.csv", row.names = F)
label <- data.frame(item=colnames(count), label=scRNA_F$Species.celltype)
write.csv(label,"07too-many cells/allcells/labels_celltype.csv", row.names = F)

#Cone
Idents(scRNA) <- "celltype"
Cone <- scRNA[,scRNA@meta.data[["celltype"]] %in% c("Cone")]
Idents(Cone) <- "Species"
#scRNA_F=Cone
scRNA_F=subset(Cone, downsample=4000)
count = as.data.frame(as.matrix(scRNA_F@assays$RNA@counts))
write.csv(count,"07too-many cells/Cone/count.csv")
label <- data.frame(item=colnames(count), label=scRNA_F$Species)
write.csv(label,"07too-many cells/Cone/labels.csv", row.names = F)

#Rod
Idents(scRNA) <- "celltype"
Rod <- scRNA[,scRNA@meta.data[["celltype"]] %in% c("Rod")]
Idents(Rod) <- "Species"
scRNA_F=Rod
count = as.data.frame(as.matrix(scRNA_F@assays$RNA@counts))
write.csv(count,"07too-many cells/Rod/count.csv")
label <- data.frame(item=colnames(count), label=scRNA_F$Species)
write.csv(label,"07too-many cells/Rod/labels.csv", row.names = F)

#Bipolar
Idents(scRNA) <- "celltype"
Bipolar <- scRNA[,scRNA@meta.data[["celltype"]] %in% c("Bipolar")]
Idents(Bipolar) <- "Species"
#scRNA_F=Bipolar
scRNA_F=subset(Bipolar, downsample=4000)
count = as.data.frame(as.matrix(scRNA_F@assays$RNA@counts))
write.csv(count,"07too-many cells/Bipolar/count.csv")
label <- data.frame(item=colnames(count), label=scRNA_F$Species)
write.csv(label,"07too-many cells/Bipolar/labels.csv", row.names = F)

#Horizontal
Idents(scRNA) <- "celltype"
Horizontal <- scRNA[,scRNA@meta.data[["celltype"]] %in% c("Horizontal")]
Idents(Horizontal) <- "Species"
scRNA_F=Horizontal
#scRNA_F=subset(Horizontal, downsample=1000)
count = as.data.frame(as.matrix(scRNA_F@assays$RNA@counts))
write.csv(count,"07too-many cells/Horizontal/count.csv")
label <- data.frame(item=colnames(count), label=scRNA_F$Species)
write.csv(label,"07too-many cells/Horizontal/labels.csv", row.names = F)

#Amacrine
Idents(scRNA) <- "celltype"
Amacrine <- scRNA[,scRNA@meta.data[["celltype"]] %in% c("Amacrine")]
Idents(Amacrine) <- "Species"
#scRNA_F=Amacrine
scRNA_F=subset(Amacrine, downsample=4000)
count = as.data.frame(as.matrix(scRNA_F@assays$RNA@counts))
write.csv(count,"07too-many cells/Amacrine/count.csv")
label <- data.frame(item=colnames(count), label=scRNA_F$Species)
write.csv(label,"07too-many cells/Amacrine/labels.csv", row.names = F)

#RGC
Idents(scRNA) <- "celltype"
RGC <- scRNA[,scRNA@meta.data[["celltype"]] %in% c("RGC")]
Idents(RGC) <- "Species"
#scRNA_F=RGC
scRNA_F=subset(RGC, downsample=4000)
count = as.data.frame(as.matrix(scRNA_F@assays$RNA@counts))
write.csv(count,"07too-many cells/RGC/count.csv")
label <- data.frame(item=colnames(count), label=scRNA_F$Species)
write.csv(label,"07too-many cells/RGC/labels.csv", row.names = F)

#Muller
Idents(scRNA) <- "celltype"
Muller <- scRNA[,scRNA@meta.data[["celltype"]] %in% c("Muller")]
Idents(Muller) <- "Species"
scRNA_F=Muller
#scRNA_F=subset(Muller, downsample=4000)
count = as.data.frame(as.matrix(scRNA_F@assays$RNA@counts))
write.csv(count,"07too-many cells/Muller/count.csv")
label <- data.frame(item=colnames(count), label=scRNA_F$Species)
write.csv(label,"07too-many cells/Muller/labels.csv", row.names = F)

################################################################################
# Identify differential expressed genes across conditions
################################################################################
# In the case of RGC, other cells are similar
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
RGC <- subset(scRNA, idents = "RGC")
Idents(RGC) <- "Species"
avg.RGC <- as.data.frame(log1p(AverageExpression(RGC, verbose = FALSE)$RNA))
avg.RGC$gene <- rownames(avg.RGC)

genes.to.label = c("SYT2", "KRT7", "ID2", "RXRG", "PDLIM5", "SEMA6A", "STMN1", "FSTL5", "ARPP21","RPH3A")
p1 <- ggplot(avg.RGC, aes(TS, Human)) + geom_point() + ggtitle("RGC")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)


