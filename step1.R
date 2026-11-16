#R4.0.4 on PC
##R version 4.1.3 (2022-03-10) -- "One Push-Up"
rm(list = ls())
library(Seurat)
load("seuratObjList_step0.RData")

cellNum<-lapply(X = seuratObjList, FUN = function(x) {
  x <- dim(x@meta.data)[1]
})


seuratObjList <- lapply(X = seuratObjList, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seuratObjList)

seuratObjList <- lapply(X = seuratObjList, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE, npcs = 50)
})
####Use the HD as reference 1
CD34.anchors <- FindIntegrationAnchors(object.list = seuratObjList, reference = 1, reduction = "cca",dims = 1:50)
# this command creates an 'integrated' data assay
CD34.combined <- IntegrateData(anchorset = CD34.anchors, dims = 1:50)
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(CD34.combined) <- "integrated"
save(CD34.combined, file = "seuratObjList_CD34_step1.RData")
CD34.combined <- ScaleData(CD34.combined, verbose = FALSE)
CD34.combined <- RunPCA(CD34.combined, verbose = FALSE)
CD34.combined <- RunUMAP(CD34.combined, dims = 1:50)
CD34.combined <- FindNeighbors(CD34.combined, reduction = "pca", dims = 1:50)
CD34.combined <- FindClusters(CD34.combined, resolution = c(0.5, 0.8, 1, 1.5,2))
CD34.combined <- RunTSNE(CD34.combined, dims = 1:50)
save(CD34.combined, file = "seuratObjList_CD34_step2.RData")
metaData<-CD34.combined@meta.data
save(metaData, file = "metaData.RData")


