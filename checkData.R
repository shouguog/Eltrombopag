rm(list = ls())
library(Seurat)

seuratObjList<-list()
seuratObj<-CreateSeuratObject(Read10X_h5("GSM8030597_DMSO_filtered_feature_bc_matrix.h5"), min.features = 100)
seuratObj@meta.data$libraryID<-"DMSO"
seuratObj@meta.data$sampleID<-"DMSO"
seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^MT-")
seuratObj <- subset(seuratObj, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
seuratObjList[["DMSO"]]<-seuratObj

seuratObj<-CreateSeuratObject(Read10X_h5("GSM8030598_TPO_filtered_feature_bc_matrix.h5"), min.features = 100)
seuratObj@meta.data$libraryID<-"TPO"
seuratObj@meta.data$sampleID<-"TPO"
seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^MT-")
seuratObj <- subset(seuratObj, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
seuratObjList[["TPO"]]<-seuratObj

seuratObj<-CreateSeuratObject(Read10X_h5("GSM8030599_AMG531_filtered_feature_bc_matrix.h5"), min.features = 100)
seuratObj@meta.data$libraryID<-"AMG531"
seuratObj@meta.data$sampleID<-"AMG531"
seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^MT-")
seuratObj <- subset(seuratObj, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
seuratObjList[["AMG531"]]<-seuratObj

seuratObj<-CreateSeuratObject(Read10X_h5("GSM8030600_Epag10_filtered_feature_bc_matrix.h5"), min.features = 100)
seuratObj@meta.data$libraryID<-"Epag10"
seuratObj@meta.data$sampleID<-"Epag10"
seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^MT-")
seuratObj <- subset(seuratObj, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
seuratObjList[["Epag10"]]<-seuratObj

seuratObj<-CreateSeuratObject(Read10X_h5("GSM8030601_Epag1_filtered_feature_bc_matrix.h5"), min.features = 100)
seuratObj@meta.data$libraryID<-"Epag1"
seuratObj@meta.data$sampleID<-"Epag1"
seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^MT-")
seuratObj <- subset(seuratObj, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
seuratObjList[["Epag1"]]<-seuratObj

save(seuratObjList, file = "seuratObjList_step0.RData")
