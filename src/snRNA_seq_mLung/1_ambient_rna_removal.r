suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(mclust))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(SoupX))
suppressPackageStartupMessages(library(scDblFinder))
options(warn=-1)
set.seed(1)


### Old 1 ###
data <- Read10X_h5("/data/amdqiao/5/data/O1_26mo_mLung_scGEX_filtered_feature_bc_matrix.h5")
mLung_O1 = CreateSeuratObject(counts = data, project = "mLung_O1", min.cells = 5)
dim(mLung_O1)
rm(data)

VlnPlot(object = mLung_O1, features = c("nCount_RNA", "nFeature_RNA"), ncol = 2, pt.size = 0)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_O1_QC_raw.pdf", width = 6, height = 6, units = "in")

mLung_O1 <- NormalizeData(mLung_O1, normalization.method = "LogNormalize")
mLung_O1 <- FindVariableFeatures(object = mLung_O1)
mLung_O1 <- ScaleData(mLung_O1)
mLung_O1 <- RunPCA(mLung_O1, features = VariableFeatures(object = mLung_O1), npcs = 100)
mLung_O1 <- mLung_O1 %>% 
  RunUMAP(reduction = "pca", dims = 1:100) %>% 
  FindNeighbors(reduction = "pca", dims = 1:100) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

toc <- Read10X_h5("/data/amdqiao/5/data/O1_26mo_mLung_scGEX_filtered_feature_bc_matrix.h5")
tod <- Read10X_h5("/data/amdqiao/5/data/O1_26mo_mLung_scGEX_raw_feature_bc_matrix.h5")

sc = SoupChannel(tod, toc)
sc = setClusters(sc, mLung_O1@meta.data[["seurat_clusters"]])
sc = autoEstCont(sc)
out = adjustCounts(sc, roundToInt=TRUE)

sce <- scDblFinder(out)

mLung_O1 = CreateSeuratObject(counts = out, project = "mLung_O1", min.cells = 5)
mLung_O1@meta.data$scDblFinder <- NULL
mLung_O1@meta.data$scDblFinder <- sce$scDblFinder.class

dim(mLung_O1)
save(mLung_O1, file = "/data/amdqiao/5/intermediate/O1_scDblFinder.RData")


### Old 2 ###
data <- Read10X_h5("/data/amdqiao/5/data/O2_26mo_mLung_scGEX_filtered_feature_bc_matrix.h5")
mLung_O2 = CreateSeuratObject(counts = data, project = "mLung_O2", min.cells = 5)
dim(mLung_O2)
rm(data)

VlnPlot(object = mLung_O2, features = c("nCount_RNA", "nFeature_RNA"), ncol = 2, pt.size = 0)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_O2_QC_raw.pdf", width = 6, height = 6, units = "in")

mLung_O2 <- NormalizeData(mLung_O2, normalization.method = "LogNormalize")
mLung_O2 <- FindVariableFeatures(object = mLung_O2)
mLung_O2 <- ScaleData(mLung_O2)
mLung_O2 <- RunPCA(mLung_O2, features = VariableFeatures(object = mLung_O2), npcs = 100)
mLung_O2 <- mLung_O2 %>% 
  RunUMAP(reduction = "pca", dims = 1:100) %>% 
  FindNeighbors(reduction = "pca", dims = 1:100) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

toc <- Read10X_h5("/data/amdqiao/5/data/O2_26mo_mLung_scGEX_filtered_feature_bc_matrix.h5")
tod <- Read10X_h5("/data/amdqiao/5/data/O2_26mo_mLung_scGEX_raw_feature_bc_matrix.h5")

sc = SoupChannel(tod, toc)
sc = setClusters(sc, mLung_O2@meta.data[["seurat_clusters"]])
sc = autoEstCont(sc)
out = adjustCounts(sc, roundToInt=TRUE)

sce <- scDblFinder(out)

mLung_O2 = CreateSeuratObject(counts = out, project = "mLung_O2", min.cells = 5)
mLung_O2@meta.data$scDblFinder <- NULL
mLung_O2@meta.data$scDblFinder <- sce$scDblFinder.class

dim(mLung_O2)
save(mLung_O2, file = "/data/amdqiao/5/intermediate/O2_scDblFinder.RData")


### Old 3 ###
data <- Read10X_h5("/data/amdqiao/5/data/O3_26mo_mLung_scGEX_filtered_feature_bc_matrix.h5")
mLung_O3 = CreateSeuratObject(counts = data, project = "mLung_O3", min.cells = 5)
dim(mLung_O3)
rm(data)

VlnPlot(object = mLung_O3, features = c("nCount_RNA", "nFeature_RNA"), ncol = 2, pt.size = 0)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_O3_QC_raw.pdf", width = 6, height = 6, units = "in")

mLung_O3 <- NormalizeData(mLung_O3, normalization.method = "LogNormalize")
mLung_O3 <- FindVariableFeatures(object = mLung_O3)
mLung_O3 <- ScaleData(mLung_O3)
mLung_O3 <- RunPCA(mLung_O3, features = VariableFeatures(object = mLung_O3), npcs = 100)
mLung_O3 <- mLung_O3 %>% 
  RunUMAP(reduction = "pca", dims = 1:100) %>% 
  FindNeighbors(reduction = "pca", dims = 1:100) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

toc <- Read10X_h5("/data/amdqiao/5/data/O3_26mo_mLung_scGEX_filtered_feature_bc_matrix.h5")
tod <- Read10X_h5("/data/amdqiao/5/data/O3_26mo_mLung_scGEX_raw_feature_bc_matrix.h5")

sc = SoupChannel(tod, toc)
sc = setClusters(sc, mLung_O3@meta.data[["seurat_clusters"]])
sc = autoEstCont(sc)
out = adjustCounts(sc, roundToInt=TRUE)

sce <- scDblFinder(out)

mLung_O3 = CreateSeuratObject(counts = out, project = "mLung_O3", min.cells = 5)
mLung_O3@meta.data$scDblFinder <- NULL
mLung_O3@meta.data$scDblFinder <- sce$scDblFinder.class

dim(mLung_O3)
save(mLung_O3, file = "/data/amdqiao/5/intermediate/O3_scDblFinder.RData")


### Young 1 ###
data <- Read10X_h5("/data/amdqiao/5/data/Y1_2mo_mLung_scGEX_filtered_feature_bc_matrix.h5")
mLung_Y1 = CreateSeuratObject(counts = data, project = "mLung_Y1", min.cells = 5)
dim(mLung_Y1)
rm(data)

VlnPlot(object = mLung_Y1, features = c("nCount_RNA", "nFeature_RNA"), ncol = 2, pt.size = 0)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_Y1_QC_raw.pdf", width = 6, height = 6, units = "in")

mLung_Y1 <- NormalizeData(mLung_Y1, normalization.method = "LogNormalize")
mLung_Y1 <- FindVariableFeatures(object = mLung_Y1)
mLung_Y1 <- ScaleData(mLung_Y1)
mLung_Y1 <- RunPCA(mLung_Y1, features = VariableFeatures(object = mLung_Y1), npcs = 100)
mLung_Y1 <- mLung_Y1 %>% 
  RunUMAP(reduction = "pca", dims = 1:100) %>% 
  FindNeighbors(reduction = "pca", dims = 1:100) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

toc <- Read10X_h5("/data/amdqiao/5/data/Y1_2mo_mLung_scGEX_filtered_feature_bc_matrix.h5")
tod <- Read10X_h5("/data/amdqiao/5/data/Y1_2mo_mLung_scGEX_raw_feature_bc_matrix.h5")

sc = SoupChannel(tod, toc)
sc = setClusters(sc, mLung_Y1@meta.data[["seurat_clusters"]])
sc = autoEstCont(sc)
out = adjustCounts(sc, roundToInt=TRUE)

sce <- scDblFinder(out)

mLung_Y1 = CreateSeuratObject(counts = out, project = "mLung_Y1", min.cells = 5)
mLung_Y1@meta.data$scDblFinder <- NULL
mLung_Y1@meta.data$scDblFinder <- sce$scDblFinder.class

dim(mLung_Y1)
save(mLung_Y1, file = "/data/amdqiao/5/intermediate/Y1_scDblFinder.RData")


### Young 2 ###
data <- Read10X_h5("/data/amdqiao/5/data/Y2_2mo_mLung_scGEX_filtered_feature_bc_matrix.h5")
mLung_Y2 = CreateSeuratObject(counts = data, project = "mLung_Y2", min.cells = 5)
dim(mLung_Y2)
rm(data)

VlnPlot(object = mLung_Y2, features = c("nCount_RNA", "nFeature_RNA"), ncol = 2, pt.size = 0)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_Y2_QC_raw.pdf", width = 6, height = 6, units = "in")

mLung_Y2 <- NormalizeData(mLung_Y2, normalization.method = "LogNormalize")
mLung_Y2 <- FindVariableFeatures(object = mLung_Y2)
mLung_Y2 <- ScaleData(mLung_Y2)
mLung_Y2 <- RunPCA(mLung_Y2, features = VariableFeatures(object = mLung_Y2), npcs = 100)
mLung_Y2 <- mLung_Y2 %>% 
  RunUMAP(reduction = "pca", dims = 1:100) %>% 
  FindNeighbors(reduction = "pca", dims = 1:100) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

toc <- Read10X_h5("/data/amdqiao/5/data/Y2_2mo_mLung_scGEX_filtered_feature_bc_matrix.h5")
tod <- Read10X_h5("/data/amdqiao/5/data/Y2_2mo_mLung_scGEX_raw_feature_bc_matrix.h5")

sc = SoupChannel(tod, toc)
sc = setClusters(sc, mLung_Y2@meta.data[["seurat_clusters"]])
sc = autoEstCont(sc)
out = adjustCounts(sc, roundToInt=TRUE)

sce <- scDblFinder(out)

mLung_Y2 = CreateSeuratObject(counts = out, project = "mLung_Y2", min.cells = 5)
mLung_Y2@meta.data$scDblFinder <- NULL
mLung_Y2@meta.data$scDblFinder <- sce$scDblFinder.class

dim(mLung_Y2)
save(mLung_Y2, file = "/data/amdqiao/5/intermediate/Y2_scDblFinder.RData")


### Young 3 ###
data <- Read10X_h5("/data/amdqiao/5/data/Y3_2mo_mLung_scGEX_filtered_feature_bc_matrix.h5")
mLung_Y3 = CreateSeuratObject(counts = data, project = "mLung_Y3", min.cells = 5)
dim(mLung_Y3)
rm(data)

VlnPlot(object = mLung_Y3, features = c("nCount_RNA", "nFeature_RNA"), ncol = 2, pt.size = 0)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_Y3_QC_raw.pdf", width = 6, height = 6, units = "in")

mLung_Y3 <- NormalizeData(mLung_Y3, normalization.method = "LogNormalize")
mLung_Y3 <- FindVariableFeatures(object = mLung_Y3)
mLung_Y3 <- ScaleData(mLung_Y3)
mLung_Y3 <- RunPCA(mLung_Y3, features = VariableFeatures(object = mLung_Y3), npcs = 100)
mLung_Y3 <- mLung_Y3 %>% 
  RunUMAP(reduction = "pca", dims = 1:100) %>% 
  FindNeighbors(reduction = "pca", dims = 1:100) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

toc <- Read10X_h5("/data/amdqiao/5/data/Y3_2mo_mLung_scGEX_filtered_feature_bc_matrix.h5")
tod <- Read10X_h5("/data/amdqiao/5/data/Y3_2mo_mLung_scGEX_raw_feature_bc_matrix.h5")

sc = SoupChannel(tod, toc)
sc = setClusters(sc, mLung_Y3@meta.data[["seurat_clusters"]])
sc = autoEstCont(sc)
out = adjustCounts(sc, roundToInt=TRUE)

sce <- scDblFinder(out)

mLung_Y3 = CreateSeuratObject(counts = out, project = "mLung_Y3", min.cells = 5)
mLung_Y3@meta.data$scDblFinder <- NULL
mLung_Y3@meta.data$scDblFinder <- sce$scDblFinder.class

dim(mLung_Y3)
save(mLung_Y3, file = "/data/amdqiao/5/intermediate/Y3_scDblFinder.RData")
