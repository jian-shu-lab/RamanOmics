library(Seurat)
library(tidyverse)
set.seed(1234)

############mLung_O1#############
O1_new_cell_gene_matrix <- read.csv("~/Desktop/O1_new_cell_gene_matrix.csv", row.names=1)
sample.raw_mLung <- t(O1_new_cell_gene_matrix)
mLung_O1 <- CreateSeuratObject(counts = sample.raw_mLung, project = "mLung_O1_ST", min.cells = 1)
mLung_O1@assays[["RNA"]]@counts@Dimnames[[1]] <- gene_O1$x
mLung_O1@assays[["RNA"]]@data@Dimnames[[1]] <- gene_O1$x
mLung_O1@meta.data[["orig.ident"]] <- "mLung_O1_ST"
dim(mLung_O1)
counts_per_cell <- Matrix::colSums(mLung_O1)
counts_per_gene <- Matrix::rowSums(mLung_O1)
mean(counts_per_cell)
mean(counts_per_gene)
mean(mLung_O1@meta.data$nFeature_RNA)
mean(mLung_O1@meta.data$nCount_RNA)
VlnPlot(mLung_O1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

mLung_O1 <- subset(mLung_O1,nCount_RNA >= 5 & nFeature_RNA >= 5)

counts_matrix <- mLung_O1@assays[["RNA"]]@counts

# Create a new Seurat object
mLung_O1 <- CreateSeuratObject(counts = counts_matrix)
mLung_O1@meta.data[["orig.ident"]] <- "mLung_O1_ST"

dim(mLung_O1)
counts_per_cell <- Matrix::colSums(mLung_O1)
counts_per_gene <- Matrix::rowSums(mLung_O1)
mean(counts_per_cell)
mean(counts_per_gene)
mean(mLung_O1@meta.data$nFeature_RNA)
mean(mLung_O1@meta.data$nCount_RNA)

mLung_O1 <- SCTransform(mLung_O1, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:50, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:50, verbose = FALSE) %>%
  FindClusters(resolution = 0.3, verbose = FALSE)

DimPlot(mLung_O1, reduction = "umap")
############mLung_O1#############
############mLung_O2#############
O2_new_cell_gene_matrix <- read.csv("~/Desktop/O2_new_cell_gene_matrix.csv", row.names=1)
sample.raw_mLung_O2 <- t(O2_new_cell_gene_matrix)
mLung_O2 <- CreateSeuratObject(counts = sample.raw_mLung_O2, project = "mLung_O2_ST", min.cells = 1)
mLung_O2@assays[["RNA"]]@counts@Dimnames[[1]] <- gene_O2$x
mLung_O2@assays[["RNA"]]@data@Dimnames[[1]] <- gene_O2$x
mLung_O2@meta.data[["orig.ident"]] <- "mLung_O2_ST"
dim(mLung_O2)
counts_per_cell <- Matrix::colSums(mLung_O2)
counts_per_gene <- Matrix::rowSums(mLung_O2)
mean(counts_per_cell)
mean(counts_per_gene)
mean(mLung_O2@meta.data$nFeature_RNA)
mean(mLung_O2@meta.data$nCount_RNA)
VlnPlot(mLung_O2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

mLung_O2 <- subset(mLung_O2,nCount_RNA >= 5 & nFeature_RNA >= 5)

dim(mLung_O2)
counts_per_cell <- Matrix::colSums(mLung_O2)
counts_per_gene <- Matrix::rowSums(mLung_O2)
mean(counts_per_cell)
mean(counts_per_gene)
mean(mLung_O2@meta.data$nFeature_RNA)
mean(mLung_O2@meta.data$nCount_RNA)

counts_matrix <- mLung_O2@assays[["RNA"]]@counts

# Create a new Seurat object
mLung_O2 <- CreateSeuratObject(counts = counts_matrix)
mLung_O2@meta.data[["orig.ident"]] <- "mLung_O2_ST"

mLung_O2 <- SCTransform(mLung_O2, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)
############mLung_O3#############
O3_new_cell_gene_matrix <- read.csv("~/Desktop/O3_new_cell_gene_matrix.csv", row.names=1)
sample.raw_mLung <- t(O3_new_cell_gene_matrix)
mLung_O3 <- CreateSeuratObject(counts = sample.raw_mLung, project = "mLung_O3_ST", min.cells = 1)
mLung_O3@assays[["RNA"]]@counts@Dimnames[[1]] <- gene_O3$x
mLung_O3@assays[["RNA"]]@data@Dimnames[[1]] <- gene_O3$x
mLung_O3@meta.data[["orig.ident"]] <- "mLung_O3_ST"
dim(mLung_O3)
counts_per_cell <- Matrix::colSums(mLung_O3)
counts_per_gene <- Matrix::rowSums(mLung_O3)
mean(counts_per_cell)
mean(counts_per_gene)
mean(mLung_O3@meta.data$nFeature_RNA)
mean(mLung_O3@meta.data$nCount_RNA)
VlnPlot(mLung_O3, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

mLung_O3 <- subset(mLung_O3,nCount_RNA >= 5 & nFeature_RNA >= 5)

dim(mLung_O3)
counts_per_cell <- Matrix::colSums(mLung_O3)
counts_per_gene <- Matrix::rowSums(mLung_O3)
mean(counts_per_cell)
mean(counts_per_gene)
mean(mLung_O3@meta.data$nFeature_RNA)
mean(mLung_O3@meta.data$nCount_RNA)

counts_matrix <- mLung_O3@assays[["RNA"]]@counts

# Create a new Seurat object
mLung_O3 <- CreateSeuratObject(counts = counts_matrix)
mLung_O3@meta.data[["orig.ident"]] <- "mLung_O3_ST"

mLung_O3 <- SCTransform(mLung_O3, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)

DimPlot(mLung_O3, reduction = "umap")
############mLung_O3#############
#################Integration#############
mLung_O.list <- c(mLung_O1,mLung_O2, mLung_O3)
features <- SelectIntegrationFeatures(object.list = mLung_O.list, nfeatures = 890)
mLung_O.list <- PrepSCTIntegration(object.list = mLung_O.list, anchor.features = features)
mLung_O.anchors <- FindIntegrationAnchors(object.list = mLung_O.list, normalization.method = "SCT",
                                          anchor.features = features)
mLung_O.combined.sct <- IntegrateData(anchorset = mLung_O.anchors, normalization.method = "SCT")

mLung_O.combined.sct <- RunPCA(mLung_O.combined.sct, verbose = FALSE)
ElbowPlot(mLung_O.combined.sct)
mLung_O.combined.sct <- RunUMAP(mLung_O.combined.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
mLung_O.combined.sct <- FindNeighbors(mLung_O.combined.sct, reduction = "pca", dims = 1:30)
mLung_O.combined.sct <- FindClusters(mLung_O.combined.sct, resolution = 0.1)

DimPlot(mLung_O.combined.sct, reduction = "umap", group.by = "orig.ident")
DimPlot(mLung_O.combined.sct, reduction = "umap")
save(mLung_O.combined.sct, file='~/Desktop/mLung_O.combined.sct.RData')
