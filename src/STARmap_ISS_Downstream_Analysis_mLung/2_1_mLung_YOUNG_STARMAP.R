library(Seurat)
library(tidyverse)
set.seed(1234)

############mLung_Y1#############
Y1_new_cell_gene_matrix <- read.csv("~/Desktop/Y1_new_cell_gene_matrix.csv", row.names=1)
sample.raw_mLung <- t(Y1_new_cell_gene_matrix)
mLung_Y1 <- CreateSeuratObject(counts = sample.raw_mLung, project = "mLung_Y1_ST", min.cells = 1)
mLung_Y1@assays[["RNA"]]@counts@Dimnames[[1]] <- gene_Y1$x
mLung_Y1@assays[["RNA"]]@data@Dimnames[[1]] <- gene_Y1$x
mLung_Y1@meta.data[["orig.ident"]] <- "mLung_Y1_ST"
dim(mLung_Y1)
counts_per_cell <- Matrix::colSums(mLung_Y1)
counts_per_gene <- Matrix::rowSums(mLung_Y1)
mean(counts_per_cell)
mean(counts_per_gene)
mean(mLung_Y1@meta.data$nFeature_RNA)
mean(mLung_Y1@meta.data$nCount_RNA)
VlnPlot(mLung_Y1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

mLung_Y1 <- subset(mLung_Y1,nCount_RNA >= 5 & nFeature_RNA >= 5)

counts_matrix <- mLung_Y1@assays[["RNA"]]@counts

# Create a new Seurat object
mLung_Y1 <- CreateSeuratObject(counts = counts_matrix)
mLung_Y1@meta.data[["orig.ident"]] <- "mLung_Y1_ST"

dim(mLung_Y1)
counts_per_cell <- Matrix::colSums(mLung_Y1)
counts_per_gene <- Matrix::rowSums(mLung_Y1)
mean(counts_per_cell)
mean(counts_per_gene)
mean(mLung_Y1@meta.data$nFeature_RNA)
mean(mLung_Y1@meta.data$nCount_RNA)

mLung_Y1 <- SCTransform(mLung_Y1, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:50, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:50, verbose = FALSE) %>%
  FindClusters(resolution = 0.3, verbose = FALSE)

DimPlot(mLung_Y1, reduction = "umap")
############mLung_Y1#############
############mLung_Y2#############
Y2_new_cell_gene_matrix <- read.csv("~/Desktop/Y2_new_cell_gene_matrix.csv", row.names=1)
sample.raw_mLung_Y2 <- t(Y2_new_cell_gene_matrix)
mLung_Y2 <- CreateSeuratObject(counts = sample.raw_mLung_Y2, project = "mLung_Y2_ST", min.cells = 1)
mLung_Y2@assays[["RNA"]]@counts@Dimnames[[1]] <- gene_Y2$x
mLung_Y2@assays[["RNA"]]@data@Dimnames[[1]] <- gene_Y2$x
mLung_Y2@meta.data[["orig.ident"]] <- "mLung_Y2_ST"
dim(mLung_Y2)
counts_per_cell <- Matrix::colSums(mLung_Y2)
counts_per_gene <- Matrix::rowSums(mLung_Y2)
mean(counts_per_cell)
mean(counts_per_gene)
mean(mLung_Y2@meta.data$nFeature_RNA)
mean(mLung_Y2@meta.data$nCount_RNA)
VlnPlot(mLung_Y2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

mLung_Y2 <- subset(mLung_Y2,nCount_RNA >= 5 & nFeature_RNA >= 5)

dim(mLung_Y2)
counts_per_cell <- Matrix::colSums(mLung_Y2)
counts_per_gene <- Matrix::rowSums(mLung_Y2)
mean(counts_per_cell)
mean(counts_per_gene)
mean(mLung_Y2@meta.data$nFeature_RNA)
mean(mLung_Y2@meta.data$nCount_RNA)

counts_matrix <- mLung_Y2@assays[["RNA"]]@counts

# Create a new Seurat object
mLung_Y2 <- CreateSeuratObject(counts = counts_matrix)
mLung_Y2@meta.data[["orig.ident"]] <- "mLung_Y2_ST"

mLung_Y2 <- SCTransform(mLung_Y2, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)
############mLung_Y3#############
Y3_new_cell_gene_matrix <- read.csv("~/Desktop/Y3_new_cell_gene_matrix.csv", row.names=1)
sample.raw_mLung <- t(Y3_new_cell_gene_matrix)
mLung_Y3 <- CreateSeuratObject(counts = sample.raw_mLung, project = "mLung_Y3_ST", min.cells = 1)
mLung_Y3@assays[["RNA"]]@counts@Dimnames[[1]] <- gene_Y3$x
mLung_Y3@assays[["RNA"]]@data@Dimnames[[1]] <- gene_Y3$x
mLung_Y3@meta.data[["orig.ident"]] <- "mLung_Y3_ST"
dim(mLung_Y3)
counts_per_cell <- Matrix::colSums(mLung_Y3)
counts_per_gene <- Matrix::rowSums(mLung_Y3)
mean(counts_per_cell)
mean(counts_per_gene)
mean(mLung_Y3@meta.data$nFeature_RNA)
mean(mLung_Y3@meta.data$nCount_RNA)
VlnPlot(mLung_Y3, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

mLung_Y3 <- subset(mLung_Y3,nCount_RNA >= 5 & nFeature_RNA >= 5)

dim(mLung_Y3)
counts_per_cell <- Matrix::colSums(mLung_Y3)
counts_per_gene <- Matrix::rowSums(mLung_Y3)
mean(counts_per_cell)
mean(counts_per_gene)
mean(mLung_Y3@meta.data$nFeature_RNA)
mean(mLung_Y3@meta.data$nCount_RNA)

counts_matrix <- mLung_Y3@assays[["RNA"]]@counts

# Create a new Seurat object
mLung_Y3 <- CreateSeuratObject(counts = counts_matrix)
mLung_Y3@meta.data[["orig.ident"]] <- "mLung_Y3_ST"

mLung_Y3 <- SCTransform(mLung_Y3, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)

DimPlot(mLung_Y3, reduction = "umap")
############mLung_Y3#############
#################Integration#############
mLung_Y.list <- c(mLung_Y1,mLung_Y2, mLung_Y3)
features <- SelectIntegrationFeatures(object.list = mLung_Y.list, nfeatures = 890)
mLung_Y.list <- PrepSCTIntegration(object.list = mLung_Y.list, anchor.features = features)
mLung_Y.anchors <- FindIntegrationAnchors(object.list = mLung_Y.list, normalization.method = "SCT",
                                          anchor.features = features)
mLung_Y.combined.sct <- IntegrateData(anchorset = mLung_Y.anchors, normalization.method = "SCT")

mLung_Y.combined.sct <- RunPCA(mLung_Y.combined.sct, verbose = FALSE)
ElbowPlot(mLung_Y.combined.sct)
mLung_Y.combined.sct <- RunUMAP(mLung_Y.combined.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
mLung_Y.combined.sct <- FindNeighbors(mLung_Y.combined.sct, reduction = "pca", dims = 1:30)
mLung_Y.combined.sct <- FindClusters(mLung_Y.combined.sct, resolution = 0.1)

DimPlot(mLung_Y.combined.sct, reduction = "umap", group.by = "orig.ident")
DimPlot(mLung_Y.combined.sct, reduction = "umap")
save(mLung_Y.combined.sct, file='~/Desktop/mLung_Y.combined.sct.RData')

