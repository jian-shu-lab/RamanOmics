suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(mclust))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(tidyverse))
options(warn=-1)
set.seed(1)


### Old 1 ###
load("/data/amdqiao/5/intermediate/O1_DoubletFinder.RData")

counts_per_cell <- Matrix::colSums(mLung_O1)
counts_per_gene <- Matrix::rowSums(mLung_O1)
mean(counts_per_cell)
mean(counts_per_gene)
mean(mLung_O1@meta.data$nFeature_RNA)

mLung_O1[["percent.mt"]] <- PercentageFeatureSet(mLung_O1, pattern = "^mt-")
VlnPlot(mLung_O1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_O1_QC_doublet_removed.pdf", width = 9, height = 6, units = "in")

feat.max <- round(mean(mLung_O1$nFeature_RNA) + 2 * sd(mLung_O1$nFeature_RNA), digits = -2)
mLung_O1 <- subset(mLung_O1, subset = nFeature_RNA > 300 & nFeature_RNA < feat.max & nCount_RNA > 300 & percent.mt < 10)
VlnPlot(mLung_O1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_O1_QC_filtered.pdf", width = 9, height = 6, units = "in")
dim(mLung_O1)

save(mLung_O1, file = "/data/amdqiao/5/intermediate/O1_filtered.RData")


### Old 2 ###
load("/data/amdqiao/5/intermediate/O2_DoubletFinder.RData")

counts_per_cell <- Matrix::colSums(mLung_O2)
counts_per_gene <- Matrix::rowSums(mLung_O2)
mean(counts_per_cell)
mean(counts_per_gene)
mean(mLung_O2@meta.data$nFeature_RNA)

mLung_O2[["percent.mt"]] <- PercentageFeatureSet(mLung_O2, pattern = "^mt-")
VlnPlot(mLung_O2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_O2_QC_doublet_removed.pdf", width = 9, height = 6, units = "in")

feat.max <- round(mean(mLung_O2$nFeature_RNA) + 2 * sd(mLung_O2$nFeature_RNA), digits = -2)
mLung_O2 <- subset(mLung_O2, subset = nFeature_RNA > 300 & nFeature_RNA < feat.max & nCount_RNA > 300 & percent.mt < 10)
VlnPlot(mLung_O2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_O2_QC_filtered.pdf", width = 9, height = 6, units = "in")
dim(mLung_O2)

save(mLung_O2, file = "/data/amdqiao/5/intermediate/O2_filtered.RData")


### Old 3 ###
load("/data/amdqiao/5/intermediate/O3_DoubletFinder.RData")

counts_per_cell <- Matrix::colSums(mLung_O3)
counts_per_gene <- Matrix::rowSums(mLung_O3)
mean(counts_per_cell)
mean(counts_per_gene)
mean(mLung_O3@meta.data$nFeature_RNA)

mLung_O3[["percent.mt"]] <- PercentageFeatureSet(mLung_O3, pattern = "^mt-")
VlnPlot(mLung_O3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_O3_QC_doublet_removed.pdf", width = 9, height = 6, units = "in")

feat.max <- round(mean(mLung_O3$nFeature_RNA) + 2 * sd(mLung_O3$nFeature_RNA), digits = -2)
mLung_O3 <- subset(mLung_O3, subset = nFeature_RNA > 300 & nFeature_RNA < feat.max & nCount_RNA > 300 & percent.mt < 10)
VlnPlot(mLung_O3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_O3_QC_filtered.pdf", width = 9, height = 6, units = "in")
dim(mLung_O3)

save(mLung_O3, file = "/data/amdqiao/5/intermediate/O3_filtered.RData")


### Young 1 ###
load("/data/amdqiao/5/intermediate/Y1_DoubletFinder.RData")

counts_per_cell <- Matrix::colSums(mLung_Y1)
counts_per_gene <- Matrix::rowSums(mLung_Y1)
mean(counts_per_cell)
mean(counts_per_gene)
mean(mLung_Y1@meta.data$nFeature_RNA)

mLung_Y1[["percent.mt"]] <- PercentageFeatureSet(mLung_Y1, pattern = "^mt-")
VlnPlot(mLung_Y1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_Y1_QC_doublet_removed.pdf", width = 9, height = 6, units = "in")

feat.max <- round(mean(mLung_Y1$nFeature_RNA) + 2 * sd(mLung_Y1$nFeature_RNA), digits = -2)
mLung_Y1 <- subset(mLung_Y1, subset = nFeature_RNA > 300 & nFeature_RNA < feat.max & nCount_RNA > 300 & percent.mt < 10)
VlnPlot(mLung_Y1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_Y1_QC_filtered.pdf", width = 9, height = 6, units = "in")
dim(mLung_Y1)

save(mLung_Y1, file = "/data/amdqiao/5/intermediate/Y1_filtered.RData")


### Young 2 ###
load("/data/amdqiao/5/intermediate/Y2_DoubletFinder.RData")

counts_per_cell <- Matrix::colSums(mLung_Y2)
counts_per_gene <- Matrix::rowSums(mLung_Y2)
mean(counts_per_cell)
mean(counts_per_gene)
mean(mLung_Y2@meta.data$nFeature_RNA)

mLung_Y2[["percent.mt"]] <- PercentageFeatureSet(mLung_Y2, pattern = "^mt-")
VlnPlot(mLung_Y2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_Y2_QC_doublet_removed.pdf", width = 9, height = 6, units = "in")

feat.max <- round(mean(mLung_Y2$nFeature_RNA) + 2 * sd(mLung_Y2$nFeature_RNA), digits = -2)
mLung_Y2 <- subset(mLung_Y2, subset = nFeature_RNA > 300 & nFeature_RNA < feat.max & nCount_RNA > 300 & percent.mt < 10)
VlnPlot(mLung_Y2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_Y2_QC_filtered.pdf", width = 9, height = 6, units = "in")
dim(mLung_Y2)

save(mLung_Y2, file = "/data/amdqiao/5/intermediate/Y2_filtered.RData")


### Young 3 ###
load("/data/amdqiao/5/intermediate/Y3_DoubletFinder.RData")

counts_per_cell <- Matrix::colSums(mLung_Y3)
counts_per_gene <- Matrix::rowSums(mLung_Y3)
mean(counts_per_cell)
mean(counts_per_gene)
mean(mLung_Y3@meta.data$nFeature_RNA)

mLung_Y3[["percent.mt"]] <- PercentageFeatureSet(mLung_Y3, pattern = "^mt-")
VlnPlot(mLung_Y3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_Y3_QC_doublet_removed.pdf", width = 9, height = 6, units = "in")

feat.max <- round(mean(mLung_Y3$nFeature_RNA) + 2 * sd(mLung_Y3$nFeature_RNA), digits = -2)
mLung_Y3 <- subset(mLung_Y3, subset = nFeature_RNA > 300 & nFeature_RNA < feat.max & nCount_RNA > 300 & percent.mt < 10)
VlnPlot(mLung_Y3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_Y3_QC_filtered.pdf", width = 9, height = 6, units = "in")
dim(mLung_Y3)

save(mLung_Y3, file = "/data/amdqiao/5/intermediate/Y3_filtered.RData")
