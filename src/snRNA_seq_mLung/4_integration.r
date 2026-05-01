suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(mclust))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(tidyverse))
options(warn=-1)
set.seed(1)


load("/data/amdqiao/5/intermediate/O1_filtered.RData")
load("/data/amdqiao/5/intermediate/O2_filtered.RData")
load("/data/amdqiao/5/intermediate/O3_filtered.RData")
load("/data/amdqiao/5/intermediate/Y1_filtered.RData")
load("/data/amdqiao/5/intermediate/Y2_filtered.RData")
load("/data/amdqiao/5/intermediate/Y3_filtered.RData")

mLung_O1 <- NormalizeData(mLung_O1, normalization.method = "LogNormalize")
mLung_O1 <- FindVariableFeatures(object = mLung_O1, selection.method="vst")
mLung_O2 <- NormalizeData(mLung_O2, normalization.method = "LogNormalize")
mLung_O2 <- FindVariableFeatures(object = mLung_O2, selection.method="vst")
mLung_O3 <- NormalizeData(mLung_O3, normalization.method = "LogNormalize")
mLung_O3 <- FindVariableFeatures(object = mLung_O3, selection.method="vst")
mLung_Y1 <- NormalizeData(mLung_Y1, normalization.method = "LogNormalize")
mLung_Y1 <- FindVariableFeatures(object = mLung_Y1, selection.method="vst")
mLung_Y2 <- NormalizeData(mLung_Y2, normalization.method = "LogNormalize")
mLung_Y2 <- FindVariableFeatures(object = mLung_Y2, selection.method="vst")
mLung_Y3 <- NormalizeData(mLung_Y3, normalization.method = "LogNormalize")
mLung_Y3 <- FindVariableFeatures(object = mLung_Y3, selection.method="vst")

features <- SelectIntegrationFeatures(object.list = c(mLung_O1, mLung_O2, mLung_O3, mLung_Y1, mLung_Y2, mLung_Y3), nfeatures = 3000)
mLung_anchors <- FindIntegrationAnchors(object.list = c(mLung_O1, mLung_O2, mLung_O3, mLung_Y1, mLung_Y2, mLung_Y3), anchor.features = features)
mLung_combined <- IntegrateData(anchorset = mLung_anchors)

rm(list = c("mLung_O1", "mLung_O2", "mLung_O3", "mLung_Y1", "mLung_Y2", "mLung_Y3"))

DefaultAssay(mLung_combined) <- "RNA"
mLung_combined[["percent.ribo"]] <- PercentageFeatureSet(mLung_combined, pattern = "Rp[ls]")
mLung_combined[["percent.mt"]] <- PercentageFeatureSet(mLung_combined, pattern = "^mt-")

DefaultAssay(mLung_combined) <- "integrated"
mLung_combined <- ScaleData(mLung_combined, do.center = TRUE, do.scale = TRUE, features = rownames(mLung_combined))
mLung_combined <- RunPCA(mLung_combined, ndims.print = 1:10, nfeatures.print = 5)

save(mLung_combined, file = "/data/amdqiao/5/intermediate/mLung_integrated.RData")

stdev <- mLung_combined@reductions$pca@stdev
var <- stdev^2
EndVar = 0
for(i in 1:length(var)){
  total <- sum(var)
  numerator <- sum(var[1:i])
  expvar <- numerator/total
  if(EndVar == 0){
    if(expvar > 0.9){
      EndVar <- EndVar + 1
      PCNum <- i
    }
  }
}
# Confirm PC's determined explain > 90% of variance
sum(var[1:PCNum])/ sum(var)

mLung_combined <- mLung_combined %>% 
  RunUMAP(reduction = "pca", dims = 1:PCNum) %>% 
  FindNeighbors(reduction = "pca", dims = 1:PCNum) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

DimPlot(mLung_combined, reduction = "umap", pt.size = 0.1)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_combined_umap_res_0.8.pdf", width = 8.5, height = 7, units = "in")
DimPlot(mLung_combined, reduction = "umap", group.by = "orig.ident", pt.size = 0.1)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_combined_umap_group_by_sample_res_0.8.pdf", width = 8.5, height = 7, units = "in")
DimPlot(mLung_combined, reduction = "umap", split.by = "orig.ident", ncol = 3, pt.size = 0.1)
ggsave("/home/amdqiao/5_sennet_raman/img/mLung_combined_umap_split_by_sample_res_0.8.pdf", width = 11, height = 7, units = "in")

save(mLung_combined, file = "/data/amdqiao/5/intermediate/mLung_clustered.RData")

markers <- FindAllMarkers(mLung_combined, min.pct = 0.25, only.pos = TRUE, assay = "RNA")
write.csv(markers, file = "/home/amdqiao/5_sennet_raman/results/mLung_markers_res_0.8.csv")
