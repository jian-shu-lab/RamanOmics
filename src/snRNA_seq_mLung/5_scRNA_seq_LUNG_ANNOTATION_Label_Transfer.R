library(Seurat)
library(tidyverse)
set.seed(1234)

load("/Users/js1009/Desktop/GSE124872_raw_counts_single_cell.RData")
mLungGSE124872 <- CreateSeuratObject(counts = raw_counts, project = "mLung", min.cells = 1)

mLungGSE124872 <- NormalizeData(mLungGSE124872) 
mLungGSE124872 <- FindVariableFeatures(mLungGSE124872, nfeatures = 3000)
mLungGSE124872 <- ScaleData(mLungGSE124872)
mLungGSE124872 <- RunPCA(mLungGSE124872, verbose = FALSE)

stdev <- mLungGSE124872@reductions$pca@stdev
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
#Confirm #PC's determined explain > 90% of variance
sum(var[1:PCNum])/ sum(var)

mLungGSE124872 <- RunUMAP(mLungGSE124872, dims = 1:PCNum)
mLungGSE124872 <- FindNeighbors(mLungGSE124872, dims = 1:PCNum) 
mLungGSE124872 <- FindClusters(mLungGSE124872, resolution = 2)

Theis_cell <- read.delim("~/Desktop/Theis_cell.txt")
mLungGSE124872@meta.data[["celltype"]] <- Theis_cell$celltype
Idents(object = mLungGSE124872) <- "celltype"
levels(mLungGSE124872)

Theis_cell <- read.delim("~/Desktop/Theis_age.txt")
mLungGSE124872@meta.data[["Age"]] <- Theis_cell$Age

DimPlot(mLungGSE124872, group.by = c("celltype"), raster = FALSE)
DimPlot(mLungGSE124872, group.by = c("Age"), raster = FALSE)
DimPlot(mLungGSE124872, raster = FALSE)


#Label Trasferring
DefaultAssay(mLung_combined) <- "RNA"

transfer.anchors <- FindTransferAnchors(
  reference = mLungGSE124872,
  normalization.method = "LogNormalize",
  query = mLung_combined,
  dims = 1:30, 
  reference.reduction = "pca",
  features = mLung_combined@assays[["RNA"]]@counts@Dimnames[[1]],
  k.anchor = 20
)

predictions <- TransferData(
  anchorset = transfer.anchors,
  anchors, 
  refdata = mLungGSE124872$celltype,
  dims = 1:30
)

### label transfer
mLung_combined <- AddMetaData(mLung_combined, metadata = predictions)

### umap integration
mLungGSE124872 <- RunUMAP(
  mLungGSE124872, 
  dims = 1:30, 
  reduction = "pca", 
  return.model = TRUE
)

mLung_combined <- MapQuery(
  anchorset = transfer.anchors, 
  reference = mLungGSE124872, 
  query = mLung_combined,
  refdata = list(celltype = "celltype"), 
  reference.reduction = "pca", 
  reduction.model = "umap"
)

Idents(object = mLung_combined) <- "predicted.celltype"
DimPlot(mLung_combined, reduction = 'ref.umap', group.by = 'predicted.celltype')

DimPlot(mLung_combined, reduction = 'umap', group.by = 'predicted.celltype')

saveRDS(mLung_combined, '~/Desktop/mLung_label_transfer.rds')


