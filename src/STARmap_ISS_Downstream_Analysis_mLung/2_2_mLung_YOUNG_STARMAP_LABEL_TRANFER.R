#entire dataset###
DefaultAssay(mLung_Y.combined.sct) <- "integrated"
DefaultAssay(mLung_Y_sc) <- "RNA"
mLung_Y_sc <- SCTransform(mLung_Y_sc, vst.flavor = "v2", verbose = FALSE, variable.features.n= 3000, min_cells=1, return.only.var.genes = FALSE)
######Integrate scRNA and ST ############

### transfer embedding
transfer.anchors <- FindTransferAnchors(
  reference = mLung_Y_sc,
  query = mLung_Y.combined.sct,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:30,
  features = mLung_Y.combined.sct@assays[["RNA"]]@counts@Dimnames[[1]]
)

predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata = mLung_Y_sc$cell_annotated,
  dims = 1:30
)

### label transfer
mLung_Y.combined.sct <- AddMetaData(mLung_Y.combined.sct, metadata = predictions)

### umap integration
mLung_Y_sc <- RunUMAP(
  mLung_Y_sc, 
  dims = 1:30, 
  return.model = TRUE
)

mLung_Y.combined.sct <- MapQuery(
  anchorset = transfer.anchors, 
  reference = mLung_Y_sc, 
  query = mLung_Y.combined.sct,
  refdata = list(celltype = "cell_annotated"), 
  reduction.model = "umap"
)

Idents(object = mLung_Y.combined.sct) <- "predicted.celltype"
DimPlot(mLung_Y.combined.sct, reduction= "ref.umap")

table(mLung_Y.combined.sct$predicted.celltype)
save(mLung_Y.combined.sct,file="~/Desktop/mLung_Y_cells_ST_ANNOTATED.RData")
write.csv(mLung_Y.combined.sct@meta.data, paste0('~/Desktop/Young_label_transferred_annotations_all.txt'))

DefaultAssay(mLung_Y.combined.sct) <- "RNA"
mLung_Y.combined.sct <- NormalizeData(mLung_Y.combined.sct)
marke <- FindAllMarkers(mLung_Y.combined.sct, only.pos = TRUE, logfc.threshYoung=0.1, min.pct = 0.01)
write.csv(marke, '~/Desktop/marke_starmap_label_transfer1.txt')


