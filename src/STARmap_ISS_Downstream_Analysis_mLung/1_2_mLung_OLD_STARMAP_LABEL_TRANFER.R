#SCT normalized the scRNA-seq data###
DefaultAssay(mLung_O.combined.sct) <- "integrated"
DefaultAssay(mLung_O_sc) <- "RNA"
mLung_O_sc <- SCTransform(mLung_O_sc, vst.flavor = "v2", verbose = FALSE, variable.features.n= 3000, min_cells=1, return.only.var.genes = FALSE)
######Integrate scRNA and ST ############

### transfer embedding
transfer.anchors <- FindTransferAnchors(
  reference = mLung_O_sc,
  query = mLung_O.combined.sct,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:30,
  features = mLung_O.combined.sct@assays[["RNA"]]@counts@Dimnames[[1]]
)

predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata = mLung_O_sc$cell_annotated,
  dims = 1:30
)

### label transfer
mLung_O.combined.sct <- AddMetaData(mLung_O.combined.sct, metadata = predictions)

### umap integration
mLung_O_sc <- RunUMAP(
  mLung_O_sc, 
  dims = 1:30, 
  return.model = TRUE
)

mLung_O.combined.sct <- MapQuery(
  anchorset = transfer.anchors, 
  reference = mLung_O_sc, 
  query = mLung_O.combined.sct,
  refdata = list(celltype = "cell_annotated"), 
  reduction.model = "umap"
)

Idents(object = mLung_O.combined.sct) <- "predicted.celltype"
DimPlot(mLung_O.combined.sct, reduction= "ref.umap")

table(mLung_O.combined.sct$predicted.celltype)
save(mLung_O.combined.sct,file="~/Desktop/mLung_O_cells_ST_ANNOTATED.RData")
