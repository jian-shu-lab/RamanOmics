suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(mclust))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DoubletFinder))
options(warn=-1)
set.seed(1)

my_paramSweep <- function(seu, PCs=1:10, sct = FALSE, num.cores=1) {
  require(Seurat); require(fields); require(parallel)
  ## Set pN-pK param sweep ranges
  pK <- c(0.0005, 0.001, 0.005, seq(0.01,0.3,by=0.01))
  pN <- seq(0.05,0.3,by=0.05)

  ## Remove pK values with too few cells
  min.cells <- round(nrow(seu@meta.data)/(1-0.05) - nrow(seu@meta.data))
  pK.test <- round(pK*min.cells)
  pK <- pK[which(pK.test >= 1)]

  ## Extract pre-processing parameters from original data analysis workflow
  orig.commands <- seu@commands

  ## Down-sample cells to 10000 (when applicable) for computational effiency
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 10000, replace=FALSE)]
    data <- seu@assays$RNA@counts[ , real.cells]
    n.real.cells <- ncol(data)
  }

  if (nrow(seu@meta.data) <= 10000){
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts
    n.real.cells <- ncol(data)
  }

  ## Iterate through pN, computing pANN vectors at varying pK
  #no_cores <- detectCores()-1
  if(num.cores>1){
    require(parallel)
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)),
                        FUN = parallel_paramSweep,
                        n.real.cells,
                        real.cells,
                        pK,
                        pN,
                        data,
                        orig.commands,
                        PCs,
                        sct,mc.cores=num.cores)
    stopCluster(cl)
  }else{
    output2 <- lapply(as.list(1:length(pN)),
                      FUN = parallel_paramSweep,
                      n.real.cells,
                      real.cells,
                      pK,
                      pN,
                      data,
                      orig.commands,
                      PCs,
                      sct)
  }

  ## Write parallelized output into list
  sweep.res.list <- list()
  list.ind <- 0
  for(i in 1:length(output2)){
    for(j in 1:length(output2[[i]])){
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }

  ## Assign names to list of results
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_" ))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)

}

my_doubletFinder <- function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, sct = FALSE, annotations = NULL) {
  require(Seurat); require(fields); require(KernSmooth)

  ## Generate new list of doublet classificatons from existing pANN vector to save time
  if (reuse.pANN != FALSE ) {
    pANN.old <- seu@meta.data[ , reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    return(seu)
  }

  if (reuse.pANN == FALSE) {
    ## Make merged real-artifical data
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating",n_doublets,"artificial doublets...",sep=" "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    # Keep track of the types of the simulated doublets
    if(!is.null(annotations)){
      stopifnot(typeof(annotations)=="character")
      stopifnot(length(annotations)==length(Cells(seu)))
      stopifnot(!any(is.na(annotations)))
      annotations <- factor(annotations)
      names(annotations) <- Cells(seu)
      doublet_types1 <- annotations[real.cells1]
      doublet_types2 <- annotations[real.cells2]
    }
    ## Store important pre-processing information
    orig.commands <- seu@commands

    ## Pre-process Seurat object
    if (sct == FALSE) {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)

      print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets,
                                     normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
                                     scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
                                     margin = orig.commands$NormalizeData.RNA@params$margin)

      print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(seu_wdoublets,
                                            selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
                                            loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
                                            clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
                                            mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
                                            dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
                                            num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
                                            binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
                                            nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
                                            mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
                                            dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)

      print("Scaling data...")
      seu_wdoublets <- ScaleData(seu_wdoublets,
                                 features = orig.commands$ScaleData.RNA$features,
                                 model.use = orig.commands$ScaleData.RNA$model.use,
                                 do.scale = orig.commands$ScaleData.RNA$do.scale,
                                 do.center = orig.commands$ScaleData.RNA$do.center,
                                 scale.max = orig.commands$ScaleData.RNA$scale.max,
                                 block.size = orig.commands$ScaleData.RNA$block.size,
                                 min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)

      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets,
                              features = orig.commands$ScaleData.RNA$features,
                              npcs = length(PCs),
                              rev.pca =  orig.commands$RunPCA.RNA$rev.pca,
                              weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
                              verbose=FALSE)
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc() # Free up memory
    }

    if (sct == TRUE) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)

      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)

      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc()
    }

    ## Compute PC distance matrix
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)

    ## Compute pANN
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    if(!is.null(annotations)){
      neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = length(levels(doublet_types1))))
    }
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
      if(!is.null(annotations)){
        for(ct in unique(annotations)){
          neighbors_that_are_doublets = neighbors[neighbors>n_real.cells]
          if(length(neighbors_that_are_doublets) > 0){
            neighbor_types[i,] <-
              table( doublet_types1[neighbors_that_are_doublets - n_real.cells] ) +
              table( doublet_types2[neighbors_that_are_doublets - n_real.cells] )
            neighbor_types[i,] <- neighbor_types[i,] / sum( neighbor_types[i,] )
          } else {
            neighbor_types[i,] <- NA
          }
        }
      }
    }
    print("Classifying doublets..")
    classifications <- rep("Singlet",n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN",pN,pK,nExp,sep="_")] <- pANN[rownames(seu@meta.data), 1]
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    if(!is.null(annotations)){
      colnames(neighbor_types) = levels(doublet_types1)
      for(ct in levels(doublet_types1)){
        seu@meta.data[, paste("DF.doublet.contributors",pN,pK,nExp,ct,sep="_")] <- neighbor_types[,ct]
      }
    }
    return(seu)
  }
}

### Old 1 ###
load("/data/amdqiao/5/intermediate/O1_scDblFinder.RData")
dim(mLung_O1)
mLung_O1 = FindVariableFeatures(mLung_O1, verbose = F)
mLung_O1 = ScaleData(mLung_O1, verbose = F)
mLung_O1 = RunPCA(mLung_O1, verbose = F, npcs = 20)
mLung_O1 = RunUMAP(mLung_O1, dims = 1:10, verbose = F)
mLung_O1 <- NormalizeData(mLung_O1, normalization.method = "LogNormalize")

sweep.res <- my_paramSweep(mLung_O1) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
# pdf("/home/amdqiao/5_sennet_raman/img/temp_pK_mLung_O1.pdf")
# barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las = 2)
# dev.off()
nExp <- round(ncol(mLung_O1) * 0.06)
mLung_O1 <- my_doubletFinder(mLung_O1, pN = 0.25, pK = 0.25, nExp = nExp, PCs = 1:10)
DF.name = colnames(mLung_O1@meta.data)[grepl("DF.classification", colnames(mLung_O1@meta.data))]
mLung_O1 = mLung_O1[, mLung_O1@meta.data[, DF.name] == "Singlet"]
mLung_O1 = mLung_O1[, mLung_O1@meta.data[, "scDblFinder"] == "singlet"]
dim(mLung_O1)
save(mLung_O1, file = "/data/amdqiao/5/intermediate/O1_DoubletFinder.RData")


### Old 2 ###
load("/data/amdqiao/5/intermediate/O2_scDblFinder.RData")
dim(mLung_O2)
mLung_O2 = FindVariableFeatures(mLung_O2, verbose = F)
mLung_O2 = ScaleData(mLung_O2, verbose = F)
mLung_O2 = RunPCA(mLung_O2, verbose = F, npcs = 20)
mLung_O2 = RunUMAP(mLung_O2, dims = 1:10, verbose = F)
mLung_O2 <- NormalizeData(mLung_O2, normalization.method = "LogNormalize")

sweep.res <- my_paramSweep(mLung_O2) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
# pdf("/home/amdqiao/5_sennet_raman/img/temp_pK_mLung_O2.pdf")
# barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las = 2)
# dev.off()
nExp <- round(ncol(mLung_O2) * 0.06)
mLung_O2 <- my_doubletFinder(mLung_O2, pN = 0.25, pK = 0.26, nExp = nExp, PCs = 1:10)
DF.name = colnames(mLung_O2@meta.data)[grepl("DF.classification", colnames(mLung_O2@meta.data))]
mLung_O2 = mLung_O2[, mLung_O2@meta.data[, DF.name] == "Singlet"]
mLung_O2 = mLung_O2[, mLung_O2@meta.data[, "scDblFinder"] == "singlet"]
dim(mLung_O2)
save(mLung_O2, file = "/data/amdqiao/5/intermediate/O2_DoubletFinder.RData")


### Old 3 ###
load("/data/amdqiao/5/intermediate/O3_scDblFinder.RData")
dim(mLung_O3)
mLung_O3 = FindVariableFeatures(mLung_O3, verbose = F)
mLung_O3 = ScaleData(mLung_O3, verbose = F)
mLung_O3 = RunPCA(mLung_O3, verbose = F, npcs = 20)
mLung_O3 = RunUMAP(mLung_O3, dims = 1:10, verbose = F)
mLung_O3 <- NormalizeData(mLung_O3, normalization.method = "LogNormalize")

sweep.res <- my_paramSweep(mLung_O3) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
# pdf("/home/amdqiao/5_sennet_raman/img/temp_pK_mLung_O3.pdf")
# barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las = 2)
# dev.off()
nExp <- round(ncol(mLung_O3) * 0.06)
mLung_O3 <- my_doubletFinder(mLung_O3, pN = 0.25, pK = 0.19, nExp = nExp, PCs = 1:10)
DF.name = colnames(mLung_O3@meta.data)[grepl("DF.classification", colnames(mLung_O3@meta.data))]
mLung_O3 = mLung_O3[, mLung_O3@meta.data[, DF.name] == "Singlet"]
mLung_O3 = mLung_O3[, mLung_O3@meta.data[, "scDblFinder"] == "singlet"]
dim(mLung_O3)
save(mLung_O3, file = "/data/amdqiao/5/intermediate/O3_DoubletFinder.RData")


### Young 1 ###
load("/data/amdqiao/5/intermediate/Y1_scDblFinder.RData")
dim(mLung_Y1)
mLung_Y1 = FindVariableFeatures(mLung_Y1, verbose = F)
mLung_Y1 = ScaleData(mLung_Y1, verbose = F)
mLung_Y1 = RunPCA(mLung_Y1, verbose = F, npcs = 20)
mLung_Y1 = RunUMAP(mLung_Y1, dims = 1:10, verbose = F)
mLung_Y1 <- NormalizeData(mLung_Y1, normalization.method = "LogNormalize")

sweep.res <- my_paramSweep(mLung_Y1) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
# pdf("/home/amdqiao/5_sennet_raman/img/temp_pK_mLung_Y1.pdf")
# barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las = 2)
# dev.off()
nExp <- round(ncol(mLung_Y1) * 0.06)
mLung_Y1 <- my_doubletFinder(mLung_Y1, pN = 0.25, pK = 0.12, nExp = nExp, PCs = 1:10)
DF.name = colnames(mLung_Y1@meta.data)[grepl("DF.classification", colnames(mLung_Y1@meta.data))]
mLung_Y1 = mLung_Y1[, mLung_Y1@meta.data[, DF.name] == "Singlet"]
mLung_Y1 = mLung_Y1[, mLung_Y1@meta.data[, "scDblFinder"] == "singlet"]
dim(mLung_Y1)
save(mLung_Y1, file = "/data/amdqiao/5/intermediate/Y1_DoubletFinder.RData")


### Young 2 ###
load("/data/amdqiao/5/intermediate/Y2_scDblFinder.RData")
dim(mLung_Y2)
mLung_Y2 = FindVariableFeatures(mLung_Y2, verbose = F)
mLung_Y2 = ScaleData(mLung_Y2, verbose = F)
mLung_Y2 = RunPCA(mLung_Y2, verbose = F, npcs = 20)
mLung_Y2 = RunUMAP(mLung_Y2, dims = 1:10, verbose = F)
mLung_Y2 <- NormalizeData(mLung_Y2, normalization.method = "LogNormalize")

sweep.res <- my_paramSweep(mLung_Y2) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
# pdf("/home/amdqiao/5_sennet_raman/img/temp_pK_mLung_Y2.pdf")
# barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las = 2)
# dev.off()
nExp <- round(ncol(mLung_Y2) * 0.06)
mLung_Y2 <- my_doubletFinder(mLung_Y2, pN = 0.25, pK = 0.28, nExp = nExp, PCs = 1:10)
DF.name = colnames(mLung_Y2@meta.data)[grepl("DF.classification", colnames(mLung_Y2@meta.data))]
mLung_Y2 = mLung_Y2[, mLung_Y2@meta.data[, DF.name] == "Singlet"]
mLung_Y2 = mLung_Y2[, mLung_Y2@meta.data[, "scDblFinder"] == "singlet"]
dim(mLung_Y2)
save(mLung_Y2, file = "/data/amdqiao/5/intermediate/Y2_DoubletFinder.RData")


### Young 3 ###
load("/data/amdqiao/5/intermediate/Y3_scDblFinder.RData")
dim(mLung_Y3)
mLung_Y3 = FindVariableFeatures(mLung_Y3, verbose = F)
mLung_Y3 = ScaleData(mLung_Y3, verbose = F)
mLung_Y3 = RunPCA(mLung_Y3, verbose = F, npcs = 20)
mLung_Y3 = RunUMAP(mLung_Y3, dims = 1:10, verbose = F)
mLung_Y3 <- NormalizeData(mLung_Y3, normalization.method = "LogNormalize")

sweep.res <- my_paramSweep(mLung_Y3) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
# pdf("/home/amdqiao/5_sennet_raman/img/temp_pK_mLung_Y3.pdf")
# barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las = 2)
# dev.off()
nExp <- round(ncol(mLung_Y3) * 0.06)
mLung_Y3 <- my_doubletFinder(mLung_Y3, pN = 0.25, pK = 0.27, nExp = nExp, PCs = 1:10)
DF.name = colnames(mLung_Y3@meta.data)[grepl("DF.classification", colnames(mLung_Y3@meta.data))]
mLung_Y3 = mLung_Y3[, mLung_Y3@meta.data[, DF.name] == "Singlet"]
mLung_Y3 = mLung_Y3[, mLung_Y3@meta.data[, "scDblFinder"] == "singlet"]
dim(mLung_Y3)
save(mLung_Y3, file = "/data/amdqiao/5/intermediate/Y3_DoubletFinder.RData")
