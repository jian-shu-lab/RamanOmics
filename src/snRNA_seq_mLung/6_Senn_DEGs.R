library(mclust)
library(Matrix)
library(SeuratDisk)
library(SeuratData)
library(tidyverse)
library(scDblFinder)
library(Seurat)
set.seed(1234)



write.csv(mSkin@meta.data[["orig.ident"]], '~/Desktop/age.txt')
age <- read.delim("~/Desktop/age.txt")
mSkin@meta.data[['Aging']] <- age$x

mSkin$Aging_Senn <- paste(mSkin$Aging, mSkin$p21.expression, sep = "_")
Idents(mSkin) <- "Aging_Senn"
DimPlot(mSkin, reduction = 'umap')
x <- FindMarkers(mSkin, ident.1 = 'Old_p21+', ident.2 = 'Old_p21-' ,only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(x, '~/Desktop/Old_sen_genes.txt')

y <- FindMarkers(mSkin, ident.1 = 'Young_p21+', ident.2 = 'Young_p21-' ,only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(y, '~/Desktop/Young_sen_genes.txt')

Idents(mSkin) <- "p21.expression"
z <- FindMarkers(mSkin, ident.1 = 'p21+', ident.2 = 'p21-' ,only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(z, '~/Desktop/ALL_sen_genes.txt')

mSkin$Aging_Senn_celltype <- paste(mSkin$Aging_Senn, mSkin$final.annotation, sep = "_")
Idents(mSkin) <- "Aging_Senn_celltype"
x1 <- FindMarkers(mSkin, ident.1 = 'Old_p21+_Interfollicular epidermis cell 5', ident.2 = 'Old_p21-_Interfollicular epidermis cell 5' ,only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(x1, '~/Desktop/Old_sen_IFE5_genes.txt')

y1 <- FindMarkers(mSkin, ident.1 = 'Young_p21+_Interfollicular epidermis cell 5', ident.2 = 'Young_p21-_Interfollicular epidermis cell 5' ,only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(y1, '~/Desktop/Young_sen_IFE5_genes.txt')

z1 <- FindMarkers(mSkin, ident.1 = 'Old_p21+_Interfollicular epidermis cell 4', ident.2 = 'Old_p21-_Interfollicular epidermis cell 4' ,only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(z1, '~/Desktop/Old_sen_IFE4_genes.txt')

y2 <- FindMarkers(mSkin, ident.1 = 'Young_p21+_Interfollicular epidermis cell 4', ident.2 = 'Young_p21-_Interfollicular epidermis cell 4' ,only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(y2, '~/Desktop/Young_sen_IFE4_genes.txt')

z2 <- FindMarkers(mSkin, ident.1 = 'Old_p21+_Interfollicular epidermis cell 3', ident.2 = 'Old_p21-_Interfollicular epidermis cell 3' ,only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(z2, '~/Desktop/Old_sen_IFE3_genes.txt')

y3 <- FindMarkers(mSkin, ident.1 = 'Young_p21+_Interfollicular epidermis cell 3', ident.2 = 'Young_p21-_Interfollicular epidermis cell 3' ,only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(y3, '~/Desktop/Young_sen_IFE3_genes.txt')


age <- read.delim("~/Desktop/age.txt")
mLung@meta.data[['Aging']] <- age$x

mLung$Aging_Senn <- paste(mLung$Aging, mLung$p21.expression, sep = "_")
Idents(mLung) <- "Aging_Senn"
DimPlot(mLung, reduction = 'umap')
x <- FindMarkers(mLung, ident.1 = 'Old_p21+', ident.2 = 'Old_p21-' ,only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(x, '~/Desktop/Old_sen_genes.txt')

y <- FindMarkers(mLung, ident.1 = 'Young_p21+', ident.2 = 'Young_p21-' ,only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(y, '~/Desktop/Young_sen_genes.txt')

mLung$Aging_Senn_celltype <- paste(mLung$Aging_Senn, mLung$final.annotation, sep = "_")
Idents(mLung) <- "Aging_Senn_celltype"
x1 <- FindMarkers(mLung, ident.1 = 'Old_p21+_Alveolar type 2 cell 2', ident.2 = 'Old_p21-_Alveolar type 2 cell 2' ,only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(x1, '~/Desktop/Old_sen_AT2_2_genes.txt')

y1 <- FindMarkers(mLung, ident.1 = 'Young_p21+_Alveolar type 2 cell 2', ident.2 = 'Young_p21-_Alveolar type 2 cell 2' ,only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(y1, '~/Desktop/Young_sen_AT2_2_genes.txt')

z1 <- FindMarkers(mLung, ident.1 = 'Old_p21+_Alveolar type 2 cell 1', ident.2 = 'Old_p21-_Alveolar type 2 cell 1' ,only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(z1, '~/Desktop/Old_sen_AT2_1_genes.txt')

y2 <- FindMarkers(mLung, ident.1 = 'Young_p21+_Alveolar type 2 cell 1', ident.2 = 'Young_p21-_Alveolar type 2 cell 1' ,only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(y2, '~/Desktop/Young_sen_AT2_1_genes.txt')

mLung$Anno_clean <- gsub(" \\d+$", "", mLung$final.annotation)
unique(mLung$Anno_clean)
mLung$Aging_Senn_celltype_clean <- paste(mLung$Aging_Senn, mLung$Anno_clean, sep = "_")
Idents(mLung) <- "Aging_Senn_celltype_clean"
x1 <- FindMarkers(mLung, ident.1 = 'Old_p21+_Alveolar type 2 cell', ident.2 = 'Old_p21-_Alveolar type 2 cell' ,only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(x1, '~/Desktop/Old_sen_AT2_genes.txt')





