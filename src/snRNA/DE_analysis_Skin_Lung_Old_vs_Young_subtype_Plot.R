
##==========================================<<<Skin>>>==========================================
library(Seurat)
library(dplyr)
library(ggplot2)
# install.packages("svglite")
setwd('/data/omicseq/Haochun/js_server/SenNet/')
skin<-readRDS('FIgure_Results/Lung/1001_mLung.rds')


###Prepare data####
skin$condition <- sapply(skin$orig.ident, function(x) {
  if (strsplit(x, "_")[[1]][2] %in% c("O1","O2","O3")) {
    "Old"
  } else {
    "Young"
  }})

##Merge subtypes for IFE and FIB
skin$final.annotation<-as.character(skin$final.annotation)
skin$annotation_merged<-ifelse(skin$final.annotation %in%c('Endothelial cell 1','Endothelial cell 2','Endothelial cell 3','Endothelial cell 4','Endothelial cell 5','Endothelial cell 6'),'Endothelial cell',skin$final.annotation)
skin$annotation_merged<-ifelse(skin$annotation_merged %in%c(),'Fibroblast',skin$annotation_merged)
table(skin$annotation_merged)


skin$celltype.samples <- paste(skin$annotation_merged, skin$condition, sep = "_")
Idents(skin) <- "celltype.samples"
DefaultAssay(skin) <- "RNA"
skin <- NormalizeData(skin, normalization.method = "LogNormalize", scale.factor = 10000)


DE_file_list<-list.files(path='FIgure_Results/mSKIN/CellType_Comparisons/',pattern = ".txt$", recursive = TRUE, full.names = TRUE)
DE_file_list
###Crispr####
# Mesenchyme_Kidney_combined <- FindMarkers(Kidney_combined, assay = "RNA", ident.1 = "Mesenchyme_Kidney_Crispr", ident.2 = "Mesenchyme_Kidney_control", test.use ="MAST", min.pct = 0.01, verbose = FALSE , min.cells.group = 10)
# head(Mesenchyme_Kidney_combined, n = 15)
# write.csv(Mesenchyme_Kidney_combined, "~/Desktop/Mesenchyme_Kidney_combined_Crispr_vs_Control.csv")
DE<-read.csv('FIgure_Results/mSKIN/CellType_Comparisons/Adipocyte 1_Old_vs_Young.txt')
DE<-DE_df[DE_df$celltype%in%c('Interfollicular epidermis cell 4','Interfollicular epidermis cell 5','Interfollicular epidermis cell 3','Interfollicular epidermis cell 2','Interfollicular epidermis cell 1'), ]

skin_sub <- subset(skin, idents = c("Interfollicular epidermis cell_Young", "Interfollicular epidermis cell_Old"))

DE %>%
  dplyr::filter((avg_log2FC > 1 | avg_log2FC < -1) & p_val_adj < 0.05) %>%
  ungroup() -> DEgenes


skin_sub <- ScaleData(skin_sub, assay="RNA", verbose = FALSE)
skin_sub$celltype.samples <- factor(skin_sub$celltype.samples, levels = c("Interfollicular epidermis cell_Young", "Interfollicular epidermis cell_Old"))
Idents(skin_sub) <- skin_sub$celltype.samples
# Create a customized heatmap
# Define colors for the heatmap
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)
sorted_genes <- DEgenes %>%
  arrange((avg_log2FC)) %>%
  pull(X)
# Plot the heatmap
p <- DoHeatmap(skin_sub, features = sorted_genes, size = 5) +
  scale_fill_gradientn(colors = heatmap_colors) +
  theme_void() +
  NoLegend()

# Customize further with ggplot2
p <- p + 
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(), # Remove x-axis ticks
    axis.text.y = element_text(size = 10)
  )

# Display the heatmap
print(p)
ggsave("FIgure_Results/Ke_ordered/heatmap_young_old_Interfollicular_epidermis_cell.svg", plot = p, width = 10, height = 8)





output_dir <- "FIgure_Results/DE_Heatmap_Skin"

DE_file_list <- list.files(path = "FIgure_Results/mSKIN/CellType_Comparisons/",
                           pattern = ".txt$", recursive = TRUE, full.names = TRUE)

DE_file_list
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)


for (de_file in DE_file_list) {
  file_base <- basename(de_file)
  celltype <- stringr::str_remove(file_base, "_Old_vs_Young.txt")
  
  idents_use <- c(paste0(celltype, "_Young"), paste0(celltype, "_Old"))
  
  if (!all(idents_use %in% Idents(skin))) {
    message("Skipping: ", celltype, " (not found in Seurat object)")
    next
  }
  
  DE <- read.csv(de_file)
  
  DEgenes <- DE %>%
    filter((avg_log2FC > 1 | avg_log2FC < -1) & p_val_adj < 0.05) %>%
    ungroup()
  
  if (nrow(DEgenes) == 0) {
    message("Skipping: ", celltype, " (no significant genes)")
    next
  }
  
  gene_col <- if ("X" %in% colnames(DEgenes)) "X" else "gene"
  sorted_genes <- DEgenes %>%
    arrange(avg_log2FC) %>%
    pull(gene_col) %>%
    unique()
  
  skin_sub <- subset(skin, idents = idents_use)
  table(Idents(skin_sub))
  skin_sub <- ScaleData(skin_sub, assay = "RNA", verbose = FALSE)
  
  skin_sub$celltype.samples <- factor(skin_sub$celltype.samples, levels = idents_use)
  Idents(skin_sub) <- skin_sub$celltype.samples
  
  p <- DoHeatmap(skin_sub, features = sorted_genes, size = 5) +
    scale_fill_gradientn(colors = heatmap_colors) +
    theme_void() +
    NoLegend() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 10)
    )
  output_path <- file.path(output_dir, paste0("heatmap_young_old_", stringr::str_replace_all(celltype, " ", "_"), ".svg"))
  ggsave(output_path, plot = p, width = 10, height = 8)
  message("Saved: ", output_path)
}




##==========================================<<<Lung>>>==========================================
library(Seurat)
library(dplyr)
library(ggplot2)
# install.packages("svglite")
setwd('/data/Haochun/SenNet/')
lung<-readRDS('FIgure_Results/Lung/1001_mLung.rds')


###Prepare data####
lung$condition <- sapply(lung$orig.ident, function(x) {
  if (strsplit(x, "_")[[1]][2] %in% c("O1","O2","O3")) {
    "Old"
  } else {
    "Young"
  }})

lung$final.annotation<-as.character(lung$final.annotation)
lung$final.annotation<-ifelse(lung$final.annotation=='Krt4/Krt13 Epithelial cell state 1','Krt13 Epithelial cell state 1',lung$final.annotation)
lung$celltype.samples <- paste(lung$final.annotation, lung$condition, sep = "_")
Idents(lung) <- "celltype.samples"
DefaultAssay(lung) <- "RNA"
lung <- NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)

output_dir <- "FIgure_Results/DE_Heatmap_Lung/"

DE_file_list <- list.files(path = "FIgure_Results/mLUNG/CellType_Comparisons",
                           pattern = ".txt$", recursive = TRUE, full.names = TRUE)

DE_file_list
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)


for (de_file in DE_file_list) {
  file_base <- basename(de_file)
  celltype <- stringr::str_remove(file_base, "_Old_vs_Young.txt")
  
  idents_use <- c(paste0(celltype, "_Young"), paste0(celltype, "_Old"))
  
  if (!all(idents_use %in% Idents(lung))) {
    message("Skipping: ", celltype, " (not found in Seurat object)")
    next
  }
  
  DE <- read.csv(de_file)
  
  DEgenes <- DE %>%
    filter((avg_log2FC > 1 | avg_log2FC < -1) & p_val_adj < 0.05) %>%
    ungroup()
  
  if (nrow(DEgenes) == 0) {
    message("Skipping: ", celltype, " (no significant genes)")
    next
  }
  
  gene_col <- if ("X" %in% colnames(DEgenes)) "X" else "gene"
  sorted_genes <- DEgenes %>%
    arrange(avg_log2FC) %>%
    pull(gene_col) %>%
    unique()
  
  lung_sub <- subset(lung, idents = idents_use)
  table(Idents(lung_sub))
  lung_sub <- ScaleData(lung_sub, assay = "RNA", verbose = FALSE)
  
  lung_sub$celltype.samples <- factor(lung_sub$celltype.samples, levels = idents_use)
  Idents(lung_sub) <- lung_sub$celltype.samples
  
  p <- DoHeatmap(lung_sub, features = sorted_genes, size = 5) +
    scale_fill_gradientn(colors = heatmap_colors) +
    theme_void() +
    NoLegend() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 10)
    )
  output_path <- file.path(output_dir, paste0("heatmap_young_old_", stringr::str_replace_all(celltype, " ", "_"), ".svg"))
  ggsave(output_path, plot = p, width = 10, height = 8)
  message("Saved: ", output_path)
}





## Subtype skin

skin<-readRDS('FIgure_Results/Skin/1001_mSkin.rds')
DE_file_list <- list.files(path = "FIgure_Results/mSKIN/CellType_Comparisons/",
                           pattern = ".txt$", recursive = TRUE, full.names = TRUE)

DE_list <- list()
for (file in DE_file_list) {
  DE_sub <- read.csv(file)
  file_base <- basename(file)
  celltype <- stringr::str_remove(file_base, "_Old_vs_Young.txt")
  DE_sub$celltype <- celltype
  DE_list[[celltype]] <- DE_sub
}

DE_df <- bind_rows(DE_list)



###Prepare data####
skin$condition <- sapply(skin$orig.ident, function(x) {
  if (strsplit(x, "_")[[1]][2] %in% c("O1","O2","O3")) {
    "Old"
  } else {
    "Young"
  }})

table(skin$final.annotation)



skin$celltype.samples <- paste(skin$final.annotation, skin$condition, sep = "_")
Idents(skin) <- "celltype.samples"
DefaultAssay(skin) <- "RNA"
skin <- NormalizeData(skin, normalization.method = "LogNormalize", scale.factor = 10000)

DE<-DE_df[DE_df$celltype%in%c('Interfollicular epidermis cell 4','Interfollicular epidermis cell 5','Interfollicular epidermis cell 3','Interfollicular epidermis cell 2','Interfollicular epidermis cell 1'), ]
# DE<-DE_df[DE_df$celltype%in%c('Fibroblast 3','Fibroblast 2','Fibroblast 1'), ]

skin_sub <- subset(skin, idents = c('Interfollicular epidermis cell 4_Old','Interfollicular epidermis cell 5_Old','Interfollicular epidermis cell 3_Old','Interfollicular epidermis cell 2_Old','Interfollicular epidermis cell 1_Old',
                                    'Interfollicular epidermis cell 4_Young','Interfollicular epidermis cell 5_Young','Interfollicular epidermis cell 3_Young','Interfollicular epidermis cell 2_Young','Interfollicular epidermis cell 1_Young'))
# skin_sub <- subset(skin, idents = c('Fibroblast 1_Old','Fibroblast 2_Old','Fibroblast 3_Old',
#                                     'Fibroblast 1_Young','Fibroblast 2_Young','Fibroblast 3_Young'))

DE %>%
  dplyr::filter((avg_log2FC > 1 | avg_log2FC < -1) & p_val_adj < 0.05) %>%
  ungroup() -> DEgenes

skin_sub<-NormalizeData(skin_sub)
skin_sub <- ScaleData(skin_sub, assay="RNA", verbose = FALSE)

# sorted_genes<-c('Cdk8','Rfx2','Tsix','Xist','Gm20388','Gm15564','Tanc2','Ano6','Gpcpd1','Il31ra','Atrnl1','Ttn','Col23a1','Mill1','Lysmd4','Sult5a1','Cmss1')
# sorted_genes<-c('Cdk8','Nuak1','Gm15564','Col1a1','Slit3','Prkg1','Dlc1','Il31ra','Hexb','Cmss1','Tmem56','Kpna1','Gphn','Celf2','Gm26870')

expr_mat <- GetAssayData(skin_sub, assay = "RNA", slot = "scale.data")
avg_expr <- sapply(levels(Idents(skin_sub)), function(ct) {
  rowMeans(expr_mat[, WhichCells(skin_sub, idents = ct), drop = FALSE])
})
d <- dist(t(avg_expr)) 
hc <- hclust(d)        
new_order <- hc$labels[hc$order]  
skin_sub$celltype.samples <- factor(skin_sub$celltype.samples)
Idents(skin_sub) <- skin_sub$celltype.samples
# skin_sub$celltype.samples <- factor(skin_sub$celltype.samples, levels = c('Fibroblast 1_Old','Fibroblast 2_Old','Fibroblast 3_Old',
#                                                                           'Fibroblast 1_Young','Fibroblast 2_Young','Fibroblast 3_Young'))
skin_sub$celltype.samples <- factor(skin_sub$celltype.samples, levels = c('Interfollicular epidermis cell 1_Old','Interfollicular epidermis cell 2_Old','Interfollicular epidermis cell 3_Old','Interfollicular epidermis cell 4_Old','Interfollicular epidermis cell 5_Old',
                                                                          'Interfollicular epidermis cell 1_Young','Interfollicular epidermis cell 2_Young','Interfollicular epidermis cell 3_Young','Interfollicular epidermis cell 4_Young','Interfollicular epidermis cell 5_Young'))
# 

# Idents(skin_sub) <- skin_sub$celltype.samples
# Create a customized heatmap
# Define colors for the heatmap
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)
# sorted_genes <- DEgenes %>%
#   arrange((avg_log2FC)) %>%
#   pull(X)

heatmap_colors <- colorRampPalette(c("#2162a4", "white", "#ad1a2c"))(100)
sorted_genes <- DEgenes %>%
  arrange((avg_log2FC)) %>%
  pull(X)

mat <- GetAssayData(skin_sub, slot = "scale.data")
mat_sub <- mat[sorted_genes, ]
dist_rows <- dist(mat_sub)
clust_rows <- hclust(dist_rows)
gene_order <- clust_rows$labels[clust_rows$order]


# avg_expr <- AverageExpression(skin_sub, group.by = "celltype.samples", slot = "scale.data", assay = "RNA",features = sorted_genes)$RNA
# dist_cols <- dist(t(avg_expr))
# clust_cols <- hclust(dist_cols)
# group_order <- clust_cols$labels[clust_cols$order]
# skin_sub$celltype.samples <- factor(skin_sub$celltype.samples, levels = group_order)
# Idents(skin_sub) <- skin_sub$celltype.samples

#"#FFF9E5","#FFD84D","#E69F00","#EAF7E3","#7FC97F","#006837"
# Plot the heatmap
p <- DoHeatmap(skin_sub, features = gene_order, size = 5,group.by = 'celltype.samples',group.colors=c("#E69F00","#E69F00","#E69F00","#E69F00","#E69F00","#7FC97F","#7FC97F","#7FC97F","#7FC97F","#7FC97F") ) +
  scale_fill_gradientn(colors = heatmap_colors) 



# Customize further with ggplot2
p <- p + 
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(), # Remove x-axis ticks
    axis.text.y = element_text(size = 10)
  )

# Display the heatmap
print(p)
ggsave("FIgure_Results/Ke_ordered_updated_2/heatmap_young_old_Skin_Interfollicular_epidermis_cell_merged_color_bar.svg", plot = p, width = 12, height = 10)


## Subtype lung

lung<-readRDS('FIgure_Results/Lung/1001_mLung.rds')
DE_file_list <- list.files(path = "FIgure_Results/mLUNG/CellType_Comparisons/",
                           pattern = ".txt$", recursive = TRUE, full.names = TRUE)

DE_list <- list()
for (file in DE_file_list) {
  DE_sub <- read.csv(file)
  file_base <- basename(file)
  celltype <- stringr::str_remove(file_base, "_Old_vs_Young.txt")
  DE_sub$celltype <- celltype
  DE_list[[celltype]] <- DE_sub
}

DE_df <- bind_rows(DE_list)



###Prepare data####
lung$condition <- sapply(lung$orig.ident, function(x) {
  if (strsplit(x, "_")[[1]][2] %in% c("O1","O2","O3")) {
    "Old"
  } else {
    "Young"
  }})

table(lung$final.annotation)



lung$celltype.samples <- paste(lung$final.annotation, lung$condition, sep = "_")
Idents(lung) <- "celltype.samples"
DefaultAssay(lung) <- "RNA"
lung <- NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)

DE<-DE_df[DE_df$celltype%in%c('Endothelial cell 1','Endothelial cell 2','Endothelial cell 3','Endothelial cell 4','Endothelial cell 5','Endothelial cell 6'), ]

lung_sub <- subset(lung, idents = c('Endothelial cell 1_Old','Endothelial cell 2_Old','Endothelial cell 3_Old','Endothelial cell 4_Old','Endothelial cell 5_Old','Endothelial cell 6_Old',
                                    'Endothelial cell 1_Young','Endothelial cell 2_Young','Endothelial cell 3_Young','Endothelial cell 4_Young','Endothelial cell 5_Young','Endothelial cell 6_Young'))
# DE<-DE_df[DE_df$celltype%in%c('T cell 1'), ]
# 
# lung_sub <- subset(lung, idents = c('T cell 1_Old',
#                                     'T cell 1_Young'))

DE %>%
  dplyr::filter((avg_log2FC > 1 | avg_log2FC < -1) & p_val_adj < 0.05) %>%
  ungroup() -> DEgenes

lung_sub<-NormalizeData(lung_sub)
lung_sub <- ScaleData(lung_sub, assay="RNA", verbose = FALSE)

# sorted_genes<-c('Klf2','Klf4','Auts2','Tnik','Setbp1','Sox5','Dusp1','Hdac9','Tnip1','Tmsb4x','Vwf','Edn1','Timp3','Prkg1','Tmem100',
#                 'Cdh13','Timp3','Thsd4','Sema5a','Sema6d','Vmp1','Gda','Pcolce2','Serpine1','Mmp16','Slc6a6','Pde4b','Pde4d',
#                 'Pde1a','Gm20388','Gm26883','Acox1','H2-Aa','H2-Eb1','Ciita','Cd74','B2m','Ly6a','Ly6c2','Xist',
#                 'Hbb-bs','Hbb-bt','Hba-a1','Hba-a2','Setbp1','Hdac9','Dusp1','Fndc1','Zbtb7c','Gm43291','AY036118','Gm20388','Gm26883','C130026I21Rik')
# sorted_genes<-c('Tnfaip3','Il6ra','Ccl5','S100a6','Hspa8','Tacr1','Foxp1','Lef1','Tcf7','Satb1','Aff3','Pecam1','St6gal1','Klf6',
#                 'Maml3','Arap2','Runx2','Denn4a','Setbp1','Tmsb4x','Hip1','Hbb-bs','Sftpc','Fm43291','Gm20917','AY038118',
#                 'C130011G23Rik','Gm43291')

# expr_mat <- GetAssayData(lung_sub, assay = "RNA", slot = "scale.data")[sorted_genes, ]
# avg_expr <- sapply(levels(Idents(lung_sub)), function(ct) {
#   rowMeans(expr_mat[, WhichCells(lung_sub, idents = ct), drop = FALSE])
# })
# d <- dist(t(avg_expr)) 
# hc <- hclust(d)        
# new_order <- hc$labels[hc$order]  
# lung_sub$celltype.samples <- factor(lung_sub$celltype.samples, levels = c('Endothelial cell 4_Old','Endothelial cell 5_Old','Endothelial cell 6_Old','Endothelial cell 3_Old','Endothelial cell 2_Old','Endothelial cell 1_Old'))
# Idents(lung_sub) <- lung_sub$celltype.samples
# 
# 
# 
# DE<-DE_df[DE_df$celltype%in%c('T cell 1'), ]
# 
# lung_sub <- subset(lung, idents = c('T cell 1_Old'))
# 
# DE %>%
#   dplyr::filter((avg_log2FC > 1 | avg_log2FC < -1) & p_val_adj < 0.05) %>%
#   ungroup() -> DEgenes
# 
# 
# lung_sub <- ScaleData(lung_sub, assay="RNA", verbose = FALSE)
# lung_sub$celltype.samples <- factor(lung_sub$celltype.samples, levels = c('Endothelial cell 2_Old','Endothelial cell 1_Old','Endothelial cell 4_Old','Endothelial cell 5_Old','Endothelial cell 3_Old','Endothelial cell 6_Old',
#                                                                           'Endothelial cell 2_Young','Endothelial cell 1_Young','Endothelial cell 4_Young','Endothelial cell 5_Young','Endothelial cell 3_Young','Endothelial cell 6_Young'))
lung_sub$celltype.samples <- factor(lung_sub$celltype.samples, levels = c('Endothelial cell 1_Old','Endothelial cell 2_Old','Endothelial cell 3_Old','Endothelial cell 4_Old','Endothelial cell 5_Old','Endothelial cell 6_Old',
                                                                          'Endothelial cell 1_Young','Endothelial cell 2_Young','Endothelial cell 3_Young','Endothelial cell 4_Young','Endothelial cell 5_Young','Endothelial cell 6_Young'))


Idents(lung_sub) <- lung_sub$celltype.samples
# Create a customized heatmap
# Define colors for the heatmap
heatmap_colors <- colorRampPalette(c("#2162a4", "white", "#ad1a2c"))(100)
sorted_genes <- DEgenes %>%
  arrange((avg_log2FC)) %>%
  pull(X)

mat <- GetAssayData(lung_sub, slot = "scale.data")
mat_sub <- mat[sorted_genes, ]
dist_rows <- dist(mat_sub)
clust_rows <- hclust(dist_rows)
gene_order <- clust_rows$labels[clust_rows$order]


# avg_expr <- AverageExpression(lung_sub, group.by = "celltype.samples", slot = "scale.data", assay = "RNA",features = sorted_genes)$RNA
# dist_cols <- dist(t(avg_expr))
# clust_cols <- hclust(dist_cols)
# group_order <- clust_cols$labels[clust_cols$order]
# lung_sub$celltype.samples <- factor(lung_sub$celltype.samples, levels = group_order)
# Idents(lung_sub) <- lung_sub$celltype.samples
# Plot the heatmap
p <- DoHeatmap(lung_sub, features =gene_order, size = 5,group.by = 'celltype.samples',group.colors=c("#E69F00","#E69F00","#E69F00","#E69F00","#E69F00","#E69F00","#7FC97F","#7FC97F","#7FC97F","#7FC97F","#7FC97F","#7FC97F")) +
  scale_fill_gradientn(colors = heatmap_colors) +
  theme_void() +
  NoLegend()

# Customize further with ggplot2
p <- p + 
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(), # Remove x-axis ticks
    axis.text.y = element_text(size = 10)
  )

# Display the heatmap
print(p)
ggsave("FIgure_Results/Ke_ordered_updated/heatmap_young_old_lung_Endothelial_merged.pdf", plot = p, width = 10, height = 16)

