library(Seurat)
library(CellChat)
library(tidyverse)
set.seed(1234)
setwd('/data/omicseq/Haochun/js_server/SenNet/FIgure_Results/')
mSkin <- readRDS('Skin/1001_mSkin.rds')
Idents(mSkin) <- mSkin$final.annotation
DimPlot(mSkin,raster=FALSE,reduction='umap',label=TRUE)

mSkin$condition <- sapply(mSkin$orig.ident, function(x) {
  if (strsplit(x, "_")[[1]][2] %in% c("O1","O2","O3")) {
    "Old"
  } else {
    "Young"
  }})

mSkin_Young <- subset(x = mSkin, subset = condition =='Young')
mSkin_Old <- subset(x = mSkin, subset = condition=='Old')
rm(mSkin)

mSkin_Young$samples<-mSkin_Young$orig.ident
mSkin_Old$samples<-mSkin_Old$orig.ident

mSkin_Young$samples<-as.factor(mSkin_Young$samples)
mSkin_Old$samples<-as.factor(mSkin_Old$samples)

cellchatYoung <- createCellChat(object = mSkin_Young, group.by = "final.annotation", assay = "RNA")
cellchatOld <- createCellChat(object = mSkin_Old, group.by = "final.annotation", assay = "RNA")

CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB)
cellchatYoung@DB <- CellChatDB.use
cellchatOld@DB <- CellChatDB.use

cellchatYoung <- subsetData(cellchatYoung)
cellchatOld <- subsetData(cellchatOld)

cellchatYoung <- identifyOverExpressedGenes(cellchatYoung)
cellchatYoung <- identifyOverExpressedInteractions(cellchatYoung)

cellchatOld <- identifyOverExpressedGenes(cellchatOld)
cellchatOld <- identifyOverExpressedInteractions(cellchatOld)

cellchatYoung <- computeCommunProb(cellchatYoung, type = "triMean", population.size = TRUE)
cellchatOld <- computeCommunProb(cellchatOld, type = "triMean", population.size = TRUE)

cellchatYoung <- filterCommunication(cellchatYoung, min.cells = 10)
cellchatOld <- filterCommunication(cellchatOld, min.cells = 10)

cellchatYoung <- computeCommunProbPathway(cellchatYoung)
cellchatOld <- computeCommunProbPathway(cellchatOld)

cellchatYoung <- aggregateNet(cellchatYoung)
cellchatOld <- aggregateNet(cellchatOld)

load('Cellchat/Skin/cellchatYoung.RData')
load('Cellchat/Skin/cellchatOld.RData')


celltypes <- levels(cellchatYoung@idents)
celltypes <- rownames(mat)

# group.colors <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(length(celltypes)), celltypes)
group.colors <- c(
  "Interfollicular epidermis cell 1" = "#A6CEE3",
  "Interfollicular epidermis cell 2" = "#6AA8CE",
  "Interfollicular epidermis cell 3" = "#2F82B9",
  "Merkel cell 1" = "#4E98A6",
  "Fibroblast 1" = "#8EC694",
  "Inner root sheath cell 1" = "#98D277",
  "Early granular keratinocyte 1" = "#60B64D",
  "Interfollicular epidermis cell 4" = "#439F34",
  "Basal keratinocyte 1" = "#9B9C64",
  "Interfollicular epidermis cell 5" = "#F29A94",
  "Muscle cell 1" = "#F16667",
  "Interfollicular keratynocite 1" = "#E62E30",
  "Hair follicle stem cell 1" = "#EA4833",
  "Adipocyte 1" = "#F59057",
  "Lagerhans cell 1" = "#FDB45D",
  "Outer bulge cell 1" = "#FE982C",
  "Fibroblast 2" = "#FC8108",
  "Sebocyte 1" = "#E59766",
  "Granular layer keratinocyte 1" = "#CEADC4",
  "Endothelial/Smooth muscle cell 1" = "#A787C0",
  "Dermal papilla cells 1" = "#7D54A5",
  "Cycling basal cell 1" = "#8D6B99",
  "T cell 1" = "#CFC099",
  "Fibroblast 3" = "#F5EB8B",
  "Macrophage 1" = "#D3A259",
  "Melanocytes 1" = "#B15928"
)

groupSize <- as.numeric(table(cellchatYoung@idents))
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchatYoung@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use= group.colors ,vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)))
netVisual_circle(cellchatYoung@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use= group.colors ,vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)))


library(RColorBrewer)

mat <- cellchatYoung@net$weight
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  safe_name <- gsub("/", "_", rownames(mat)[i])
  
  svg(filename = paste0("Cellchat/Skin/split_circle_plot_Young/circle_plot_Young_withLegend_", safe_name, ".svg"), width = 8, height = 8)
  
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat) * 10,
    title.name = rownames(mat)[i],
     color.use= group.colors  
  )
  
  dev.off()
}

plot.new()

legend(
  "center",  
  legend = names(group.colors),
  col = group.colors,
  pch = 19,
  cex = 0.9,
  pt.cex = 1,
  x.intersp = 0.6,
  y.intersp = 0.9,
  bty = "n",
  title = "Cell types",
  ncol = 2
)

groupSize <- as.numeric(table(cellchatOld@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchatOld@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use= group.colors ,vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)))
netVisual_circle(cellchatOld@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use= group.colors ,vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)))


mat <- cellchatOld@net$weight

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  safe_name <- gsub("/", "_", rownames(mat)[i])
  
  svg(filename = paste0("Cellchat/Skin/split_circle_plot_Old/circle_plot_Old_withLegend_", safe_name, ".svg"), width = 8, height = 8)
  
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat) * 10,
    title.name = rownames(mat)[i],
    # vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)),  
    color.use= group.colors  
  )
  
  dev.off()
}

plot.new()

legend(
  "center",  
  legend = names(group.colors),
  col = group.colors,
  pch = 19,
  cex = 0.9,
  pt.cex = 1,
  x.intersp = 0.6,
  y.intersp = 0.9,
  bty = "n",
  title = "Cell types",
  ncol = 2
)


load('Cellchat/Skin/cellchatYoung.RData')
load('Cellchat/Skin/cellchatOld.RData')

object.list <- list(Old = cellchatOld, Young = cellchatYoung)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat<- liftCellChat(cellchat, group.new = levels(cellchat@idents$joint))

# rm(cellchatOld)
# rm(cellchatYoung)

save(object.list, file = "~/Desktop/Hear/cellCHAT/cellchat_object.list_mouse_HEART_Old_CTL.RData")
save(cellchat, file = "~/Desktop/Hear/cellCHAT/cellchat_merged_mouse_HEART_Old_CTL.RData")


gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
library(reticulate)
py_install("umap-learn")
py_module_available("umap")  
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

rankSimilarity(cellchat, type = "functional")
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
gg1 + gg2


df.rank <- rankNet(cellchat, mode = "comparison", measure = "weight", stacked = FALSE, do.stat = FALSE)$data
library(tidyr)
df.wide <- df.rank %>%
  select(name,contribution,group)%>%
  tidyr::pivot_wider(names_from = group, values_from = contribution)
# library(ggplot2)
#   
# ggplot(df.wide, aes(x = Young, y = Old)) +
#   geom_point(color = "steelblue", size = 3, alpha = 0.7) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
#   scale_x_log10() +
#   scale_y_log10() +
#   theme_classic(base_size = 14) +
#   labs(x = "Pathway strength (Young, log10)",
#        y = "Pathway strength (Old, log10)",
#        title = "Pathway Communication Strength: Young vs Old (log10 scale)") +
#   ggrepel::geom_text_repel(
#     data = df.wide %>%
#       mutate(diff = abs(Old - Young)) %>%
#       arrange(desc(diff)) %>%
#       head(100),
#     aes(label = name),
#     size = 3
#   )

library(dplyr)
library(ggplot2)
library(ggforce)

pseudo <- 1e-5
df.wide <- df.wide %>%
  mutate(
    Old = Old + pseudo,
    Young = Young + pseudo,
    total = Old + Young,
    pct_old = Old / (Old + Young)
  )

df.pies <- df.wide %>%
  mutate(x = Young, y = Old) %>%
  rowwise() %>%
  do({
    data.frame(
      pathway = .$name,
      x = .$x,
      y = .$y,
      start = c(0, 2*pi*.$pct_old),
      end   = c(2*pi*.$pct_old, 2*pi),
      fill  = c("Old", "Young")
    )
  })

p<-ggplot(df.pies) +
  geom_arc_bar(aes(x0 = x, y0 = y, r0 = 0, r = 0.03, start = start, end = end, fill = fill),
               color = "black") +
  scale_fill_manual(values = c("Old" = "#ee918b", "Young" = "#b3cde6")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_log10() + scale_y_log10() +
  coord_equal() +
  theme_classic() +
  labs(x = "Young strength (log10+pseudo)", y = "Old strength (log10+pseudo)",
       title = "Pathway strength with pie-chart markers")+
  ggrepel::geom_text_repel(
    data = df.pies %>% distinct(pathway, x, y),
    aes(x = x, y = y, label = pathway),
    size = 3,
    color = "black",
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = 50
  )
ggsave("Cellchat/Skin/pathway_scatter_pie_labeled.tiff", plot = p, width = 10, height = 8, dpi = 300, compression = "lzw")
ggsave("Cellchat/Skin/pathway_scatter_pie_labeled.pdf",  plot = p, width = 10, height = 8)   
ggsave("Cellchat/Skin/pathway_scatter_pie_labeled.svg",  plot = p, width = 10, height = 8)


df_old.net <- subsetCommunication(cellchatOld)
df_young.net <- subsetCommunication(cellchatYoung)

write.csv(df_old.net,'Cellchat/Skin/old_net.csv')
write.csv(df_young.net,'Cellchat/Skin/young_net.csv')


unique(df_young.net$pathway_name)
pathways.show<-'EGF'
svg("Cellchat/Skin/test.svg", width = 18, height = 12)
circlize::circos.clear()
netVisual_aggregate(cellchatYoung, signaling = pathways.show, layout = "chord",color.use = group.colors,show.legend = T,signaling.name = pathways.show)
dev.off()

# load('Cellchat/Skin/cellchatYoung.RData')
# load('Cellchat/Skin/cellchatOld.RData')
# 
# df_young.net<-read.csv("Cellchat/Skin/young_net.csv")
# df_old.net<-read.csv("Cellchat/Skin/old_net.csv")
pathways.all <- unique(df_young.net$pathway_name)

for (pw in pathways.all) {
  safe_pw <- gsub("[/\\:*?\"<>| ]", "_", pw)
  
  svg_file <- paste0("Cellchat/Skin/split_pathway_plot_Young/aggregate_chord_", safe_pw, ".pdf")
  # tiff(svg_file, width = 18, height = 12, units = "in", res = 300, compression = "lzw")
  pdf(svg_file, width = 18, height = 12) 

  # svg(svg_file, width = 18, height = 12)
  circlize::circos.clear()
  
  tryCatch({
    netVisual_aggregate(
      object = cellchatYoung,
      signaling = pw,
      layout = "chord",
      color.use = group.colors,
      show.legend = TRUE,
      signaling.name = pw
    )
  }, error = function(e) {
    message("Failed to draw ", pw, ": ", e$message)
  })
  
  dev.off()
}



pathways.all <- unique(df_old.net$pathway_name)

for (pw in pathways.all) {
  safe_pw <- gsub("[/\\:*?\"<>| ]", "_", pw)
  
  svg_file <- paste0("Cellchat/Skin/split_pathway_plot_Old/aggregate_chord_", safe_pw, ".pdf")
  # tiff(svg_file, width = 18, height = 12, units = "in", res = 300, compression = "lzw")
  pdf(svg_file, width = 18, height = 12) 
  # svg(svg_file, width = 18, height = 12)
  circlize::circos.clear()
  
  tryCatch({
    netVisual_aggregate(
      object = cellchatOld,
      signaling = pw,
      layout = "chord",
      color.use = group.colors,
      show.legend = TRUE,
      signaling.name = pw
    )
  }, error = function(e) {
    message("Failed to draw ", pw, ": ", e$message)
  })
  
  dev.off()
}
netVisual_bubble(cellchatOld, remove.isolate = FALSE)




##Lung
library(Seurat)
library(CellChat)
library(tidyverse)
set.seed(1234)
setwd('/data/omicseq/Haochun/js_server/SenNet/FIgure_Results/')
mlung <- readRDS('Lung/1001_mLung.rds')
Idents(mlung) <- mlung$final.annotation
DimPlot(mlung,raster=FALSE,reduction='umap',label=TRUE)

mlung$condition <- sapply(mlung$orig.ident, function(x) {
  if (strsplit(x, "_")[[1]][2] %in% c("O1","O2","O3")) {
    "Old"
  } else {
    "Young"
  }})

mlung_Young <- subset(x = mlung, subset = condition =='Young')
mlung_Old <- subset(x = mlung, subset = condition=='Old')
rm(mlung)

mlung_Young$samples<-mlung_Young$orig.ident
mlung_Old$samples<-mlung_Old$orig.ident

mlung_Young$samples<-as.factor(mlung_Young$samples)
mlung_Old$samples<-as.factor(mlung_Old$samples)

cellchatYoung <- createCellChat(object = mlung_Young, group.by = "final.annotation", assay = "RNA")
cellchatOld <- createCellChat(object = mlung_Old, group.by = "final.annotation", assay = "RNA")

CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB)
cellchatYoung@DB <- CellChatDB.use
cellchatOld@DB <- CellChatDB.use

cellchatYoung <- subsetData(cellchatYoung)
cellchatOld <- subsetData(cellchatOld)

cellchatYoung <- identifyOverExpressedGenes(cellchatYoung)
cellchatYoung <- identifyOverExpressedInteractions(cellchatYoung)

cellchatOld <- identifyOverExpressedGenes(cellchatOld)
cellchatOld <- identifyOverExpressedInteractions(cellchatOld)

cellchatYoung <- computeCommunProb(cellchatYoung, type = "triMean", population.size = TRUE)
cellchatOld <- computeCommunProb(cellchatOld, type = "triMean", population.size = TRUE)

cellchatYoung <- filterCommunication(cellchatYoung, min.cells = 10)
cellchatOld <- filterCommunication(cellchatOld, min.cells = 10)

cellchatYoung <- computeCommunProbPathway(cellchatYoung)
cellchatOld <- computeCommunProbPathway(cellchatOld)

cellchatYoung <- aggregateNet(cellchatYoung)
cellchatOld <- aggregateNet(cellchatOld)




library(RColorBrewer)
mat <- cellchatYoung@net$weight
celltypes <- levels(cellchatYoung@idents)
celltypes <- rownames(mat)
group.colors <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(length(celltypes)), celltypes)


groupSize <- as.numeric(table(cellchatYoung@idents))
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchatYoung@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use= group.colors ,vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)))
netVisual_circle(cellchatYoung@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use= group.colors ,vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)))


mat <- cellchatYoung@net$weight
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  safe_name <- gsub("/", "_", rownames(mat)[i])
  
  svg(filename = paste0("Cellchat/Lung/split_circle_plot_Young/circle_plot_Young_", safe_name, ".svg"), width = 8, height = 8)
  
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat) * 10,
    title.name = rownames(mat)[i],
    vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)),
    color.use= group.colors  
  )
  
  dev.off()
}

plot.new()

legend(
  "center",  
  legend = names(group.colors),
  col = group.colors,
  pch = 19,
  cex = 0.9,
  pt.cex = 1,
  x.intersp = 0.6,
  y.intersp = 0.9,
  bty = "n",
  title = "Cell types",
  ncol = 2
)

groupSize <- as.numeric(table(cellchatOld@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchatOld@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use= group.colors ,vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)))
netVisual_circle(cellchatOld@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use= group.colors ,vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)))


mat <- cellchatOld@net$weight

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  safe_name <- gsub("/", "_", rownames(mat)[i])
  
  svg(filename = paste0("Cellchat/Lung/split_circle_plot_Old/circle_plot_Old_", safe_name, ".svg"), width = 8, height = 8)
  
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat) * 10,
    title.name = rownames(mat)[i],
    vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)),
    color.use= group.colors  
  )
  
  dev.off()
}

plot.new()

legend(
  "center",  
  legend = names(group.colors),
  col = group.colors,
  pch = 19,
  cex = 0.9,
  pt.cex = 1,
  x.intersp = 0.6,
  y.intersp = 0.9,
  bty = "n",
  title = "Cell types",
  ncol = 2
)

save(cellchatYoung,file = 'Cellchat/Skin/cellchatYoung.RData')
save(cellchatOld,file = 'Cellchat/Skin/cellchatOld.RData')


load('Cellchat/Lung/cellchatYoung.RData')
load('Cellchat/Lung/cellchatOld.RData')


object.list <- list(Old = cellchatOld, Young= cellchatYoung)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat<- liftCellChat(cellchat, group.new = levels(cellchat@idents$joint))
# 
# rm(cellchatOld)
# rm(cellchatYoung)

save(object.list, file = "~/Desktop/Hear/cellCHAT/cellchat_object.list_mouse_HEART_Old_CTL.RData")
save(cellchat, file = "~/Desktop/Hear/cellCHAT/cellchat_merged_mouse_HEART_Old_CTL.RData")


gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
library(reticulate)
py_install("umap-learn")
py_module_available("umap")  
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)


rankSimilarity(cellchat, type = "functional")
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2


df.rank <- rankNet(cellchat, mode = "comparison", measure = "weight", stacked = FALSE, do.stat = FALSE)$data
library(tidyr)
df.wide <- df.rank %>%
  select(name,contribution,group)%>%
  tidyr::pivot_wider(names_from = group, values_from = contribution)
# library(ggplot2)
#   
# ggplot(df.wide, aes(x = Young, y = Old)) +
#   geom_point(color = "steelblue", size = 3, alpha = 0.7) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
#   scale_x_log10() +
#   scale_y_log10() +
#   theme_classic(base_size = 14) +
#   labs(x = "Pathway strength (Young, log10)",
#        y = "Pathway strength (Old, log10)",
#        title = "Pathway Communication Strength: Young vs Old (log10 scale)") +
#   ggrepel::geom_text_repel(
#     data = df.wide %>%
#       mutate(diff = abs(Old - Young)) %>%
#       arrange(desc(diff)) %>%
#       head(100),
#     aes(label = name),
#     size = 3
#   )

library(dplyr)
library(ggplot2)
library(ggforce)

pseudo <- 1e-5
df.wide <- df.wide %>%
  mutate(
    Old = Old + pseudo,
    Young = Young + pseudo,
    total = Old + Young,
    pct_old = Old / (Old + Young)
  )

df.pies <- df.wide %>%
  mutate(x = Young, y = Old) %>%
  rowwise() %>%
  do({
    data.frame(
      pathway = .$name,
      x = .$x,
      y = .$y,
      start = c(0, 2*pi*.$pct_old),
      end   = c(2*pi*.$pct_old, 2*pi),
      fill  = c("Old", "Young")
    )
  })

p<-ggplot(df.pies) +
  geom_arc_bar(aes(x0 = x, y0 = y, r0 = 0, r = 0.05, start = start, end = end, fill = fill),
               color = "black") +
  scale_fill_manual(values = c("Old" = "#ee918b", "Young" = "#b3cde6")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_log10() + scale_y_log10() +
  coord_equal() +
  theme_classic() +
  labs(x = "Young strength (log10+pseudo)", y = "Old strength (log10+pseudo)",
       title = "Pathway strength with pie-chart markers")
  # ggrepel::geom_text_repel(
  #   data = df.pies %>% distinct(pathway, x, y),
  #   aes(x = x, y = y, label = pathway),
  #   size = 3,
  #   color = "black",
  #   box.padding = 0.3,
  #   point.padding = 0.2,
  #   max.overlaps = 50
  # )
ggsave("Cellchat/Lung/pathway_scatter_pie.tiff", plot = p, width = 10, height = 8, dpi = 300, compression = "lzw")
ggsave("Cellchat/Lung/pathway_scatter_pie.pdf",  plot = p, width = 10, height = 8)   
ggsave("Cellchat/Lung/pathway_scatter_pie.svg",  plot = p, width = 10, height = 8)


df_old.net <- subsetCommunication(cellchatOld)
df_young.net <- subsetCommunication(cellchatYoung)

write.csv(df_old.net,'Cellchat/Lung/old_net.csv')
write.csv(df_young.net,'Cellchat/Lung/young_net.csv')


unique(df_young.net$pathway_name)
pathways.show<-'EGF'
svg("Cellchat/lung/test.svg", width = 18, height = 12)
circlize::circos.clear()
netVisual_aggregate(cellchatYoung, signaling = pathways.show, layout = "chord",color.use = group.colors,show.legend = T,signaling.name = pathways.show)
dev.off()


df_young.net<-read.csv('Cellchat/Lung/young_net.csv')




# 
# load('Cellchat/Lung/cellchatYoung.RData')
# load('Cellchat/Lung/cellchatOld.RData')
# 
# df_young.net<-read.csv("Cellchat/Lung/young_net.csv")
# df_old.net<-read.csv("Cellchat/Lung/old_net.csv")


pathways.all <- unique(df_young.net$pathway_name)

for (pw in pathways.all) {
  safe_pw <- gsub("[/\\:*?\"<>| ]", "_", pw)
  
  svg_file <- paste0("Cellchat/Lung/split_pathway_plot_Young/aggregate_chord_", safe_pw, ".pdf")
  # tiff(svg_file, width = 18, height = 12, units = "in", res = 300, compression = "lzw")
  pdf(svg_file, width = 18, height = 12)
  # svg(svg_file, width = 18, height = 12)
  circlize::circos.clear()
  
  tryCatch({
    netVisual_aggregate(
      object = cellchatYoung,
      signaling = pw,
      layout = "chord",
      color.use = group.colors,
      show.legend = TRUE,
      signaling.name = pw
    )
  }, error = function(e) {
    message("Failed to draw ", pw, ": ", e$message)
  })
  
  dev.off()
}



pathways.all <- unique(df_old.net$pathway_name)

for (pw in pathways.all) {
  safe_pw <- gsub("[/\\:*?\"<>| ]", "_", pw)
  
  svg_file <- paste0("Cellchat/Lung/split_pathway_plot_Old/aggregate_chord_", safe_pw, ".pdf")
  # tiff(svg_file, width = 18, height = 12, units = "in", res = 300, compression = "lzw")
  pdf(svg_file, width = 18, height = 12)
  # svg(svg_file, width = 18, height = 12)
  circlize::circos.clear()
  
  tryCatch({
    netVisual_aggregate(
      object = cellchatOld,
      signaling = pw,
      layout = "chord",
      color.use = group.colors,
      show.legend = TRUE,
      signaling.name = pw
    )
  }, error = function(e) {
    message("Failed to draw ", pw, ": ", e$message)
  })
  
  dev.off()
}


save(cellchatYoung,file = 'Cellchat/Lung/cellchatYoung.RData')
save(cellchatOld,file = 'Cellchat/Lung/cellchatOld.RData')


netVisual_bubble(cellchatOld, remove.isolate = FALSE)

##Specific
#Skin
load('Cellchat/Skin/cellchatYoung.RData')
load('Cellchat/Skin/cellchatOld.RData')


celltypes <- levels(cellchatYoung@idents)
# celltypes <- rownames(mat)

# group.colors <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(length(celltypes)), celltypes)
group.colors <- c(
  "Interfollicular epidermis cell 1" = "#00BFFF",   # Deep Sky Blue
  "Interfollicular epidermis cell 2" = "#0096C7",   # Strong cyan-blue
  "Interfollicular epidermis cell 3" = "#0077B6",   # Rich blue
  "Merkel cell 1" = "#0081A7",                      # Vibrant teal
  "Fibroblast 1" = "#00C49A",                       # Bright turquoise
  "Inner root sheath cell 1" = "#72DD78",           # Saturated light green
  "Early granular keratinocyte 1" = "#32CD32",      # Lime Green
  "Interfollicular epidermis cell 4" = "#228B22",   # Forest Green
  "Basal keratinocyte 1" = "#BFAF00",               # Strong olive/yellow-green
  "Interfollicular epidermis cell 5" = "#FF6B6B",   # Watermelon red
  "Muscle cell 1" = "#FF3B30",                      # Apple Red
  "Interfollicular keratynocite 1" = "#FF0000",     # Pure red
  "Hair follicle stem cell 1" = "#FF4500",          # Orange Red
  "Adipocyte 1" = "#FF7F00",                        # Bright Orange
  "Lagerhans cell 1" = "#FFAA00",                   # Vivid Yellow-Orange
  "Outer bulge cell 1" = "#FF8C00",                 # Dark Orange
  "Fibroblast 2" = "#FFA500",                      # Orange
  "Sebocyte 1" = "#E67300",                         # Burnt orange
  "Granular layer keratinocyte 1" = "#DA70D6",      # Orchid
  "Endothelial/Smooth muscle cell 1" = "#A020F0",   # Purple
  "Dermal papilla cells 1" = "#8A2BE2",             # Blue Violet
  "Cycling basal cell 1" = "#9932CC",               # Dark Orchid
  "T cell 1" = "#FFD700",                           # Gold
  "Fibroblast 3" = "#FFEA00",                       # Bright yellow
  "Macrophage 1" = "#FFB000",                       # Golden yellow
  "Melanocytes 1" = "#A0522D"                       # Saddle brown
)

plot.new()

legend(
  "center",  
  legend = names(group.colors),
  col = group.colors,
  pch = 19,
  cex = 0.9,
  pt.cex = 1,
  x.intersp = 0.6,
  y.intersp = 0.9,
  bty = "n",
  title = "Cell types",
  ncol = 2
)

groupSize <- as.numeric(table(cellchatYoung@idents))


par(mfcol = c(5, 2), mar = c(0.5, 1, 1, 1), oma = c(0, 0, 0, 0), xpd = TRUE)
target_cells <- c(
  "Interfollicular epidermis cell 1",
  "Interfollicular epidermis cell 2",
  "Interfollicular epidermis cell 3",
  "Interfollicular epidermis cell 4",
  "Interfollicular epidermis cell 5"
)





mat_young <- cellchatYoung@net$weight
for (celltype in target_cells) {
  i <- which(rownames(mat_young) == celltype)
  mat2 <- matrix(0, nrow = nrow(mat_young), ncol = ncol(mat_young), dimnames = dimnames(mat_young))
  mat2[i, ] <- mat_young[i, ]
  
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat_young),
    edge.width.max = 10,
    title.name = paste0(celltype, " (Young)"),
    color.use = group.colors,
    vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)),
    arrow.size = 0.1,
    alpha.edge = 10
  )
}

mat_old <- cellchatOld@net$weight
for (celltype in target_cells) {
  i <- which(rownames(mat_old) == celltype)
  mat2 <- matrix(0, nrow = nrow(mat_old), ncol = ncol(mat_old), dimnames = dimnames(mat_old))
  mat2[i, ] <- mat_old[i, ]
  
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat_old),
    edge.width.max = 10,
    title.name = paste0(celltype, " (Old)"),
    color.use = group.colors,
    vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)),
    arrow.size = 0.1,
    alpha.edge = 10
  )
}



par(mfcol = c(3, 2), mar = c(0.5, 1, 1, 1), oma = c(0, 0, 0, 0), xpd = TRUE)
target_cells <- c(
  "Fibroblast 1",
  "Fibroblast 2",
  "Fibroblast 3"
)

mat_young <- cellchatYoung@net$weight
for (celltype in target_cells) {
  i <- which(rownames(mat_young) == celltype)
  mat2 <- matrix(0, nrow = nrow(mat_young), ncol = ncol(mat_young), dimnames = dimnames(mat_young))
  mat2[i, ] <- mat_young[i, ]
  
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat_young),
    edge.width.max = 5,
    title.name = paste0(celltype, " (Young)"),
    color.use = group.colors,
    vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)),
    arrow.size = 0.1,
    alpha.edge = 10
  )
}

mat_old <- cellchatOld@net$weight
for (celltype in target_cells) {
  i <- which(rownames(mat_old) == celltype)
  mat2 <- matrix(0, nrow = nrow(mat_old), ncol = ncol(mat_old), dimnames = dimnames(mat_old))
  mat2[i, ] <- mat_old[i, ]
  
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat_old),
    edge.width.max = 5,
    title.name = paste0(celltype, " (Old)"),
    color.use = group.colors,
    vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)),
    arrow.size = 0.1,
    alpha.edge = 10
  )
}




par(mfcol = c(1, 2), mar = c(0.5, 1, 1, 1), oma = c(0, 0, 0, 0), xpd = TRUE)
target_cells <- c(
  "Interfollicular epidermis cell 1",
  "Interfollicular epidermis cell 2",
  "Interfollicular epidermis cell 3",
  "Interfollicular epidermis cell 4",
  "Interfollicular epidermis cell 5"
)
mat_young <- cellchatYoung@net$weight
i <- which(rownames(mat_young) %in%target_cells )
mat2 <- matrix(0, nrow = nrow(mat_young), ncol = ncol(mat_young), dimnames = dimnames(mat_young))
mat2[i, ] <- mat_young[i, ]

netVisual_circle(
  mat2,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  edge.weight.max = max(mat_young),
  edge.width.max = 15,
  title.name = paste0("IFE", " (Young)"),
  color.use = group.colors,
  # vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)),
  arrow.size = 0.3,
  alpha.edge = 100
)
mat_old <- cellchatOld@net$weight
i <- which(rownames(mat_old) %in%target_cells)
mat2 <- matrix(0, nrow = nrow(mat_old), ncol = ncol(mat_old), dimnames = dimnames(mat_old))
mat2[i, ] <- mat_old[i, ]

netVisual_circle(
  mat2,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  edge.weight.max = max(mat_old),
  edge.width.max = 15,
  title.name = paste0('IFE', " (Old)"),
  color.use = group.colors,
  # vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)),
  arrow.size = 0.3,
  alpha.edge = 100
)


par(mfcol = c(1, 2), mar = c(0.5, 1, 1, 1), oma = c(0, 0, 0, 0), xpd = TRUE)
target_cells <- target_cells <- c(
  "Fibroblast 1",
  "Fibroblast 2",
  "Fibroblast 3"
)
mat_young <- cellchatYoung@net$weight
i <- which(rownames(mat_young) %in%target_cells )
mat2 <- matrix(0, nrow = nrow(mat_young), ncol = ncol(mat_young), dimnames = dimnames(mat_young))
mat2[i, ] <- mat_young[i, ]

netVisual_circle(
  mat2,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  edge.weight.max = max(mat_young),
  edge.width.max = 15,
  title.name = paste0("FIB", " (Young)"),
  color.use = group.colors,
  # vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)),
  arrow.size = 0.3,
  alpha.edge = 100
)
mat_old <- cellchatOld@net$weight
i <- which(rownames(mat_old) %in%target_cells)
mat2 <- matrix(0, nrow = nrow(mat_old), ncol = ncol(mat_old), dimnames = dimnames(mat_old))
mat2[i, ] <- mat_old[i, ]

netVisual_circle(
  mat2,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  edge.weight.max = max(mat_old),
  edge.width.max = 15,
  title.name = paste0('FIB', " (Old)"),
  color.use = group.colors,
  # vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)),
  arrow.size = 0.3,
  alpha.edge = 100
)







##Specific
#Lung
load('Cellchat/Lung/cellchatOld.RData')
load('Cellchat/Lung/cellchatYoung.RData')


celltypes <- levels(cellchatYoung@idents)
groupSize <- as.numeric(table(cellchatYoung@idents))


# group.colors <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(length(celltypes)), celltypes)
group.colors <- c(
  "Alveolar type 2 cell 1"             = "#00BFFF",  # DeepSkyBlue
  "Endothelial cell 1"                 = "#1E90FF",  # DodgerBlue
  "Club cell 1"                        = "#0077B6",  # Rich Blue
  "Macrophage 1"                       = "#FF3030",  # Bright Red
  "Alveolar type 2 cell 2"             = "#00CED1",  # Dark Turquoise
  "Alveolar type 1 cell 1"             = "#00B2EE",  # Vivid Sky Blue
  "Endothelial cell 2"                 = "#4682B4",  # Steel Blue
  "Endothelial cell 3"                 = "#5B5EA6",  # Indigo Blue
  "Fibroblast 1"                       = "#3CB371",  # Medium Sea Green
  "Club cell 2"                        = "#1CA9C9",  # Cyan Process
  "Endothelial cell 4"                 = "#0066CC",  # Strong Blue
  "Fibroblast 2"                       = "#2E8B57",  # Sea Green
  "T cell 1"                           = "#FF4500",  # Orange Red
  "Monocyte 1"                         = "#FF6347",  # Tomato
  "Red blood cell 1"                   = "#DC143C",  # Crimson
  "Mesenchymal cell 1"                 = "#228B22",  # Forest Green
  "B cell 1"                           = "#FF1493",  # Deep Pink
  "Endothelial cell 5"                = "#4169E1",  # Royal Blue
  "Monocyte 2"                         = "#FF6A6A",  # Light Coral Red
  "Dendritic cell 1"                   = "#FF0000",  # Red
  "Smooth muscle cell 1"              = "#A52A2A",  # Brown
  "Alveolar type 1 cell 2"            = "#20B2AA",  # LightSeaGreen
  "Ciliated cell 1"                   = "#FF8C00",  # Dark Orange
  "Goblet cell 1"                     = "#FFA500",  # Orange
  "Lipofibroblast 1"                  = "#32CD32",  # LimeGreen
  "Smooth muscle cell 2"              = "#8B4513",  # SaddleBrown
  "Krt4/Krt13 Epithelial cell state 1" = "#800080", # Purple
  "Endothelial cell 6"                = "#6495ED"   # Cornflower Blue
)
plot.new()

legend(
  "center",  
  legend = names(group.colors),
  col = group.colors,
  pch = 19,
  cex = 0.9,
  pt.cex = 1,
  x.intersp = 0.6,
  y.intersp = 0.9,
  bty = "n",
  title = "Cell types",
  ncol = 2
)


par(mfcol = c(6, 2), mar = c(0.5, 1, 1, 1), oma = c(0, 0, 0, 0), xpd = TRUE)
target_cells <- c(
  "Endothelial cell 1" ,
  "Endothelial cell 2" ,
  "Endothelial cell 3" ,
  "Endothelial cell 4" ,
  "Endothelial cell 5" ,
  "Endothelial cell 6" 
)

mat_young <- cellchatYoung@net$weight
for (celltype in target_cells) {
  i <- which(rownames(mat_young) == celltype)
  mat2 <- matrix(0, nrow = nrow(mat_young), ncol = ncol(mat_young), dimnames = dimnames(mat_young))
  mat2[i, ] <- mat_young[i, ]
  
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat_young),
    edge.width.max = 10,
    title.name = paste0(celltype, " (Young)"),
    color.use = group.colors,
    vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)),
    arrow.size = 0.1,
    alpha.edge = 10
  )
}

mat_old <- cellchatOld@net$weight
for (celltype in target_cells) {
  i <- which(rownames(mat_old) == celltype)
  mat2 <- matrix(0, nrow = nrow(mat_old), ncol = ncol(mat_old), dimnames = dimnames(mat_old))
  mat2[i, ] <- mat_old[i, ]
  
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat_old),
    edge.width.max = 10,
    title.name = paste0(celltype, " (Old)"),
    color.use = group.colors,
    vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)),
    arrow.size = 0.1,
    alpha.edge = 10
  )
}



par(mfcol = c(1, 2), mar = c(0.5, 1, 1, 1), oma = c(0, 0, 0, 0), xpd = TRUE)
target_cells <- c(
  "T cell 1" 
)

mat_young <- cellchatYoung@net$weight
for (celltype in target_cells) {
  i <- which(rownames(mat_young) == celltype)
  mat2 <- matrix(0, nrow = nrow(mat_young), ncol = ncol(mat_young), dimnames = dimnames(mat_young))
  mat2[i, ] <- mat_young[i, ]
  
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat_young),
    edge.width.max = 20,
    title.name = paste0(celltype, " (Young)"),
    color.use = group.colors,
    vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)),
    arrow.size = 0.1,
    alpha.edge = 10
  )
}

mat_old <- cellchatOld@net$weight
for (celltype in target_cells) {
  i <- which(rownames(mat_old) == celltype)
  mat2 <- matrix(0, nrow = nrow(mat_old), ncol = ncol(mat_old), dimnames = dimnames(mat_old))
  mat2[i, ] <- mat_old[i, ]
  
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat_old),
    edge.width.max = 20,
    title.name = paste0(celltype, " (Old)"),
    color.use = group.colors,
    vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)),
    arrow.size = 0.1,
    alpha.edge = 10
  )
}




par(mfcol = c(1, 2), mar = c(0.5, 1, 1, 1), oma = c(0, 0, 0, 0), xpd = TRUE)
target_cells <- target_cells <- c(
  "Endothelial cell 1" ,
  "Endothelial cell 2" ,
  "Endothelial cell 3" ,
  "Endothelial cell 4" ,
  "Endothelial cell 5" ,
  "Endothelial cell 6" 
)
mat_young <- cellchatYoung@net$weight
i <- which(rownames(mat_young) %in%target_cells )
mat2 <- matrix(0, nrow = nrow(mat_young), ncol = ncol(mat_young), dimnames = dimnames(mat_young))
mat2[i, ] <- mat_young[i, ]

netVisual_circle(
  mat2,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  edge.weight.max = max(mat_young),
  edge.width.max = 15,
  title.name = paste0("Endo", " (Young)"),
  color.use = group.colors,
  # vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)),
  arrow.size = 0.3,
  alpha.edge = 100
)
mat_old <- cellchatOld@net$weight
i <- which(rownames(mat_old) %in%target_cells)
mat2 <- matrix(0, nrow = nrow(mat_old), ncol = ncol(mat_old), dimnames = dimnames(mat_old))
mat2[i, ] <- mat_old[i, ]

netVisual_circle(
  mat2,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  edge.weight.max = max(mat_old),
  edge.width.max = 15,
  title.name = paste0('Endo', " (Old)"),
  color.use = group.colors,
  # vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)),
  arrow.size = 0.3,
  alpha.edge = 100
)


par(mfcol = c(1, 2), mar = c(0.5, 1, 1, 1), oma = c(0, 0, 0, 0), xpd = TRUE)
target_cells <- target_cells <- c(
  'T cell 1'
)
mat_young <- cellchatYoung@net$weight
i <- which(rownames(mat_young) %in%target_cells )
mat2 <- matrix(0, nrow = nrow(mat_young), ncol = ncol(mat_young), dimnames = dimnames(mat_young))
mat2[i, ] <- mat_young[i, ]

netVisual_circle(
  mat2,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  edge.weight.max = max(mat_young),
  edge.width.max = 50,
  title.name = paste0("T Cell 1", " (Young)"),
  color.use = group.colors,
  # vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)),
  arrow.size = 0.3,
  alpha.edge = 100
)
mat_old <- cellchatOld@net$weight
i <- which(rownames(mat_old) %in%target_cells)
mat2 <- matrix(0, nrow = nrow(mat_old), ncol = ncol(mat_old), dimnames = dimnames(mat_old))
mat2[i, ] <- mat_old[i, ]

netVisual_circle(
  mat2,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  edge.weight.max = max(mat_old),
  edge.width.max = 50,
  title.name = paste0('T cell 1', " (Old)"),
  color.use = group.colors,
  # vertex.label.color = rep(rgb(1, 1, 1, alpha = 0), length(group.colors)),
  arrow.size = 0.3,
  alpha.edge = 100
)

###SenNet Plot

cellchat<-readRDS('amanda_Sennet/0829_cellchat_mLung_old_p21.rds')
# load('Cellchat/Skin/cellchatOld.RData')

object.list <- list(Old = cellchatOld, Young = cellchatYoung)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat<- liftCellChat(cellchat, group.new = levels(cellchat@idents$joint))

# rm(cellchatOld)
# rm(cellchatYoung)

save(object.list, file = "~/Desktop/Hear/cellCHAT/cellchat_object.list_mouse_HEART_Old_CTL.RData")
save(cellchat, file = "~/Desktop/Hear/cellCHAT/cellchat_merged_mouse_HEART_Old_CTL.RData")


gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
library(reticulate)
py_install("umap-learn")
py_module_available("umap")  
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

rankSimilarity(cellchat, type = "functional")
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
gg1 + gg2


df.rank <- rankNet(cellchat, mode = "comparison", measure = "weight", stacked = FALSE, do.stat = FALSE)$data
library(tidyr)
df.wide <- df.rank %>%
  select(name,contribution,group)%>%
  tidyr::pivot_wider(names_from = group, values_from = contribution)
# library(ggplot2)
#   
# ggplot(df.wide, aes(x = Young, y = Old)) +
#   geom_point(color = "steelblue", size = 3, alpha = 0.7) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
#   scale_x_log10() +
#   scale_y_log10() +
#   theme_classic(base_size = 14) +
#   labs(x = "Pathway strength (Young, log10)",
#        y = "Pathway strength (Old, log10)",
#        title = "Pathway Communication Strength: Young vs Old (log10 scale)") +
#   ggrepel::geom_text_repel(
#     data = df.wide %>%
#       mutate(diff = abs(Old - Young)) %>%
#       arrange(desc(diff)) %>%
#       head(100),
#     aes(label = name),
#     size = 3
#   )

library(dplyr)
library(ggplot2)
library(ggforce)

pseudo <- 1e-5
df.wide <- df.wide %>%
  mutate(
    Old = Old + pseudo,
    Young = Young + pseudo,
    total = Old + Young,
    pct_old = Old / (Old + Young)
  )

df.pies <- df.wide %>%
  mutate(x = Young, y = Old) %>%
  rowwise() %>%
  do({
    data.frame(
      pathway = .$name,
      x = .$x,
      y = .$y,
      start = c(0, 2*pi*.$pct_old),
      end   = c(2*pi*.$pct_old, 2*pi),
      fill  = c("Old", "Young")
    )
  })

p<-ggplot(df.pies) +
  geom_arc_bar(aes(x0 = x, y0 = y, r0 = 0, r = 0.03, start = start, end = end, fill = fill),
               color = "black") +
  scale_fill_manual(values = c("Old" = "#ee918b", "Young" = "#b3cde6")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_log10() + scale_y_log10() +
  coord_equal() +
  theme_classic() +
  labs(x = "Young strength (log10+pseudo)", y = "Old strength (log10+pseudo)",
       title = "Pathway strength with pie-chart markers")+
  ggrepel::geom_text_repel(
    data = df.pies %>% distinct(pathway, x, y),
    aes(x = x, y = y, label = pathway),
    size = 3,
    color = "black",
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = 50
  )
ggsave("Cellchat/Skin/pathway_scatter_pie_labeled.tiff", plot = p, width = 10, height = 8, dpi = 300, compression = "lzw")
ggsave("Cellchat/Skin/pathway_scatter_pie_labeled.pdf",  plot = p, width = 10, height = 8)   
ggsave("Cellchat/Skin/pathway_scatter_pie_labeled.svg",  plot = p, width = 10, height = 8)
