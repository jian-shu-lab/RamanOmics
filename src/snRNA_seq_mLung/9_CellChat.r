suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(mclust))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(CellChat))
options(stringsAsFactors = FALSE)
set.seed(1)

### Old ###
mLung_combined <- readRDS("/data/amdqiao/5/intermediate/0719_mLung.rds")
mLung_combined <- subset(mLung_combined, subset = orig.ident %in% c("mLung_O1", "mLung_O2", "mLung_O3"))
dim(mLung_combined)

# Add p21 metadata
p21_positive_cells <- subset(mLung_combined, subset = Cdkn1a > 0)
p21_negative_cells <- subset(mLung_combined, subset = Cdkn1a <= 0)
p21_positive_cells[["p21.expression"]] <- "p21_Positive"
p21_negative_cells[["p21.expression"]] <- "p21_Negative"
metadata_positive <- data.frame(p21_positive_cells[["p21.expression"]])
metadata_negative <- data.frame(p21_negative_cells[["p21.expression"]])
metadata <- rbind(metadata_positive, metadata_negative)
mLung_combined <- AddMetaData(mLung_combined, metadata = metadata)

# Combine two conditions into a new column
metadata <- data.frame(mLung_combined[[c("p21.expression", "seurat_clusters")]])
head(metadata)
metadata$combined.condition <- paste(metadata$p21.expression, metadata$seurat_clusters, sep = "_c")
metadata$p21.expression <- NULL
metadata$seurat_clusters <- NULL
head(metadata)
mLung_combined <- AddMetaData(mLung_combined, metadata = metadata)

Idents(mLung_combined) <- "combined.condition"
Idents(mLung_combined) <- factor(
    x = Idents(mLung_combined),
    levels = sort(levels(mLung_combined))
)
table(Idents(mLung_combined))

# Get data from Seurat
mLung_data <- GetAssayData(mLung_combined, assay = "RNA", slot = "data")
mLung_labels <- Idents(mLung_combined)
mLung_meta <- data.frame(group = mLung_labels, row.names = names(mLung_labels))

# Create CellChat object
cellchatmLung <- createCellChat(object = mLung_data, meta = mLung_meta, group.by = "group")

# Add cell information
levels(cellchatmLung@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchatmLung@idents)) # number of cells in each cell group

# Select L-R database
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show database structure
dplyr::glimpse(CellChatDB$interaction)

# Use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # Subset database

# Set the database to use in the object
cellchatmLung@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchatmLung <- subsetData(cellchatmLung)

cellchatmLung <- identifyOverExpressedGenes(cellchatmLung)
cellchatmLung <- identifyOverExpressedInteractions(cellchatmLung)

# Compute the communication probability and infer cellular communication network
cellchatmLung <- computeCommunProb(cellchatmLung, population.size = TRUE)
cellchatmLung <- filterCommunication(cellchatmLung, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchatmLung <- computeCommunProbPathway(cellchatmLung)

# Calculate the aggregated cell-cell communication network
cellchatmLung <- aggregateNet(cellchatmLung)

groupSize <- as.numeric(table(cellchatmLung@idents))
par(mar = c(1,1,1,1), mfrow = c(1,2), xpd = TRUE)

# Number of interactions in mLung data
pdf("./0828_mLung_old_p21_#_interactions.pdf", width = 15, height = 15)
netVisual_circle(cellchatmLung@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F)
dev.off()

# Interaction weights/strength in mLung data
pdf("./0828_mLung_old_p21_interaction_weights.pdf", width = 15, height = 15)
netVisual_circle(cellchatmLung@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F)
dev.off()

saveRDS(cellchatmLung, file = "/data/amdqiao/5/intermediate/0828_cellchat_mLung_old_p21.rds")

mat <- cellchatmLung@net$weight
groupSize <- as.numeric(table(cellchatmLung@idents))

for (i in 1:nrow(mat)) {
	filename <- paste("./interaction_weights_each_cell_type/Old_", rownames(mat)[i], ".pdf", sep = "")

    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]

    pdf(filename, width = 25, height = 25)
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    dev.off()
}

cellchatmLung@netP$pathways


### Young ###
mLung_combined <- readRDS("/data/amdqiao/5/intermediate/0719_mLung.rds")
mLung_combined <- subset(mLung_combined, subset = orig.ident %in% c("mLung_Y1", "mLung_Y2", "mLung_Y3"))
dim(mLung_combined)

# Add p21 metadata
p21_positive_cells <- subset(mLung_combined, subset = Cdkn1a > 0)
p21_negative_cells <- subset(mLung_combined, subset = Cdkn1a <= 0)
p21_positive_cells[["p21.expression"]] <- "p21_Positive"
p21_negative_cells[["p21.expression"]] <- "p21_Negative"
metadata_positive <- data.frame(p21_positive_cells[["p21.expression"]])
metadata_negative <- data.frame(p21_negative_cells[["p21.expression"]])
metadata <- rbind(metadata_positive, metadata_negative)
mLung_combined <- AddMetaData(mLung_combined, metadata = metadata)

# Combine two conditions into a new column
metadata <- data.frame(mLung_combined[[c("p21.expression", "seurat_clusters")]])
head(metadata)
metadata$combined.condition <- paste(metadata$p21.expression, metadata$seurat_clusters, sep = "_c")
metadata$p21.expression <- NULL
metadata$seurat_clusters <- NULL
head(metadata)
mLung_combined <- AddMetaData(mLung_combined, metadata = metadata)

Idents(mLung_combined) <- "combined.condition"
Idents(mLung_combined) <- factor(
    x = Idents(mLung_combined),
    levels = sort(levels(mLung_combined))
)
table(Idents(mLung_combined))

# Get data from Seurat
mLung_data <- GetAssayData(mLung_combined, assay = "RNA", slot = "data")
mLung_labels <- Idents(mLung_combined)
mLung_meta <- data.frame(group = mLung_labels, row.names = names(mLung_labels))

# Create CellChat object
cellchatmLung <- createCellChat(object = mLung_data, meta = mLung_meta, group.by = "group")

# Add cell information
levels(cellchatmLung@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchatmLung@idents)) # number of cells in each cell group

# Select L-R database
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show database structure
dplyr::glimpse(CellChatDB$interaction)

# Use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # Subset database

# Set the database to use in the object
cellchatmLung@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchatmLung <- subsetData(cellchatmLung)

cellchatmLung <- identifyOverExpressedGenes(cellchatmLung)
cellchatmLung <- identifyOverExpressedInteractions(cellchatmLung)

# Compute the communication probability and infer cellular communication network
cellchatmLung <- computeCommunProb(cellchatmLung, population.size = TRUE)
cellchatmLung <- filterCommunication(cellchatmLung, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchatmLung <- computeCommunProbPathway(cellchatmLung)

# Calculate the aggregated cell-cell communication network
cellchatmLung <- aggregateNet(cellchatmLung)

groupSize <- as.numeric(table(cellchatmLung@idents))
par(mar = c(1,1,1,1), mfrow = c(1,2), xpd = TRUE)

# Number of interactions in mLung data
pdf("./0828_mLung_young_p21_#_interactions.pdf", width = 15, height = 15)
netVisual_circle(cellchatmLung@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F)
dev.off()

# Interaction weights/strength in mLung data
pdf("./0828_mLung_young_p21_interaction_weights.pdf", width = 15, height = 15)
netVisual_circle(cellchatmLung@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F)
dev.off()

saveRDS(cellchatmLung, file = "/data/amdqiao/5/intermediate/0828_cellchat_mLung_young_p21.rds")

mat <- cellchatmLung@net$weight
groupSize <- as.numeric(table(cellchatmLung@idents))

for (i in 1:nrow(mat)) {
	filename <- paste("./interaction_weights_each_cell_type/Young_", rownames(mat)[i], ".pdf", sep = "")

    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]

    pdf(filename, width = 25, height = 25)
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    dev.off()
}

cellchatmLung@netP$pathways

sessionInfo()
