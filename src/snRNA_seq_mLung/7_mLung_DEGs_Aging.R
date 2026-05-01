library(mclust)
library(Matrix)
library(SeuratDisk)
library(SeuratData)
library(tidyverse)
library(scDblFinder)
library(Seurat)
set.seed(1234)


write.csv(mLung@meta.data[["orig.ident"]], '~/Desktop/age.txt')
age <- read.delim("~/Desktop/age.txt")
mLung@meta.data[['Aging']] <- age$x


mLung$Aging_celltype <- paste(mLung$final.annotation, mLung$Aging, sep = "_")
Idents(mLung) <- "Aging_celltype"

# Get unique cell types (without "Old" or "Young")
cell_types <- unique(gsub("_Old|_Young", "", Idents(mLung)))

# Create a folder to save results
output_dir <- "~/Desktop/CellType_Comparisons"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Loop through each cell type
for (cell_type in cell_types) {
  # Define identities for "Old" and "Young"
  old_identity <- paste0(cell_type, "_Old")
  young_identity <- paste0(cell_type, "_Young")
  
  # Check if both identities exist in the dataset
  if (!(old_identity %in% Idents(mLung) & young_identity %in% Idents(mLung))) {
    cat(paste("Skipping", cell_type, "- identities not found\n"))
    next
  }
  
  # Perform differential expression analysis
  markers <- FindMarkers(
    mLung,
    ident.1 = old_identity,
    ident.2 = young_identity,
    only.pos = FALSE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
  
  # Save the results
  output_file <- file.path(output_dir, paste0(cell_type, "_Old_vs_Young.txt"))
  write.csv(markers, file = output_file)
  
  cat(paste("Results saved for", cell_type, "in", output_file, "\n"))
}

