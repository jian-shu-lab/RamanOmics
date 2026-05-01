library(clusterProfiler)
library(org.Mm.eg.db)
library(EnhancedVolcano)
library(ggplot2)

# Define the directory containing the DEG CSV files
deg_folder <- '~/Desktop/CellType_Comparisons/'

# List all CSV files in the directory
deg_files <- list.files(deg_folder, pattern = "*.txt", full.names = TRUE)

# Create a directory for saving GO results if it doesn't exist
go_output_dir <- '~/Desktop/GO_mLUNG/'
dir.create(go_output_dir, showWarnings = FALSE)

# Loop through each DEG file
for (deg_file in deg_files) {
  
  # Read the CSV file
  deg_data <- read.csv(deg_file, row.names = 1)
  
  # Extract the cell type name from the file name
  celltype <- gsub(".*CellType_Comparisons/(.*)_Old_vs_Young.txt", "\\1", deg_file)
  
  ### GO Analysis for Upregulated Genes ###
  filtered_gene_list_up <- rownames(deg_data[deg_data$avg_log2FC >= 1 & deg_data$p_val_adj < 0.05, ])
  
  # Perform GO analysis only if there are more than 8 upregulated genes
  if (length(filtered_gene_list_up) > 0) {
    
    # Convert gene symbols to Entrez IDs
    genelist_up <- bitr(filtered_gene_list_up, fromType = "SYMBOL", 
                        toType = c("ENTREZID", "GENENAME"), 
                        OrgDb = org.Mm.eg.db)
    
    if (!is.null(genelist_up)) {
      # Perform GO enrichment analysis
      ego_up <- enrichGO(gene = genelist_up$ENTREZID,
                         OrgDb = org.Mm.eg.db,
                         ont = 'BP',
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH")
      
      # Create and save dotplot for upregulated genes
      if (nrow(ego_up@result) > 0) {
        p_up <- dotplot(ego_up, showCategory = 30)
        ggsave(paste0(Sys.Date(), '_UP_', celltype, '_Old_vs_Young_go-all.png'), 
               plot = p_up, width = 8, height = 10, units = 'in', dpi = 320)
        
        # Save GO results to file
        write.csv(ego_up@result, file = paste0(go_output_dir, 'UP_', celltype, '_Old_vs_Young_go-all.txt'))
        print(paste("GO analysis for upregulated genes in", celltype, "saved."))
      }
    }
  } else {
    print(paste("Skipping upregulated GO analysis for", celltype, ": fewer than 9 genes"))
  }
  
  ### GO Analysis for Downregulated Genes ###
  filtered_gene_list_down <- rownames(deg_data[deg_data$avg_log2FC <= -1 & deg_data$p_val_adj < 0.05, ])
  
  # Perform GO analysis only if there are more than 8 downregulated genes
  if (length(filtered_gene_list_down) > 0) {
    
    # Convert gene symbols to Entrez IDs
    genelist_down <- bitr(filtered_gene_list_down, fromType = "SYMBOL", 
                          toType = c("ENTREZID", "GENENAME"), 
                          OrgDb = org.Mm.eg.db)
    
    if (!is.null(genelist_down)) {
      # Perform GO enrichment analysis
      ego_down <- enrichGO(gene = genelist_down$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = 'BP',
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH")
      
      # Create and save dotplot for downregulated genes
      if (nrow(ego_down@result) > 0) {
        p_down <- dotplot(ego_down, showCategory = 30)
        ggsave(paste0(Sys.Date(), '_DOWN_', celltype, '_Old_vs_Young_go-all.png'), 
               plot = p_down, width = 8, height = 10, units = 'in', dpi = 320)
        
        # Save GO results to file
        write.csv(ego_down@result, file = paste0(go_output_dir, 'DOWN_', celltype, '_Old_vs_Young_go-all.txt'))
        print(paste("GO analysis for downregulated genes in", celltype, "saved."))
      }
    }
  } else {
    print(paste("Skipping downregulated GO analysis for", celltype, ": fewer than 9 genes"))
  }
}
