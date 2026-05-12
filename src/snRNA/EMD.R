# BiocManager::install("EMDomics")
library(EMDomics)
library(Seurat)
seu_lung<-readRDS('/data/framont/Codes_for_Haochun/1001_mLung.rds')
seu_lung$age_group <- ifelse(grepl("O", seu_lung$orig.ident), "Old", "Young")
cell_types <- unique(seu_lung$final.annotation)


pwd<-'/data/Haochun/SenNet/emd_output_Lung'
files<-list.files(path = pwd, pattern = "\\.txt$", full.names = TRUE)
files<-gsub("^log_|\\.txt$", "", basename(files))
files
files <- gsub("_", " ",files)
setdiff(as.character(cell_types),files)


pwd<-'/data/Haochun/SenNet/emd_output_Lung/'
files<-list.files(path = pwd, pattern = "\\.csv$", full.names = TRUE)
files<-gsub("^emd_|\\.csv$", "", basename(files))
files
files <- gsub("_", " ",files)
setdiff(as.character(cell_types),files)
cell_types<-c("Early granular keratinocyte 1","Melanocytes 1"     )

# for (ct in cell_types) {
#   ct='Ciliated cell 1'
#   subset_ct <- subset(seu_lung, subset = final.annotation == ct)
#   avg_expr <- AverageExpression(subset_ct, group.by = "orig.ident", assays = "RNA", slot = "data")$RNA
#   avg_expr<-as.data.frame(avg_expr)
#   avg_expr <- avg_expr[rowSums(avg_expr != 0) > 1, ]
#   label<-c('Old','Old','Old',"Young",'Young','Young')
#   names(label)<-c("mLung_O1", "mLung_O2", "mLung_O3", "mLung_Y1", "mLung_Y2" ,"mLung_Y3")
#   results<-calculate_emd(avg_expr[1:100,],label,nperm=1000,parallel=FALSE)
#   
# }

log <- function(msg) {
  cat(paste0(Sys.time(), " - ", msg, "\n"))
}

library(Seurat)
library(EMDomics)
library(doParallel)
library(foreach)

ncores <- 2
registerDoParallel(cores = ncores)

output_dir <- "/data/Haochun/SenNet/emd_output_Lung/"
dir.create(output_dir, showWarnings = FALSE)
cell_types<-'Krt4/Krt13 Epithelial cell state 1'
foreach(ct = cell_types, .packages = c("Seurat", "EMDomics")) %dopar% {
  
  if( ct=='Krt4/Krt13 Epithelial cell state 1'){
    ct_fix='Krt4_Krt13 Epithelial cell state 1'
    log_file <- file.path(output_dir, paste0("log_", gsub(" ", "_", ct_fix), ".txt"))
  }else{
    log_file <- file.path(output_dir, paste0("log_", gsub(" ", "_", ct), ".txt"))
  }

  log <- function(msg) {
    cat(paste0(Sys.time(), " - ", msg, "\n"), file = log_file, append = TRUE)
  }
  
  log(paste0("Start cell type: ", ct))
  
  subset_ct <- subset(seu_lung, subset = final.annotation == ct)
  avg_expr <- AverageExpression(subset_ct, group.by = "orig.ident", assays = "RNA", slot = "data")$RNA
  avg_expr <- as.data.frame(avg_expr)
  avg_expr <- avg_expr[rowSums(avg_expr != 0) > 1, ]
  
  label <- c('Old','Old','Old',"Young",'Young','Young')
  names(label) <- c("mLung_O1", "mLung_O2", "mLung_O3", "mLung_Y1", "mLung_Y2", "mLung_Y3")
  
  gene_chunks <- split(rownames(avg_expr), ceiling(seq_along(rownames(avg_expr)) / 100))
  print(ct)
  emd_results_list <- list()
  
  for (i in seq_along(gene_chunks)) {
    print(i)
    genes <- gene_chunks[[i]]
    log(paste0("Processing chunk ", i, " / ", length(gene_chunks)))
    expr_sub<-avg_expr[genes, ]
    label_sub<-label[colnames(avg_expr)]
    res <- calculate_emd(expr_sub, label_sub, nperm = 1000, parallel = FALSE, verbose = FALSE)
    emd_results_list[[i]] <- res$emd
  }
  
  emd_combined <- do.call(rbind, emd_results_list)
  if( ct=='Krt4/Krt13 Epithelial cell state 1'){
    ct_fix<-'Krt4_Krt13 Epithelial cell state 1'
    write.csv(emd_combined, file = file.path(output_dir, paste0("emd_", gsub(" ", "_", ct_fix), ".csv")), row.names = TRUE)
    
  }else{
    write.csv(emd_combined, file = file.path(output_dir, paste0("emd_", gsub(" ", "_", ct), ".csv")), row.names = TRUE)
    
  }

  log("Finished and wrote output.")
}


##Draw plot
csv_dir <- "/data/Haochun/SenNet/emd_output_Skin"

files <- list.files(csv_dir, pattern = "^emd_.*\\.csv$", full.names = TRUE)
files<-files[c(1,3:length(files))]
sig_counts <- c()

for (file in files) {
  df <- read.csv(file, row.names = 1)
  cell_type <- gsub("^emd_|\\.csv$", "", basename(file))
  cell_type <- gsub("_", " ", cell_type)

  sig_counts[cell_type] <- sum(df$q.value < 0.05, na.rm = TRUE)
}


sig_df <- data.frame(
  CellType = names(sig_counts),
  SigGeneCount = as.integer(sig_counts)
)

sig_df <- sig_df[order(sig_df$SigGeneCount, decreasing = FALSE), ]

library(ggplot2)
#Barplot
library(ggbreak)
ggplot(sig_df, aes(x = reorder(CellType, -SigGeneCount), y = SigGeneCount)) +
  geom_bar(stat = "identity", fill = "#a4a0a0", color = "black", size = 0.4) + 
  geom_text(aes(label = SigGeneCount), 
            vjust = -0.5, size = 3) +  
  labs(x = "Cell Type", y = "Genes with q < 0.05",
       title = "Number of Significant Genes per Cell Type") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),                
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.line = element_line(color = "black"),   
    axis.ticks = element_line(color = "black")
  )


library(scales)
squish_transform <- function(break_point = 150, lower_scale = 1, upper_scale = 3) {
  trans_new(
    name = "squish",
    transform = function(y) ifelse(y <= break_point,
                                   y / lower_scale,
                                   break_point / lower_scale + (y - break_point) / upper_scale),
    inverse = function(z) ifelse(z <= break_point / lower_scale,
                                 z * lower_scale,
                                 break_point + (z - break_point / lower_scale) * upper_scale)
  )
}

p<-ggplot(sig_df, aes(x = reorder(CellType, -SigGeneCount), y = SigGeneCount)) +
  geom_bar(stat = "identity", fill = "#a4a0a0", color = "black", size = 0.4) + 
  geom_text(aes(label = SigGeneCount), vjust = -0.5, size = 3) +  
  labs(x = "Cell Type", y = "Genes with q < 0.05",
       title = "Number of Significant Genes per Cell Type") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),                
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.line = element_line(color = "black"),   
    axis.ticks = element_line(color = "black")
  ) +
  scale_y_continuous(trans = squish_transform(break_point = 150,
                                              lower_scale = 2,   
                                              upper_scale = 0.3  
  ))

ggsave("/data/Haochun/SenNet/FIgure_Results/Ke_ordered_updated/EMD_Barplot_Skin_scaled.svg", p, width = 8, height = 6, dpi = 300)

#lollipop plot

p<-ggplot(sig_df, aes(x = reorder(CellType, SigGeneCount,decreasing=TRUE), y = SigGeneCount)) +
  geom_segment(aes(xend = CellType, y = 0, yend = SigGeneCount),
               color = "#fea443", size = 1.2,alpha = 0.8) +  
  geom_text(aes(label = SigGeneCount), 
            vjust = -1, size = 3) +  
  geom_point(size = 4, color = "#45c4b2",alpha=1) +  
  labs(x = "Cell Type", y = "Genes with q < 0.05",
       title = "Number of Significant Genes per Cell Type") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),                
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.line = element_line(color = "black"),   
    axis.ticks = element_line(color = "black")
  )
ggsave("/data/Haochun/SenNet/FIgure_Results/Ke_ordered_updated/EMD_lollipop_Skin.svg", p, width = 8, height = 6, dpi = 300)




