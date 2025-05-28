################
#### Function to generate a heatmap for one tissue type, multiple genes, and multiple cell types. 
################


generate_heatmaps_many_t_g_single_c <- function(tissues, genestoplot, celltypestoplot, all_fmi, output_file, title = "Heatmap") {
  
  # Initialize an empty dataframe to store z-scores
  z_scores_all <- data.frame(matrix(ncol = length(genestoplot) + 1, nrow = 0))
  colnames(z_scores_all) <- c(genestoplot, "tissue")
  
  # Iterate through the specified tissues
  for (tissue in tissues) {
    # Subset data for the current tissue
    subtissue_fmi <- subset(all_fmi, Tissue == tissue)
    cells_subtissue_fmi <- subset(subtissue_fmi, CellTypeManual.l3 %in% celltypestoplot)
    
    # Extract expression data for genes of interest
    expr_data <- FetchData(cells_subtissue_fmi, vars = genestoplot)
    expr_data$condition <- cells_subtissue_fmi@meta.data$GA_Condition
    expr_data$celltype <- cells_subtissue_fmi@meta.data$CellTypeManual.l3
    
    # Subset the expression data for the desired cell type
    sub_expr_data <- expr_data[expr_data$celltype %in% celltypestoplot, ]
    if (nrow(sub_expr_data) == 0) next # Skip if no data for the current cell type
    
    # Calculate mean expression per condition
    result <- sub_expr_data %>%
      group_by(condition) %>%
      summarise(across(genestoplot, mean))
    
    # Convert to matrix and calculate z-scores
    data_matrix <- result %>%
      column_to_rownames("condition") %>%
      as.matrix()
    z_scores <- as.data.frame(scale(data_matrix))
    z_scores$tissue <- tissue
    
    # Append the z-scores to the master dataframe
    z_scores_all <- rbind(z_scores_all, z_scores)
  }
  
  # save heatmap
  png(filename = output_file, width = 3000, height = 2000, res = 300)
  pheatmap(z_scores_all[,-ncol(z_scores_all)], 
           main = title,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           legend = TRUE,
           angle_col = 45,
           treeheight_col = 1,
           cellwidth = 30,
           cellheight = 7,
           fontsize_col = 7.5,
           fontsize_row = 6,
           border_color = NA,
           gaps_row = seq(length(genestoplot), nrow(z_scores_all), by = length(genestoplot)),
           gaps_col = c(3))  # Add gaps after column 3
  dev.off()
}


################ 
#### Example Usage
################

library(dplyr)
library(Seurat)
library(pheatmap)
library(gridExtra)
library(tibble)

# Read seurat object. Remove the donor FJJ

all_fmi<- readRDS(file = "/Users/ysanchez/Documents/Projects-analysis/FMI-sc-all/objects/FMI-all-20patients-20240827.rds")
nojj<-unique(all_fmi$sample_id)[!grepl("F2044JJ",unique(all_fmi$sample_id))]
all_fmi <- subset(all_fmi, sample_id %in% nojj)

# Example Usage
tissues <- c("PVBP", "Myometrium", "CAM")  # List of tissues to include
genestoplot <- c( "GSDMA", "RMDN3","ERAP2")  # Genes to plot
celltypestoplot <- "Macrophage"  # Cell type to analyze
output_file <- "/Users/ysanchez/Documents/Projects-analysis/FMI-sc-all/plots/20250401_Macrophages_proportions/example.png"

# Call the function
generate_heatmaps_many_t_g_single_c(tissues, genestoplot, celltypestoplot, all_fmi, output_file, title = "Macrophages example")

