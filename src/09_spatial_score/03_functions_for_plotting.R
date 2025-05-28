
################
#### Function to add the hypoxia score into the seurat object
################

add_score <- function(seurat_obj, file, barcode_col, new_metadata_col) {
  # Check the order of the barcodes in the Seurat object and the file
  cellnames <- Cells(seurat_obj)
  idx <- match(cellnames, file[[barcode_col]])  # Match Seurat barcodes to file barcodes
  
  # Reorder the file based on the Seurat object cell order
  file <- file[idx,]
  
  # Check for any unmatched barcodes
  if (any(is.na(idx))) {
    stop("Some barcodes in the Seurat object do not match the barcodes in the file.")
  }
  
  # Add the new metadata column to the Seurat object
  seurat_obj[[new_metadata_col]] <- file[[new_metadata_col]]
  
  # Return the updated Seurat object
  return(seurat_obj)
}

################
#### Function to make a boxplot per gene and donor facet by condition
################

# Function that takes 'gene_column' as an argument, filters it and returns the boxplot
boxplot_gene_donor <- function(metadata, gene_column,colors_conditions,region) {
  
  # Use enquo to capture the column name and evaluate it later
  gene_column <- enquo(gene_column)
  
  # Filter rows where the selected 'gene' column is > 0
  filtered_data <- metadata %>%
    filter(!!gene_column > 0)
  
  filtered_data <- filtered_data %>%
    filter(regions %in% region)
  
  # reorder by donor gestational age
  filtered_data$donor2 <- factor(filtered_data$donor2, levels = donors_to_plot, ordered = TRUE)
  
  # Create the bar plot using ggplot2
  
  barplot <-  ggplot(filtered_data, aes(x=donor2, y=!!gene_column, fill = GA_Condition)) +
    geom_boxplot(color="black") +  
    scale_fill_manual(values = colors_conditions)  + 
    theme(
      axis.text.x = element_text(angle = 90, size = 12),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      strip.text = element_text(size = 12, color = "black"),
      legend.position = "right",
      legend.text = element_text(size = 22),
      legend.title= element_blank()) +
    labs(x= "", y = gene_column) + NoLegend() +
    ggtitle(paste0(gene_column, " in spatial domain: ", region)) +
    # scale_x_discrete(labels=c("Early_Control" = "Control", "Early_PE" = "PE",
    #                           "Late_Control" = "Control", "Late_PE" = "PE")) +
    facet_wrap(.~ GA_Condition, scales = "free_x", ncol = 4)  # Facet by Condition 
  return(barplot)
}

################
#### Function to make a boxplot per gene facet by condition
################

plot_box_plot <- function(object, region,var_to_plot) {
  
  meta <- object@meta.data %>% filter(regions %in% region) 
  
  p1 <- meta %>% ggplot(aes(x=GA_Condition, y={{ var_to_plot }}, fill = GA_Condition)) +
    geom_boxplot(color="black")+  
    scale_fill_manual(values = colors_conditions2) + 
    theme(
      panel.background = element_blank(),
      # Add axis line
      axis.line = element_line(colour = "grey"),
      axis.text.x = element_text(angle = 0, size = 12),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 14),
      strip.text = element_text(size = 12, color = "black"),
      legend.position = "top",
      legend.text = element_text(size = 18),
      legend.title= element_blank()) +
    
    labs(x= "", y = "hypoxia_score")  + NoLegend() +
    scale_x_discrete(labels=c("Early_Control" = "Control", "Early_PE" = "PE",
                              "Late_Control" = "Control", "Late_PE" = "PE")) +
    facet_wrap(~ GA_Cat, scales = "free_x") +
    theme(strip.background =element_rect(fill="white"))
  
  p1  <-  p1 + facet_grid(. ~ GA_Cat, scales = "free", space = "free") + NoLegend() 
  
  return(p1)
}

