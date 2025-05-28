
################
##### Function CV. Calculate the CV per capture area. seuratobj will have the deconvoluted density as a column. df_metadata will be a data frame with the information of donor, condition, etc.
################

cv_spatial <- function(seuobj, region, column, df_metadata) {
  
  # Initialize a list to store Moran's I results for each slice
  cv_all_slices <- list()
  
  # coefficient of variation
  cv <- function(x) { sd(x) / mean(x) }
  
  # Subset by the region
  seuobj <- subset(seuobj, subset = regions == region)
  # print(unique(seuobj@meta.data$regions))
  
  # Loop through each slice in the Seurat object
  for (slice_name in unique(seuobj@meta.data$Image_name)) {
    cat("Processing slice:", slice_name, "\n")
    
    temp <- seuobj
    
    # subset the object for the relevant image
    temp <- subset(temp, Image_name == slice_name)
    
    
    ad <- temp@meta.data[[column]]
    
    
    # Calculate CV per capture area. This is a vector, each spot will have only one value
    cv_slice  <- cv(ad)
    
    print(cv_slice)
    
    # Store results in the list with slice name as the key
    cv_all_slices [[slice_name]] <- cv_slice
  }
  
  # Combine the results into a data frame for easy comparison
  cv_df <- data.frame(
    Image_name = names(cv_all_slices),
    CV = unlist(cv_all_slices)
  )
  
  # Drop rows with NA values
  cv_df <- cv_df %>% drop_na()
  
  # Add other relevant metadata
  cv_df <- inner_join(cv_df, df_metadata, by = "Image_name")
  
  return(cv_df)
}

################
##### Function to make stats. Takes as an input the df generated above
################

# stats
stats_cv <- function(df) {
  data_early_p <- df %>% filter(GA_Cat == "Early")
  data_late_p <- df %>% filter(GA_Cat == "Late")
  
  
  early <- wilcox.test(data_early_p$CV ~ data_early_p$GA_Condition, data = data_early_p)
  print(early)
  # Wilcoxon Rank-Sum Test (Mann-Whitney U Test)
  # Assuming 'data' is a data frame with columns 'abundance_per_spot'  and 'condition'
  # data_late_pl <- data_late_pl %>% filter(GA_Condition == "Late")
  late <- wilcox.test(data_late_p$CV ~ data_late_p$GA_Condition, data = data_late_p)
  print(late)
}

################
##### Function that takes 'gene_column' as an argument, filters it and returns the plot
################

boxplot_gene <- function(metadata, gene_column,colors_conditions) {
  
  # Use enquo to capture the column name and evaluate it later
  gene_column <- enquo(gene_column)
  
  # Filter rows where the selected 'gene' column is > 0
  filtered_data <- metadata %>%
    filter(!!gene_column > 0)
  
  # filtered_data <- metadata
  # filtered_data <- filtered_data %>%
  #   filter(CellTypeManual.l3 %in% celltype)
  
  # Create the bar plot using ggplot2
  
  barplot <-  ggplot(filtered_data, aes(x=GA_Condition, y=!!gene_column, fill = GA_Condition)) +
    geom_boxplot(color="black") +  
    scale_fill_manual(values = colors_conditions)  + 
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "grey"),
      strip.background =element_rect(fill="white"),
      axis.text.x = element_text(angle = 0, size = 12),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 14),
      strip.text = element_text(size = 12, color = "black"),
      legend.position = "left",
      legend.text = element_text(size = 18),
      legend.title= element_blank()) +
    labs(x= "", y = gene_column)  + NoLegend() +
    scale_x_discrete(labels=c("Early_Control" = "Control", "Early_PE" = "PE",
                              "Late_Control" = "Control", "Late_PE" = "PE")) +
    facet_wrap(~ GA_Cat, scales = "free_x")  # Facet by Condition 
  return(barplot)
}

################
##### Function that takes 'gene_column' as an argument, filters it and returns the plot. Draws a significance bar. Modify as necessary.
################

# # Early_PE, Early_Control, Late_PE, Late_Control
# colors_conditions <- c("#0072B2", "#E69F00", "#56B4E9", "#F0E442")

boxplot_sig <- function(metadata, gene_column,colors_conditions) {
  
  # Use enquo to capture the column name and evaluate it later
  gene_column <- enquo(gene_column)
  
  # Filter rows where the selected 'gene' column is > 0
  filtered_data <- metadata %>%
    filter(!!gene_column > 0)
  
  barplot <-  ggplot(filtered_data, aes(x=GA_Condition, y=!!gene_column, fill = GA_Condition)) +
    geom_boxplot(color="black") +  
    scale_fill_manual(values = colors_conditions)  + 
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "grey"),
      strip.background =element_rect(fill="white"),
      axis.text.x = element_text(angle = 0, size = 12),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 14),
      strip.text = element_text(size = 12, color = "black"),
      legend.position = "left",
      legend.text = element_text(size = 18),
      legend.title= element_blank()) +
    labs(x= "", y = gene_column)  + NoLegend() +
    scale_x_discrete(labels=c("Early_Control" = "Control", "Early_PE" = "PE",
                              "Late_Control" = "Control", "Late_PE" = "PE")) +
    facet_wrap(~ GA_Cat, scales = "free_x") + 
    # Add significance bracket with manual asterisks
    geom_signif(comparisons = list(c("Early_Control", "Early_PE")), 
                annotations = "*",  # Change this to "**" or "*" depending on your p-value
                y_position = max(data$coly) * 1.05,  # Adjust the position
                tip_length = 0.02,  # Length of the bracket tips
                textsize = 6) + # Size of the asterisks +
    
    scale_x_discrete(labels=c("Early_Control" = "Control", "Early_PE" = "PE",
                              "Late_Control" = "Control", "Late_PE" = "PE")) # Facet by Condition 
  return(barplot)
}
