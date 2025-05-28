################
#### FUNCTION - this function creates a dataframe to plot the results of cell2loc (already in the metadata of an object, grouping by donor or GA_Condtion).
#### df_spots is a daframe containing the information about number of spots per region per patient
################
 
df_c2loc_region_donor <- function(seuratobj, df_spots_patient, donor_patient, spatial_region) {
  
  meta <- seuratobj@meta.data 
  # subset by var and region
  meta <- meta %>% filter(donor %in% donor_patient) 
  meta <- meta %>% filter(regions %in% spatial_region)
  
  # Select the results from the metatada and sum the cols - this gives a vector. The abundance per cell type is collapsed into one vector 
  meta  <- meta[, grepl("^c2loc_", colnames(meta))]  %>% summarise(across(everything(), sum)) 
  
  # The abundance has to be divided by the total number of spots in a given region.
  # append the number of spots
  df_spots_patient <- df_spots_patient %>% filter(regions %in% spatial_region)
  # print("debug: regions")
  # print(unique(df_spots_patient$regions)) 
  
  # make it compatible with ggplot
  meta  <- meta  %>% pivot_longer(cols = everything(), names_to = "CellTypeManual.l2", values_to = "predicted_abundance")
  
  # add the metadata of donor
  meta$donor <- donor_patient
  
  # add the number of spots per region
  meta <- meta %>% inner_join(df_spots_patient)
  
  # calculate the abundance per spot
  meta <-  meta %>% mutate(abundance_per_spot = predicted_abundance/total_spots)
  
  # calculate the proportion per cell type
  meta <-   meta %>% mutate(proportion_region = abundance_per_spot / sum(abundance_per_spot))
  
  # Change the name of CellTypeManual.l2
  meta$CellTypeManual.l2 <- gsub("c2loc_", "", meta$CellTypeManual.l2)
  
  # Make a cumulative abundance
  meta <- meta %>% arrange(desc(proportion_region)) %>% mutate(Cumulative_normalised_abundance = cumsum(proportion_region))
  
  # Add the cols of cell groups
  meta <- meta %>% inner_join(cellgroups)
  
  # Add the info about 
  meta <- meta %>% inner_join(donor.order.condition)
  # order by most contribution to less contribution
  # meta3$CellTypeManual.l2 = factor(meta3$CellTypeManual.l2 , levels = meta3$CellTypeManual.l2, ordered = TRUE)
  
  return(meta)
}

################
#### # Function that takes 'gene_column' as an argument, filters it and returns the plot. Plot per donor.Condition
################

boxplot_gene_donor <- function(metadata, column, region, colors_conditions) {
  
  # Use enquo to capture the column name and evaluate it later
  column <- enquo(column)
  
  filtered_data <- metadata %>%
    filter(regions %in% region)
  
  # Filter rows where the selected 'gene' column is > 0
  filtered_data <- filtered_data %>%
    select(donor.Condition2,column, donor2, donor,GA_Condition) 
  
  # Add the info about 
  # filtered_data <- filtered_data %>% inner_join(donor.order.condition)
  
  # reorder by donor gestational age
  filtered_data$donor.Condition2 <- factor(filtered_data$donor.Condition2, levels =  donor.order.condition$donor2, ordered = TRUE)
  
  # Create the bar plot using ggplot2
  
  plot <-  ggplot(filtered_data, aes(x=donor.Condition2, y=!!column, fill = GA_Condition)) +
    geom_boxplot(color="black") +  
    scale_fill_manual(values = colors_conditions)  + 
    theme(
      axis.text.x = element_text(angle = 90, size = 8),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      strip.text = element_text(size = 12, color = "black"),
      legend.position = "right",
      legend.text = element_text(size = 22),
      legend.title= element_blank()) +
    labs(x= "", y = column)  + NoLegend() +
    # scale_x_discrete(labels=c("Early_Control" = "Control", "Early_PE" = "PE",
    #                           "Late_Control" = "Control", "Late_PE" = "PE")) +
    facet_grid(.~ GA_Condition, scales = "free_x")  # Facet by Condition 
  
  return(plot)
}