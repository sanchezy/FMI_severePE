# Some functions for plotting the results from EdgeR. 

###########################################################################
# Function to create a dataframe that is used to make plots

df_results_per_all_contrast_tissue <- function(tissue) {
  
  # this is for one one contrast, then loop through the rest
  contrast = "Early_disease" 
  
  # string to be removed
  string1 <- paste("_",tissue,sep="")
  
  # This is the full table of all cells: ~0 + GA_Condition: Early
  model2 <- read.csv(paste0(dir, "model2-fulltable/Alldata-TopTags-EdgeR-pseudobulk-allcelltypes-",tissue,"-factor-",exp_var,"-contrast-",contrast,"-",annotation_level,"-20240828.csv"))
  
  # make an extra column with the effect
  model2$effect <- ifelse(model2$logFC<0, "down_regulated", "up_regulated") 
  
  # Early
  number_of_de = data.frame()
  temp = data.frame()
  number_of_de <- model2 %>% filter(padj < p_thr)  %>% filter(abs(logFC) > logFC_threshold) %>%
    group_by(cell_type) %>%
    dplyr::count() %>%
    # filter(n>1) %>%
    arrange(desc(n))
  
  # Add a column of Tissue
  number_of_de$Tissue <- tissue
  number_of_de$Contrast <- contrast
  
  # Replace the name tissue in the celltype column 
  number_of_de$CellTypeManual.l2 <- number_of_de$cell_type 
  number_of_de$CellTypeManual.l2 <- str_replace(number_of_de$CellTypeManual.l2, string1,"")
  
  # Add the information on cell groups
  number_of_de <- number_of_de %>% inner_join(cellgroups)
  number_of_de$cellgroup = factor(number_of_de$cellgroup, levels = celltype.order.l1, ordered = TRUE)
  
  # Change the sign of the negative
  # number_of_de$n <- ifelse(number_of_de$effect == "down_regulated",-number_of_de$n, number_of_de$n)
  
  # now loop through the rest of the tissues and append the list number_of_de_model2_early
  
  contrasts <- list("Late_disease" , "AveragePE", "AverageGA")
  
  for (comparison in contrasts) {
    print(paste("Analysis for contrast: ", comparison))
    
    # This is the full table of all cells: ~0 + GA_Condition: Early
    model2 <- read.csv(paste0(dir, "model2-fulltable/Alldata-TopTags-EdgeR-pseudobulk-allcelltypes-",tissue,"-factor-",exp_var,"-contrast-",contrast,"-",annotation_level,"-20240828.csv"))
    
    # make an extra column with the effect
    model2$effect <- ifelse(model2$logFC<0, "down_regulated", "up_regulated") 
    
    # head(model2_early)
    
    temp <- model2 %>% filter(padj < p_thr)  %>% filter(abs(logFC) > logFC_threshold) %>%
      group_by(cell_type) %>%
      dplyr::count() %>%
      arrange(desc(n))
    
    # Add a column of Tissue
    # print(tissue)
    temp$Tissue <- tissue
    temp$Contrast <- comparison
    
    # Replace the name tissue in the celltype column 
    temp$CellTypeManual.l2 <- temp$cell_type 
    temp$CellTypeManual.l2 <- str_replace(temp$CellTypeManual.l2, string1,"")
    
    # Add the information on cell groups
    temp <- temp %>% inner_join(cellgroups)
    temp$cellgroup = factor(temp$cellgroup, levels = celltype.order.l1, ordered = TRUE)
    
    # Change the sign of the negative
    # temp$n <- ifelse(temp$effect == "down_regulated",-temp$n, temp$n)
    
    number_of_de <- rbind(temp,number_of_de)
    # rm(temp)
  }
  
  # This is the dataframe that is used for plotting
  contrast.order <- c("Early_disease","Late_disease","AveragePE", "AverageGA")
  number_of_de$Contrast = factor(number_of_de$Contrast, levels = contrast.order, ordered = TRUE)
  number_of_de$CellTypeManual.l2 = factor(number_of_de$CellTypeManual.l2, levels = celltype_order.l2, ordered = TRUE)
  
  return(number_of_de)
}


###########################################################################
# Same as before but remove the constrast averageGA
df_results_per_contrast_tissue <- function(tissue) {
  
  # this is for one one contrast, then loop through the rest
  contrast = "Early_disease" 
  
  # string to be removed
  string1 <- paste("_",tissue,sep="")
  
  # This is the full table of all cells: ~0 + GA_Condition: Early
  model2 <- read.csv(paste0(dir, "model2-fulltable/Alldata-TopTags-EdgeR-pseudobulk-allcelltypes-",tissue,"-factor-",exp_var,"-contrast-",contrast,"-",annotation_level,"-20240828.csv"))
  
  # # make an extra column with the effect
  # model2$effect <- ifelse(model2$logFC<0, "down_regulated", "up_regulated") 
  
  # Early
  number_of_de = data.frame()
  temp = data.frame()
  number_of_de <- model2 %>% filter(padj < p_thr)  %>% filter(abs(logFC) > logFC_threshold) %>%
    group_by(cell_type) %>%
    dplyr::count() %>%
    # filter(n>1) %>%
    arrange(desc(n))
  
  # Add a column of Tissue
  number_of_de$Tissue <- tissue
  number_of_de$Contrast <- contrast
  
  # Replace the name tissue in the celltype column 
  number_of_de$CellTypeManual.l3 <- number_of_de$cell_type 
  number_of_de$CellTypeManual.l3 <- str_replace(number_of_de$CellTypeManual.l3, string1,"")
  
  # Add the information on cell groups
  number_of_de <- number_of_de %>% inner_join(cellgroups)
  number_of_de$cellgroup = factor(number_of_de$cellgroup, levels = celltype.order.l1, ordered = TRUE)
  
  # Change the sign of the negative
  # number_of_de$n <- ifelse(number_of_de$effect == "down_regulated",-number_of_de$n, number_of_de$n)
  
  # now loop through the rest of the tissues and append the list number_of_de_model2_early
  
  contrasts <- list("Late_disease" , "AveragePE")
  
  for (comparison in contrasts) {
    print(paste("Analysis for contrast: ", comparison))
    
    # This is the full table of all cells: ~0 + GA_Condition: Early
    model2 <- read.csv(paste0(dir, "model2-fulltable/Alldata-TopTags-EdgeR-pseudobulk-allcelltypes-",tissue,"-factor-",exp_var,"-contrast-",comparison,"-",annotation_level,"-20240828.csv"))
    
    # make an extra column with the effect
    model2$effect <- ifelse(model2$logFC<0, "down_regulated", "up_regulated") 
    
    # head(model2_early)
    
    temp <- model2 %>% filter(padj < p_thr)  %>% filter(abs(logFC) > logFC_threshold) %>%
      group_by(cell_type) %>%
      dplyr::count() %>%
      arrange(desc(n))
    
    # Add a column of Tissue
    # print(tissue)
    temp$Tissue <- tissue
    temp$Contrast <- comparison
    
    # Replace the name tissue in the celltype column 
    temp$CellTypeManual.l3 <- temp$cell_type 
    temp$CellTypeManual.l3 <- str_replace(temp$CellTypeManual.l3, string1,"")
    
    # Add the information on cell groups
    temp <- temp %>% inner_join(cellgroups)
    temp$cellgroup = factor(temp$cellgroup, levels = celltype.order.l1, ordered = TRUE)
    
    # Change the sign of the negative
    # temp$n <- ifelse(temp$effect == "down_regulated",-temp$n, temp$n)
    
    number_of_de <- rbind(temp,number_of_de)
    # rm(temp)
  }
  
  # This is the dataframe that is used for plotting
  contrast.order <- c("Early_disease","Late_disease","AveragePE")
  number_of_de$Contrast = factor(number_of_de$Contrast, levels = contrast.order, ordered = TRUE)
  number_of_de$CellTypeManual.l3 = factor(number_of_de$CellTypeManual.l3, levels = celltype_order.l3, ordered = TRUE)
  
  return(number_of_de)
}