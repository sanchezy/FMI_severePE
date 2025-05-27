### Function to create a dataframe from the results of propeller. This function is used to make relevant visualizations

df_results_per_contrast_tissue <- function(tissue) {
  
  # this is for one one contrast, then loop through the rest
  contrast = "Early_disease" 
  resultsprop <- data.frame()
  
  resultsprop <- read.csv(paste(dir, "propeller-tstatistic-contrast-",contrast,"-",annotation_level_l3,"-", tissue,"-20240914.csv",sep =""))
  resultsprop <- resultsprop %>% filter(P.Value < p_thr)
  
  
  # make an extra column with the effect
  resultsprop$effect <- ifelse(resultsprop$PropRatio < 1, "decreased", "increased") 
  # add extra columns of metadata
  resultsprop$Tissue <- tissue
  resultsprop$contrast <- contrast
  
  # Change the name of the cell type 
  colnames(resultsprop)[1] <- "CellTypeManual.l3"
  
  # Remove the columns 2 and 3 called "PropMean.GA_ConditionEarlyControl" and "PropMean.GA_ConditionEarly_PE". This is neccesary to append the results of all the constrast together
  resultsprop <- subset(resultsprop, select = c(-2,-3))
  
  # Add the information on cell groups
  resultsprop <- resultsprop %>% inner_join(cellgroups)
  
  resultsprop$cellgroup = factor(resultsprop$cellgroup, levels = celltype.order.l1, ordered = TRUE)
  
  # Remove cell types that come from one donor
  celltypes_one_donor <- celltypes_one_donor_all %>% filter(Tissue == tissue) %>% filter(ga == ga) 
  # head(celltypes_one_donor)
  
  resultsprop <- resultsprop %>% filter(!CellTypeManual.l3 %in% celltypes_one_donor$CellTypeManual.l3)
  
  # now loop through the rest of the tissues and append the list to results prop
  contrasts <- list("Late_disease" , "averagePE")
  
  for (comparison in contrasts) {
    # print(paste("Analysis for contrast: ", comparison))  
    
    temp <- data.frame()
    temp <- read.csv(paste(dir, "propeller-tstatistic-contrast-",comparison,"-",annotation_level_l3,"-", tissue,"-20240914.csv",sep =""))
    temp <- temp %>% filter(P.Value < p_thr)
    
    # make an extra column with the effect
    temp$effect <- ifelse(temp$PropRatio < 1, "decreased", "increased") 
    # add extra columns of metadata
    temp$Tissue <- tissue
    temp$contrast <- comparison
    
    # Change the name of the columns to rbind it together
    colnames(temp)[1] <- "CellTypeManual.l3"
    
    # remove the cols of "PropContrast". Note that averagePE will have 4 cols, one per contrast 
    if (comparison == "averagePE"){
      temp <- subset(temp, select = c(-2:-5))  
    }
    else {
      temp <- subset(temp, select = c(-2:-3))  
    }
    
    # Add the information on cell groups
    temp <- temp %>% inner_join(cellgroups)
    
    temp$cellgroup = factor(temp$cellgroup, levels = celltype.order.l1, ordered = TRUE)
    
    # Remove cell types that come from one donor
    celltypes_one_donor <- celltypes_one_donor_all %>% filter(Tissue == tissue) %>% filter(ga == ga) 
    temp <- temp %>% filter(!CellTypeManual.l3 %in% celltypes_one_donor$CellTypeManual.l3)
    
    resultsprop <- rbind(temp,resultsprop)
    head(resultsprop)
    
    # print(paste("Done - contrast: ", comparison))
    rm(temp)
  }
  
  # This is the dataframe that is used for plotting
  contrast.order <- c("Early_disease","Late_disease","averagePE")
  resultsprop$contrast = factor(resultsprop$contrast, levels = contrast.order, ordered = TRUE)
  resultsprop$CellTypeManual.l3 = factor(resultsprop$CellTypeManual.l3, levels = celltype_order.l3, ordered = TRUE)
  
  return(resultsprop)
}