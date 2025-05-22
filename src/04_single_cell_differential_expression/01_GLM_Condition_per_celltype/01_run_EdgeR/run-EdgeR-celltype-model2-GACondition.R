

# This script loops through tissues and celltypes running a DE using EdgeR. 

# This script implements the model: design <- model.matrix(~0 + GA_Condition, y$samples).
# The contrast of interest are constructed using the makeContrasts function.
# my.contrast <- makeContrasts(
#   Early_disease = GA_ConditionEarly_PE - GA_ConditionEarly_Control, 
#   Late_disease = GA_ConditionLate_PE - GA_ConditionLate_Control,
#   averageGA = (GA_ConditionEarly_PE + GA_ConditionEarly_Control)/2 - (GA_ConditionLate_PE + GA_ConditionLate_Control)/2,
#   averagePE = (GA_ConditionEarly_PE + GA_ConditionLate_PE)/2 - (GA_ConditionEarly_Control+ GA_ConditionLate_Control)/2,
#   levels=design)


library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(stringr)
library(dplyr)
library(SingleCellExperiment)
library(scater)
library(edgeR)
library(scran)

exp_var <- "GA_Condition" # Explanatory variable: GA_Condition: 4 categories

# directory to save the results
dir <- "/home/ssd/ysanchez/Projects/FMI-all-singlecell-20230308/outputs/EdgeR-single-cell-20240828/EdgeR-single-cell-model2-GA_Condition-20240828-level2/"

# Choose the tissues and explanatory variable
tissues <- list("CAM", "Myometrium", "PBMC", "Placenta")
annotation_level <- "CellTypeManual.l3"

# Read the data
data.1 <- readRDS(file = "~/Projects/FMI-all-singlecell-20230308/outputs/objects/FMI-all-20patients-20240827.rds")
data.1 <- subset(data.1, subset = donor %in% "F2044JJ", invert = TRUE)

#Change the default assay
DefaultAssay(data.1) <- "RNA"

#Remove the variable features
VariableFeatures(data.1, assay = "RNA") <- NULL

for (tissue in tissues) {
  print(paste("Analysis for tissue: ", tissue))
  
  # Read the data again
  data.1 <- readRDS(file = "~/Projects/FMI-all-singlecell-20230308/outputs/objects/FMI-all-20patients-20240827.rds")
  data.1 <- subset(data.1, subset = donor %in% "F2044JJ", invert = TRUE)
  
  #Change the default assay
  DefaultAssay(data.1) <- "RNA"
  
  #Remove the variable features
  VariableFeatures(data.1, assay = "RNA") <- NULL
  
  # Subset by location 
  data.1 <- subset(data.1, subset = Tissue_old %in% tissue)
  
  
  ### Prepare the data
  
  # Create extra columns in the metadata with the information about celltype and tissue and condition
  data.1@meta.data$CellTypeManual.l3.Tissue <- paste(data.1@meta.data$CellTypeManual.l3, data.1@meta.data$Tissue_old,sep= "_")
  
  # Replace symbol "/" by "_" as the paths gets confused with the paths.
  # data.1@meta.data$CellTypeManual.l3.Tissue <- str_replace_all(data.1@meta.data$CellTypeManual.l3.Tissue_old,'/','_')
  # unique(data.1@meta.data$CellTypeManual.l3.Tissue)
  
  
  # LIST 1: remove clusters that come from <2 donor per cluster
  cellstofilterout1 <- data.1@meta.data  %>% group_by(CellTypeManual.l3.Tissue,GA_Condition)  %>% 
    summarize(n_donor = n_distinct(donor))  %>%
    filter(n_donor < 2)
  
  # LIST 2: remove early or late cells that do not have enough cells in Early_PE or Early_Control or Late_PE or Late_Control. Less than 25 per category
  cellstofilterout2 <- data.1@meta.data %>% group_by(CellTypeManual.l3.Tissue,GA_Condition) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    filter(n < 25)
  
  cellstofilterout <- c(cellstofilterout1$CellTypeManual.l3.Tissue, cellstofilterout2$CellTypeManual.l3.Tissue)
  
  cellstofilterout <- unique(cellstofilterout)
  
  print(paste("Filtering n number of cells: ", length(cellstofilterout)))
  print(cellstofilterout)  
  
  
  # Filter out these cells
  Idents(data.1) = data.1$CellTypeManual.l3.Tissue
  
  # follow this sintaxt: subset(x = pbmc, idents = c("CD4 T cells", "CD8 T cells"), invert = TRUE)
  data.1 <- subset(x = data.1, idents = cellstofilterout, invert = TRUE)
  
  
  # Some dataframes for storing the information per celltype
  DE.Early_disease = data.frame()
  DE.Late_disease = data.frame()
  DE.AverageGA = data.frame()
  DE.AveragePE = data.frame()
  
  DE.Early_disease_filtered = data.frame()
  DE.Late_disease_filtered = data.frame()
  DE.AverageGA_filtered = data.frame()
  DE.AveragePE_filtered = data.frame()
  
  # Make a loop per cell type
  cell_types = unique(data.1$CellTypeManual.l3.Tissue)
  for (cell_type in cell_types) {
    
    print("Subsetting celltype:")
    print(cell_type)
    print("")
    
    # Subset the Seurat object to the relevant cell type column
    Idents(data.1) = data.1$CellTypeManual.l3.Tissue
    sub = data.1 %>% subset(idents = cell_type)
    
    # here set up a function per Tissue at the moment subset
    # Transform the single-cell object to 
    data.1.sc <- as.SingleCellExperiment(sub)
    
    # Agregate across cells of the same donor and same celltype, create a pseudobulk per donor
    summed <- aggregateAcrossCells(data.1.sc, id=colData(data.1.sc)[,c("donor")])
    
    # Create a DGEList object 
    group1 <- summed$GA_Condition  
    
    # Group considered to make a the DGElist
    print("Group considered to make a the DGElist:")
    print(group1)
    print("")
    
    y <- DGEList(counts(summed), samples=colData(summed), group=group1, genes=rownames(summed))
    # filter out genes with low counts
    print("Dimensions before subsetting:")
    print(dim(y))
    print("")
    
    # Remove genes that are low expressed. Group represent the variable to compare (here Condition)
    keep <- filterByExpr(y, group=group1, min.count=10, min.total.count=20)
    y <- y[keep, , keep.lib.sizes=FALSE]
    
    print("Dimensions after subsetting:")
    print(dim(y))
    print("")
    
    # Calculate effective library sizes - TMM normalization procedure
    y <- calcNormFactors(y)
    
    design <- model.matrix(~0 + GA_Condition, y$samples)
    
    y <- estimateDisp(y, design)
    
    fit <- glmFit(y, design, robust=TRUE)
    
    # The contrast of interest can be constructed using the makeContrasts function.
    
    my.contrast <- makeContrasts(
      Early_disease = GA_ConditionEarly_PE - GA_ConditionEarly_Control,
      Late_disease = GA_ConditionLate_PE - GA_ConditionLate_Control,
      averageGA = (GA_ConditionEarly_PE + GA_ConditionEarly_Control)/2 - (GA_ConditionLate_PE + GA_ConditionLate_Control)/2,
      averagePE = (GA_ConditionEarly_PE + GA_ConditionLate_PE)/2 - (GA_ConditionEarly_Control+ GA_ConditionLate_Control)/2,
      interactionPE =  (GA_ConditionEarly_PE + GA_ConditionEarly_Control) - (GA_ConditionLate_PE - GA_ConditionLate_Control),
      levels=design)
    
    print("My contrasts:")
    print(my.contrast)
    
    # Fit one model per contrast 
    # lrt <-  glmLRT(fit, contrast=my.contrast[,"group"])
    lrt.Early_disease <- glmLRT(fit, contrast=my.contrast[,"Early_disease"])
    lrt.Late_disease <- glmLRT(fit, contrast=my.contrast[,"Late_disease"])
    lrt.averageGA <- glmLRT(fit, contrast=my.contrast[,"averageGA"])
    lrt.averagePE <- glmLRT(fit, contrast=my.contrast[,"averagePE"])
    
    # The top set of most significant differentially expressed transcripts can be retrieved using topTags. This will be a table to check the fold changes.
    tt.Early_disease <- as.data.frame(topTags(lrt.Early_disease, n = Inf))
    tt.Late_disease <- as.data.frame(topTags(lrt.Late_disease, n = Inf))
    tt.averageGA <- as.data.frame(topTags(lrt.averageGA, n = Inf))
    tt.averagePE <- as.data.frame(topTags(lrt.averagePE, n = Inf))
    
    # Add a column of padjusted values using a bonferroni correction. This table should have the same number of DE after filtering by adjusting "tt.Early_disease$padj < 0.1"
    tt.Early_disease$padj <- p.adjust(tt.Early_disease$PValue, method="BH")
    tt.Late_disease$padj <- p.adjust(tt.Late_disease$PValue, method="BH")
    tt.averageGA$padj <- p.adjust(tt.averageGA$PValue, method="BH")
    tt.averagePE$padj  <- p.adjust(tt.averagePE$PValue, method="BH")
    
   
    #### Save the raw data - one table per contrast
    # this is the fulltable: concatenate all celltypes into one file. One per contrast
    # Early_disease
    tt.Early_disease <- tt.Early_disease %>%  mutate(cell_type =!!cell_type)
    DE.Early_disease %<>% bind_rows(tt.Early_disease)
    write.csv(DE.Early_disease,file=paste(dir, "model2-fulltable/Alldata-TopTags-EdgeR-pseudobulk-allcelltypes-",tissue,"-factor-",exp_var,"-contrast-Early_disease-",annotation_level,"-20240828.csv", sep=""))
    
    print("final nrow of all raw results - Contrast Early_disease:")
    print(nrow(DE.Early_disease))
    
    # Late_disease
    tt.Late_disease <- tt.Late_disease %>%  mutate(cell_type =!!cell_type)
    DE.Late_disease %<>% bind_rows(tt.Late_disease)
    write.csv(DE.Late_disease,file=paste(dir, "model2-fulltable/Alldata-TopTags-EdgeR-pseudobulk-allcelltypes-",tissue,"-factor-",exp_var,"-contrast-Late_disease-",annotation_level,"-20240828.csv", sep=""))
    
    print("final nrow of all raw results - Contrast Late_disease:")
    print(nrow(DE.Late_disease))
    
    # AverageGA
    tt.averageGA <- tt.averageGA %>%  mutate(cell_type =!!cell_type)
    DE.AverageGA %<>% bind_rows(tt.averageGA)
    write.csv(DE.AverageGA,file=paste(dir, "model2-fulltable/Alldata-TopTags-EdgeR-pseudobulk-allcelltypes-",tissue,"-factor-",exp_var,"-contrast-AverageGA-",annotation_level,"-20240828.csv", sep=""))
    
    print("final nrow of all raw results - Contrast AverageGA:")
    print(nrow(DE.AverageGA))
    
    # AveragePE
    tt.averagePE <- tt.averagePE %>%  mutate(cell_type =!!cell_type)
    DE.AveragePE %<>% bind_rows(tt.averagePE)
    write.csv(DE.AveragePE,file=paste(dir, "model2-fulltable/Alldata-TopTags-EdgeR-pseudobulk-allcelltypes-",tissue,"-factor-",exp_var,"-contrast-AveragePE-",annotation_level,"-20240828.csv", sep=""))
    
    print("final nrow of all raw results - Contrast AveragePE:")
    print(nrow(DE.AveragePE))
    
  }
  
  print(paste("All done for tissue: ", tissue))
}
