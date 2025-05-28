# This script runs dawnn a DA tool for single-cell analysisby Dr George Hall: https://www.biorxiv.org/content/10.1101/2023.05.05.539427v1.full.pdf+html

# This script quantify the OVERLAP between CONDITION and GESTATIONAL AGE
# This script runs DA using (1) *Condition* (label_1 = "Condition", label_2 = "PE") and (2) *GA_Category* (label_1 = "Early",label_2 = "Late") and save the output as a csv file.

# This analysis was done by tissue (loops through tissues). Filtering out patient FJJ.

# load libraries
library(dplyr)
library(Seurat)
library(reticulate)
library(dawnn)
library(ggplot2)
library(harmony)
library(cowplot)
library(ggvenn)

dirsave <- "/Users/ysanchez/Documents/Projects-analysis/FMI-sc-all/outputs/dawnn-results-20250227/dawnn_overlap_20250227/"


########## 
# Read the seurat object:
data.obj <- readRDS(file = "/Users/ysanchez/Documents/Projects-analysis/FMI-sc-all/objects/FMI-all-20patients-20240827.rds")

# Filter out JJ 
data.obj <- subset(data.obj, subset = donor %in% "F2044JJ", invert = TRUE)

data.1 <- data.obj

########## Loop through tissues
# Create a list of unique tissues and loop 
tissues <-unique(data.1@meta.data$Tissue)

print(table(data.1@meta.data$Tissue, data.1@meta.data$Condition))

for (i in 1:length(tissues)) {
  
  # Read it again. Otherwise, the loop does not work.
  data.1 <- data.obj  

  data.1 <- subset(data.1, subset = Tissue %in% tissues[i])
  
  print(paste ("Doing tissue: ", tissues[i]))
  
  # dawnn will run only if the object has more than 1001 cells.
  if (nrow(data.1@meta.data) < 1001) {
    
    print(paste ("Skipping tissue: ", tissues[i], "not enough cells: n_cells ",nrow(data.1@meta.data)))
    
  } else {
    
    ########## CLUSTERING PIPELINE
    # # #This is global-scaling normalization method "LogNormalize" that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result
    data.1 <- NormalizeData(data.1,normalization.method = "LogNormalize", scale.factor = 10000)
    # 
    # #Calculate a subset of features that exhibit high cell-to-cell variation in the dataset
    data.1 <- FindVariableFeatures(data.1, selection.method = "vst", nfeatures = 2000)
    # 
    # #Scale the data so that the mean expression across cells is 0 and the variance across cells is 1
    data.1 <- ScaleData(data.1,vars.regress = "percent.mt")
    
    # PCA calculation
    data.1 <- RunPCA(data.1, verbose = FALSE)
    
    #Methods available. If length(g$Var1) >1 (equal=2), then do the harmony workflow. Otherwise, do the normal workflow because there is data of only one modality.
    g <- table(data.1@meta.data$Method)
    g <-as.data.frame(g)
    
    if (length(g$Var1) > 1) {
      
      reductionname ="harmony"
      #Harmony integration - only when there are cells from two methods: single-cell and single-nuclei. If there is only one, then follow the normal analysis
      data.1 <- data.1 %>%
        RunHarmony(group.by.vars="Method",reduction.save=reductionname)
      
      # calculate umap
      data.1 <- data.1 %>%
        RunUMAP(reduction = reductionname, dims = 1:25) %>%
        FindNeighbors(reduction = reductionname, dims = 1:25) 
      
    } else {
      reductionname ="pca"
      # calculate umap
      data.1 <- data.1 %>%
        RunUMAP(reduction = reductionname, dims = 1:25) %>%
        FindNeighbors(reduction = reductionname, dims = 1:25) 
    }
    print(paste("reduction name:", reductionname))
    
    ########## Run dawnn on Condition: Control vs PE 
    data.1 <- run_dawnn(data.1, label_names = "Condition", label_1 = "Control",
                        label_2 = "PE", reduced_dim =  reductionname, alpha = 0.05)
    
    
    # Save the output in a text file
    output_da_condition <- data.1@meta.data %>% dplyr::select(dawnn_scores, dawnn_lfc, dawnn_p_vals, dawnn_da_verdict,CellTypeManual.l1,CellTypeManual.l2,CellTypeManual.l3, Condition,Tissue)
    
    # Save the output in a text file
    print(paste("Done dawnn - condition for tissue:",tissues[i],"-",reductionname))
    
    print(table(data.1@meta.data$dawnn_da_verdict))

    #Rename the columns 
    data.1@meta.data$dawnn_scores_condition <- data.1@meta.data$dawnn_scores
    data.1@meta.data$dawnn_lfc_condition <- data.1@meta.data$dawnn_lfc
    data.1@meta.data$dawnn_p_vals_condition <- data.1@meta.data$dawnn_p_vals
    data.1@meta.data$dawnn_da_verdict_condition <- data.1@meta.data$dawnn_da_verdict
    
    ########## Run dawnn on Gestational Age: Early vs Late
    # Run dawnn
    data.1 <- run_dawnn(data.1, label_names = "GA_Category", label_1 = "Early",
                        label_2 = "Late", reduced_dim =  reductionname, alpha = 0.05)
    
    # Save the output in a text file
    output_da_ga <- data.1@meta.data  %>% dplyr::select(dawnn_scores, dawnn_lfc, dawnn_p_vals, dawnn_da_verdict,CellTypeManual.l1,CellTypeManual.l2,CellTypeManual.l3, Condition,Tissue)
    
    # Save the output in a text file
    print(paste("done dawnn - ga for tissue:",tissues[i],"-",reductionname))
    print(table(data.1@meta.data$dawnn_da_verdict))
    
    
    #Rename the columns 
    data.1@meta.data$dawnn_scores_ga <- data.1@meta.data$dawnn_scores
    data.1@meta.data$dawnn_lfc_ga <- data.1@meta.data$dawnn_lfc
    data.1@meta.data$dawnn_p_vals_ga <- data.1@meta.data$dawnn_p_vals
    data.1@meta.data$dawnn_da_verdict_ga <- data.1@meta.data$dawnn_da_verdict
    
   
    ########## Calculate OVERLAP in any of the comparisons (condition or ga)
    
    # Create a column that checks the output of DA.
    data.1@meta.data$joinDA <- ifelse(data.1@meta.data$dawnn_da_verdict_condition == 'TRUE' & data.1@meta.data$dawnn_da_verdict_ga == "TRUE",'DA_overlap',  
                                      ifelse(data.1@meta.data$dawnn_da_verdict_condition == 'TRUE', 'DA_Condition',
                                             ifelse(data.1@meta.data$dawnn_da_verdict_ga == 'TRUE', 'DA_GA','No_DA')))
    
    print(paste("This the results of overlap for tissue:", tissues[i], sep =" "))
    table(data.1@meta.data$joinDA)
    
    output_join_da <- data.1@meta.data %>% dplyr::select(dawnn_scores_condition, dawnn_lfc_condition, dawnn_p_vals_condition, dawnn_da_verdict_condition,dawnn_scores_ga, dawnn_lfc_ga, dawnn_p_vals_ga, dawnn_da_verdict_ga,joinDA,CellTypeManual.l1,CellTypeManual.l2,CellTypeManual.l3, Condition,GA_Category,Tissue)
    
    write.csv(output_join_da, paste0(dirsave,"dawnn-overlap-joinDA-",tissues[i], ".csv"), row.names = TRUE)
    
  }
}