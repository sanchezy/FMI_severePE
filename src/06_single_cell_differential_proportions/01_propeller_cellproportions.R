# This script run propeller to quantify changes in cell proportions. https://phipsonlab.github.io/propeller-paper-analysis/pbmcJP.html
# Obtained from: https://github.com/Oshlack/speckle
# Note: This analysis removes the donor "FJJ". Tissues: Placenta, Myometrium, PBMC are done with this script. 

library(Seurat)
library(speckle)
library(limma)
library(ggplot2)
library(edgeR)
library(patchwork)
library(cowplot)
library(gridGraphics)
library(dplyr)

# directory where to save the results
dir <- "/Users/ysanchez/Documents/Projects-analysis/FMI-sc-all/outputs/propeller_results_20240914/"

# Open the seurat data
data.1 <- readRDS(file = "/Users/ysanchez/Documents/Projects-analysis/FMI-sc-all/objects/FMI-all-20patients-20240827.rds")
data.1 <- subset(data.1, subset = donor %in% "F2044JJ", invert = TRUE)

# Parameters: Placenta, Myometrium, PBMC. Change as needed.
tissue <- "PBMC"
# Parameters: "CellTypeManual.l2" or "CellTypeManual.l3"
annotation_level <- "CellTypeManual.l2" 

# subset by tissue
data.1 <- subset(data.1, subset = Tissue_old %in% tissue)

# Create a dataframe with the counts, proportions and transformed proportions using the 'asin' transformation
props <- getTransformedProps(clusters = data.1$CellTypeManual.l3, sample = data.1$donor, transform = "asin")

print("Proportions done")

# Gather the sampleinfor for matrix construction
sampleinfo <- data.1@meta.data %>% group_by(donor, GA_Condition) %>% summarise(n = n())
sampleinfo <- as.data.frame(sampleinfo)
sampleinfo$GA_Condition <- as.factor(sampleinfo$GA_Condition)
sampleinfo$donor <- as.factor(sampleinfo$donor)


print("Sample_info done") 
print(sampleinfo)

# Matrix design the same as EdgeR for differential expression
 design <- model.matrix(~0 + sampleinfo$GA_Condition, sampleinfo$donor)


# Change the colnames to match the EdgeR pipeline
 colnames(design) <- c("GA_ConditionEarly_Control","GA_ConditionEarly_PE","GA_ConditionLate_Control","GA_ConditionLate_PE")
 
 print("This is the design matrix for all the tissues")
 print(design)

# The contrast of interest can be constructed using the makeContrasts function.
my.contrast <- makeContrasts(
  Early_disease = GA_ConditionEarly_PE - GA_ConditionEarly_Control,
  Late_disease = GA_ConditionLate_PE - GA_ConditionLate_Control,
  averageGA = (GA_ConditionEarly_PE + GA_ConditionEarly_Control)/2 - (GA_ConditionLate_PE + GA_ConditionLate_Control)/2,
  averagePE = (GA_ConditionEarly_PE + GA_ConditionLate_PE)/2 - (GA_ConditionEarly_Control+ GA_ConditionLate_Control)/2,
  interactionPE =  (GA_ConditionEarly_PE + GA_ConditionEarly_Control) - (GA_ConditionLate_PE - GA_ConditionLate_Control),
  levels=design)

#  Early_disease
propeller_early <- propeller.ttest(props, design = design, contrasts = my.contrast[,"Early_disease"], robust=TRUE,trend=FALSE,sort=TRUE)
write.csv(propeller_early, file = paste(dir, "propeller-tstatistic-contrast-Early_disease-",annotation_level,"-", tissue,"-20240914.csv",sep =""))

#  Late_disease
propeller_late <- propeller.ttest(props, design = design, contrasts = my.contrast[,"Late_disease"], robust=TRUE,trend=FALSE,sort=TRUE)
write.csv(propeller_late, file = paste(dir, "propeller-tstatistic-contrast-Late_disease-",annotation_level,"-", tissue,"-20240914.csv",sep =""))

#  averagePE
propeller_averagePE <- propeller.ttest(props, design = design, contrasts = my.contrast[,"averagePE"], robust=TRUE,trend=FALSE,sort=TRUE)
write.csv(propeller_averagePE, file = paste(dir, "propeller-tstatistic-contrast-averagePE-",annotation_level,"-", tissue,"-20240914.csv",sep =""))

# averageGA
propeller_averageGA <- propeller.ttest(props, design = design, contrasts = my.contrast[,"averageGA"], robust=TRUE,trend=FALSE,sort=TRUE)
write.csv(propeller_averageGA, file = paste(dir, "propeller-tstatistic-contrast-averageGA-",annotation_level,"-", tissue,"-20240914.csv",sep =""))

print(paste("All done for tissue", tissue, sep = " "))
