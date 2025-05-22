# This script runs EdgeR to quantify Macrophages differences **across tissues only in PE** in either Early or Late disease.
# Here I performed a paired differential expression analysis across tissues within the PE group, accounting for donor-level effects.
# 
# This answers: What are the genes that different depending on the tissue? 
#   
# The matrix is:  *design <- model.matrix(~0 + Condition.Tissue, y$samples)*

# Author: Yara E.Sanchez Corrales

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
library(ggvenn)

# Call the function to run the results per contrast.
source("get_de_results.R")

dir <- "/home/ssd/ysanchez/Projects/FMI-all-singlecell-20230308/outputs/EdgeR-single-cell-Macrophages-20240408/"

#########  Read the data and subset accordingly
data.1 <- readRDS(file = "~/Projects/FMI-all-singlecell-20230308/outputs/objects/FMI-all-20patients-20240827.rds")
data.1 <- subset(data.1, subset = donor %in% "F2044JJ", invert = TRUE)

#Change the default assay
DefaultAssay(data.1) <- "RNA"

#Remove the variable features
VariableFeatures(data.1, assay = "RNA") <- NULL

# Subset by Macrophages
data.1 <- subset(data.1, subset = CellTypeManual.l3 %in% "Macrophage")

# Subset by Condition - only PE
data.1 <- subset(data.1, subset = Condition %in% "PE")


# Make an extra metadata column
data.1@meta.data$CellTypeManual.l3.Tissue <- paste(data.1@meta.data$CellTypeManual.l3, data.1@meta.data$Tissue,sep= "_")
data.1@meta.data$Condition.Tissue <- paste(data.1@meta.data$GA_Condition, data.1@meta.data$Tissue,sep= "_")
unique(data.1@meta.data$Condition.Tissue)

######### Transform the single-cell object to SingleCellExperiment and make a DGEList object
data.1.sc <- as.SingleCellExperiment(data.1)

# Agregate across cells of the same donor and same celltype, create a pseudobulk per donor
summed <- aggregateAcrossCells(data.1.sc, id=colData(data.1.sc)[,c("donor","Condition.Tissue")])


# Create a DGEList object. Instead of GA.Condition
group1 <- summed$Condition.Tissue

# Group considered to make a the DGElist
print("Group considered to make a the DGElist:")
print(group1)
print("")

y <- DGEList(counts(summed), samples=colData(summed), group=group1, genes=rownames(summed))
# filter out genes with low counts
print("Dimensions before subsetting:")
print(dim(y))
print("")

######### Filter low expressed gene, calculate norm factors and build a design matrix.

# Remove genes that are low expressed. Group represent the variable to compare (here Condition)
keep <- filterByExpr(y, group=group1, min.count=10, min.total.count=20)
y <- y[keep, , keep.lib.sizes=FALSE]

print("Dimensions after subsetting:")
print(dim(y))
print("")

# Calculate effective library sizes - TMM normalization procedure
y <- calcNormFactors(y)

design <- model.matrix(~0 + Condition.Tissue, y$samples)
# design <- model.matrix(~0 + Condition.Tissue + donor, y$samples)

# colnames(design) <- levels(y$samples$group)
#print(design)

#print(colnames(design))

y <- estimateDisp(y, design)

fit <- glmFit(y, design, robust=TRUE)

######### Establish the contrasts 
# The contrast of interest can be constructed using the makeContrasts function.
my.contrast <- makeContrasts(
  PVBP_Myometrium_Early = Condition.TissueEarly_PE_PVBP - Condition.TissueEarly_PE_Myometrium,
  PVBP_CAM_Early = Condition.TissueEarly_PE_PVBP - Condition.TissueEarly_PE_CAM,
  Myometrium_CAM_Early =  Condition.TissueEarly_PE_Myometrium - Condition.TissueEarly_PE_CAM,
  
  PVBP_Myometrium_Late = Condition.TissueLate_PE_PVBP - Condition.TissueLate_PE_Myometrium,
  PVBP_CAM_Late = Condition.TissueLate_PE_PVBP - Condition.TissueLate_PE_CAM,
  Myometrium_CAM_Late =  Condition.TissueLate_PE_Myometrium - Condition.TissueLate_PE_CAM,
  
  levels=design)

######### Call the function called 'get_de_results' per contrasts
# PVBP_Myometrium
PVBP_Myometrium_Early <- get_de_results(fit,  my.contrast, "PVBP_Myometrium_Early",output_dir = dir)
PVBP_Myometrium_Late <- get_de_results(fit,  my.contrast, "PVBP_Myometrium_Late",output_dir = dir)

# PVBP_CAM
PVBP_CAM_Early <- get_de_results(fit,  my.contrast, "PVBP_CAM_Early",output_dir = dir)
PVBP_CAM_Late <- get_de_results(fit,  my.contrast, "PVBP_CAM_Late",output_dir = dir)

# Myometrium_CAM
Myometrium_CAM_Early <- get_de_results(fit,  my.contrast, "Myometrium_CAM_Early",output_dir = dir)
Myometrium_CAM_Late <- get_de_results(fit,  my.contrast, "Myometrium_CAM_Late",output_dir = dir)

print("All done")
