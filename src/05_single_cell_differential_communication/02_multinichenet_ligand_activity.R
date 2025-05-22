
# This script create the plot from the ligand activity calculated using cell-to-cell differential communication on multinichenetr.

# Author: Yara E. Sanchez Corrales

library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr)
library(Seurat)
library(stringr) 
library(tidyr)
library(ggpubr)

tissue= "Placenta"

# path to the results
path_results <- "~/Projects/FMI-all-singlecell-20230614/outputs/multinichenet_results_20240910/"

# path to save the plots
pathsave = paste0("~/Projects/FMI-all-singlecell-20230614/outputs/multinichenet_results_20240910/",tissue,"/plots_",tissue,"/")

# read the results
multinichenet_output <- readRDS(paste0(path_results,tissue,"/",tissue, "_multinichenet_output.rds"))

###### Comparisons

contrast_tbl = tibble(contrast =
                        c("Early_PE-Early_Control","Late_PE-Late_Control"),
                      group = c("Early_PE","Late_PE"))

###### Ligand activity of the myeloid compartment

IFN <- c("IFNL3","IFNL1","IFNG","IFNB1","IFNA16","IFNA13","IFNA1","IFITM1")
# Filter by cell types
celltypes <- c("CD14_Monocyte","CD16_Monocyte","Macrophage")

# Filter the data frames
filtered_group_prioritization_tbl <- multinichenet_output$prioritization_tables$group_prioritization_tbl[
  multinichenet_output$prioritization_tables$group_prioritization_tbl$receiver %in% celltypes, ]

filtered_sample_prioritization_tbl <- multinichenet_output$prioritization_tables$sample_prioritization_tbl[
  multinichenet_output$prioritization_tables$sample_prioritization_tbl$receiver %in% celltypes, ]

filtered_ligand_activities_target_de_tbl <- multinichenet_output$prioritization_tables$ligand_activities_target_de_tbl[
  multinichenet_output$prioritization_tables$ligand_activities_target_de_tbl$receiver %in% celltypes, ]

# Create a new list preserving the structure of multinichenet_output
filtered_multi <- multinichenet_output
filtered_multi$prioritization_tables$group_prioritization_tbl <- filtered_group_prioritization_tbl
filtered_multi$prioritization_tables$sample_prioritization_tbl <- filtered_sample_prioritization_tbl
filtered_multi$prioritization_tables$ligand_activities_target_de_tbl <- filtered_ligand_activities_target_de_tbl

# Selecting a given list of activity genes
plot_oi = make_ligand_activity_plots(
  filtered_multi$prioritization_tables, 
  ligands_oi = IFN, 
  contrast_tbl,
  widths = NULL)

# Plot only a set of ligands
plot_oi
ggsave(paste(pathsave, tissue ,"-IFN1-ligand-activity-in-Myeloid-receiver-cells_20241127.png",sep=''), plot_oi, width=12,height=3, bg = "white")