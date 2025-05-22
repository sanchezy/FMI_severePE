
# This script create plots and visualisations from the  cell-to-cell differential communication using multinichenetr.
# Author: Yara E. Sanchez Corrales

library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr)
library(Seurat)
library(stringr) 
library(tidyr)
library(ggpubr)

tissue= "Placenta_Myometrium"

# path to the results
path_results <- "~/Projects/FMI-all-singlecell-20230614/outputs/multinichenet_results_20240910/"

# path to save the plots
pathsave = paste0("~/Projects/FMI-all-singlecell-20230614/outputs/multinichenet_results_20240910/",tissue,"/plots_",tissue,"/")

# read the results
multinichenet_output <- readRDS(paste0(path_results,tissue,"/",tissue, "_multinichenet_output.rds"))

# order the donor: Early_Control, Early_PE, Late_Control, Late_PE,
donor.order <- c("FVQ", 
                 "FLR",
                 "FVB", 
                 
                 "FCM",
                 "FAM",
                 "FGS",
                 "FEP",
                 "FVS",
                 
                 "FLJ",
                 "FJD",
                 "FVT",
                 "FYI",
                 "FNS",
                 "FRK",
                 
                 "FFM",
                 "FRH",
                 "FSH",
                 "FND",
                 "FMM")




###### Comparisons

contrast_tbl = tibble(contrast =
                        c("Early_PE-Early_Control","Late_PE-Late_Control"),
                      group = c("Early_PE","Late_PE"))

###### Bubble plot - STB (sender) LEP to the immune
d <- multinichenet_output$prioritization_tables$group_prioritization_tbl %>% filter(receptor %in% "LEPR") %>% filter(sender %in% "STB") %>% filter(prioritization_score > 0.5)
d <- d %>% filter(!receiver %in% "Intermediate.macrophage") 
order <- c("STB", "Endothelial", "LED", "CD56_NK","CD14_Monocyte","CD16_Monocyte", "Macrophage")

group_oi = "Early_PE"



d$receiver <- factor(d$receiver, levels = order)


plot_oi = make_sample_lr_prod_activity_plots(
  multinichenet_output$prioritization_tables,
  d)

plot_oi

ggsave(paste0(pathsave ,tissue,"-LEPR-prioritized_tbl_score_thr_0.5_20250414.png",sep=''), plot_oi, width=17,height=6,  bg = "white",  units = 'in', dpi = 300)

###### Bubble plot - STB (sender) TGFB1 to the immune

# Get the TGFB1 that are within the first 200 prioritization interaction
h <- multinichenet_output$prioritization_tables$group_prioritization_tbl
h <- h %>% filter(sender %in% "STB")  %>% filter(ligand %in% "TGFB1") %>% filter(prioritization_score > 0.8) %>%  filter(receiver %in% cellstypes)

pplot_oi = make_sample_lr_prod_activity_plots(
  multinichenet_output$prioritization_tables,
  h)

plot_oi
ggsave(paste0(pathsave ,tissue,"STB-TGFB1-full-threshold-score-0.8-prioritized_tbl_selected_celltypes_20241127.png",sep=''), plot_oi, width=16,height=5, bg = "white")


