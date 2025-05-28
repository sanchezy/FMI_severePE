# This script makes venn diagrams color-coded by the frequency of the overlap of DA-gestational_age and DA-condition
# This scripts reads the output of the script: source("~/Documents/Projects-analysis/FMI-sc-all/scripts/dawnn_20240312/dawnn_overlap/run-dawnn-tissue-overlap.R")

library(ggVennDiagram)
library(dplyr)
library(ggplot2)

dir <- "/Users/ysanchez/Documents/Projects-analysis/FMI-sc-all/outputs/dawnn-results-20250227/dawnn_overlap_20250227/"
# Read the text file

 # tissues <- c("PVBP","Myometrium","PBMC","CAM")
 tissues <- c("CAM")
# tissues <- c("Placenta")

for (i in 1:length(tissues)) {
results <- read.csv(paste0(dir,"dawnn-overlap-joinDA-",tissues[i], ".csv"))


# Calculate the overlap of dawnn_da_verdict_condition and ga == TRUE
# How many cells are detected as da in both condition and ga?
da_cells_condition = results %>% dplyr::filter(dawnn_da_verdict_condition == TRUE)
da_cells_condition <- da_cells_condition$X

da_cells_ga = results  %>% dplyr::filter(dawnn_da_verdict_ga == TRUE)
da_cells_ga <- da_cells_ga$X

x = list(
  co = da_cells_condition,
  ga = da_cells_ga)


 plot4 <- ggVennDiagram(x, shape_id = "201",label_size=11, label = "percent") + scale_fill_gradient(low="grey90",high = "red")

# plot4
ggsave(paste(dir,"plots-overlap/venn_overlap_20250228/sc-FMI-venn-overlap-",tissues[i],".png",sep=''), plot4, width=6,height=6, bg = "white")


}