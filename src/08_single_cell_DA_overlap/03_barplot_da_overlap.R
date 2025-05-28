# This scripts plots the barplot of overlap for Figure.

library(dplyr)
library(ggplot2)
library(reshape2) 

dirsave <- "/Users/ysanchez/Documents/Projects-analysis/FMI-sc-all/outputs/dawnn-results-20250227/dawnn_overlap_20250227/"

# results as percentage of overlap
results <- data.frame( 
           tissue = c("PVBP","Myometrium","CAM","PBMC"),
           overlap = c(45, 51,61, 73))

level_order = c("PVBP","Myometrium","CAM","PBMC")

data_melted <- melt(results, id.vars = "tissue")

# Now create the actual plot
p1 <- ggplot(data_melted, aes(x = factor(tissue, level = level_order), y = value, fill = variable)) +
  geom_col() + 
  geom_hline(yintercept = 50, linetype = "dashed") +
   # scale_fill_manual(values = "navy") +
  #geom_text(aes(label = value), position = position_stack(vjust = 1.1),size=5)  + 
  # facet_grid(. ~ GA_Condition, scales = "free",space = "free") + 
  theme(
    axis.text.x = element_text(angle = 0, size = 12),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.title= element_blank()) +
  labs(y = "percentage of DA overlap", x = "") + NoLegend()

ggsave(paste0(dirsave,"barplot_percentage_of_DA_overlap_20250228.png"), p1, width=4,height=4, bg = "white",  units = 'in', dpi = 300)

print(p1)