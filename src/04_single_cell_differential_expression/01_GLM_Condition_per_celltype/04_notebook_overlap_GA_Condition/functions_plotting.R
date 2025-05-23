
####### Functions related to the figure on overlap of gestational age and condition

# Make a funtion to calculate the number of DE genes on a given tissue and in a given contrast
df_results_per_contrast_tissue <- function(tissue, contrast) {
  
  # string to be removed
  string1 <- paste("_",tissue,sep="")
  
  # This is the full table of all cells: ~0 + GA_Condition: Early
  model2 <- read.csv(paste(dir, "model2-fulltable/Alldata-TopTags-EdgeR-pseudobulk-allcelltypes-",tissue,"-factor-GA_Condition-contrast-",contrast,"-CellTypeManual.l3-20240828.csv",sep=""))
  
  
  # Number of genes DE in a given contrast
  number_of_de = data.frame()
  temp = data.frame()
  number_of_de <- model2 %>% filter(padj < p_thr)  %>% filter(abs(logFC) > logFC_threshold) %>%
    group_by(cell_type) %>%
    dplyr::count() %>%
    # filter(n>1) %>%
    arrange(desc(n))
  
  # Add extra columns
  number_of_de$Tissue <- tissue
  number_of_de$Contrast <- contrast
  
  
  # Replace the name tissue in the celltype column 
  number_of_de$CellTypeManual.l3 <- number_of_de$cell_type 
  number_of_de$CellTypeManual.l3 <- str_replace(number_of_de$CellTypeManual.l3, string1,"")
  
  # Add the information on cell groups
  number_of_de <- number_of_de %>% inner_join(cellgroups)
  number_of_de$cellgroup = factor(number_of_de$cellgroup, levels = celltype.order.l1, ordered = TRUE)
  
  return(number_of_de)
}

# Function to make the plot
barplot_number_genes_contrast <- function(tissue, df) {
  
  barplot <- ggplot(df, aes(x = CellTypeManual.l3, y= n, fill = Contrast)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.5) + 
    scale_fill_manual(values = c("#C43E96","#5AAA46")) +
    theme(axis.text.x = element_text(size=10, angle = 90),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 12),
          legend.title= element_blank()) +
    labs( x = "", y = "number of genes differentially expressed") 
  
  
  barplot <- barplot + facet_grid( ~ cellgroup, scales = "free", space = "free") +  theme(legend.position = "top") + theme(strip.text = element_text(size = 12, color = "black")) + NoLegend()
  
  # barplot <- barplot + scale_x_discrete(expand = c(0, 0))
  
  return(barplot)
}
