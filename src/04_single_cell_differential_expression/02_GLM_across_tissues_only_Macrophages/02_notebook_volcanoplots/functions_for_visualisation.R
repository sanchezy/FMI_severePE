# functions to make plots

###### Get the number of DEG of a given contrast. Asumme results from one contrasts are in one table.
number_of_deg <- function(table,p_thr) { 
  number_DEG <- table %>% filter(padj < p_thr) %>% 
    group_by(cell_type) %>%
    dplyr::count() %>%
    filter(n>1) %>%
    arrange(desc(n))
  return(number_DEG)
}

###### make a volcano plot, used Enhancedvolcano, load it beforehand.

plot_volcano <- function(table, comparison){ 
  
  results <- table
  
  plot1 <-  EnhancedVolcano(results,
                            lab = results$genes,
                            x = 'logFC',
                            y = 'PValue',
                            pCutoff = 0.05, 
                            pCutoffCol = "padj",
                            labSize = 3.0,
                            legendLabSize = 12,
                            legendIconSize = 4.0,
                            drawConnectors = FALSE,
                            widthConnectors = 0.25,
                            colConnectors = 'grey50',
                            
                            title = paste("Cell type: ", results$cell_type, sep=""),
                            subtitle = bquote(italic(comparison))
  )
  
  return(plot1)
}

###### make a summarised table

make_a_deg_table <- function(table, p_thr,logFC_thr) { 
  
  results <- table
  
  up = list()
  down = list()
  
  print("Results for celltype:")
  print(results$cell_type[1])
  
  results <- results %>% filter(padj < p_thr)  %>% select("genes","logFC","padj","cell_type")
  
  # genes up
  print("Genes upregulated in disease")
  up <- results %>% filter(logFC > logFC_thr) %>% arrange(desc(logFC))
  print(kable(up)) 
  cat("\n")
  
  # genes down
  print("Genes downregulated in disease")
  down <- results %>% filter(logFC < -logFC_thr) %>% arrange(desc(logFC))
  print(kable(down)) 
  cat("\n")
  
}