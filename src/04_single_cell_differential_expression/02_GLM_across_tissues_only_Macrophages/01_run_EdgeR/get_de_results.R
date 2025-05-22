
# functio to run EdgeR to test for differences across tissues for a specific cell type.
get_de_results <- function(fit, contrast_matrix, contrast_name, n = Inf,output_dir = ".", filename = NULL) {
  
  # Extract the contrast vector
  contrast_vector <- contrast_matrix[, contrast_name]
  
  # Perform LRT
  lrt <- edgeR::glmLRT(fit, contrast = contrast_vector)
  
  # Get results
  df <- as.data.frame(edgeR::topTags(lrt, n = n))
  
  # Adjust p-values (FDR)
  df$padj <- p.adjust(df$PValue, method = "BH")
  
  # Build filename if not provided
  if (is.null(filename)) {
    filename <- paste0("contrast_", contrast_name, ".csv")
  }
  
  # Ensure output directory ends with slash
  if (!grepl("/$", output_dir)) {
    output_dir <- paste0(output_dir, "/")
  }
  
  # Write to CSV
  write.csv(df, file = paste0(output_dir, "Across-Tissues-PE-TopTags-EdgeR-pseudobulk-Macrophages-factor-contrast-", contrast_name,"-20250410.csv", sep=""), row.names = TRUE)
  
  
  return(df)
}