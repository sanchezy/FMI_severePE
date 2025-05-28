
################
##### Function that loops through samples(images) and returns a table called nearest_neighbour_region with nearest neigbour region per spot
################

nebs_per_image <- function(seurat_obj,ImageName) {
  
  # Ensure Image_name is valid
  if (!ImageName %in% names(seurat_obj@images)) {
    stop(paste("Error: Image_name", ImageName, "not found in seurat_obj@images"))
  }
  
  # Read the metadata, information about regions is stored here
  data <- seurat_obj@meta.data %>% filter(Image_name %in% ImageName)
  
  # Access coordinates for the specified Image_name
  coordinates <- seurat_obj@images[[ImageName]]@coordinates
  
  # Ensure coordinates exist
  if (is.null(coordinates)) {
    stop(paste("Error: Coordinates not found for Image_name", ImageName))
  }
  
  # Use nn2 to find the nearest neighbours
  nebs <- nn2(coordinates[, 2:3], k = 2)[[1]]  # Extract nearest neighbour indices
  
  # discard column 1 (these will almost always be the cell itself)
  nebs = nebs[,2] 
  
  # supply these row indices to the dataset to return the phenotypes of the neighbours stored in the column called "regions"
  # Get the regions of the nearest neighbours
  neighbour_regions <- data[nebs, "regions"]
  
  # Add it as a column to the data
  data$nearest_neighbour_region <- neighbour_regions
  
  # Select the required columns
  result <- data %>% select(Barcode, nearest_neighbour_region)
  
  return(result)
  
}

################
##### Function that takes the neigh information and returns a matrix of frequency. It also filters by a parameter such as: GA_Condition
################
 
neigh_matrix <- function(seurat_obj, condition) {
  
  # filter by condition 
  data <- seurat_obj@meta.data %>% filter(GA_Condition %in% condition)
  
  # filter by the spots that have a max SC.CTB
  data <- data %>% filter(SC.CTB_Abundance %in% "SC.CTB")
  
  # get it as a dataframe
  nebs <- data$nearest_neighbour_region
  
  # return a matrix
  nebs = as.matrix(table(data$regions, nebs))
  
  return(nebs)
  
}
