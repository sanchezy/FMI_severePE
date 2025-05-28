
# This script calculate a score using "AddModule" function in Seurat and produce a csv file.
# Calculation happens per capture area and outputs a file with the score per donor.
# Update on 20250210: This uses an updated hypoxia score gene list:
# ("HTRA1", "HTRA4", "FSTL3", "EGLN3", "TMEM45A")

library(dplyr)
library(Seurat)
library(patchwork)
library(knitr)
library(ggplot2)
library(stringr)
library(tidyverse)
library(magrittr)
library(paletteer)

# path_save
dir <- "~/Projects/FMI-all-Spatial-20240312/outputs/Spatial_scores/hypoxia_score/"

# order the donor: Early_Control, Early_PE, Late_Control, Late_PE,
tissue = "Placenta"

# Open the object
spatial.1 <- readRDS(file = "~/Projects/FMI-all-Spatial-20240312/objects/FMI-Placenta-Spatial-20241120.rds")

# Subset the object, select only the villi
spatial.1 <- subset(spatial.1, subset = regions %in% c("villi"))

# Create a column of a simplified name
spatial.1@meta.data$donor2 <- paste0("F", str_sub(spatial.1@meta.data$donor, start= -2))

# list of libraries or capture area
libraries <- unique(spatial.1@meta.data$library_id)

donors_to_plot <- c("FVQ", 
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
                    "FRK",
                    "FNS",
                    
                    "FFM",
                    "FRH",
                    "FND",
                    "FSH",
                    "FMM")

# Selected genes from above to create a score:
 hypoxia <- list(c("HTRA1", "HTRA4", "FSTL3", "EGLN3", "TMEM45A"))

# Loop through donors and calculate the ifn score
all_scores <- list()  # Initialize an empty list

for (i in 1:length(libraries)) {
  
 
  # Subset the object per donor
  spatial.temp <- subset(spatial.1, subset = library_id %in% libraries[i])
  
  # Check if the subsetting is ok, if not complain
  if (length(unique(spatial.temp@meta.data$library_id)) == 1) { 
    print(paste0("subsetting library: ", unique(spatial.temp@meta.data$library_id)))
  } else {
    print("subsetting did not work")
    next
  }
  
  # normalise the 
  spatial.temp <- NormalizeData(spatial.temp, normalization.method = "LogNormalize", scale.factor = 1000000)
  
  # calculate the ifn score
   spatial.temp <- AddModuleScore(
    object = spatial.temp,
    features = hypoxia,
    ctrl = 30,
    name = 'hypoxia'
  )
   
   # Extract the score and add to the list
   score <- spatial.temp@meta.data %>% select(hypoxia1, Tissue, donor2, library_id)
   all_scores[[i]] <- score
  
}

# Combine all the scores into a single data frame
final_scores <- do.call(rbind, all_scores)

# Write to a CSV file
write.csv(final_scores, paste0(dir, tissue, "-hypoxia-score-per-capture-area-five-genes-20250320.csv"), row.names = TRUE)

print("All done")
