# This scripts loads all the outputs from CellBender - one h5 file per library - and merged into a single object using a loop.
# for multiple libraries, each barcode is appended with a number (to avoid the same barcode in multiple cells). This makes the object compatible with the output of cellranger aggregate that can be loaded into cloupe.

library(dplyr)
library(Seurat)
library(stringr)

# Create a list of files to append. The order of the files is very important
list <- c(
  "F1676VQ-AB2-Sc-4-1",
  "F1676VQ-AB2-Sc-4-2",
  "F1676VQ-CD-Sc-4",
  "F1676VQ-E-Sc-4")

nsamples <- length(list)

# Path to the CellBender Output
pathtoCellBenderOutput <-'../Projects/FMI-all-singlecell-20230308/CellBender-temporal/'

#Open the first file: for instance: CellBender-output-20230228/F1676VQ-AB2-Sc-4-1-CellBender-out_filtered.h5
data.file <- paste0(pathtoCellBenderOutput,list[1],"-CellBender-out_filtered.h5")
data.data <- Read10X_h5(filename = data.file, use.names = TRUE)

# create Seurat object called cellbender
cellbender <- CreateSeuratObject(counts = data.data, project = "FMI", assay = "RNA")

# cellbender$sampleCellBender <- list[1]
cellbender$sample_id <- list[1]
rm(data.data)
rm(data.file)

# # Create a list of files to append
for (i in 2:nsamples) {
  
  data.file <- paste0(pathtoCellBenderOutput,list[i],"-CellBender-out_filtered.h5")
  data.data <- Read10X_h5(filename = data.file, use.names = TRUE)
 
   # create Seurat object
  temp <- CreateSeuratObject(counts = data.data, project = "FMI", assay = "RNA")
  
  temp$sample_id <- list[i]
  
  #Change the name of Barcodes appending a dash (the barcodes of the first library have a -1, the second library -2, etc. This makes the object compatible with the output of cellranger aggregate and cloupe)
  cellsnames <- Cells(temp)
  cellsnames <- as.data.frame(cellsnames)
  
  newstring <-paste0("-",i)
  
  cellsnames$cellsnames<-str_replace(cellsnames$cellsnames, "-1", newstring)
  temp <- RenameCells(temp, new.names = cellsnames$cellsnames)
  print(newstring)
  
  # Merge two objects
  temp2 <- merge(x = cellbender, y=temp, project = "FMI")
  cellbender <- temp2
  rm(data.data)
  rm(data.file)
  rm(cellsnames)
  rm(temp)
  rm(temp2)
}

## Add metadata 

# Depending on your libraries, you could add metadata afterwards. For instance, I add: 1)condition, 2)  age, etc.

# Save the object
saveRDS(cellbender, file = "../../outputs/FMI.allcells.cellbender_20230327_temp.rds")
