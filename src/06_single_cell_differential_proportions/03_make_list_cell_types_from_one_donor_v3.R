# Make a list of cells types per tissue per gestational age that come from n_donor < 2 and only exists in one condition (for instance, only PE)
# These are the filters used for propeller analysis in Feb/March 2024.


library(dplyr)
library(patchwork)
library(Seurat)

data.obj <- readRDS(file = "/Users/ysanchez/Documents/Projects-analysis/FMI-sc-all/objects/FMI-all-20patients-20240827.rds")
data.obj <- subset(data.obj, subset = donor %in% "F2044JJ", invert = TRUE)

data.1 <- data.obj@meta.data

print(unique(data.1$donor))

filteredcells_early = data.frame()
filteredcells_late = data.frame()

temp = data.frame()

######## GESTATIONAL AGE = EARLY
ga = "Early"

#First do it for the placenta 
tissue = "Placenta"

# Subset by location 
data.1 <- subset(data.1, subset = Tissue_old %in% tissue)

# LIST 1: remove clusters that come from <2 donor per cluster in either Early_PE, Early_Control (this is clusters that came from one donor)
cellstofilterout1  <- data.1 %>% filter(GA_Category == ga ) %>% group_by(CellTypeManual.l3, GA_Condition)  %>% 
  summarize(n_donor = n_distinct(donor))  %>%
  filter(n_donor < 2) 

# LIST 2: remove clusters that come only from one condition
cellstofilterout2  <- data.1 %>% filter(GA_Category == ga ) %>% group_by(CellTypeManual.l3)  %>% 
  summarize(n_condition = n_distinct(GA_Condition))  %>%
  filter(n_condition < 2) 

cellstofilterout <- c(cellstofilterout1$CellTypeManual.l3, cellstofilterout2$CellTypeManual.l3)
cellstofilterout <- unique(cellstofilterout)

filteredcells_early <- data.frame(cellstofilterout, tissue, ga)

print("Early")
print((filteredcells_early))


# Now, loop through the rest of the tissues and append the dataframe
tissues <- list("Myometrium","CAM","PBMC")

for (tissue in tissues) {
  print(paste("Analysis for tissue: ", tissue))
  
  data.1 <- data.obj@meta.data
  
  # Subset by location 
  data.1 <- subset(data.1, subset = Tissue_old %in% tissue)
  
  # LIST 1: remove clusters that come from <2 donor per cluster in either Early_PE, Early_Control (this is clusters that came from one donor)
  cellstofilterout1  <- data.1 %>% filter(GA_Category == ga ) %>% group_by(CellTypeManual.l3, GA_Condition)  %>% 
    summarize(n_donor = n_distinct(donor))  %>%
    filter(n_donor < 2) 
  
  # LIST 2: remove clusters that come only from one condition
  cellstofilterout2  <- data.1 %>% filter(GA_Category == ga ) %>% group_by(CellTypeManual.l3)  %>% 
    summarize(n_condition = n_distinct(GA_Condition))  %>%
    filter(n_condition < 2) 
  
  cellstofilterout <- c(cellstofilterout1$CellTypeManual.l3, cellstofilterout2$CellTypeManual.l3)
  cellstofilterout <- unique(cellstofilterout)

  temp <- data.frame(cellstofilterout, tissue, ga)
  
  # Append the dataframe
  filteredcells_early <- rbind(filteredcells_early,temp)
  print("Done early")
}

######## GESTATIONAL AGE = LATE 

ga = "Late"

#First do it for the placenta 
tissue = "Placenta"

data.1 <- data.obj@meta.data

# Subset by location 
data.1 <- subset(data.1, subset = Tissue_old %in% tissue)

# LIST 1: remove clusters that come from <2 donor per cluster in either Early_PE, Early_Control (this is clusters that came from one donor)
cellstofilterout1  <- data.1 %>% filter(GA_Category == ga ) %>% group_by(CellTypeManual.l3, GA_Condition)  %>% 
  summarize(n_donor = n_distinct(donor))  %>%
  filter(n_donor < 2) 

# LIST 2: remove clusters that come only from one condition
cellstofilterout2  <- data.1 %>% filter(GA_Category == ga ) %>% group_by(CellTypeManual.l3)  %>% 
  summarize(n_condition = n_distinct(GA_Condition))  %>%
  filter(n_condition < 2) 

cellstofilterout <- c(cellstofilterout1$CellTypeManual.l3, cellstofilterout2$CellTypeManual.l3)
cellstofilterout <- unique(cellstofilterout)

filteredcells_late <- data.frame(cellstofilterout, tissue, ga)

print("Placenta_Late")
print(head(filteredcells_late))

# loop through tissues and ga and append the dataframe

# Now, loop through the rest of the tissues and append the dataframe
tissues <- list("Myometrium","CAM","PBMC")

for (tissue in tissues) {
  print(paste("Analysis for tissue: ", tissue))
  
  data.1 <- data.obj@meta.data
  
  # Subset by location 
  data.1 <- subset(data.1, subset = Tissue_old %in% tissue)
  
  # LIST 1: remove clusters that come from <2 donor per cluster in either Early_PE, Early_Control (this is clusters that came from one donor)
  cellstofilterout1  <- data.1 %>% filter(GA_Category == ga ) %>% group_by(CellTypeManual.l3, GA_Condition)  %>% 
    summarize(n_donor = n_distinct(donor))  %>%
    filter(n_donor < 2) 
  
  # LIST 2: remove clusters that come only from one condition
  cellstofilterout2  <- data.1 %>% filter(GA_Category == ga ) %>% group_by(CellTypeManual.l3)  %>% 
    summarize(n_condition = n_distinct(GA_Condition))  %>%
    filter(n_condition < 2) 
  
  cellstofilterout <- c(cellstofilterout1$CellTypeManual.l3, cellstofilterout2$CellTypeManual.l3)
  cellstofilterout <- unique(cellstofilterout)
  
  temp <- data.frame(cellstofilterout, tissue, ga)
  
  # Append the dataframe
  filteredcells_late <- rbind(filteredcells_late,temp)
}

# Append both early and late and save the dataframe
filteredcells <- rbind(filteredcells_early,filteredcells_late)
print("Head all filtered cells")
colnames(filteredcells)[1] <- "CellTypeManual.l3"
colnames(filteredcells)[2] <- "Tissue"

print(head(filteredcells))

# save it
write.csv(filteredcells,"/Users/ysanchez/Documents/Projects-analysis/FMI-sc-all/outputs/metadata/sc-FMI-all-single-cell-celltypes-l3-from-one-donor_20240914.csv")

print("All done")