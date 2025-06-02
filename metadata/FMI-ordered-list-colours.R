# Cell types ordered and colours.

########### Colours
colors_cellgroups <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2","#CC79A7","#E5CCFF")
cellgroups_transparent <- alpha(colors_cellgroups, 0.6) 

# order of cell types:
celltype.order.l1 <- c("B","T/NK","Myeloid","Trophoblast","Stromal", "Vascular")

# Early_PE, Early_Control, Late_PE, Late_Control
# colors_conditions <- c("#C43E96", "#F2B342", "#95AA46", "#925A44")
colors_conditions2 <- c("#0072B2", "#E69F00", "#56B4E9", "#F0E442")

# colors tissue: Placenta, Myometrium, CAM, PBMC 
colours_tissue <- c("#66CC99", "#F0E442", "#D55E00","#9999CC")

tissue.order <- c("PVBP","Myometrium","CAM","PBMC")

########### Order of cell types

# order the donor: Early_Control, Early_PE, Late_Control, Late_PE,
donor.order <- c("FVQ", "FVB", "FLR","FCM","FGS","FVS","FEP","FAM","FRK","FLJ","FJD","FVT","FYI","FNS","FRH","FND","FFM","FSH","FMM")
condition.order <- c("Early_Control","Early_Control","Early_Control","Early_PE","Early_PE","Early_PE","Early_PE","Early_PE","Late_Control","Late_Control","Late_Control","Late_Control","Late_Control","Late_Control","Late_PE","Late_PE","Late_PE","Late_PE","Late_PE")
ga_cat <- c("Early","Early","Early","Early","Early","Early","Early","Early","Late","Late","Late","Late","Late","Late","Late","Late","Late","Late","Late")

donor.order.condition <- data.frame(donor = donor.order, condition = condition.order, GA_cat = ga_cat )

# Including FJJ
donor.order2 <- c("F1676VQ", "F1828VB", "F1836LR","F2044JJ","F1678CM","F1686GS","F1756VS","F1906EP","F1918AM","F1668RK","F1758LJ","F1844JD","F1892VT","F1950PYI","F1958NS","F1682RH","F1788ND","F1888FM","F1974SH","F2026MM")
condition.order2 <- c("Early_Control","Early_Control","Early_Control","Early_Control","Early_PE","Early_PE","Early_PE","Early_PE","Early_PE","Late_Control","Late_Control","Late_Control","Late_Control","Late_Control","Late_Control","Late_PE","Late_PE","Late_PE","Late_PE","Late_PE")
ga_cat2 <- c("Early","Early","Early","Early","Early","Early","Early","Early","Early","Late","Late","Late","Late","Late","Late","Late","Late","Late","Late","Late")

donor.order.conditionJJ <- data.frame(donor = donor.order2, condition = condition.order2)


# Note the change of tissue = "Placenta" for "PVBP".
tissue.order <- c("PVBP","Myometrium","CAM","PBMC")

#  Cell type order
# Order all : c("B","Myeloid","Stromal","T/NK","Trophoblast")

celltype_order.l2 <- c("B","plasma/plasmablast","CD4","CD8","MAIT","gdT","NKT","NK","ILC","Monocyte","Macrophage","Hofbauer","Dendritic",  "pDC","Mast","Neutrophils","Platelet","LED", "Endothelial","Decidual","Fibroblast","Muscle","Epithelium","STB" ,   "CTB",    "EVT" ,   "SC-CTB" ,"SC-EVT")

#  Annotation level 2
celltype_order_B.l2 <- c("B","plasma/plasmablast")
celltype_order_T.l2 <- c("CD4","CD8","MAIT","gdT","NKT","NK","ILC")
celltype_order_Myeloid.l2 <- c("Monocyte","Macrophage","Hofbauer","Dendritic",  "pDC","Mast","Neutrophils","Platelet")
celltype_order_Stromal.l2 <- c("LED", "Endothelial","Decidual","Fibroblast","Muscle","Epithelium")
celltype_order_Tropho.l2 <- c("STB" ,   "CTB",    "EVT" ,   "SC-CTB" ,"SC-EVT")

# Annotation level 3
# Order all : c("B","Myeloid","Stromal","T/NK","Trophoblast")
celltype_order.l3 <- c("Immature_Naive_B",
                       "Naive_B",
                       "Memory_B",
                       "ABC",
                       "plasmablast",
                       "plasma",
                       
                       "Fetal_CD4_T_Naive_CM",
                       "CD4_T_Naive_CM",
                       "CD4_Th",
                       "FoxP3-Treg",
                       "CD4_TEM",
                       "CD4_T_CTL",
                       
                       "Fetal_CD8_T_Naive_CM",
                       "CD8_T_Naive_CM",
                       "GZMK_CD8_T",
                       "CD8_TEM",
                       "CD8_T_Exhausted",
                       
                       "MAIT",
                       "gdT",
                       "Fetal_NKT",
                       "NKT",
                       "NKT-proliferating",
                       "CD56_NK",
                       "CD16_NK",
                       "NK",
                       "ILC",
                       
                       "CD14_Monocyte",
                       "Fetal-Monocyte",
                       "CD16_Monocyte",
                       "Nonclassical-monocyte",
                       "Fetal-Nonclassical-monocyte",
                       "Intermediate-macrophage",
                       "Macrophage",
                       "Macrophage-proliferating",
                       "Hofbauer",
                       "cDC1",
                       "cDC1-proliferating",
                       "cDC2",
                       "cDC3",
                       "pDC",
                       "Mast",
                       "Neutrophils",
                       "Platelet",
                       
                       "LED",
                       "Endothelial",
                       "Decidual_stromal",
                       "Fetal_fibroblast",
                       "Smooth_muscle",
                       "SC-Epithelial",
                       "Uterine_Epithelial",
                       
                       "CTB-proliferating",
                       "CTB",
                       "SC-CTB",
                       "EVT",
                       "SC-EVT",
                       "eEVT",
                       "GC",
                       "STB"
)

# Order per cell group: c("B","Myeloid","Stromal","T/NK","Trophoblast")
celltype_order_B.l3 <- c("Immature_Naive_B",
                         "Naive_B",
                         "Memory_B",
                         "ABC",
                         "plasmablast",
                         "plasma")


celltype_order_T.l3 <- c("Fetal_CD4_T_Naive_CM",
                         "CD4_T_Naive_CM",
                         "CD4_Th",
                         "FoxP3-Treg",
                         "CD4_TEM",
                         "CD4_T_CTL",
                         
                         "Fetal_CD8_T_Naive_CM",
                         "CD8_T_Naive_CM",
                         "GZMK_CD8_T",
                         "CD8_TEM",
                         "CD8_T_Exhausted",
                         
                         "MAIT",
                         "gdT",
                         "Fetal_NKT",
                         "NKT",
                         "NKT-proliferating",
                         "CD56_NK",
                         "CD16_NK",
                         "NK",
                         "ILC")

celltype_order_Myeloid.l3 <- c("CD14_Monocyte",
                               "Fetal-Monocyte",
                               "CD16_Monocyte",
                               "Nonclassical-monocyte",
                               "Fetal-Nonclassical-monocyte",
                               "Intermediate-macrophage",
                               "Macrophage",
                               "Macrophage-proliferating",
                               "Hofbauer",
                               "cDC1",
                               "cDC1-proliferating",
                               "cDC2",
                               "cDC3",
                               "pDC",
                               "Mast",
                               "Neutrophils",
                               "Platelet")

celltype_order_Stromal.l3 <- c(
  "Decidual_stromal",
  "Fetal_fibroblast",
  "Smooth_muscle",
  "SC-Epithelial",
  "Uterine_Epithelial")

celltype_order_Vascular.l3 <- c("LED",
                                "Endothelial")

celltype_order_Tropho.l3 <- c("CTB-proliferating",
                              "CTB",
                              "SC-CTB",
                              "EVT",
                              "SC-EVT",
                              "eEVT",
                              "GC",
                              "STB" )



# Dataframes per group
B.group <- data.frame(CellTypeManual.l3 = celltype_order_B.l3, cellgroup = "B")
T.group <- data.frame(CellTypeManual.l3 = celltype_order_T.l3, cellgroup = "T/NK")
Myeloid.group <- data.frame(CellTypeManual.l3 = celltype_order_Myeloid.l3, cellgroup = "Myeloid")
Trophoblast.group <- data.frame(CellTypeManual.l3 = celltype_order_Tropho.l3, cellgroup = "Trophoblast")
Stromal.group <- data.frame(CellTypeManual.l3 = celltype_order_Stromal.l3, cellgroup = "Stromal")
Vascular.group <- data.frame(CellTypeManual.l3 = celltype_order_Vascular.l3, cellgroup = "Vascular")
# concatenate the dataframes
cellgroups <- rbind(B.group,T.group,Myeloid.group,Trophoblast.group,Stromal.group, Vascular.group)

################ 
# Define consistent colors for spatial_domains
spatial_domains_colors <- c(
  "villi" = "#E69F00", # Orange
  "decidua" = "#009E73", #  Green
  "vasculature_intermediate_villi" = "#56B4E9", # Sky Blue
  "fibrin" = "#F0E442", # Yellow
  "muscle" = "#0072B2", # Blue
  "chorion" = "#D55E00", # Red
  "stromal_layer" = "#56B4E9", # Sky Blue
  "background" = "#000000" # Black or grey
)