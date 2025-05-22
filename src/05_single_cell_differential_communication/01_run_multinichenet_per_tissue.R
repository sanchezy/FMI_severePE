# This script runs the cell-to-cell differential communication on multinichenetr. It follows closely the vigette of the tool.

# Author: Yara E. Sanchez Corrales

library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr)
library(Seurat)

####### setup
# path to save the output
path_network = "~/Projects/FMI-all-singlecell-20230614/outputs/multinichenet_results/"

# Subset by Tissue: Placenta, Myometrium, CAM, PBMC
tissue = "Placenta_Myometrium"
annotation_level = "CellTypeManual.l3"

path = paste0("~/Projects/FMI-all-singlecell-20230614/outputs/multinichenet_results_20240910/",tissue,"/")

# organism = "human"
# if(organism == "human"){
#   lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
#   lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
#   ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
#   colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
#   rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
# } else if(organism == "mouse"){
#   lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
#   lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
#   ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
#   colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
#   rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
# }
# 
# # saveRDS(lr_network, paste0(path,"lr_network_multinichenet.rds"))
# saveRDS(ligand_target_matrix, paste0(path,"ligand_target_matrix_multinichenet.rds"))

#  The database for organism = "human" was saved once, now it is just loaded
lr_network <- readRDS(paste0(path_network,"lr_network_multinichenet.rds"))
ligand_target_matrix <- readRDS(paste0(path_network,"ligand_target_matrix_multinichenet.rds"))

####### Read the single-cell object and subset by tissue. Remove donor FJJ.

# read the latest object and subset by a Tissue
data.1 <- readRDS(file = "/home/ysanchez/mnt/vrtx-directory-singlecell4/outputs/objects/FMI-all-20patients-20240827.rds")
# Remove the patient JJ
data.1 <- subset(data.1, subset = donor %in% "F2044JJ", invert = TRUE)

# subset by tissue
data.1 <- subset(data.1, subset = Tissue_old %in% c("Placenta", "Myometrium"))

####### Filter cells types coming from one donor or with less than 25 cells.

# This filtering was also done in EdgeR for DGE 

# Create extra columns in the metadata with the information about celltype and tissue and condition
data.1@meta.data$CellTypeManual.l3.Tissue <- paste(data.1@meta.data$CellTypeManual.l3, data.1@meta.data$Tissue_old,sep= "_")

# Replace symbol "/" by "_" as the paths gets confused with the paths.
# data.1@meta.data$CellTypeManual.l3.Tissue <- str_replace_all(data.1@meta.data$CellTypeManual.l3.Tissue_old,'/','_')
# unique(data.1@meta.data$CellTypeManual.l3.Tissue)

# LIST 1: remove clusters that come from <2 donor per cluster
cellstofilterout1 <- data.1@meta.data  %>% group_by(CellTypeManual.l3.Tissue,GA_Condition)  %>% 
  summarize(n_donor = n_distinct(donor))  %>%
  filter(n_donor < 2)

# LIST 2: remove early or late cells that do not have enough cells in Early_PE or Early_Control or Late_PE or Late_Control. Less than 20 per category
cellstofilterout2 <- data.1@meta.data %>% group_by(CellTypeManual.l3.Tissue,GA_Condition) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  filter(n < 25)

cellstofilterout <- c(cellstofilterout1$CellTypeManual.l3.Tissue, cellstofilterout2$CellTypeManual.l3.Tissue)

cellstofilterout <- unique(cellstofilterout)

print(paste("Filtering n number of cells: ", length(cellstofilterout)))
print(cellstofilterout)  

# Filter out these cells
Idents(data.1) = data.1$CellTypeManual.l3.Tissue

# follow this sintaxt: subset(x = pbmc, idents = c("CD4 T cells", "CD8 T cells"), invert = TRUE)
data.1 <- subset(x = data.1, idents = cellstofilterout, invert = TRUE)

# Track cell types that were filtered out
write.csv(cellstofilterout, file=paste0(path, tissue,"-cellsfilteredout-multinichenet-",annotation_level,"-20240910.csv"))

####### convert to a sce object
# If you start from a Seurat object, you can convert it easily to a SingleCellExperiment via sce = Seurat::as.SingleCellExperiment(seurat_obj, assay = "RNA")
sce = Seurat::as.SingleCellExperiment(data.1, assay = "RNA")

sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()

# Make sure the names are syntactically valid
SummarizedExperiment::colData(sce)$donor = SummarizedExperiment::colData(sce)$donor %>% make.names()
SummarizedExperiment::colData(sce)$GA_Condition = SummarizedExperiment::colData(sce)$GA_Condition %>% make.names()
SummarizedExperiment::colData(sce)$CellTypeManual.l3 = SummarizedExperiment::colData(sce)$CellTypeManual.l3 %>% make.names()
SummarizedExperiment::colData(sce)$CellTypeManual.l2 = SummarizedExperiment::colData(sce)$CellTypeManual.l2 %>% make.names()
SummarizedExperiment::colData(sce)$CellTypeManual.l1 = SummarizedExperiment::colData(sce)$CellTypeManual.l1 %>% make.names()

#######  Step 1: Extract cell type abundance and expression information from receiver and sender cell types, and link this expression information for ligands of the sender cell types to the corresponding receptors of the receiver cell types

# Define important metadata. 
sample_id = "donor"
group_id = "GA_Condition"
celltype_id = "CellTypeManual.l3"
covariates = NA
batches = NA

# Define set of senders and receiver cells types. Here we want to compare All-vs-All 
senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()

min_cells = 10

abundance_expression_info = get_abundance_expression_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, min_cells = min_cells, senders_oi = senders_oi, receivers_oi = receivers_oi, lr_network = lr_network, batches = batches)

####### Step 2: Perform genome-wide differential expression analysis of receiver and sender cell types to define DE genes between the conditions of interest. Based on this analysis, we can define the logFC/p-value of ligands in senders and receptors in receivers, and define the set of affected target genes in the receiver.

contrasts_oi = c("'Early_PE-Early_Control','Late_PE-Late_Control'")

contrast_tbl = tibble(contrast =
                        c("Early_PE-Early_Control","Late_PE-Late_Control"),
                      group = c("Early_PE","Late_PE"))

## Perform the DE analysis for each cell type.
DE_info = get_DE_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, batches = batches, covariates = covariates, contrasts_oi = contrasts_oi, min_cells = min_cells)

celltype_de = DE_info$celltype_de$de_output_tidy

# Combine DE information for ligand-senders and receptors-receivers (similar to step1 - abundance_expression_info$sender_receiver_info)
sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network)
  
####### Step 3: Predict NicheNet ligand activities and NicheNet ligand-target links based on these differential expression results
logFC_threshold = 0.50
p_val_threshold = 0.05
fraction_cutoff = 0.05

p_val_adj = FALSE 

top_n_target = 150

# Run the NicheNet ligand activity analysis
ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(get_ligand_activities_targets_DEgenes(
  receiver_de = celltype_de,
  receivers_oi = receivers_oi,
  ligand_target_matrix = ligand_target_matrix,
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = p_val_adj,
  top_n_target = top_n_target
)))

#######  Step 4: Use the information collected above to prioritize all sender-ligand—receiver-receptor pairs.

# In the 3 previous steps, we calculated expression, differential expression and NicheNet activity information. Now we will combine these different types of information in one prioritization scheme.
# We will set our preference for this dataset as follows - and recommend the user to use the same weights by default:

prioritizing_weights_DE = c("de_ligand" = 1,
                            "de_receptor" = 1)
prioritizing_weights_activity = c("activity_scaled" = 2)

prioritizing_weights_expression_specificity = c("exprs_ligand" = 2,
                                                "exprs_receptor" = 2)

prioritizing_weights_expression_sufficiency = c("frac_exprs_ligand_receptor" = 1)

prioritizing_weights_relative_abundance = c( "abund_sender" = 0,
                                             "abund_receiver" = 0)


prioritizing_weights = c(prioritizing_weights_DE, 
                         prioritizing_weights_activity, 
                         prioritizing_weights_expression_specificity,
                         prioritizing_weights_expression_sufficiency, 
                         prioritizing_weights_relative_abundance)

# Make necessary grouping data frame

sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group")
}

# Run the prioritization

prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  prioritizing_weights = prioritizing_weights,
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender
))

# In multi-sample datasets, we have the opportunity to look whether expression of ligand-receptor across all samples is correlated with the expression of their by NicheNet predicted target genes. This is what we will do with the following line of code:
lr_target_prior_cor = lr_target_prior_cor_inference(prioritization_tables$group_prioritization_tbl$receiver %>% unique(), abundance_expression_info, celltype_de, grouping_tbl, prioritization_tables, ligand_target_matrix, logFC_threshold = logFC_threshold, p_val_threshold = p_val_threshold, p_val_adj = p_val_adj)

#######  Step 5: Save the results. 
multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
) 
multinichenet_output_test = make_lite_output(multinichenet_output)

save = TRUE
if(save == TRUE){
  saveRDS(multinichenet_output, paste0(path, tissue ,"_multinichenet_output.rds"))
}

print("All done")