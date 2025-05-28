################
#### Script to export seurat object to h5ad (compatible with scanpy and other python tools). Use this h5ad to train cell2location model at level 2 and 3. 
#### Similar to prepare a reference to Xenium.
#### Author: Yara E. Sanchez Corrales
################

# Load libraries and Seurat object.
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggvenn)
library(stringr)
library(knitr)
library(zellkonverter)
library(SingleCellExperiment)
library(scater)
library(scran)
library(SeuratDisk)

# Read the single cell seurat object
data.1 <- readRDS(file = "~/Projects/FMI-all-singlecell-20230308/outputs/objects/FMI-all-20patients-20240827.rds")
data.1

# This is a list of the whole genome - number of genes are 36601
q <- read.csv("/home/ssd/ysanchez/Projects/FMI-all-singlecell-20230308/scripts/reference_Xenium/gene_id_and_symbol_10X_whole_genome.csv")
print(head(q))

# check rownames(data.1) are in the same order as q$gene_name
# Add an extra column to check the order of barcodes is correct

check1 <- ifelse(rownames(data.1) == q$gene_name, "OK","check")
print(table(check1))
# setdiff(rownames(data.1),q$V2)

##########  Create a new object
# # get only the raw counts
# rename the rownames to ensemblID
counts <- GetAssayData(data.1,assay = "RNA",slot = "counts") 
# rownames(counts) <- q$gene_id

data.1 <- CreateSeuratObject(counts = counts, meta.data = data.1@meta.data)
print(head(rownames(data.1)))

#  Add information of gene names
data.1@assays$RNA@meta.features$gene_name <- q$gene_name

# Add Ensembl ID as part of the meta.features
data.1@assays$RNA@meta.features$ensembleID <- q$gene_id
# data.1@assays$RNA@meta.features$feature_type < q$feature_type


# Check meta features
print(dplyr::glimpse(GetAssay(data.1)@meta.features))

# Output
## Rows: 58,604
## Columns: 12
## $ feature_type          <fct> Gene Expression, Gene Expression, Gene Expressio…
## $ highly.variable       <lgl> FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,…
## $ mvp.mean              <dbl> 6.398244e-05, 2.274395e-03, 6.175251e-05, 1.3728…
## $ mvp.dispersion        <dbl> 0.8350443, 2.4422800, 1.2953346, 2.6563521, NaN,…
## $ mvp.dispersion.scaled <dbl> -0.5739472, 0.5332035, -0.2568744, 0.6806679, 0.…
## $ mean                  <dbl> 3.881082e-05, 1.079995e-03, 3.267138e-05, 4.7735…
## $ std                   <dbl> 5.573522e-03, 3.173095e-02, 5.634017e-03, 8.0407…
## $ ensembl_version       <chr> "ENSG00000223972.5", "ENSG00000227232.5", "ENSG0…
## $ feature_is_filtered   <lgl> FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,…
## $ feature_name          <fct> DDX11L1, WASH7P, MIR6859-1, MIR1302-2HG, MIR1302…
## $ feature_reference     <fct> NCBITaxon:9606, NCBITaxon:9606, NCBITaxon:9606, …
## $ feature_biotype       <fct> gene, gene, gene, gene, gene, gene, gene, gene, …

# save the object
# saveRDS(data.1, file = "~/Projects/FMI-all-singlecell-20230308/outputs/objects/FMI-all-20patients-20240827-reduced.rds")

data.1 %>% dplyr::glimpse()

# Check the slots
print(data.1@assays$RNA@counts[23:6,67:75])


print(data.1@assays$RNA@data[23:6,67:75])

# data.1@assays$RNA@scale.data[23:6,7:10]

##### Confirm that all the values stored in the matrix are integers. If your counts contain non-integer values, they are potentially normalized and not true UMI counts. 
# Please make sure your data have not been normalized. Otherwise, downstream analysis will be incorrect.

all.equal(GetAssayData(object=data.1, assay="RNA", slot="counts")@x, as.integer(GetAssayData(object=data.1, assay="RNA", slot="counts")@x))

# Output
## [1] TRUE

# Take a look at the first set of values, they should all be integer 

head(GetAssayData(object=data.1, assay="RNA", slot="counts")@x)

# Output
## [1]  1  2 18  1  2  1

#####  Do a standard transformation of counts using Seurat

SaveH5Seurat(object = data.1, filename = "~/Projects/FMI-all-singlecell-20230308/outputs/objects/FMI-all-20patients-20240827-reduced.h5Seurat")
Convert(source = "~/Projects/FMI-all-singlecell-20230308/outputs/objects/FMI-all-20patients-20240827-reduced.h5Seurat", dest = "h5ad")

print("All done")
