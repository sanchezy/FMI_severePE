This folder contains scripts for running cellbender and create a seurat object from the outputs. These steps use the HPC cluster at UCL (myriad)

1) Run cellbender per library. Use the script: `01_CellBender-FMI-loop-myriad.sh`
2) Concatenate the outputs of cellbender into a Seurat object for further processing. Example: `02_script-make-a-Seurat-object-from-cellbender.R` 