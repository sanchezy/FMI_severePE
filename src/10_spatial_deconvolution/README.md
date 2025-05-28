This folder contains notebooks and scripts related to Visium spatial deconvolution using [cell2location] ( https://colab.research.google.com/github/BayraktarLab/cell2location/blob/master/docs/notebooks/cell2location_tutorial.ipynb) 


1) `01_export_h5ad_object.R`: script to export an h5ad object from a Seurat object. Important to train the single-cell model in part 1.
2) `02_cell2location_train_singlecell_vrtx_part1.py`: part 1, train a single-cell model.
3) `03_cell2location_check_training_singlecell_vrtx_part2.py`: part 2, check the training of the single-cell model is reasonable.
4) `04_cell2location_map_loop_part3.py`: part 3, cell2location deconvolution of Visium. This script loops through capture areas.
5) `05_cell2location_plot_spatial_vrtx_part4.py`: part 4, visualisation of the deconvolution. Plot one or several cell types.
6) `06_cell2loc_estimate_gene_expression_celltype_part5.py`: part 5, estimate the gene expression from the deconvolution.
7) `07_plot_genes_per_cell_type.py`: function to estimate gene expression per cell type. Fine tune this function.
8) `08_bash_scripts`: scripts to prepare the data and run the above scripts in a HPC cluster. 

