This folder contains notebooks and scripts related to differential expression analysis using Generalised linear models (GLM) per cell type and comparing conditions. 

The model is given by the following matrix: `design <- model.matrix(~0 + GA_Condition, y$samples)`

1) `01_run_EdgeR`: scripts to run EdgeR. We use a HPC at UCL.
2) `02_notebook_number_of_DEG`: notebook plotting the number of DEG per contrast and cell type.
3) `03_notebook_volcano_plots`: notebook about ploting volcano plots. Note the volcano plots are not rendered here (code is commented) but tables are visible.
