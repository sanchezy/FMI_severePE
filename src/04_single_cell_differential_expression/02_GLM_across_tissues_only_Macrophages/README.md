This folder contains notebooks and scripts related to differential expression analysis using Generalised linear models (GLM) on Macrophages across tissues in disease. 

The model is given by the following matrix: `design <- model.matrix(~0 + Condition.Tissue, y$samples)`

1) `01_run_EdgeR`: scripts to run EdgeR. We use a HPC at UCL.
2) `02_notebook_volcanoplots `: notebook plotting the number of DEG and volcano plots per contrast and cell type.

