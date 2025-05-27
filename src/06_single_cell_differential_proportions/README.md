This folder contains scripts and notebook related to differential proportion analysis using Generalised linear models (GLM) per cell type and comparing conditions implemented in [`propeller`](https://phipsonlab.github.io/propeller-paper-analysis/pbmcJP.html)

This folder contains the following:

1) `01_propeller_cellproportions.R`: scripts to run cell proportion analysis for tissues PBMC, PVBP, Myometrium.
2) `02_propeller_cellproportions_CAM.R`: scripts to run cell proportion analysis for tissue = CAM. This tissue was profiled using two methods, thus 'Method' was included in the matrix design.
3) `03_make_list_cell_types_from_one_donor_v3.R`: script to make a list of cell types that come from one donor or only one condition.
4) `04_notebook_propeller_results_PropRatio_per_contrast_v3-2025-05-27.md`: notebook visualising the results from the analysis per contrast.
5) `05_function_create_df.R`: function used in the notebook to create visualisations of results. 