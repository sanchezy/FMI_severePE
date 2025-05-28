This folder contains bash scripts to prepare the Spatialranger output and to run [cell2location](https://colab.research.google.com/github/BayraktarLab/cell2location/blob/master/docs/notebooks/cell2location_tutorial.ipynb) in a HPC cluster (vrtx)

## Contents

1) `01_script_change_name_tissue_position_v3.sh`: script to change the name of a file called 'tissue_position.csv' needed to run cell2location.
2) `02_script_rm_tissue_positions_list.sh`: script to remove the file above, need it to load the data in Seurat.
3) `03_run_cell2loc_vrtx.sh`: script to run cell2location script. Modify the script accordingly. 
4) `list_Visium_samples_20240904.txt`: list of Visium capture areas. 

