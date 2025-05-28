# Cell2location part 4: the spatial map has happened. Here just load it to make some plots - loop through samples and save the plots.
# Load the manual annotations and subset accordingly.

#  First activate the enviroment: conda activate cell2loc_env
# run as: python /home/ssd/ysanchez/Projects/cell2location/scripts/cell2location_check_map_spatial_vrtx.py

# Load the packages
import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
import gc


import cell2location 

import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns

# silence scanpy that prints a lot of warnings
import warnings
warnings.filterwarnings('ignore')

# Set paths to data and results used through the document:
sp_data_folder = '/home/ssd/ysanchez/Projects//Projects/FMI-Spatial-transcriptomics/Spaceranger-output/'
results_folder = '/home/ssd/ysanchez/Projects/cell2location/single_cell_training/results_CellTypeManual_l2_20240830'

print("Checking training from:")
print(results_folder)

samples_list = ["Visium-F1668RK-CAM-Sp-4-GEX", "Visium-F1676VQ-CAM-1-sp","Visium-F1678CM-CAM-Sp-4-GEX",
                "Visium-F1682RH-CAM-Sp-4-GEX","Visium-F1686GS-Sp-CAM","Visium-F1756VS-Sp-CAM","Visium-F1758LJ-Sp-CAM",
                "Visium-F1788ND-Sp-CAM","Visium-F1828VB-Sp-CAM","Visium-F1836LR-Sp-CAM","Visium-F1844JD-CAM-sp",
                "Visium-F1888FM-Sp-CAM", "Visium-F1892VT-Sp-CAM","Visium-F1906EP-CAM-sp","Visium-F1918AM-CAM-sp","Visium-F1950PYI-CAM-sp",
                "Visium-F1958NS-CAM-sp","Visium-F1974SH-CAM-sp","Visium-F2026MM-CAM-sp"]



# This is result folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'

# This is the directory where manual annotations are:
manual = '/home/ssd/ysanchez/Projects/FMI-all-Spatial-20240312/outputs/manual_annotations_cam/manual_annotations_cam_v2'


for sample_name in samples_list:
    print("###########################")
    print("Processing", sample_name)

    run_name = f'{results_folder}/cell2location_map/cell2location_trained_model_{sample_name}'

    # Open the Visium file
    #sample_name = "Visium-F1888FM-Sp-CAM"
    adata_file = f"{results_folder}/cell2location_map/cell2location_results_{sample_name}.h5ad"
    adata = sc.read_h5ad(adata_file)

    mod = cell2location.models.Cell2location.load(f"{run_name}", adata)

    print("This is the model:")
    print(mod)                                              
    #add 5% quantile, representing confident cell abundance, 'at least this amount is present', 
    #to adata.obs with nice names for plotting. This line adds the values of abundace per cell type to "obs"
    adata.obs[adata.uns['mod']['factor_names']] = adata.obsm['q05_cell_abundance_w_sf']


    print("check spatial object")
    print(adata)

    #Open the manual annotations and subset the spots barcodes by the region called chorion 
    #HERE!!!!
    # Use an f-string to dynamically insert sample_name into the file path
    file_path = f"{manual}/{sample_name}_regions_YSC1.csv"
    
    # Read the CSV file with the manual annotations
    manual_regions = pd.read_csv(file_path)
    
    # Plot only the chorion region
    manual_regions = manual_regions[manual_regions['regions'] == 'chorion']

    # Print the Barcodes
    #print("Print barcodes corresponding to a chorion region")
    #print(manual_regions.Barcode)

    print("Print Obs slot, is it here the spots barcodes?")
    print(adata.obs)

    # Filter by the chorion region
    adata  = adata[manual_regions.Barcode].copy()

    print("check spatial object after subsetting")
    print(adata)

    # # select one slide
    from cell2location.utils import select_slide
    slide = select_slide(adata, sample_name)

    print("print the abundances estimated by the tool")
    print("Plotting celltypes")

    # select up to 6 clusters

    clust_labels = ['SC-CTB','EVT','CTB']
    ##Define your cluster colors
    clust_colors = {'SC-CTB': 'green', 'EVT': 'red', 'CTB':'blue'}

    #clust_col = [clust_colors[label] for label in clust_labels]

    # print("cluster labels")
    # print(clust_labels)

    # #print("cluster colours")
    # #print(clust_col)

    # print("Plotting celltypes")
    # # # plot in spatial coordinates
    # with mpl.rc_context({'axes.facecolor':  'black',
    #                  'figure.figsize': [4.5, 5]}):

    #     sc.pl.spatial(slide, cmap='magma',
    #               # show first 8 cell types
    #               color= clust_labels,
    #               ncols=2, size=1.3,
    #               img_key='hires',
    #               # limit color scale at 99.2% quantile of cell abundance
    #               vmin=0, vmax='p99.2'
    #              )

    #plt.savefig(f"{results_folder}/Spatial_plots_SC.EVT_EVT/plot_celltypes_EVT_SC-CTB_magma_{sample_name}.png",
    #            bbox_inches='tight')
    #plt.close()



    print("Now we use cell2location plotter that allows showing multiple cell types in one panel")
    # # plot in spatial coordinates
    from cell2location.plt import plot_spatial

    with mpl.rc_context({"figure.figsize": (15, 15)}):
        fig = plot_spatial(
            adata=slide,
            labels=clust_labels,
            reorder_cmap= [3,4,0],
            color=list({'SC-CTB': 'green', 'EVT': 'red', 'CTB':'blue'}),
            max_color_quantile=0.992,
            circle_diameter=6,
            show_img=True,
            colorbar_position='bottom',
            colorbar_shape={"horizontal_gaps": 0.2},
        )

    plt.savefig(f"{results_folder}/Spatial_plots_SC.EVT_EVT/plot_{sample_name}_merged_EVT_SC.CTB_CTB_newcolours.png",
                bbox_inches='tight')
    plt.close()
    
    print("Done", sample_name)
    print("###########################")

print("All done")