# Cell2location part 5: Estimate cell-type specific expression of every gene in the spatial data. Used to visualise the results. Loop throug donors
# Uses some custom made functions from plot_genes_per_cell_type.py
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


# samples: LR, PYI and NS got no decidua
samples_list = ["Visium-F1676VQ-Sp-CD-2",
                #"Visium-F1836LR-Sp-CD",
                #"Visium-F1836LR-Sp-CD-2",
                "Visium-F1828VB-Sp-CD",
                "Visium-F1678CM-CD-Sp-4-GEX", 
                "Visium-F1918AM-CD-sp",
                "Visium-F1686GS-Sp-CD",
                "Visium-F1906EP-CD-sp",
                "Visium-F1756VS-Sp-CD",
                "Visium-F1758LJ-Sp-CD",
                "Visium-F1844JD-CD-sp",
                "Visium-F1892VT-Sp-CD",
                #"Visium-F1950PYI-CD-sp",
                #"Visium-F1958NS-CD-sp",
                "Visium-F1668RK-CD-Sp-4-GEX",    
                "Visium-F1888FM-Sp-CD",
                "Visium-F1682RH-CD-Sp-4-GEX",
                "Visium-F1974SH-CD-sp",
                "Visium-F1788ND-Sp-CD",
                "Visium-F2026MM-CD-sp"]


# This is result folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'

# This is the directory where manual annotations are:
manual = '/home/ssd/ysanchez/Projects/FMI-all-Spatial-20240312/outputs/manual_annotations_myometrium/manual_annotations_myometrium_v2/'


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

    
    # Filter by the chorion region
     #adata  = adata[manual_regions.Barcode].copy()

     #print("check spatial object after subsetting")
     #print(adata)

    # run export_posterior
    adata = mod.export_posterior(
    adata, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': False}
    )

    # Compute expected expression per cell type
    expected_dict = mod.module.model.compute_expected_per_cell_type(
    mod.samples["post_sample_q05"], mod.adata_manager
    )

    # Add to anndata layers - one per cell type
    for i, n in enumerate(mod.factor_names_):
        adata.layers[n] = expected_dict['mu'][i]

    # Save anndata object with results
    #adata_file = f"{run_name}/sp.h5ad"
    #adata.write(adata_file)
    print("Check spatial object after computing expected expression per cell type")
    print(adata)

    #Open the manual annotations and subset the spots barcodes by the region called decidua 
    #HERE!!!!
    # # Use an f-string to dynamically insert sample_name into the file path
    file_path = f"{manual}/{sample_name}_regions_YSC1.csv"
    
    # Read the CSV file with the manual annotations
    manual_regions = pd.read_csv(file_path)
    
    # Plot only the chorion region
    manual_regions = manual_regions[manual_regions['regions'] == 'decidua']

    # Filter by the chorion region
    adata  = adata[manual_regions.Barcode].copy()

    # print("check spatial object after subsetting")
    # print(adata)

    image_key = "lowres"  # Replace with the correct key for your image (e.g., "hires" or "lowres")
    image = adata.uns["spatial"][sample_name]["images"][image_key]

    
    width = image.shape[1]
    height = image.shape[0]

    # Get the image size (height and width in pixels)
    height, width = image.shape[:2]

    print(f"Image size: {width} x {height} pixels")

    # print("Print Obs slot, is it here the spots barcodes?")
    # print(adata.obs)
    #print(adata.uns["spatial"])

    #print("Print layer CTB")
    #print(adata.layers["CTB"])
    
    # list cell types and genes for plotting
    #ctypes = ['Decidual', 'EVT']
    ctypes = ['Decidual']
    genes = ['TIMP3']

    with mpl.rc_context({'axes.facecolor':  'black'}):
    
    # # select one slide
        from cell2location.utils import select_slide
    slide = select_slide(adata, sample_name)

    # export the estimated expression 
    from plot_genes_per_cell_type import export_genes_per_cell_type
    export_genes_per_cell_type(slide, genes=['TIMP3'], ctypes=['Decidual'],sample_name = sample_name, output_dir=f"{results_folder}/output_csvs_decidua")

    
    # Import only the function. Plot the gene
    from plot_genes_per_cell_type import plot_genes_per_cell_type_custom
    plot_genes_per_cell_type_custom(slide, genes, ctypes,30000,30000);


    plt.savefig(f"{results_folder}/Spatial_plots_Myo_decidua_filtered_resized/{sample_name}_estimated_Gene_Expression_Decidual.png",
                 bbox_inches='tight')
    plt.close()
    
    print("Done", sample_name)
    print("###########################")

print("All done")