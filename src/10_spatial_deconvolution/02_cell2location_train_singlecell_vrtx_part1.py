# Cell2location Part 1: Train the single_cell model. In this script, I train the single-cell model using the latest single-cell reference and level2 of annotation.
# Update: the single-cell reference was updated in Aug 2024. I used the script export_h5ad_object_2024083.Rmd to export the single-cell object.
# First activate the conda env in vrtx-1: conda activate cell2loc_env

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
results_folder = '/home/ssd/ysanchez/Projects/cell2location/single_cell_training/results_CellTypeManual_l2_20240830'

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'


######## Single-cell data
# Read data
adata_ref = sc.read_h5ad('/home/ssd/ysanchez/Projects/FMI-all-singlecell-20230308/outputs/objects/FMI-all-20patients-20240827-reduced.h5ad')
print("Reading single cell reference")
print(adata_ref.shape)

print("%%%%%%%%%%%%%%%%%%")
print("singlecell object:")
print(adata_ref)

adata_ref.var['SYMBOL'] = adata_ref.var.index

# rename 'GeneID-2' as necessary for your data
adata_ref.var.set_index('ensembleID', drop=True, inplace=True)

#print(adata_ref.__dict__['_raw'].__dict__['_var'])
#adata_ref.__dict__['_raw'].__dict__['_var'] = adata_ref.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})

#del(adata_ref.var['_index']) 

print("%%%%%%%%%%%%%%%%%%")
print("singlecell object shape")
print(adata_ref)

from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_ref = adata_ref[:, selected].copy()

# Check again the object
print("%%%%%%%%%%%%%%%%%%")
print(adata_ref)
print("adata shape here:")
print(adata_ref.shape)
print("adata cellID")
print(adata_ref.var['SYMBOL'])

# Estimation of reference cell type signatures (NB regression)
# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        #batch_key='Individual',
                        # cell type, covariate used for constructing signatures
                        labels_key='CellTypeManual.l2',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        #categorical_covariate_keys=['Method']
                       )

# create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

# view anndata_setup as a sanity check
print("view model setup")
mod.view_anndata_setup()

# TRain the model
mod.train(max_epochs=350, use_gpu=False)

print("model finished training")

mod.plot_history(20)

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': False}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/FMI-all-20patients-20240827-reduced-cell2loc-CellTypeManual-l2.h5ad"
print(adata_file)

adata_ref.write(adata_file)
adata_file

adata_ref = mod.export_posterior(
    adata_ref, use_quantiles=True,
    # choose quantiles
    #add_to_obsm=["q95", "q0001"],
    sample_kwargs={'batch_size': 2500, 'use_gpu': False}
)

# Check adata_ref
print("%%%%%%%%%%%%%%%%%%")
print("singlecell object:")
print(adata_ref)
print("adata cellID")
print(adata_ref.var['SYMBOL'])
print("Model:")
print(mod)
print("Finished")