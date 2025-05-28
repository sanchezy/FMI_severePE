# Cell2location Part 2: check the trained single-cell model. This script will plot the diagnostic plots for checking the single-cell training is correct.
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

print("Checking training from:")
print(results_folder)
# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

#load the pre-trained model
adata_file = f"{ref_run_name}/FMI-all-20patients-20240827-reduced-cell2loc-CellTypeManual-l2.h5ad"
adata_ref = sc.read_h5ad(adata_file)

mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)

# view anndata_setup as a sanity check
mod.view_anndata_setup()

# check the elbow plot loss history. It should go down and level off
print("check mod atributes:")
print(mod)


# If more training is needed, one needs to start the trainging again: mod.train(max_epochs=250, use_gpu=False), export the posteriors and save the model and sc object again.
# In this section, we export the estimated cell abundance (summary of the posterior distribution).
#adata_ref = mod.export_posterior(
#    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': False}
#)


# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': False}
)

#Examine QC plots. 
print("QC plot shold be a diagonal")
#print(mod.plot_QC())

mod.plot_QC()
plt.savefig(f"{results_folder}/reconstruction_accuracy_histogram_CellTypeManual_l2.png",
                bbox_inches='tight')
plt.close()

mod.plot_history(0)
plt.savefig(f"{results_folder}/training_ELBO_history_CellTypeManual_l2.png",
                   bbox_inches='tight')
plt.close()

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' 
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}' 
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']

print("Check the means_per_cluster")
print(inf_aver.iloc[0:5, 0:5])

