# Cell2location Part 3: Cell2location_map - deconvolute the Visium, loop through samples and produce and output. This script follows very closely the colab: https://colab.research.google.com/github/BayraktarLab/cell2location/blob/master/docs/notebooks/cell2location_tutorial.ipynb

# First activate the conda env in vrtx-1: conda activate cell2loc_env (use script run_cell2loc_map_loop_vrtx.sh)


import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib as mpl

import cell2location
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
from cell2location.utils import select_slide
from cell2location.plt import plot_spatial

# Set paths to data and results used through the document:
sp_data_folder = '/home/ssd/ysanchez/Projects/FMI-Spatial-transcriptomics/Spaceranger-output/'
results_folder = '/home/ssd/ysanchez/Projects/cell2location/single_cell_training/results_CellTypeManual_l2_20240226'

print("Checking training from:")
print(results_folder)
# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'


#load the pre-trained model
adata_file = f"{ref_run_name}/FMI-all-20patients-20240110-reduced-cell2loc-CellTypeManual-l2.h5ad"
adata_ref = sc.read_h5ad(adata_file)

mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)

# view anndata_setup as a sanity check
print("This is the single-cell trained model:")
mod.view_anndata_setup()

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
print("Check inf_aver.index")
print(inf_aver.index)


# IMPORTANT: Because of spaceranger making breaking changes regarding the
# structure and names of its output files (!!!), there is now a bug by default
# with scanpy that only becomes apparent late in the pipeline. To prevent this,
# the file $SAMPLE/outs/spatial/tissue_positions.csv needs to be copied to
# $SAMPLE/outs/spatial/tissue_positions_list.csv and the header line deleted
# (see: https://github.com/scverse/scanpy/issues/2499#issuecomment-1607268186)


# This list was run on the 20240520
samples_list = ["Visium-F1668RK-A-Sp-4-GEX", "Visium-F1668RK-B-Sp-4-GEX","Visium-F1668RK-CD-Sp-4-GEX",
                "Visium-F1676VQ-A-Sp-4-GEX", "Visium-F1676VQ-B-Sp-4-GEX","Visium-F1676VQ-CD-Sp-4-GEX","Visium-F1676VQ-Sp-CD-2",
                "Visium-F1678CM-A-Sp-4-GEX","Visium-F1678CM-B-Sp-4-GEX","Visium-F1678CM-CD-Sp-4-GEX",
                "Visium-F1682RH-A-Sp-4-GEX","Visium-F1682RH-B-Sp-4-GEX","Visium-F1682RH-CD-Sp-4-GEX"]

             

for sample_name in samples_list:

    print("Processing", sample_name)

    adata_vis = sc.read_visium(f"/home/ssd/ysanchez/Projects/FMI-Spatial-transcriptomics/Spaceranger-output/{sample_name}/outs")
    
    #adata_vis.obs['sample'] = list(adata_vis.uns['spatial'].keys())[0]
    adata_vis.obs['sample'] = sample_name

    print("%%%%%%%%%%%%%%%%%%")
    print("spatial object after reading it from spaceranger")
    print(adata_vis)

  # This code would change symbols to ENSEMBL names, a best practice probably
    adata_vis.var['SYMBOL'] = adata_vis.var_names
    #adata.var['ensembleID'] = adata.var['genes_ids']
    
    adata_vis.var.set_index('gene_ids', drop=True, inplace=True)
    print("spatial object after renaming vars")
    print(adata_vis)

    # Uncommenting this line breaks things later
    # adata_vis.var.set_index('gene_ids', drop=True, inplace=True)
    adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]
    # remove MT genes for spatial mapping (keeping their counts in the object)
    adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
    adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]

    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    # prepare anndata for cell2location model
    cell2location.models.Cell2location.setup_anndata(adata = adata_vis)

    # create and train the model
    # Need to choose N_cells_per_location and detection_alpha (see
    # https://github.com/BayraktarLab/cell2location/blob/master/docs/images/Note_on_selecting_hyperparameters.pdf))
    mod = cell2location.models.Cell2location(
        adata_vis, cell_state_df=inf_aver, 
        # the expected average cell abundance: tissue-dependent 
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location=8,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection:
        detection_alpha=20
    ) 

    # Train cell2location (default max_epochs=30000)
    mod.train(max_epochs=3000,
              # train using full data (batch_size=None)
              batch_size=None,
              # use all data points in training because
              # we need to estimate cell abundance at all locations
              train_size=1,
              use_gpu=False,
             )
    
    mod.plot_history(100)
    plt.savefig(f"{run_name}/training_spatial_plot_history_CellTypeManual_l2_{sample_name}.png",
                   bbox_inches='tight')
    plt.close()


    # Visualise ELBO plot
    # mod.plot_history(100)

    adata_vis = mod.export_posterior(
        adata_vis,
        sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': False}
    )

    # Save model
    trained_model_name = f"{run_name}/cell2location_trained_model_{sample_name}"
    mod.save(trained_model_name, overwrite=True)
    print("model saved")
    
    # Save anndata object with results
    adata_file = f"{run_name}/cell2location_results_{sample_name}.h5ad"

    # Save a csv file with the results
    adata_vis.write(adata_file)
    print("This is adata_file")
    print(adata_file)

    # Save the csv of the  5% quantile of the posterior distribution, representing the value of cell abundance that the model has high confidence in (aka ‘at least this amount is present’)
    # to adata.obs with nice names for plotting
    adata_export_5 =  f"{run_name}/cell2location_results_q05_cell_abundance_w_sf_{sample_name}.csv"

    adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
    print("This is the q05_cell_abundance_w_sf")
    print(adata_vis.obsm['q05_cell_abundance_w_sf'])

    adata_vis.obsm['q05_cell_abundance_w_sf'].to_csv(adata_export_5)
    print("q05_cell_abundance saved")
    
    print("All done")