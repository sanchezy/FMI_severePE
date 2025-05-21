#!/usr/bin/env python

#This script Runs scrublet per library given as a list
# Outputs a text file with doblets scores and barcodes
# Edit: 1) samplesID, input_dir, output_dir and expected_doublet_rate
# It takes the output from cellbender as input, this requiered an extra function (from https://github.com/broadinstitute/CellBender/issues/128)

# Author: Yara E. Sanchez Corrales

import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import h5py

import tables
import scipy.sparse as sp
from typing import Dict, Optional

# Give a list of samples ID (these are the directory names from cellranger or cellbender). Note the order must be the same as the order of aggregate from cellranger.

samplesID = [
"F1758LJ-CAM-sn", 
"F1676VQ-CAM-sn",
"F1678CM-CAM-sn",
"F1682RH-CAM-sn", 
"F1686GS-CAM-sn",
"F1668RK-CAM-sn"]


# Initialise a variable that will add a string to the barcodes.
i = 1 

# Create an outputs into an existing directory:
output_dir = '/lustre/scratch/scratch/sejjys1/Projects-Output/FMI-all-CellRanger7-20230222/scrublet-using-CellBender-input-5thRun/'

# This function was copied from: https://github.com/broadinstitute/CellBender/issues/128
def anndata_from_h5(file: str,
                    analyzed_barcodes_only: bool = True) -> 'anndata.AnnData':
    """Load an output h5 file into an AnnData object for downstream work.

    Args:
        file: The h5 file
        analyzed_barcodes_only: False to load all barcodes, so that the size of
            the AnnData object will match the size of the input raw count matrix.
            True to load a limited set of barcodes: only those analyzed by the
            algorithm. This allows relevant latent variables to be loaded
            properly into adata.obs and adata.obsm, rather than adata.uns.

    Returns:
        adata: The anndata object, populated with inferred latent variables
            and metadata.

    """

    d = dict_from_h5(file)
    X = sp.csc_matrix((d.pop('data'), d.pop('indices'), d.pop('indptr')),
                      shape=d.pop('shape')).transpose().tocsr()

    # check and see if we have barcode index annotations, and if the file is filtered
    barcode_key = [k for k in d.keys() if (('barcode' in k) and ('ind' in k))]
    if len(barcode_key) > 0:
        max_barcode_ind = d[barcode_key[0]].max()
        filtered_file = (max_barcode_ind >= X.shape[0])
    else:
        filtered_file = True

    if analyzed_barcodes_only:
        if filtered_file:
            # filtered file being read, so we don't need to subset
            print('Assuming we are loading a "filtered" file that contains only cells.')
            pass
        elif 'barcode_indices_for_latents' in d.keys():
            X = X[d['barcode_indices_for_latents'], :]
            d['barcodes'] = d['barcodes'][d['barcode_indices_for_latents']]
        elif 'barcodes_analyzed_inds' in d.keys():
            X = X[d['barcodes_analyzed_inds'], :]
            d['barcodes'] = d['barcodes'][d['barcodes_analyzed_inds']]
        else:
            print('Warning: analyzed_barcodes_only=True, but the key '
                  '"barcodes_analyzed_inds" or "barcode_indices_for_latents" '
                  'is missing from the h5 file. '
                  'Will output all barcodes, and proceed as if '
                  'analyzed_barcodes_only=False')

    # Construct the anndata object.
    adata = anndata.AnnData(X=X,
                            obs={'barcode': d.pop('barcodes').astype(str)},
                            var={'gene_name': (d.pop('gene_names') if 'gene_names' in d.keys()
                                               else d.pop('name')).astype(str)},
                            dtype=X.dtype)
    adata.obs.set_index('barcode', inplace=True)
    adata.var.set_index('gene_name', inplace=True)

    # For CellRanger v2 legacy format, "gene_ids" was called "genes"... rename this
    if 'genes' in d.keys():
        d['id'] = d.pop('genes')

    # For purely aesthetic purposes, rename "id" to "gene_id"
    if 'id' in d.keys():
        d['gene_id'] = d.pop('id')

    # If genomes are empty, try to guess them based on gene_id
    if 'genome' in d.keys():
        if np.array([s.decode() == '' for s in d['genome']]).all():
            if '_' in d['gene_id'][0].decode():
                print('Genome field blank, so attempting to guess genomes based on gene_id prefixes')
                d['genome'] = np.array([s.decode().split('_')[0] for s in d['gene_id']], dtype=str)

    # Add other information to the anndata object in the appropriate slot.
    _fill_adata_slots_automatically(adata, d)

    # Add a special additional field to .var if it exists.
    if 'features_analyzed_inds' in adata.uns.keys():
        adata.var['cellbender_analyzed'] = [True if (i in adata.uns['features_analyzed_inds'])
                                            else False for i in range(adata.shape[1])]

    if analyzed_barcodes_only:
        for col in adata.obs.columns[adata.obs.columns.str.startswith('barcodes_analyzed')
                                     | adata.obs.columns.str.startswith('barcode_indices')]:
            try:
                del adata.obs[col]
            except Exception:
                pass
    else:
        # Add a special additional field to .obs if all barcodes are included.
        if 'barcodes_analyzed_inds' in adata.uns.keys():
            adata.obs['cellbender_analyzed'] = [True if (i in adata.uns['barcodes_analyzed_inds'])
                                                else False for i in range(adata.shape[0])]

    return adata


def dict_from_h5(file: str) -> Dict[str, np.ndarray]:
    """Read in everything from an h5 file and put into a dictionary."""
    d = {}
    with tables.open_file(file) as f:
        # read in everything
        for array in f.walk_nodes("/", "Array"):
            d[array.name] = array.read()
    return d


def _fill_adata_slots_automatically(adata, d):
    """Add other information to the adata object in the appropriate slot."""

    for key, value in d.items():
        try:
            if value is None:
                continue
            value = np.asarray(value)
            if len(value.shape) == 0:
                adata.uns[key] = value
            elif value.shape[0] == adata.shape[0]:
                if (len(value.shape) < 2) or (value.shape[1] < 2):
                    adata.obs[key] = value
                else:
                    adata.obsm[key] = value
            elif value.shape[0] == adata.shape[1]:
                if value.dtype.name.startswith('bytes'):
                    adata.var[key] = value.astype(str)
                else:
                    adata.var[key] = value
            else:
                adata.uns[key] = value
        except Exception:
            print('Unable to load data into AnnData: ', key, value, type(value))

############ This is scrublet code ############


for ID in samplesID:
    #Load counts matrix and gene list
    #Load the raw counts matrix as a scipy sparse matrix with cells as rows and genes as columns.
    # CellBender-output-20220912 was a test using  --total-droplets-included 20000. CellBender-output-20230228 took the output of cellranger7
    input_dir = '/lustre/scratch/scratch/sejjys1/Projects-Output/FMI-all-CellRanger7-20230222/CellBender-output-20230902/'
    
    adata = anndata_from_h5(input_dir + ID + '-CellBender-out_filtered.h5')
    adata.var_names_make_unique()

    #Initialize Scrublet object. Specify the expected doublet rate.
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.076)

    #Run the default pipeline, which includes:
    #Doublet simulation
    #Normalization, gene filtering, rescaling, PCA
    #Doublet score calculation
    #Doublet score threshold detection and doublet calling
    # Calling the main function
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                            min_cells=3, 
                                                            min_gene_variability_pctl=85, 
                                                            n_prin_comps=30)

    #print(help(scrub.scrub_doublets))                                                         

    # set manually the threshold in case the automatic threshold did not work
    scrub.call_doublets(threshold=0.25)

    # this histogram should be bimodal
    scrub.plot_histogram()
    #ptl.plot = scrub.plot_histogram();
    plt.savefig(output_dir + ID + '_scrublet_hist' +  '.png')

    print('Running UMAP...')
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    scrub.plot_embedding('UMAP', order_points=True)
    plt.savefig(output_dir + ID + '_scrublet_umap.png')
    #print('Done.')

    # Export the data into a csv file
    df = pd.DataFrame({
        'doublet_score': scrub.doublet_scores_obs_,
        'predicted_doublet': scrub.predicted_doublets_
    })

    #print(df)
  
    # Print number of doublets detected
    n_doublets = df.groupby('predicted_doublet').count()
    print(n_doublets)

    # Add the barcodes as a column in the last position
    df['barcodes'] = adata.obs_names
    # This step is not done yet. Get rid off the last string. It is always 1.
    #df.barcodes = df['barcodes'].str[:-1] 
    #print(barcodes)

    # Done in a different script. Add a string in the barcodes, this will allow to avoid repetions when concatenating textfiles from different libraries
    #df.barcodes = df.barcodes + str(i)

    df.to_csv(output_dir + ID + '_scrublet_output_table.csv', index=False)
    i = i+1 
    #print(df)

    # save the ambient RNA from cellbender 
    pd.DataFrame(adata.var).to_csv(output_dir + ID + '-ambient-expression-from-cellbender.csv')

    print('Done for sample ' + ID)


