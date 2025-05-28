# Function used to estimate gene expression per cell type. Needs fine tuning depending on the tissue.
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np

import os


def plot_genes_per_cell_type(slide, genes, ctypes):
    n_genes = len(genes)
    n_ctypes = len(ctypes)
    fig, axs = plt.subplots(
        nrows=n_genes, ncols=n_ctypes + 1, figsize=(4.5 * (n_ctypes + 1) + 2, 5 * n_genes + 1), squeeze=False
    )
    # axs = axs.reshape((n_genes, n_ctypes+1))

    # plots of every gene
    for j in range(n_genes):
        # limit color scale at 99.2% quantile of gene expression (computed across cell types)
        quantile_across_ct = np.array(
            [
                np.quantile(slide.layers[n][:, slide.var["SYMBOL"] == genes[j]].toarray(), 0.992)
                for n in slide.uns["mod"]["factor_names"]
            ]
        )
        quantile_across_ct = np.partition(quantile_across_ct.flatten(), -2)[-2]
        sc.pl.spatial(
            slide,
            cmap="magma",
            color=genes[j],
            # layer=ctypes[i],
            gene_symbols="SYMBOL",
            ncols=4,
            size=1.3,
            img_key="hires",
            # limit color scale at 99.2% quantile of gene expression
            vmin=0,
            vmax="p99.2",
            ax=axs[j, 0],
            show=False,
        )

        # plots of every cell type
        for i in range(n_ctypes):
            sc.pl.spatial(
                slide,
                cmap="magma",
                color=genes[j],
                layer=ctypes[i],
                gene_symbols="SYMBOL",
                ncols=4,
                size=1.3,
                img_key="hires",
                # limit color scale at 99.2% quantile of gene expression
                vmin=0,
                vmax=quantile_across_ct,
                ax=axs[j, i + 1],
                show=False,
            )
            axs[j, i + 1].set_title(f"{genes[j]} {ctypes[i]}")

    return fig, axs

# changed function to be more custom made 
def plot_genes_per_cell_type_custom(slide, genes, ctypes, im_width,im_height):
    n_genes = len(genes)
    n_ctypes = len(ctypes)
    fig, axs = plt.subplots(
        nrows=n_genes, ncols=n_ctypes + 1, figsize=(4.5 * (n_ctypes + 1) + 2, 5 * n_genes + 1), squeeze=False
    )
    # axs = axs.reshape((n_genes, n_ctypes+1))

    # plots of every gene (in total)
    for j in range(n_genes):
        # limit color scale at 99.2% quantile of gene expression (computed across cell types)
        quantile_across_ct = np.array(
            [
                np.quantile(slide.layers[n][:, slide.var["SYMBOL"] == genes[j]].toarray(), 0.992)
                for n in slide.uns["mod"]["factor_names"]
            ]
        )
        quantile_across_ct = np.partition(quantile_across_ct.flatten(), -2)[-2]
        sc.pl.spatial(
            slide,
            cmap="magma",
            color=genes[j],
            # layer=ctypes[i],
            gene_symbols="SYMBOL",
            ncols=4,
            size=1.3,
            img_key="hires",
            # limit color scale at 99.2% quantile of gene expression
            vmin=0,
            vmax="p99.2",
            ax=axs[j, 0],
            show=False,
        )

        # plots of every gene decomposed by cell type. No cropping
        for i in range(n_ctypes):
            sc.pl.spatial(
                slide,
                cmap="magma",
                color=genes[j],
                layer=ctypes[i],
                gene_symbols="SYMBOL",
                ncols=4,
                size=1.3,
                img_key="hires",
                # limit color scale at 99.2% quantile of gene expression
                vmin=0,
                vmax=120,
                alpha_img = 0.7,
                ax=axs[j, i + 1],
                show=False,
                crop_coord=[0, im_width, 0, im_height],
                #alpha=0.5,
            )
            axs[j, i + 1].set_title(f"{genes[j]} {ctypes[i]}")

    return fig, axs


def export_genes_per_cell_type(slide, genes, ctypes,sample_name, output_dir="output_csvs"):
    """
    Export a CSV file for each unique combination of gene, slice, and cell type.
    
    Parameters:
    - slide: AnnData object containing spatial transcriptomics data.
    - genes: List of genes to export data for.
    - ctypes: List of cell types to export data for.
    - output_dir: Directory to save the exported CSV files (default: "output_csvs").
    """
   # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Retrieve a concise name for the slice
    slice_name = slide.uns.get("name", "unknown_slice")
    slice_name = str(slice_name).replace(" ", "_")  # Sanitize the name

    for gene in genes:
        for ctype in ctypes:
            # Extract data for this gene and cell type
            data = slide.obs.copy()
            data[f"{gene}_expression"] = slide.layers[ctype][:, slide.var["SYMBOL"] == gene].toarray().flatten()

            # Subset to only include the desired columns
            selected_columns = [f"{gene}_expression", "sample", "array_row", "array_col"]
            subset_data = data[selected_columns]

           # Construct the output file path
            output_file = os.path.join(output_dir, f"{gene}_{sample_name}_{ctype}.csv")

            # Save to CSV
            try:
                subset_data.to_csv(output_file, index=True)
                print(f"Exported: {output_file}")
            except Exception as e:
                print(f"Failed to export {output_file}: {e}")
           
