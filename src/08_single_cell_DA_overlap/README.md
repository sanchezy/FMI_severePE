We used [dawnn] (https://github.com/george-hall-ucl/dawnn), a DA tool for single-cell analysis by Dr George Hall to measure the overlap between *Condition* and *Gestational age*. Preprint can be found here: [dawnn-preprint] (https://www.biorxiv.org/content/10.1101/2023.05.05.539427v1.full.pdf+html)

This folder contains the following scripts:

1) `01_run_dawnn_tissue_overlap.R`: scripts to run DA taking as a main variable *Condition* or *Gestational age* and measuring the overlap. 
2) `02_plot_da-overlap-venn-diagrams.R`: scripts to visualise the results by plotting a venn diagram.
3) `03_barplot_da_overlap.R`: script to make a barplot of the overlap. 