#!/bin/bash

# This script deletes the new file called "tissue_positions_list.csv" this file was created to be opened by cell2location but it prevents Seurat to open it.
# IMPORTANT: Because of spaceranger making breaking changes regarding the
# structure and names of its output files (!!!), there is now a bug by default
# with scanpy that only becomes apparent late in the pipeline. To prevent this,
# the file $SAMPLE/outs/spatial/tissue_positions.csv needs to be copied to
# $SAMPLE/outs/spatial/tissue_positions_list.csv and the header line deleted
# (see: https://github.com/scverse/scanpy/issues/2499#issuecomment-1607268186)
# Find this file at: ysanchez@194.82.242.226:Projects/cell2location/scripts

# The text file contains a list of libraries
cat list_Visium_samples_20240904.txt | while read ID; do

##remove
rm /home/ssd/ysanchez/Projects/FMI-Spatial-transcriptomics/Spaceranger-output/$ID/outs/spatial/tissue_positions_list.csv 
done

