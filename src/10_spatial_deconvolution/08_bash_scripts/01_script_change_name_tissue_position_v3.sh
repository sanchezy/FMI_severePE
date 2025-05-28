#!/bin/bash

# This script makes the new file called "tissue_positions_list.csv" and deletes the first line.
# IMPORTANT: Because of spaceranger making breaking changes regarding the
# structure and names of its output files (!!!), there is now a bug by default
# with scanpy that only becomes apparent late in the pipeline. To prevent this,
# the file $SAMPLE/outs/spatial/tissue_positions.csv needs to be copied to
# $SAMPLE/outs/spatial/tissue_positions_list.csv and the header line deleted
# (see: https://github.com/scverse/scanpy/issues/2499#issuecomment-1607268186)

# The text file contains a list of libraries
cat list_Visium_samples_20240904.txt | while read ID; do

printf $ID
cp /home/ssd/ysanchez/Projects/FMI-Spatial-transcriptomics/Spaceranger-output/$ID/outs/spatial/tissue_positions.csv /home/ssd/ysanchez/Projects/FMI-Spatial-transcriptomics/Spaceranger-output/$ID/outs/spatial/tissue_positions_list.csv

#remove the header line
sed -i 1d /home/ssd/ysanchez/Projects/FMI-Spatial-transcriptomics/Spaceranger-output/$ID/outs/spatial/tissue_positions_list.csv
done

