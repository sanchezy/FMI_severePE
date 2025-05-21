# This is a sample script to create the input file to run cellranger multi. For each sample name in the list, the loop below creates a csv file that will be the input to cellranger multi. Give a list of files here or within a file.
# This script was run in a HPC Myriad at UCL

# This is a list.
for ID in "F1758LJ-CAM-sn" "F1676VQ-CAM-sn" "F1678CM-CAM-sn" "F1682RH-CAM-sn" "F1686GS-CAM-sn" "F1668RK-CAM-sn"; do

cat << EOF > input-multi-csv-files-forcecells-Cellranger7/cellrangermulti-csv-input-$ID.csv

# this is the csv file that will be the input to cellranger multi 
[gene-expression]
reference,/shared/ucl/apps/spaceranger/references/refdata-gex-GRCh38-2020-A
force-cells, 10000
[vdj]
reference,/lustre/scratch/scratch/sejjys1/references/vdj_reference
[libraries]
fastq_id,fastqs,feature_types,subsample_rate
${ID}_GEX,/lustre/scratch/scratch/sejjys1/Data/FMI-S2-snGEX-VDJ-fastq/snGEX_VDJ/,gene expression,
${ID}_VDJ,/lustre/scratch/scratch/sejjys1/Data/FMI-S2-snGEX-VDJ-fastq/snGEX_VDJ/,vdj-t,
EOF
echo "Finishing cvs file for sample $ID..."
done