#!/bin/bash
echo "This script creates a spatialranger script..."

ID=("F1682RH-B-Sp-4-GEX" "F1682RH-A-Sp-4-GEX" "F1682RH-CD-Sp-4-GEX" "F1682RH-CAM-Sp-4-GEX")
area=("A1" "B1" "C1" "D1")
image=("F1682RH_V10T06-054_A1.jpg" "F1682RH_V10T06-054_B1.jpg"  "F1682RH_V10T06-054_C1.jpg" "F1682RH_V10T06-054_D1.jpg")


for i in ${!ID[@]}; do

    cat << EOF > SpaceRanger-inputs/runSpacerangerSlurm-fmi-3rdRun-${ID[$i]}.sh
#!/usr/bin/env bash
#
#!!! This script runs spaceranger in vrtx-1 cluster
#
# =============================================================================
# Job Script
# =============================================================================
#
#SBATCH -J spa${ID[$i]}
#SBATCH --export=ALL
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --signal=2
#SBATCH --no-requeue
### Alternatively: --ntasks=1 --cpus-per-task={NUM_THREADS}
###   Consult with your cluster administrators to find the combination that
###   works best for single-node, multi-threaded applications on your system.
#SBATCH --mem=64G
#SBATCH -o ../Projects/FMI-Spatial-transcriptomics/logs/slurm-Spaceranger-fmi-${ID[$i]}.out
#SBATCH -e ../Projects/FMI-Spatial-transcriptomics/logs/slurm-Spaceranger-fmi-${ID[$i]}.err

cd ../Projects/FMI-Spatial-transcriptomics/
module load spaceranger
spaceranger count --id=Visium-${ID[$i]} --fastqs=../../fastq-files/FMI-pilot-fastq/F1682RH/ --slide=V10T06-054 --area=${area[$i]} --image=FMI-pilot-Visum-images/${image[$i]} --transcriptome=../../references/refdata-gex-GRCh38-2020-A  --jobmode=local  --localcores=8  --localmem=57 --sample=${ID[$i]}


EOF
echo "Finishing sh file runSpacerangerSlurm-fmi-3rdRun-${ID[$i]}.sh ..."
done