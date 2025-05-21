#!/usr/bin/env bash
#
#!!! This script runs spaceranger in vrtx-1 cluster
#
# =============================================================================
# Job Script
# =============================================================================
#
#SBATCH -J spaceranger
#SBATCH --export=ALL
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --signal=2
#SBATCH --no-requeue
### Alternatively: --ntasks=1 --cpus-per-task={NUM_THREADS}
###   Consult with your cluster administrators to find the combination that
###   works best for single-node, multi-threaded applications on your system.
#SBATCH --mem=64G
#SBATCH -o slurm-Spaceranger-fmi-F1668RK-B.out
#SBATCH -e slurm-Spaceranger-fmi-error-F1668RK-B.out

export PATH=/opt/spaceranger-1.3.0:$PATH
spaceranger count --id=Visium-F1668RK-B-Sp-4-GEX --fastqs=FMI-pilot-fastq/F1668RK-B-Sp-4-GEX --slide=V10T06-057 --area=A1 --image=FMI-pilot-Visum-images/F1668RK-jpg/F1668RK_V10T06-057_A1.jpg --transcriptome=references/refdata-gex-GRCh38-2020-A  --jobmode=local  --localcores=8  --localmem=57 --sample=F1668RK-B-Sp-4-GEX

