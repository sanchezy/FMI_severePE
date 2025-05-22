#!/usr/bin/env bash

#SBATCH -J DE-EdgeR
#SBATCH --export=ALL
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --signal=2
#SBATCH --no-requeue
### Alternatively: --ntasks=1 --cpus-per-task={NUM_THREADS}
###   Consult with your cluster administrators to find the combination that
###   works best for single-node, multi-threaded applications on your system.
#SBATCH --mem=72G
#SBATCH -o /home/ssd/ysanchez/Projects/FMI-all-singlecell-20230308/outputs/EdgeR-single-cell-20240828/EdgeR-single-cell-model2-GA_Condition-20240828-level2/logs/DE-EdgeR-cluster-model2-GA_Condition-CellTypeManual.l2-20240828.out
#SBATCH -e /home/ssd/ysanchez/Projects/FMI-all-singlecell-20230308/outputs/EdgeR-single-cell-20240828/EdgeR-single-cell-model2-GA_Condition-20240828-level2/logs/DE-EdgeR-cluster-model2-GA_Condition-CellTypeManual.l2-20240828.err

module load R
Rscript /home/ssd/ysanchez/Projects/FMI-all-singlecell-20230308/scripts/Differential_Expression_20240828/run-EdgeR-celltype-model2-GACondition.R