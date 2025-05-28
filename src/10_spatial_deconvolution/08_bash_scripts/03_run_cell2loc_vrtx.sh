#!/usr/bin/env bash
#
#
#
# =============================================================================
# Job Script
# =============================================================================
#
#SBATCH -J cell2lc2
#SBATCH --export=ALL
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --signal=2
#SBATCH --no-requeue
### Alternatively: --ntasks=1 --cpus-per-task={NUM_THREADS}
###   Consult with your cluster administrators to find the combination that
###   works best for single-node, multi-threaded applications on your system.
#SBATCH --mem=64G
#SBATCH --open-mode=append
#SBATCH -o /home/ssd/ysanchez/Projects/cell2location/logs/estimate_gene_expression-loop_CellTypeManual2-20250123.out
#SBATCH -e /home/ssd/ysanchez/Projects/cell2location/logs/estimate_gene_expression-loop_CellTypeManual2-20250123.err

source /home/ssd/ysanchez/miniconda3/etc/profile.d/conda.sh 

conda activate cell2loc_env

# Change the name of the script accordingly to the step.
python /home/ssd/ysanchez/Projects/cell2location/scripts/cell2loc_estimate_gene_expression_celltype_loop_20250114.py 
