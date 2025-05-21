#!/bin/bash -l
#$ -l h_rt=42:00:00
#$ -l mem=64G
#$ -N scrubloop4
#$ -V
#$ -pe smp 8
#$ -cwd
#$ -o /lustre/scratch/scratch/sejjys1/Projects-Output/FMI-all-CellRanger7-20230222/scrublet-using-CellBender-input-5thRun/logs/$JOB_ID-scrublet-out
#$ -e /lustre/scratch/scratch/sejjys1/Projects-Output/FMI-all-CellRanger7-20230222/scrublet-using-CellBender-input-5thRun/logs/$JOB_ID-scrublet-error

module load python3/recommended
python3 /lustre/scratch/scratch/sejjys1/Scripts/scrublet-scRNA-FMI-loop-myriad-taking-cellbender-as-input.py 