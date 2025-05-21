# This script runs CellBender in UCL myriad cluster. # For each sample name in the list, the loop below creates a file called
# CellBender_single_job_submission.sh and submits this using qsub. This file is
# then overwritten in the next iteration of the loop. This script runs in Job Submission mode.

ffor ID in "F1758LJ-CAM-sn" "F1676VQ-CAM-sn" "F1678CM-CAM-sn" "F1682RH-CAM-sn" "F1686GS-CAM-sn" "F1668RK-CAM-sn"; do

    cat << EOF > CellBender_single_job_submission.sh
#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -l mem=64G
#$ -N cBen$ID
#$ -V
#$ -pe smp 8
#$ -cwd
#$ -o /lustre/scratch/scratch/sejjys1/Projects-Output/FMI-all-CellRanger7-20230222/CellBender-output-20230902/logs/CellBender-OUT-${ID}-$JOB_ID
#$ -e /lustre/scratch/scratch/sejjys1/Projects-Output/FMI-all-CellRanger7-20230222/CellBender-output-20230902/logs/CellBender-ERROR-${ID}-$JOB_ID

module load python3/recommended
cellbender remove-background --input /lustre/scratch/scratch/sejjys1/Projects-Output/FMI-all-CellRanger7-20230222/cellranger7-output-20230826/$ID/outs/multi/count/raw_feature_bc_matrix --output /lustre/scratch/scratch/sejjys1/Projects-Output/FMI-all-CellRanger7-20230222/CellBender-output-20230902/$ID-CellBender-out.h5  --expected-cells 10000  --total-droplets-included 20000 

EOF

    qsub CellBender_single_job_submission.sh
    echo "Submitted job for id=$ID"

    # Wait for sometime between jobs to avoid the cluster thinking you are a robot
    sleep 800
done






