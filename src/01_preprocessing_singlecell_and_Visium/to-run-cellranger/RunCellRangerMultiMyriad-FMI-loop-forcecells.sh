#!/bin/bash

# This script takes the csv file created before and runs cellranger multi in a loop for samples given in a list (ID)
# This script was run in a HPC Myriad at UCL

# This an example of script.
for ID in "F1758LJ-CAM-sn" "F1676VQ-CAM-sn" "F1678CM-CAM-sn" "F1682RH-CAM-sn" "F1686GS-CAM-sn" "F1668RK-CAM-sn"; do

# this is part of the script create and sh file and sends a job per $ID library to the cluster
    cat << EOF > cellrangermulti_single_job_submission.sh
#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -l mem=64G
#$ -N FM$ID
#$ -V
#$ -pe smp 8
#$ -cwd
#$ -o /lustre/scratch/scratch/sejjys1/Projects-Output/FMI-all-CellRanger7-20230222/logs/FMI-5thRun-default-$ID-$JOB_ID.out
#$ -e /lustre/scratch/scratch/sejjys1/Projects-Output/FMI-all-CellRanger7-20230222/logs/FMI-5thRun-default-$ID-$JOB_ID.err

export PATH=/lustre/scratch/scratch/sejjys1/cellranger7/cellranger-7.1.0:$PATH

cd Scratch/Projects-Output/FMI-all-CellRanger7-20230222/cellranger7-output-20230826/

cellranger multi --id=$ID --csv=/lustre/scratch/scratch/sejjys1/Projects-Output/FMI-all-CellRanger7-20230222/cellranger7-output-20230826/input-csv-files/cellrangermulti-csv-input-default-$ID.csv
EOF

    qsub cellrangermulti_single_job_submission.sh
    echo "Submitted job for id=$ID"

    # Wait for X min between jobs to avoid the cluster thinking you are a robot. Adjust accordingly
    sleep 400
done


