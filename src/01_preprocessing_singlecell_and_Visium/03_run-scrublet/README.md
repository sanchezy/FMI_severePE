This folder contains scripts for running scrublet, a tool to identfy doublets. The output will be concatenated and added as a metadata into the Seurat object. These scripts use the HPC cluster at UCL (myriad)

1) Script to run scrublet per library as a loop: `01_scrublet-scRNA-FMI-loop-myriad-taking-cellbender-as-input.py`
2) Run scrublet from a HPC cluster: `02_run-scrublet-loop-HPC-myriad.sh` 
3) The outputs can be concatenated into a single file that can be integrated as metadata in the Seurat object: `03_scrublet-mergeoutput-file.py`