This folder contains scripts for the preprocess of single-cell data and spatial data.

For **single-cell analysis**:

* `01_run-cellranger` - Scripts to go from Fastq files to count matrices. We used [cellranger7](https://www.10xgenomics.com/support/jp/software/cell-ranger/8.0/getting-started/cr-what-is-cell-ranger#).
* `02_run-cellbender` - Scripts for Background removal. We used [cellbender](https://cellbender.readthedocs.io/en/latest/introduction/index.html)
* `03_run-scrublet` - Scripts for doublet detection: [scrublet](https://github.com/swolock/scrublet)
* Classification fetal or maternal. The scripts were done by Jose Moreno and can be found here: [Freemuxlet](https://github.com/josemovi/fetal_maternal_freemuxlet)

For **spatial analysis**:
* `04_run-spaceranger` - Scripts to run [spaceranger](https://www.10xgenomics.com/support/cn/software/space-ranger/3.0/analysis/running-pipelines/poly-a-based-assay-count-spatial-gex)