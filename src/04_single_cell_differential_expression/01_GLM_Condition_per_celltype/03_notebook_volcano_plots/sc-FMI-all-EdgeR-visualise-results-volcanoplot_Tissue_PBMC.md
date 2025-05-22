Single-cell FMI: check the output of EdgeR - do a volcano plot
================
true
May 21, 2025

# Introduction

We want to visualise the results of a DE analysis. The goal is to find
the DE genes in PE compared to a control. Check for one tissue at a
time. Patient FJJ has been filtered out.

-   model2: \~0 + GA\_Condition

The contrasts for model 2 are: <br>

-   Early\_disease = GA\_ConditionEarly\_PE -
    GA\_ConditionEarly\_Control,
-   Late\_disease = GA\_ConditionLate\_PE - GA\_ConditionLate\_Control,
-   AverageGA = (GA\_ConditionEarly\_PE +
    GA\_ConditionEarly\_Control)/2 - (GA\_ConditionLate\_PE +
    GA\_ConditionLate\_Control)/2,
-   AveragePE = (GA\_ConditionEarly\_PE + GA\_ConditionLate\_PE)/2 -
    (GA\_ConditionEarly\_Control+ GA\_ConditionLate\_Control)/2

*Note*: volcano plots are not rendered in this notebook.

``` r
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(edgeR)
library(scran)
library(EnhancedVolcano)
require(data.table)
require(knitr)
library(kableExtra)

source("/home/ssd/ysanchez/Projects/FMI-all-singlecell-20230308/scripts/single-cell-annotation/FMI-celltypes-ordered-list.R")
```

> Tissue: PBMC

``` r
# This notebook is for tissue and annotation level:
tissue= "PBMC" # "CAM", "PBMC", "Myometrium"
annotation_level <- "CellTypeManual.l3"

# List of cell types
cellgroups$CellTypeManual.l3 <- paste0(cellgroups$CellTypeManual.l3,"_",tissue)
```

# Volcano plots in Early disease

\[1\] “Results for celltype:” \[1\] “CD4\_T\_Naive\_CM\_PBMC” \[1\]
“Genes upregulated in disease”

| genes    |    logFC |      padj | cell\_type              |
|:---------|---------:|----------:|:------------------------|
| TRBV20-1 | 10.94350 | 0.0067726 | CD4\_T\_Naive\_CM\_PBMC |
| TRBV5-1  | 10.46952 | 0.0090602 | CD4\_T\_Naive\_CM\_PBMC |

\[1\] “Genes downregulated in disease”

| genes      |     logFC |      padj | cell\_type              |
|:-----------|----------:|----------:|:------------------------|
| HIST2H2AA4 | -3.462886 | 0.0185173 | CD4\_T\_Naive\_CM\_PBMC |

\[1\] “Results for celltype:” \[1\] “CD4\_Th\_PBMC” \[1\] “Genes
upregulated in disease”

| genes    |    logFC |      padj | cell\_type    |
|:---------|---------:|----------:|:--------------|
| TRBV20-1 | 9.702783 | 0.0000298 | CD4\_Th\_PBMC |
| TRBV7-3  | 9.089353 | 0.0454882 | CD4\_Th\_PBMC |

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type    |
|:------|----------:|----------:|:--------------|
| ZFP36 | -1.749902 | 0.0134292 | CD4\_Th\_PBMC |
| DUSP1 | -2.386139 | 0.0001113 | CD4\_Th\_PBMC |
| JUN   | -3.831316 | 0.0000001 | CD4\_Th\_PBMC |
| FOS   | -4.175597 | 0.0000000 | CD4\_Th\_PBMC |
| EGR1  | -6.801658 | 0.0001113 | CD4\_Th\_PBMC |

\[1\] “Results for celltype:” \[1\] “FoxP3-Treg\_PBMC” \[1\] “Genes
upregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type       |
|:------|----------:|----------:|:-----------------|
| ZFP36 | -1.686127 | 0.0266505 | FoxP3-Treg\_PBMC |
| DUSP1 | -2.808374 | 0.0006381 | FoxP3-Treg\_PBMC |
| IFIT2 | -3.751140 | 0.0375352 | FoxP3-Treg\_PBMC |
| FOS   | -7.135763 | 0.0000655 | FoxP3-Treg\_PBMC |

\[1\] “Results for celltype:” \[1\] “CD16\_NK\_PBMC” \[1\] “Genes
upregulated in disease”

| genes |    logFC |      padj | cell\_type     |
|:------|---------:|----------:|:---------------|
| PLET1 | 3.046683 | 0.0312236 | CD16\_NK\_PBMC |
| ERAP2 | 2.093320 | 0.0083760 | CD16\_NK\_PBMC |

\[1\] “Genes downregulated in disease”

| genes      |     logFC |      padj | cell\_type     |
|:-----------|----------:|----------:|:---------------|
| RGPD2      | -2.299031 | 0.0312236 | CD16\_NK\_PBMC |
| SOCS3      | -2.510564 | 0.0368424 | CD16\_NK\_PBMC |
| IFNG       | -2.732798 | 0.0010634 | CD16\_NK\_PBMC |
| AC020911.2 | -3.203531 | 0.0013482 | CD16\_NK\_PBMC |
| EGR1       | -5.181526 | 0.0220055 | CD16\_NK\_PBMC |

\[1\] “Results for celltype:” \[1\] “Nonclassical-monocyte\_PBMC” \[1\]
“Genes upregulated in disease”

| genes |    logFC |      padj | cell\_type                  |
|:------|---------:|----------:|:----------------------------|
| ERAP2 | 3.054127 | 0.0009645 | Nonclassical-monocyte\_PBMC |
| STAT2 | 2.053342 | 0.0477182 | Nonclassical-monocyte\_PBMC |

\[1\] “Genes downregulated in disease”

| genes      |     logFC |      padj | cell\_type                  |
|:-----------|----------:|----------:|:----------------------------|
| ZFP36      | -1.735887 | 0.0477182 | Nonclassical-monocyte\_PBMC |
| MID1IP1    | -1.753753 | 0.0318676 | Nonclassical-monocyte\_PBMC |
| NINJ1      | -1.862440 | 0.0201832 | Nonclassical-monocyte\_PBMC |
| SERTAD1    | -1.918022 | 0.0498404 | Nonclassical-monocyte\_PBMC |
| ATF3       | -2.158960 | 0.0132716 | Nonclassical-monocyte\_PBMC |
| NXT1       | -2.406974 | 0.0028895 | Nonclassical-monocyte\_PBMC |
| FOS        | -2.427009 | 0.0104282 | Nonclassical-monocyte\_PBMC |
| AL118516.1 | -2.534099 | 0.0201832 | Nonclassical-monocyte\_PBMC |
| CCL3       | -3.640773 | 0.0201832 | Nonclassical-monocyte\_PBMC |
| G0S2       | -3.688710 | 0.0097110 | Nonclassical-monocyte\_PBMC |
| JUN        | -4.709236 | 0.0002220 | Nonclassical-monocyte\_PBMC |
| C15orf48   | -5.220348 | 0.0104282 | Nonclassical-monocyte\_PBMC |
| CXCL8      | -5.407211 | 0.0002220 | Nonclassical-monocyte\_PBMC |

\[1\] “Results for celltype:” \[1\] “cDC2\_PBMC” \[1\] “Genes
upregulated in disease”

| genes   |    logFC |      padj | cell\_type |
|:--------|---------:|----------:|:-----------|
| ZFP36L2 | 2.114112 | 0.0455299 | cDC2\_PBMC |

\[1\] “Genes downregulated in disease”

| genes    |     logFC |      padj | cell\_type |
|:---------|----------:|----------:|:-----------|
| ZFP36    | -2.047590 | 0.0056385 | cDC2\_PBMC |
| LGALS3   | -2.118716 | 0.0419571 | cDC2\_PBMC |
| H1FX     | -2.188622 | 0.0497542 | cDC2\_PBMC |
| IER2     | -2.206944 | 0.0083979 | cDC2\_PBMC |
| ARL4C    | -2.289995 | 0.0179144 | cDC2\_PBMC |
| NFKBIA   | -2.571819 | 0.0056679 | cDC2\_PBMC |
| PID1     | -2.617159 | 0.0063366 | cDC2\_PBMC |
| IRF2BP2  | -2.700343 | 0.0444585 | cDC2\_PBMC |
| NXT1     | -2.746845 | 0.0176660 | cDC2\_PBMC |
| DUSP1    | -2.751480 | 0.0011582 | cDC2\_PBMC |
| PLAUR    | -2.887472 | 0.0083979 | cDC2\_PBMC |
| USP12    | -3.022998 | 0.0424591 | cDC2\_PBMC |
| INSIG1   | -3.032639 | 0.0179144 | cDC2\_PBMC |
| BHLHE40  | -3.150120 | 0.0076866 | cDC2\_PBMC |
| ETV3     | -3.208006 | 0.0179144 | cDC2\_PBMC |
| SNHG15   | -3.320195 | 0.0056679 | cDC2\_PBMC |
| IER3     | -3.416281 | 0.0419571 | cDC2\_PBMC |
| CCL3     | -3.454520 | 0.0266178 | cDC2\_PBMC |
| KLF6     | -3.492043 | 0.0002655 | cDC2\_PBMC |
| DUSP2    | -3.739978 | 0.0036086 | cDC2\_PBMC |
| MAFF     | -3.889036 | 0.0465660 | cDC2\_PBMC |
| CXCL8    | -3.988129 | 0.0100140 | cDC2\_PBMC |
| PHLDA2   | -4.331877 | 0.0179144 | cDC2\_PBMC |
| CCL4     | -4.347891 | 0.0085720 | cDC2\_PBMC |
| TNF      | -4.770122 | 0.0323380 | cDC2\_PBMC |
| C15orf48 | -5.118193 | 0.0036086 | cDC2\_PBMC |
| ATF3     | -5.603244 | 0.0000952 | cDC2\_PBMC |
| ID1      | -7.080716 | 0.0323380 | cDC2\_PBMC |
| GEM      | -7.929383 | 0.0176660 | cDC2\_PBMC |
| JUN      | -8.275962 | 0.0000952 | cDC2\_PBMC |
| PHLDA1   | -8.286856 | 0.0108695 | cDC2\_PBMC |

# Volcano plots in Late disease

\[1\] “Results for celltype:” \[1\] “CD4\_Th\_PBMC” \[1\] “Genes
upregulated in disease”

| genes  |    logFC |      padj | cell\_type    |
|:-------|---------:|----------:|:--------------|
| EGR1   | 7.942570 | 0.0000000 | CD4\_Th\_PBMC |
| FOS    | 6.577105 | 0.0000000 | CD4\_Th\_PBMC |
| JUN    | 4.480472 | 0.0000000 | CD4\_Th\_PBMC |
| DUSP1  | 2.767956 | 0.0000193 | CD4\_Th\_PBMC |
| KLF6   | 2.288042 | 0.0006745 | CD4\_Th\_PBMC |
| TUBA1A | 1.834679 | 0.0191376 | CD4\_Th\_PBMC |
| CITED2 | 1.770992 | 0.0312434 | CD4\_Th\_PBMC |
| IER2   | 1.673163 | 0.0191376 | CD4\_Th\_PBMC |

\[1\] “Genes downregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Results for celltype:” \[1\] “CD14\_Monocyte\_PBMC” \[1\] “Genes
upregulated in disease”

| genes      |     logFC |      padj | cell\_type           |
|:-----------|----------:|----------:|:---------------------|
| AKAP12     | 10.895092 | 0.0047969 | CD14\_Monocyte\_PBMC |
| CH25H      |  6.012226 | 0.0426339 | CD14\_Monocyte\_PBMC |
| IL1R2      |  5.402509 | 0.0180157 | CD14\_Monocyte\_PBMC |
| AREG       |  5.252419 | 0.0306685 | CD14\_Monocyte\_PBMC |
| BNC2       |  5.132996 | 0.0073711 | CD14\_Monocyte\_PBMC |
| SERHL2     |  4.299025 | 0.0306685 | CD14\_Monocyte\_PBMC |
| PALD1      |  3.992115 | 0.0055548 | CD14\_Monocyte\_PBMC |
| LRP5       |  3.941946 | 0.0306685 | CD14\_Monocyte\_PBMC |
| LIMS2      |  3.712652 | 0.0377722 | CD14\_Monocyte\_PBMC |
| AC008667.1 |  3.516714 | 0.0306685 | CD14\_Monocyte\_PBMC |
| COL23A1    |  3.232298 | 0.0182490 | CD14\_Monocyte\_PBMC |

\[1\] “Genes downregulated in disease”

| genes    |     logFC |      padj | cell\_type           |
|:---------|----------:|----------:|:---------------------|
| APOBEC3B | -7.445528 | 0.0047969 | CD14\_Monocyte\_PBMC |

\[1\] “Results for celltype:” \[1\] “Platelet\_PBMC” \[1\] “Genes
upregulated in disease”

| genes |    logFC |      padj | cell\_type     |
|:------|---------:|----------:|:---------------|
| RPL10 | 4.556047 | 0.0131716 | Platelet\_PBMC |
| RPLP1 | 4.254077 | 0.0131716 | Platelet\_PBMC |

\[1\] “Genes downregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

# Volcano plots in averagePE

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes    |    logFC |      padj | cell\_type    |
|:---------|---------:|----------:|:--------------|
| TRBV20-1 | 5.475839 | 0.0000022 | CD4\_Th\_PBMC |
| MT-ATP8  | 1.613412 | 0.0113381 | CD4\_Th\_PBMC |

\[1\] “Genes downregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Genes downregulated in disease”

| genes    |     logFC |      padj | cell\_type |
|:---------|----------:|----------:|:-----------|
| ATF3     | -2.494196 | 0.0437563 | cDC2\_PBMC |
| C15orf48 | -3.014870 | 0.0368867 | cDC2\_PBMC |
| JUN      | -4.177545 | 0.0088463 | cDC2\_PBMC |
| PHLDA1   | -4.670877 | 0.0437563 | cDC2\_PBMC |

``` r
sessionInfo()
```

    ## R version 4.1.1 (2021-08-10)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Red Hat Enterprise Linux Server 7.6 (Maipo)
    ## 
    ## Matrix products: default
    ## BLAS:   /apps/R/4.1.1/lib64/R/lib/libRblas.so
    ## LAPACK: /apps/R/4.1.1/lib64/R/lib/libRlapack.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] kableExtra_1.4.0            knitr_1.45                  data.table_1.14.2           EnhancedVolcano_1.13.2     
    ##  [5] ggrepel_0.9.1               scran_1.22.1                edgeR_3.36.0                limma_3.50.3               
    ##  [9] scater_1.22.0               scuttle_1.4.0               SingleCellExperiment_1.16.0 SummarizedExperiment_1.24.0
    ## [13] Biobase_2.54.0              GenomicRanges_1.46.1        GenomeInfoDb_1.30.0         IRanges_2.28.0             
    ## [17] S4Vectors_0.32.3            BiocGenerics_0.40.0         MatrixGenerics_1.6.0        matrixStats_0.61.0         
    ## [21] forcats_0.5.1               purrr_1.0.2                 readr_2.1.1                 tidyr_1.1.4                
    ## [25] tibble_3.1.6                tidyverse_1.3.1             stringr_1.4.0               cowplot_1.1.1              
    ## [29] ggplot2_3.3.5               SeuratObject_4.1.3          Seurat_4.2.1                dplyr_1.0.7                
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] utf8_1.2.2                spatstat.explore_3.0-5    reticulate_1.24           tidyselect_1.1.1         
    ##   [5] htmlwidgets_1.5.4         grid_4.1.1                BiocParallel_1.28.3       Rtsne_0.15               
    ##   [9] munsell_0.5.0             ScaledMatrix_1.2.0        codetools_0.2-18          ica_1.0-2                
    ##  [13] statmod_1.4.36            future_1.23.0             miniUI_0.1.1.1            withr_2.5.0              
    ##  [17] spatstat.random_3.0-1     colorspace_2.0-2          progressr_0.10.0          highr_0.9                
    ##  [21] rstudioapi_0.13           ROCR_1.0-11               tensor_1.5                listenv_0.8.0            
    ##  [25] labeling_0.4.2            GenomeInfoDbData_1.2.7    polyclip_1.10-0           farver_2.1.0             
    ##  [29] parallelly_1.30.0         vctrs_0.6.5               generics_0.1.1            xfun_0.41                
    ##  [33] R6_2.5.1                  ggbeeswarm_0.6.0          rsvd_1.0.5                locfit_1.5-9.4           
    ##  [37] bitops_1.0-7              spatstat.utils_3.0-1      DelayedArray_0.20.0       assertthat_0.2.1         
    ##  [41] promises_1.2.0.1          scales_1.1.1              beeswarm_0.4.0            gtable_0.3.0             
    ##  [45] beachmat_2.10.0           globals_0.14.0            goftest_1.2-3             rlang_1.1.1              
    ##  [49] systemfonts_1.0.3         splines_4.1.1             lazyeval_0.2.2            spatstat.geom_3.0-3      
    ##  [53] broom_0.7.11              yaml_2.2.2                reshape2_1.4.4            abind_1.4-5              
    ##  [57] modelr_0.1.8              backports_1.4.1           httpuv_1.6.5              tools_4.1.1              
    ##  [61] ellipsis_0.3.2            RColorBrewer_1.1-2        ggridges_0.5.3            Rcpp_1.0.8               
    ##  [65] plyr_1.8.6                sparseMatrixStats_1.6.0   zlibbioc_1.40.0           RCurl_1.98-1.5           
    ##  [69] deldir_1.0-6              pbapply_1.5-0             viridis_0.6.2             zoo_1.8-9                
    ##  [73] haven_2.4.3               cluster_2.1.2             fs_1.5.2                  magrittr_2.0.1           
    ##  [77] scattermore_0.7           lmtest_0.9-39             reprex_2.0.1              RANN_2.6.1               
    ##  [81] fitdistrplus_1.1-6        hms_1.1.1                 patchwork_1.1.1           mime_0.12                
    ##  [85] evaluate_0.23             xtable_1.8-4              readxl_1.3.1              gridExtra_2.3            
    ##  [89] compiler_4.1.1            KernSmooth_2.23-20        crayon_1.4.2              htmltools_0.5.8.1        
    ##  [93] later_1.3.0               tzdb_0.2.0                lubridate_1.8.0           DBI_1.1.2                
    ##  [97] dbplyr_2.1.1              MASS_7.3-55               Matrix_1.5-4.1            cli_3.6.1                
    ## [101] metapod_1.2.0             parallel_4.1.1            igraph_1.2.11             pkgconfig_2.0.3          
    ## [105] sp_1.5-1                  plotly_4.10.0             spatstat.sparse_3.0-0     xml2_1.3.3               
    ## [109] svglite_2.1.0             vipor_0.4.5               dqrng_0.3.0               XVector_0.34.0           
    ## [113] rvest_1.0.2               digest_0.6.29             sctransform_0.3.5         RcppAnnoy_0.0.19         
    ## [117] spatstat.data_3.0-0       rmarkdown_2.25            cellranger_1.1.0          leiden_0.3.9             
    ## [121] uwot_0.1.14               DelayedMatrixStats_1.16.0 shiny_1.7.1               lifecycle_1.0.4          
    ## [125] nlme_3.1-155              jsonlite_1.7.3            BiocNeighbors_1.12.0      viridisLite_0.4.0        
    ## [129] fansi_1.0.2               pillar_1.6.5              lattice_0.20-45           fastmap_1.1.1            
    ## [133] httr_1.4.2                survival_3.2-13           glue_1.6.1                png_0.1-7                
    ## [137] bluster_1.4.0             stringi_1.7.6             BiocSingular_1.10.0       irlba_2.3.5              
    ## [141] future.apply_1.8.1
