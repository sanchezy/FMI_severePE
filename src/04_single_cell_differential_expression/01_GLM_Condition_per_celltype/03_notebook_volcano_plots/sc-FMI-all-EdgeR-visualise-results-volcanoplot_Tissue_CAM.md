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

> Tissue: CAM

``` r
# This notebook is for tissue and annotation level:
tissue= "CAM" # "CAM", "PBMC", "Myometrium"
annotation_level <- "CellTypeManual.l3"

# List of cell types
cellgroups$CellTypeManual.l3 <- paste0(cellgroups$CellTypeManual.l3,"_",tissue)
```

# Volcano plots in Early disease

\[1\] “Results for celltype:” \[1\] “Naive\_B\_CAM” \[1\] “Genes
upregulated in disease”

| genes    |    logFC |      padj | cell\_type    |
|:---------|---------:|----------:|:--------------|
| IGHV3-48 | 6.424203 | 0.0274415 | Naive\_B\_CAM |
| MT-ND6   | 3.280482 | 0.0035052 | Naive\_B\_CAM |

\[1\] “Genes downregulated in disease”

| genes    |     logFC |      padj | cell\_type    |
|:---------|----------:|----------:|:--------------|
| IGLV1-40 | -3.886220 | 0.0274415 | Naive\_B\_CAM |
| IGLV1-51 | -4.176320 | 0.0274415 | Naive\_B\_CAM |
| IGHV2-5  | -6.023976 | 0.0170950 | Naive\_B\_CAM |
| IGLV7-46 | -6.267133 | 0.0274415 | Naive\_B\_CAM |
| IGKV1-8  | -7.344755 | 0.0020166 | Naive\_B\_CAM |

\[1\] “Results for celltype:” \[1\] “Memory\_B\_CAM” \[1\] “Genes
upregulated in disease”

| genes    |    logFC |      padj | cell\_type     |
|:---------|---------:|----------:|:---------------|
| IGKV1-27 | 6.893779 | 0.0378087 | Memory\_B\_CAM |
| MT-ND6   | 4.367017 | 0.0010469 | Memory\_B\_CAM |
| GAS5     | 2.503783 | 0.0028277 | Memory\_B\_CAM |
| RPL37A   | 2.113222 | 0.0257637 | Memory\_B\_CAM |
| MT-ATP8  | 2.112808 | 0.0177264 | Memory\_B\_CAM |
| ARHGAP24 | 1.831176 | 0.0215038 | Memory\_B\_CAM |
| BLK      | 1.593747 | 0.0177264 | Memory\_B\_CAM |
| TPT1     | 1.584365 | 0.0010469 | Memory\_B\_CAM |
| RPL21    | 1.476183 | 0.0215038 | Memory\_B\_CAM |
| RPL37    | 1.130793 | 0.0177264 | Memory\_B\_CAM |

\[1\] “Genes downregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Results for celltype:” \[1\] “Fetal\_CD4\_T\_Naive\_CM\_CAM”
\[1\] “Genes upregulated in disease”

| genes   |    logFC |      padj | cell\_type                    |
|:--------|---------:|----------:|:------------------------------|
| TRBV5-1 | 9.019514 | 0.0242328 | Fetal\_CD4\_T\_Naive\_CM\_CAM |
| GAS5    | 3.255586 | 0.0000820 | Fetal\_CD4\_T\_Naive\_CM\_CAM |
| MT-CYB  | 3.011195 | 0.0001561 | Fetal\_CD4\_T\_Naive\_CM\_CAM |
| ERAP2   | 2.969298 | 0.0119681 | Fetal\_CD4\_T\_Naive\_CM\_CAM |
| MT-ATP8 | 2.511039 | 0.0022775 | Fetal\_CD4\_T\_Naive\_CM\_CAM |
| MT-CO3  | 2.404215 | 0.0173619 | Fetal\_CD4\_T\_Naive\_CM\_CAM |
| MT-CO1  | 2.361523 | 0.0119681 | Fetal\_CD4\_T\_Naive\_CM\_CAM |
| MT-ND3  | 2.316882 | 0.0107896 | Fetal\_CD4\_T\_Naive\_CM\_CAM |
| SNHG6   | 2.255587 | 0.0064601 | Fetal\_CD4\_T\_Naive\_CM\_CAM |
| RPL37A  | 2.149755 | 0.0173619 | Fetal\_CD4\_T\_Naive\_CM\_CAM |
| HLA-B   | 1.979276 | 0.0242328 | Fetal\_CD4\_T\_Naive\_CM\_CAM |
| HLA-A   | 1.736704 | 0.0262065 | Fetal\_CD4\_T\_Naive\_CM\_CAM |

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type                    |
|:------|----------:|----------:|:------------------------------|
| RGS1  | -3.472491 | 0.0409509 | Fetal\_CD4\_T\_Naive\_CM\_CAM |

\[1\] “Results for celltype:” \[1\] “CD4\_T\_Naive\_CM\_CAM” \[1\]
“Genes upregulated in disease”

| genes      |    logFC |      padj | cell\_type             |
|:-----------|---------:|----------:|:-----------------------|
| TRBV20-1   | 8.083717 | 0.0042356 | CD4\_T\_Naive\_CM\_CAM |
| TRBV5-1    | 7.665877 | 0.0276901 | CD4\_T\_Naive\_CM\_CAM |
| TRBV6-5    | 7.635681 | 0.0004865 | CD4\_T\_Naive\_CM\_CAM |
| TRBV12-3   | 7.306769 | 0.0096516 | CD4\_T\_Naive\_CM\_CAM |
| TRBV6-1    | 7.176391 | 0.0138419 | CD4\_T\_Naive\_CM\_CAM |
| FMN1       | 7.159738 | 0.0043999 | CD4\_T\_Naive\_CM\_CAM |
| TRBV6-2    | 7.144818 | 0.0037884 | CD4\_T\_Naive\_CM\_CAM |
| TRAV13-1   | 6.691626 | 0.0426001 | CD4\_T\_Naive\_CM\_CAM |
| TRAV12-1   | 6.139680 | 0.0276901 | CD4\_T\_Naive\_CM\_CAM |
| TRBV7-2    | 5.896750 | 0.0402911 | CD4\_T\_Naive\_CM\_CAM |
| TRBV3-1    | 5.052441 | 0.0402911 | CD4\_T\_Naive\_CM\_CAM |
| TMEM8A     | 3.638339 | 0.0051096 | CD4\_T\_Naive\_CM\_CAM |
| CD59       | 3.627907 | 0.0487492 | CD4\_T\_Naive\_CM\_CAM |
| MT-CYB     | 3.317381 | 0.0000105 | CD4\_T\_Naive\_CM\_CAM |
| MT-ND6     | 3.116637 | 0.0000925 | CD4\_T\_Naive\_CM\_CAM |
| MT-ND3     | 3.094021 | 0.0000194 | CD4\_T\_Naive\_CM\_CAM |
| MT-ATP8    | 2.645489 | 0.0001906 | CD4\_T\_Naive\_CM\_CAM |
| MT-ND5     | 2.534688 | 0.0008817 | CD4\_T\_Naive\_CM\_CAM |
| MT-CO3     | 2.409936 | 0.0179660 | CD4\_T\_Naive\_CM\_CAM |
| RPL37A     | 2.281672 | 0.0003826 | CD4\_T\_Naive\_CM\_CAM |
| MT-CO1     | 2.161608 | 0.0276901 | CD4\_T\_Naive\_CM\_CAM |
| MT-ND1     | 2.157840 | 0.0096516 | CD4\_T\_Naive\_CM\_CAM |
| GAS5       | 2.129428 | 0.0040726 | CD4\_T\_Naive\_CM\_CAM |
| RPS16      | 2.045756 | 0.0013472 | CD4\_T\_Naive\_CM\_CAM |
| ARHGEF18   | 2.032930 | 0.0381464 | CD4\_T\_Naive\_CM\_CAM |
| AC079793.1 | 1.998768 | 0.0276901 | CD4\_T\_Naive\_CM\_CAM |
| RNPEPL1    | 1.876785 | 0.0276901 | CD4\_T\_Naive\_CM\_CAM |
| MT-ATP6    | 1.847675 | 0.0276901 | CD4\_T\_Naive\_CM\_CAM |
| RPSA       | 1.827872 | 0.0027048 | CD4\_T\_Naive\_CM\_CAM |
| MT-ND4L    | 1.810929 | 0.0276901 | CD4\_T\_Naive\_CM\_CAM |
| SNHG6      | 1.707365 | 0.0197037 | CD4\_T\_Naive\_CM\_CAM |
| TPT1       | 1.701537 | 0.0018825 | CD4\_T\_Naive\_CM\_CAM |
| RPL21      | 1.592928 | 0.0047584 | CD4\_T\_Naive\_CM\_CAM |
| ANXA1      | 1.538982 | 0.0276901 | CD4\_T\_Naive\_CM\_CAM |
| RPL4       | 1.429282 | 0.0211445 | CD4\_T\_Naive\_CM\_CAM |
| CD52       | 1.394957 | 0.0447595 | CD4\_T\_Naive\_CM\_CAM |
| RPL37      | 1.336834 | 0.0413431 | CD4\_T\_Naive\_CM\_CAM |
| RPS8       | 1.326034 | 0.0381464 | CD4\_T\_Naive\_CM\_CAM |

\[1\] “Genes downregulated in disease”

| genes      |     logFC |      padj | cell\_type             |
|:-----------|----------:|----------:|:-----------------------|
| LINC00910  | -1.406553 | 0.0434022 | CD4\_T\_Naive\_CM\_CAM |
| AC138123.1 | -1.538551 | 0.0179660 | CD4\_T\_Naive\_CM\_CAM |
| HSPA1A     | -1.993267 | 0.0006231 | CD4\_T\_Naive\_CM\_CAM |
| CARD9      | -2.196552 | 0.0305331 | CD4\_T\_Naive\_CM\_CAM |
| HES4       | -4.280272 | 0.0013472 | CD4\_T\_Naive\_CM\_CAM |

\[1\] “Results for celltype:” \[1\] “CD4\_Th\_CAM” \[1\] “Genes
upregulated in disease”

| genes       |    logFC |      padj | cell\_type   |
|:------------|---------:|----------:|:-------------|
| TRBV20-1    | 8.024264 | 0.0038760 | CD4\_Th\_CAM |
| TRBV12-4    | 7.561271 | 0.0004744 | CD4\_Th\_CAM |
| TRBV7-9     | 7.215779 | 0.0004009 | CD4\_Th\_CAM |
| TRAV38-2DV8 | 7.050892 | 0.0095032 | CD4\_Th\_CAM |
| TRBV4-1     | 6.981139 | 0.0026552 | CD4\_Th\_CAM |
| TRBV5-1     | 6.688028 | 0.0023005 | CD4\_Th\_CAM |
| TRBV6-5     | 6.497534 | 0.0023356 | CD4\_Th\_CAM |
| TRAV8-4     | 6.433082 | 0.0021505 | CD4\_Th\_CAM |
| TRBV5-6     | 6.154810 | 0.0101137 | CD4\_Th\_CAM |
| TRAV29DV5   | 6.118530 | 0.0125617 | CD4\_Th\_CAM |
| TRBV6-1     | 5.957391 | 0.0230848 | CD4\_Th\_CAM |
| TRAV12-3    | 5.945035 | 0.0283766 | CD4\_Th\_CAM |
| TRAV9-2     | 5.940004 | 0.0165268 | CD4\_Th\_CAM |
| TRBV27      | 5.895121 | 0.0303112 | CD4\_Th\_CAM |
| TRAV27      | 5.574479 | 0.0229111 | CD4\_Th\_CAM |
| TRAV12-2    | 4.995929 | 0.0389919 | CD4\_Th\_CAM |
| MT-CYB      | 3.458458 | 0.0039819 | CD4\_Th\_CAM |
| MT-ND3      | 3.038751 | 0.0095032 | CD4\_Th\_CAM |
| MT-ND6      | 2.726901 | 0.0010000 | CD4\_Th\_CAM |
| MT-ATP8     | 2.626867 | 0.0011535 | CD4\_Th\_CAM |
| ERAP2       | 2.612446 | 0.0002939 | CD4\_Th\_CAM |
| MT-ND5      | 2.551548 | 0.0033988 | CD4\_Th\_CAM |
| AL136456.1  | 2.220897 | 0.0192371 | CD4\_Th\_CAM |
| GAS5        | 2.131902 | 0.0001606 | CD4\_Th\_CAM |
| PIK3R5      | 2.068889 | 0.0025676 | CD4\_Th\_CAM |
| AP001011.1  | 2.057456 | 0.0169851 | CD4\_Th\_CAM |
| FAM13A      | 1.981070 | 0.0023005 | CD4\_Th\_CAM |
| PLEKHB2     | 1.942460 | 0.0378189 | CD4\_Th\_CAM |
| AL627171.2  | 1.830068 | 0.0389919 | CD4\_Th\_CAM |
| MT-ND4L     | 1.810013 | 0.0303112 | CD4\_Th\_CAM |
| RPL37A      | 1.761748 | 0.0023005 | CD4\_Th\_CAM |
| IVNS1ABP    | 1.743659 | 0.0229111 | CD4\_Th\_CAM |
| SNHG6       | 1.718436 | 0.0025676 | CD4\_Th\_CAM |
| ATXN7L1     | 1.667425 | 0.0229040 | CD4\_Th\_CAM |
| LINC01619   | 1.641932 | 0.0061722 | CD4\_Th\_CAM |
| SCFD1       | 1.627456 | 0.0122214 | CD4\_Th\_CAM |
| DOCK11      | 1.605133 | 0.0101137 | CD4\_Th\_CAM |
| UBQLN2      | 1.593355 | 0.0051233 | CD4\_Th\_CAM |
| PATJ        | 1.588921 | 0.0055699 | CD4\_Th\_CAM |
| TTC14       | 1.571730 | 0.0331547 | CD4\_Th\_CAM |
| RFX3        | 1.533545 | 0.0258361 | CD4\_Th\_CAM |
| MXI1        | 1.434907 | 0.0253165 | CD4\_Th\_CAM |
| RPSA        | 1.434770 | 0.0043671 | CD4\_Th\_CAM |
| RPS16       | 1.347996 | 0.0224661 | CD4\_Th\_CAM |
| TXK         | 1.323351 | 0.0307420 | CD4\_Th\_CAM |
| OXA1L       | 1.309265 | 0.0320470 | CD4\_Th\_CAM |
| OGT         | 1.292450 | 0.0165268 | CD4\_Th\_CAM |
| TASOR       | 1.274466 | 0.0493291 | CD4\_Th\_CAM |
| TAF15       | 1.221354 | 0.0165268 | CD4\_Th\_CAM |
| TPT1        | 1.198009 | 0.0185339 | CD4\_Th\_CAM |
| CD52        | 1.142037 | 0.0434475 | CD4\_Th\_CAM |
| EIF3E       | 1.136896 | 0.0216335 | CD4\_Th\_CAM |
| RPL21       | 1.135717 | 0.0260889 | CD4\_Th\_CAM |
| RPL4        | 1.067056 | 0.0441168 | CD4\_Th\_CAM |
| RAD21       | 1.066416 | 0.0230848 | CD4\_Th\_CAM |
| NOP53       | 1.050720 | 0.0482185 | CD4\_Th\_CAM |

\[1\] “Genes downregulated in disease”

| genes      |      logFC |      padj | cell\_type   |
|:-----------|-----------:|----------:|:-------------|
| PPIB       | -0.9160104 | 0.0434168 | CD4\_Th\_CAM |
| HSPE1      | -0.9257682 | 0.0383412 | CD4\_Th\_CAM |
| MGST3      | -0.9738105 | 0.0493291 | CD4\_Th\_CAM |
| BANF1      | -0.9907895 | 0.0378189 | CD4\_Th\_CAM |
| VEGFB      | -1.0460950 | 0.0377336 | CD4\_Th\_CAM |
| DENR       | -1.0685341 | 0.0329706 | CD4\_Th\_CAM |
| ODC1       | -1.1771769 | 0.0300756 | CD4\_Th\_CAM |
| PLGRKT     | -1.1875923 | 0.0356945 | CD4\_Th\_CAM |
| TNF        | -1.2397834 | 0.0320470 | CD4\_Th\_CAM |
| CKLF       | -1.2788987 | 0.0434475 | CD4\_Th\_CAM |
| STK17A     | -1.2985201 | 0.0041330 | CD4\_Th\_CAM |
| UPF2       | -1.3486926 | 0.0025676 | CD4\_Th\_CAM |
| CYB5R3     | -1.3526323 | 0.0331547 | CD4\_Th\_CAM |
| DNPH1      | -1.3847425 | 0.0258361 | CD4\_Th\_CAM |
| DUSP2      | -1.4623293 | 0.0210005 | CD4\_Th\_CAM |
| EFHD2      | -1.6006444 | 0.0260889 | CD4\_Th\_CAM |
| CAPG       | -1.6282504 | 0.0260889 | CD4\_Th\_CAM |
| JOSD2      | -1.6471055 | 0.0101598 | CD4\_Th\_CAM |
| AC093512.2 | -1.7326363 | 0.0129978 | CD4\_Th\_CAM |
| AC138123.1 | -1.7550241 | 0.0002588 | CD4\_Th\_CAM |
| FKBP2      | -2.0143390 | 0.0021505 | CD4\_Th\_CAM |
| HIST2H2AA4 | -2.1254850 | 0.0173761 | CD4\_Th\_CAM |
| HSPA1A     | -2.1562206 | 0.0023005 | CD4\_Th\_CAM |
| Z93241.1   | -2.1723439 | 0.0211208 | CD4\_Th\_CAM |
| GZMK       | -2.5551369 | 0.0054950 | CD4\_Th\_CAM |
| OASL       | -2.5676747 | 0.0243750 | CD4\_Th\_CAM |
| HSPA1B     | -2.5828487 | 0.0267742 | CD4\_Th\_CAM |
| CCL5       | -3.5789771 | 0.0000001 | CD4\_Th\_CAM |
| IFNG       | -4.7513904 | 0.0131315 | CD4\_Th\_CAM |
| CCL4       | -5.8061097 | 0.0000063 | CD4\_Th\_CAM |
| HSPA6      | -6.9222406 | 0.0023005 | CD4\_Th\_CAM |

\[1\] “Results for celltype:” \[1\] “FoxP3-Treg\_CAM” \[1\] “Genes
upregulated in disease”

| genes       |    logFC |      padj | cell\_type      |
|:------------|---------:|----------:|:----------------|
| TRBV7-9     | 8.519951 | 0.0033197 | FoxP3-Treg\_CAM |
| TRBV20-1    | 8.513363 | 0.0006824 | FoxP3-Treg\_CAM |
| TRBV5-1     | 8.065942 | 0.0006824 | FoxP3-Treg\_CAM |
| TRBV4-2     | 8.041803 | 0.0358578 | FoxP3-Treg\_CAM |
| TRAV9-2     | 7.990852 | 0.0016785 | FoxP3-Treg\_CAM |
| TRBV4-1     | 7.785614 | 0.0101535 | FoxP3-Treg\_CAM |
| TRAV13-1    | 7.688351 | 0.0032515 | FoxP3-Treg\_CAM |
| TRBV6-5     | 7.357245 | 0.0016785 | FoxP3-Treg\_CAM |
| TRBV29-1    | 7.244264 | 0.0235662 | FoxP3-Treg\_CAM |
| TRAV38-2DV8 | 6.988898 | 0.0097406 | FoxP3-Treg\_CAM |
| TRAV12-2    | 6.895602 | 0.0050748 | FoxP3-Treg\_CAM |
| TRAV4       | 6.656578 | 0.0275683 | FoxP3-Treg\_CAM |
| TRBV10-3    | 6.445782 | 0.0310320 | FoxP3-Treg\_CAM |
| TRAV17      | 6.434997 | 0.0097406 | FoxP3-Treg\_CAM |
| CADM1       | 3.892394 | 0.0101535 | FoxP3-Treg\_CAM |
| TRBV19      | 3.886351 | 0.0142494 | FoxP3-Treg\_CAM |
| MT-CYB      | 3.659779 | 0.0091229 | FoxP3-Treg\_CAM |
| MT-ND6      | 2.853187 | 0.0097406 | FoxP3-Treg\_CAM |
| ERAP2       | 2.782459 | 0.0369165 | FoxP3-Treg\_CAM |
| GAS5        | 2.703480 | 0.0004630 | FoxP3-Treg\_CAM |
| MT-ATP8     | 2.571936 | 0.0453228 | FoxP3-Treg\_CAM |
| MT-ND5      | 2.412526 | 0.0142494 | FoxP3-Treg\_CAM |
| BCAS3       | 2.222741 | 0.0201298 | FoxP3-Treg\_CAM |
| PAN3        | 2.168128 | 0.0324308 | FoxP3-Treg\_CAM |
| CXCR4       | 2.118622 | 0.0017546 | FoxP3-Treg\_CAM |
| SNHG6       | 1.425414 | 0.0441392 | FoxP3-Treg\_CAM |

\[1\] “Genes downregulated in disease”

| genes      |     logFC |      padj | cell\_type      |
|:-----------|----------:|----------:|:----------------|
| PDCD1      | -1.820582 | 0.0201298 | FoxP3-Treg\_CAM |
| FKBP2      | -1.844999 | 0.0369165 | FoxP3-Treg\_CAM |
| CSF1       | -2.398870 | 0.0107343 | FoxP3-Treg\_CAM |
| AC093512.2 | -2.428793 | 0.0006824 | FoxP3-Treg\_CAM |
| ARG2       | -2.545003 | 0.0078673 | FoxP3-Treg\_CAM |
| ACTG2      | -3.062928 | 0.0177289 | FoxP3-Treg\_CAM |

\[1\] “Results for celltype:” \[1\] “CD4\_TEM\_CAM” \[1\] “Genes
upregulated in disease”

| genes    |     logFC |      padj | cell\_type    |
|:---------|----------:|----------:|:--------------|
| TRBV20-1 | 9.7813311 | 0.0151710 | CD4\_TEM\_CAM |
| TRBV5-1  | 9.3778707 | 0.0046748 | CD4\_TEM\_CAM |
| TRAV13-1 | 8.2995600 | 0.0234167 | CD4\_TEM\_CAM |
| TRAV8-6  | 7.9894483 | 0.0433877 | CD4\_TEM\_CAM |
| TRAV9-2  | 7.7938684 | 0.0386766 | CD4\_TEM\_CAM |
| TRBV6-5  | 7.7541842 | 0.0221088 | CD4\_TEM\_CAM |
| TRBV12-3 | 7.6958631 | 0.0486064 | CD4\_TEM\_CAM |
| TRAV8-4  | 7.4815541 | 0.0108630 | CD4\_TEM\_CAM |
| TRAV17   | 7.4029458 | 0.0341064 | CD4\_TEM\_CAM |
| TRBV27   | 7.3596710 | 0.0062579 | CD4\_TEM\_CAM |
| TRAV19   | 6.9003405 | 0.0128910 | CD4\_TEM\_CAM |
| TRAV12-1 | 6.3535897 | 0.0093948 | CD4\_TEM\_CAM |
| TRBV7-3  | 5.1167270 | 0.0339148 | CD4\_TEM\_CAM |
| MMS22L   | 3.5898643 | 0.0433877 | CD4\_TEM\_CAM |
| MT-ND6   | 3.2702472 | 0.0021295 | CD4\_TEM\_CAM |
| MT-ATP8  | 2.8790387 | 0.0221088 | CD4\_TEM\_CAM |
| MT-ND5   | 2.7225762 | 0.0128910 | CD4\_TEM\_CAM |
| ERAP2    | 2.5806387 | 0.0002556 | CD4\_TEM\_CAM |
| GAS5     | 2.1695771 | 0.0110287 | CD4\_TEM\_CAM |
| WDR36    | 2.1519926 | 0.0042052 | CD4\_TEM\_CAM |
| MAP2K5   | 1.7773038 | 0.0339148 | CD4\_TEM\_CAM |
| CXCR4    | 1.4347110 | 0.0267265 | CD4\_TEM\_CAM |
| RALGAPA2 | 1.4054795 | 0.0433877 | CD4\_TEM\_CAM |
| TMEM259  | 1.3892024 | 0.0042052 | CD4\_TEM\_CAM |
| RPL21    | 1.1844972 | 0.0412500 | CD4\_TEM\_CAM |
| LMF2     | 1.1507802 | 0.0369137 | CD4\_TEM\_CAM |
| CBLB     | 0.9457251 | 0.0221088 | CD4\_TEM\_CAM |

\[1\] “Genes downregulated in disease”

| genes      |      logFC |      padj | cell\_type    |
|:-----------|-----------:|----------:|:--------------|
| STK17A     | -0.8045210 | 0.0412500 | CD4\_TEM\_CAM |
| CMTM3      | -0.8868624 | 0.0184453 | CD4\_TEM\_CAM |
| DNPH1      | -1.0148103 | 0.0221088 | CD4\_TEM\_CAM |
| LDHA       | -1.0452154 | 0.0477698 | CD4\_TEM\_CAM |
| C1D        | -1.0471378 | 0.0290670 | CD4\_TEM\_CAM |
| CITED4     | -1.2750919 | 0.0002556 | CD4\_TEM\_CAM |
| HIST2H2AA4 | -1.6996177 | 0.0221088 | CD4\_TEM\_CAM |
| RGS16      | -1.9279127 | 0.0166573 | CD4\_TEM\_CAM |
| HLA-DRB5   | -2.1335109 | 0.0108630 | CD4\_TEM\_CAM |
| MXRA7      | -2.2102644 | 0.0002556 | CD4\_TEM\_CAM |
| CCL4L2     | -4.2796411 | 0.0002556 | CD4\_TEM\_CAM |

\[1\] “Results for celltype:” \[1\] “NKT\_CAM” \[1\] “Genes upregulated
in disease”

| genes    |    logFC |      padj | cell\_type |
|:---------|---------:|----------:|:-----------|
| ZSWIM3   | 7.278566 | 0.0422447 | NKT\_CAM   |
| FILIP1L  | 7.248613 | 0.0422447 | NKT\_CAM   |
| YBX3     | 5.579479 | 0.0044750 | NKT\_CAM   |
| RANBP17  | 3.647258 | 0.0093256 | NKT\_CAM   |
| MT-ND6   | 3.195500 | 0.0093256 | NKT\_CAM   |
| L3MBTL4  | 1.977897 | 0.0422447 | NKT\_CAM   |
| CPB2-AS1 | 1.692273 | 0.0422447 | NKT\_CAM   |
| VPS8     | 1.274299 | 0.0422447 | NKT\_CAM   |

\[1\] “Genes downregulated in disease”

| genes  |     logFC |      padj | cell\_type |
|:-------|----------:|----------:|:-----------|
| FBXO34 | -1.086351 | 0.0414997 | NKT\_CAM   |
| DNPH1  | -1.413652 | 0.0481016 | NKT\_CAM   |
| TYMS   | -1.573214 | 0.0093256 | NKT\_CAM   |
| CARD9  | -2.072049 | 0.0414997 | NKT\_CAM   |
| CD160  | -2.672589 | 0.0422447 | NKT\_CAM   |

\[1\] “Results for celltype:” \[1\] “CD14\_Monocyte\_CAM” \[1\] “Genes
upregulated in disease”

| genes      |    logFC |      padj | cell\_type          |
|:-----------|---------:|----------:|:--------------------|
| DIRC3      | 8.021873 | 0.0159290 | CD14\_Monocyte\_CAM |
| MFGE8      | 6.246822 | 0.0483553 | CD14\_Monocyte\_CAM |
| HFM1       | 6.059505 | 0.0159290 | CD14\_Monocyte\_CAM |
| ADAMTS2    | 4.740404 | 0.0114325 | CD14\_Monocyte\_CAM |
| VSIG4      | 3.764858 | 0.0024200 | CD14\_Monocyte\_CAM |
| MT-ND6     | 3.730815 | 0.0000012 | CD14\_Monocyte\_CAM |
| DST        | 3.063744 | 0.0114325 | CD14\_Monocyte\_CAM |
| TANC2      | 2.946984 | 0.0176967 | CD14\_Monocyte\_CAM |
| AL627171.2 | 2.866205 | 0.0159290 | CD14\_Monocyte\_CAM |
| AC068587.4 | 2.694221 | 0.0247609 | CD14\_Monocyte\_CAM |
| ERAP2      | 2.569637 | 0.0075349 | CD14\_Monocyte\_CAM |
| LARS2      | 2.182789 | 0.0483308 | CD14\_Monocyte\_CAM |
| GAS5       | 2.180559 | 0.0058726 | CD14\_Monocyte\_CAM |
| MMP17      | 2.067694 | 0.0486678 | CD14\_Monocyte\_CAM |
| TMEM245    | 2.042189 | 0.0144527 | CD14\_Monocyte\_CAM |
| IGF2BP2    | 1.941379 | 0.0247609 | CD14\_Monocyte\_CAM |
| TANGO6     | 1.941166 | 0.0073743 | CD14\_Monocyte\_CAM |
| BTRC       | 1.925525 | 0.0159290 | CD14\_Monocyte\_CAM |
| ZNF700     | 1.908010 | 0.0356226 | CD14\_Monocyte\_CAM |
| MS4A4E     | 1.869371 | 0.0247609 | CD14\_Monocyte\_CAM |
| INPP5A     | 1.748435 | 0.0331589 | CD14\_Monocyte\_CAM |
| MCU        | 1.713410 | 0.0303173 | CD14\_Monocyte\_CAM |
| SNTB1      | 1.626240 | 0.0114325 | CD14\_Monocyte\_CAM |
| FAR2       | 1.544232 | 0.0462351 | CD14\_Monocyte\_CAM |
| FBXL20     | 1.524908 | 0.0390754 | CD14\_Monocyte\_CAM |
| ZNF609     | 1.481643 | 0.0397666 | CD14\_Monocyte\_CAM |
| PDSS2      | 1.430096 | 0.0159290 | CD14\_Monocyte\_CAM |
| DOCK10     | 1.427723 | 0.0486678 | CD14\_Monocyte\_CAM |
| SLMAP      | 1.402149 | 0.0167819 | CD14\_Monocyte\_CAM |
| SYNE3      | 1.374824 | 0.0435676 | CD14\_Monocyte\_CAM |
| PPHLN1     | 1.360722 | 0.0282226 | CD14\_Monocyte\_CAM |
| ATXN7L1    | 1.344813 | 0.0159290 | CD14\_Monocyte\_CAM |
| FTX        | 1.267731 | 0.0462351 | CD14\_Monocyte\_CAM |
| GBF1       | 1.192611 | 0.0282226 | CD14\_Monocyte\_CAM |
| ADAM17     | 1.137192 | 0.0282226 | CD14\_Monocyte\_CAM |

\[1\] “Genes downregulated in disease”

| genes      |     logFC |      padj | cell\_type          |
|:-----------|----------:|----------:|:--------------------|
| SRI        | -1.120512 | 0.0404813 | CD14\_Monocyte\_CAM |
| PHF5A      | -1.220969 | 0.0267045 | CD14\_Monocyte\_CAM |
| THOC7      | -1.222439 | 0.0153057 | CD14\_Monocyte\_CAM |
| SMCO4      | -1.244651 | 0.0483308 | CD14\_Monocyte\_CAM |
| COPS2      | -1.275189 | 0.0247609 | CD14\_Monocyte\_CAM |
| CNN2       | -1.297157 | 0.0462351 | CD14\_Monocyte\_CAM |
| AL096870.8 | -1.298406 | 0.0282226 | CD14\_Monocyte\_CAM |
| ADI1       | -1.299678 | 0.0282226 | CD14\_Monocyte\_CAM |
| PFDN1      | -1.304898 | 0.0394706 | CD14\_Monocyte\_CAM |
| SNRPD3     | -1.313967 | 0.0176967 | CD14\_Monocyte\_CAM |
| TMCO1      | -1.369400 | 0.0176967 | CD14\_Monocyte\_CAM |
| ECE1       | -1.375850 | 0.0404813 | CD14\_Monocyte\_CAM |
| FKBP2      | -1.382680 | 0.0114309 | CD14\_Monocyte\_CAM |
| CTSC       | -1.394297 | 0.0056719 | CD14\_Monocyte\_CAM |
| TMEM35B    | -1.401680 | 0.0159290 | CD14\_Monocyte\_CAM |
| EXOSC3     | -1.402056 | 0.0491833 | CD14\_Monocyte\_CAM |
| CDK2AP1    | -1.434562 | 0.0114309 | CD14\_Monocyte\_CAM |
| SPINT2     | -1.465088 | 0.0073743 | CD14\_Monocyte\_CAM |
| IL1RN      | -1.474139 | 0.0278135 | CD14\_Monocyte\_CAM |
| SIGLEC7    | -1.475629 | 0.0159290 | CD14\_Monocyte\_CAM |
| OTUD1      | -1.490628 | 0.0414825 | CD14\_Monocyte\_CAM |
| AC005280.2 | -1.494312 | 0.0114325 | CD14\_Monocyte\_CAM |
| BATF       | -1.519135 | 0.0204499 | CD14\_Monocyte\_CAM |
| MRPS18C    | -1.551395 | 0.0006030 | CD14\_Monocyte\_CAM |
| CROCC      | -1.567637 | 0.0114309 | CD14\_Monocyte\_CAM |
| RAB24      | -1.648310 | 0.0247609 | CD14\_Monocyte\_CAM |
| AL021155.5 | -1.697045 | 0.0278135 | CD14\_Monocyte\_CAM |
| TMEM107    | -1.745350 | 0.0124810 | CD14\_Monocyte\_CAM |
| ALDH1A2    | -1.756363 | 0.0075824 | CD14\_Monocyte\_CAM |
| GK         | -1.783567 | 0.0247609 | CD14\_Monocyte\_CAM |
| AQP9       | -1.831884 | 0.0247609 | CD14\_Monocyte\_CAM |
| DHX40      | -1.985389 | 0.0114309 | CD14\_Monocyte\_CAM |
| TESC       | -2.176803 | 0.0278135 | CD14\_Monocyte\_CAM |
| HSPB1      | -2.203517 | 0.0004854 | CD14\_Monocyte\_CAM |
| HSPA1B     | -2.204670 | 0.0267045 | CD14\_Monocyte\_CAM |
| MIR3945HG  | -2.308292 | 0.0114309 | CD14\_Monocyte\_CAM |
| CCR5AS     | -2.322602 | 0.0288191 | CD14\_Monocyte\_CAM |
| HIST2H2AA4 | -2.425398 | 0.0332086 | CD14\_Monocyte\_CAM |
| KCNJ15     | -2.440475 | 0.0114325 | CD14\_Monocyte\_CAM |
| BAG3       | -2.440998 | 0.0267045 | CD14\_Monocyte\_CAM |
| PROK2      | -2.442404 | 0.0247609 | CD14\_Monocyte\_CAM |
| TNFRSF10C  | -2.465714 | 0.0006030 | CD14\_Monocyte\_CAM |
| Z93241.1   | -2.617085 | 0.0024200 | CD14\_Monocyte\_CAM |
| STEAP4     | -2.640956 | 0.0114309 | CD14\_Monocyte\_CAM |
| THBD       | -2.687436 | 0.0486678 | CD14\_Monocyte\_CAM |
| NOTCH2NLB  | -2.831095 | 0.0247609 | CD14\_Monocyte\_CAM |
| G0S2       | -2.896204 | 0.0024276 | CD14\_Monocyte\_CAM |
| AL049651.1 | -3.348548 | 0.0024276 | CD14\_Monocyte\_CAM |
| C15orf48   | -3.996809 | 0.0462351 | CD14\_Monocyte\_CAM |
| CCL4L2     | -4.054442 | 0.0424658 | CD14\_Monocyte\_CAM |
| RPS4Y1     | -6.050787 | 0.0278135 | CD14\_Monocyte\_CAM |
| CLU        | -7.531582 | 0.0075824 | CD14\_Monocyte\_CAM |

\[1\] “Results for celltype:” \[1\] “CD16\_Monocyte\_CAM” \[1\] “Genes
upregulated in disease”

| genes        |    logFC |      padj | cell\_type          |
|:-------------|---------:|----------:|:--------------------|
| MT-ND6       | 5.340796 | 0.0000004 | CD16\_Monocyte\_CAM |
| DIRC3        | 4.335649 | 0.0072485 | CD16\_Monocyte\_CAM |
| STARD13      | 4.000554 | 0.0060448 | CD16\_Monocyte\_CAM |
| F13A1        | 3.914534 | 0.0405098 | CD16\_Monocyte\_CAM |
| BNC2         | 3.576286 | 0.0223041 | CD16\_Monocyte\_CAM |
| TBC1D16      | 3.215864 | 0.0341744 | CD16\_Monocyte\_CAM |
| SLC1A3       | 2.855539 | 0.0169047 | CD16\_Monocyte\_CAM |
| BCAS3        | 2.397001 | 0.0116266 | CD16\_Monocyte\_CAM |
| LRBA         | 2.389666 | 0.0271414 | CD16\_Monocyte\_CAM |
| RASAL2       | 2.325275 | 0.0223041 | CD16\_Monocyte\_CAM |
| PPM1L        | 2.300756 | 0.0231221 | CD16\_Monocyte\_CAM |
| PIAS2        | 2.294450 | 0.0273424 | CD16\_Monocyte\_CAM |
| XYLT1        | 2.158864 | 0.0169047 | CD16\_Monocyte\_CAM |
| CPEB3        | 2.059415 | 0.0223041 | CD16\_Monocyte\_CAM |
| SPECC1       | 2.057290 | 0.0223041 | CD16\_Monocyte\_CAM |
| HDAC9        | 2.035031 | 0.0341245 | CD16\_Monocyte\_CAM |
| OSBPL3       | 1.978678 | 0.0169047 | CD16\_Monocyte\_CAM |
| SIPA1L3      | 1.938934 | 0.0185836 | CD16\_Monocyte\_CAM |
| TFCP2        | 1.914541 | 0.0115844 | CD16\_Monocyte\_CAM |
| RNASE6       | 1.869460 | 0.0226378 | CD16\_Monocyte\_CAM |
| TET3         | 1.799787 | 0.0058822 | CD16\_Monocyte\_CAM |
| ARHGAP22     | 1.724190 | 0.0152444 | CD16\_Monocyte\_CAM |
| SNX29        | 1.715314 | 0.0263082 | CD16\_Monocyte\_CAM |
| TANGO6       | 1.691233 | 0.0413587 | CD16\_Monocyte\_CAM |
| MAP3K5       | 1.671366 | 0.0116266 | CD16\_Monocyte\_CAM |
| ADAMTSL4-AS1 | 1.616044 | 0.0282113 | CD16\_Monocyte\_CAM |
| CDK17        | 1.586614 | 0.0496453 | CD16\_Monocyte\_CAM |
| C20orf194    | 1.522525 | 0.0337562 | CD16\_Monocyte\_CAM |
| LINC01619    | 1.508318 | 0.0378186 | CD16\_Monocyte\_CAM |
| NCOA2        | 1.503954 | 0.0340013 | CD16\_Monocyte\_CAM |
| ARIH1        | 1.387130 | 0.0340013 | CD16\_Monocyte\_CAM |
| INPP5D       | 1.351285 | 0.0307109 | CD16\_Monocyte\_CAM |

\[1\] “Genes downregulated in disease”

| genes      |     logFC |      padj | cell\_type          |
|:-----------|----------:|----------:|:--------------------|
| MRPL18     | -1.356864 | 0.0341744 | CD16\_Monocyte\_CAM |
| WDR74      | -1.441716 | 0.0413587 | CD16\_Monocyte\_CAM |
| CD151      | -1.455472 | 0.0337562 | CD16\_Monocyte\_CAM |
| CD300C     | -1.484475 | 0.0240453 | CD16\_Monocyte\_CAM |
| ARL6IP1    | -1.525218 | 0.0383545 | CD16\_Monocyte\_CAM |
| PDCD5      | -1.704562 | 0.0223041 | CD16\_Monocyte\_CAM |
| PPA1       | -1.771987 | 0.0174862 | CD16\_Monocyte\_CAM |
| TREM1      | -1.774384 | 0.0223041 | CD16\_Monocyte\_CAM |
| C1D        | -1.793854 | 0.0169047 | CD16\_Monocyte\_CAM |
| IL2RG      | -1.948533 | 0.0124197 | CD16\_Monocyte\_CAM |
| CBWD1      | -1.992183 | 0.0089984 | CD16\_Monocyte\_CAM |
| SIRPB1     | -2.001635 | 0.0223041 | CD16\_Monocyte\_CAM |
| CCL3L1     | -2.455324 | 0.0340013 | CD16\_Monocyte\_CAM |
| SLC39A8    | -2.462203 | 0.0115337 | CD16\_Monocyte\_CAM |
| HIST2H2AA4 | -2.498836 | 0.0199403 | CD16\_Monocyte\_CAM |
| DEFB1      | -2.606993 | 0.0340013 | CD16\_Monocyte\_CAM |
| CCL3       | -2.685017 | 0.0169047 | CD16\_Monocyte\_CAM |
| C15orf48   | -3.503364 | 0.0223041 | CD16\_Monocyte\_CAM |
| VMO1       | -4.475082 | 0.0017436 | CD16\_Monocyte\_CAM |
| CCL4L2     | -4.550519 | 0.0060448 | CD16\_Monocyte\_CAM |
| CCL4       | -5.740481 | 0.0007680 | CD16\_Monocyte\_CAM |
| G0S2       | -6.415936 | 0.0001216 | CD16\_Monocyte\_CAM |

\[1\] “Results for celltype:” \[1\] “Macrophage\_CAM” \[1\] “Genes
upregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type      |
|:------|----------:|----------:|:----------------|
| TRAF1 | -2.210212 | 0.0277576 | Macrophage\_CAM |
| CXCL3 | -2.362022 | 0.0277576 | Macrophage\_CAM |
| GPR84 | -2.697627 | 0.0277576 | Macrophage\_CAM |

\[1\] “Results for celltype:” \[1\] “EVT\_CAM” \[1\] “Genes upregulated
in disease”

| genes      |     logFC |      padj | cell\_type |
|:-----------|----------:|----------:|:-----------|
| RPL30      | 11.466283 | 0.0135030 | EVT\_CAM   |
| RPS23      | 10.911637 | 0.0264971 | EVT\_CAM   |
| LY6E       | 10.706933 | 0.0131427 | EVT\_CAM   |
| S100A6     | 10.706685 | 0.0131427 | EVT\_CAM   |
| FAU        | 10.500993 | 0.0131427 | EVT\_CAM   |
| ELOB       | 10.223556 | 0.0131427 | EVT\_CAM   |
| COX6C      | 10.093270 | 0.0131427 | EVT\_CAM   |
| RPS15A     |  9.973222 | 0.0131427 | EVT\_CAM   |
| COX4I1     |  9.891996 | 0.0269994 | EVT\_CAM   |
| FTL        |  9.399466 | 0.0131427 | EVT\_CAM   |
| RPS28      |  8.860598 | 0.0135030 | EVT\_CAM   |
| HSPB1      |  8.643520 | 0.0131427 | EVT\_CAM   |
| IFI6       |  7.669921 | 0.0131427 | EVT\_CAM   |
| RPL18A     |  7.418521 | 0.0413946 | EVT\_CAM   |
| RPL41      |  7.304755 | 0.0248912 | EVT\_CAM   |
| FTH1       |  7.281810 | 0.0258570 | EVT\_CAM   |
| RPL37A     |  6.391821 | 0.0383308 | EVT\_CAM   |
| GDA        |  5.917532 | 0.0141285 | EVT\_CAM   |
| AC024230.1 |  5.497460 | 0.0131427 | EVT\_CAM   |
| C1QTNF7    |  4.484411 | 0.0378817 | EVT\_CAM   |
| AC026167.1 |  3.977433 | 0.0392552 | EVT\_CAM   |

\[1\] “Genes downregulated in disease”

| genes     |     logFC |     padj | cell\_type |
|:----------|----------:|---------:|:-----------|
| LINC02109 | -7.598038 | 2.39e-05 | EVT\_CAM   |

\[1\] “Results for celltype:” \[1\] “Smooth\_muscle\_CAM” \[1\] “Genes
upregulated in disease”

| genes  |    logFC |      padj | cell\_type          |
|:-------|---------:|----------:|:--------------------|
| CORIN  | 8.109637 | 0.0020346 | Smooth\_muscle\_CAM |
| PAPPA2 | 7.717291 | 0.0062000 | Smooth\_muscle\_CAM |
| HLA-B  | 6.822006 | 0.0253130 | Smooth\_muscle\_CAM |
| GRIA1  | 4.879476 | 0.0253130 | Smooth\_muscle\_CAM |
| PDGFD  | 4.147795 | 0.0276363 | Smooth\_muscle\_CAM |

\[1\] “Genes downregulated in disease”

| genes      |     logFC |      padj | cell\_type          |
|:-----------|----------:|----------:|:--------------------|
| TNFAIP3    | -3.671241 | 0.0253130 | Smooth\_muscle\_CAM |
| PDE1C      | -3.923392 | 0.0082061 | Smooth\_muscle\_CAM |
| NAMPT      | -4.071637 | 0.0251227 | Smooth\_muscle\_CAM |
| LINC01705  | -4.562816 | 0.0052538 | Smooth\_muscle\_CAM |
| SH3BP4     | -4.589180 | 0.0133085 | Smooth\_muscle\_CAM |
| MT-ND1     | -4.630791 | 0.0179186 | Smooth\_muscle\_CAM |
| ADAMTS14   | -4.960155 | 0.0008832 | Smooth\_muscle\_CAM |
| MT-ATP6    | -5.068161 | 0.0082061 | Smooth\_muscle\_CAM |
| MT-CO2     | -5.085273 | 0.0452265 | Smooth\_muscle\_CAM |
| MT-ND4L    | -5.164130 | 0.0283945 | Smooth\_muscle\_CAM |
| MT-CO3     | -5.222929 | 0.0268386 | Smooth\_muscle\_CAM |
| MT-ND2     | -5.325770 | 0.0062000 | Smooth\_muscle\_CAM |
| MT-ND4     | -5.787150 | 0.0111814 | Smooth\_muscle\_CAM |
| FP671120.4 | -7.230133 | 0.0317879 | Smooth\_muscle\_CAM |
| PABPC1     | -7.359343 | 0.0020346 | Smooth\_muscle\_CAM |
| RPS24      | -7.960806 | 0.0020346 | Smooth\_muscle\_CAM |

# Volcano plots in Late disease

\[1\] “Results for celltype:” \[1\] “Naive\_B\_CAM” \[1\] “Genes
upregulated in disease”

| genes      |    logFC |      padj | cell\_type    |
|:-----------|---------:|----------:|:--------------|
| IGHV5-10-1 | 8.420236 | 0.0000420 | Naive\_B\_CAM |
| IGLV3-19   | 7.238801 | 0.0009346 | Naive\_B\_CAM |
| RPS4Y1     | 6.559610 | 0.0000420 | Naive\_B\_CAM |
| IGHV3-33   | 3.577425 | 0.0015784 | Naive\_B\_CAM |

\[1\] “Genes downregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Results for celltype:” \[1\] “Memory\_B\_CAM” \[1\] “Genes
upregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Genes downregulated in disease”

| genes    |     logFC |      padj | cell\_type     |
|:---------|----------:|----------:|:---------------|
| ANXA4    | -3.972298 | 0.0155245 | Memory\_B\_CAM |
| HSPB1    | -6.179776 | 0.0325828 | Memory\_B\_CAM |
| IGKV1-39 | -9.831277 | 0.0029708 | Memory\_B\_CAM |

\[1\] “Results for celltype:” \[1\] “Fetal\_CD4\_T\_Naive\_CM\_CAM”
\[1\] “Genes upregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Genes downregulated in disease”

| genes    |     logFC |      padj | cell\_type                    |
|:---------|----------:|----------:|:------------------------------|
| RGS1     | -2.298837 | 0.0421507 | Fetal\_CD4\_T\_Naive\_CM\_CAM |
| TRBV11-2 | -3.515802 | 0.0343506 | Fetal\_CD4\_T\_Naive\_CM\_CAM |
| TRBV10-3 | -4.829502 | 0.0000096 | Fetal\_CD4\_T\_Naive\_CM\_CAM |
| TRBV6-1  | -5.550681 | 0.0000005 | Fetal\_CD4\_T\_Naive\_CM\_CAM |

\[1\] “Results for celltype:” \[1\] “CD4\_T\_Naive\_CM\_CAM” \[1\]
“Genes upregulated in disease”

| genes   |    logFC |      padj | cell\_type             |
|:--------|---------:|----------:|:-----------------------|
| FMN1    | 5.511262 | 0.0055167 | CD4\_T\_Naive\_CM\_CAM |
| HBG2    | 5.246769 | 0.0077601 | CD4\_T\_Naive\_CM\_CAM |
| SOX4    | 4.361920 | 0.0114160 | CD4\_T\_Naive\_CM\_CAM |
| U2AF1   | 3.249218 | 0.0288364 | CD4\_T\_Naive\_CM\_CAM |
| TOX     | 2.839507 | 0.0360731 | CD4\_T\_Naive\_CM\_CAM |
| FAM241A | 1.748368 | 0.0170781 | CD4\_T\_Naive\_CM\_CAM |
| HLA-C   | 1.235068 | 0.0170781 | CD4\_T\_Naive\_CM\_CAM |

\[1\] “Genes downregulated in disease”

| genes     |     logFC |      padj | cell\_type             |
|:----------|----------:|----------:|:-----------------------|
| TRAC      | -1.474074 | 0.0055167 | CD4\_T\_Naive\_CM\_CAM |
| MX2       | -1.614140 | 0.0395752 | CD4\_T\_Naive\_CM\_CAM |
| HSBP1     | -1.901340 | 0.0281398 | CD4\_T\_Naive\_CM\_CAM |
| PPP1CB    | -1.919295 | 0.0281398 | CD4\_T\_Naive\_CM\_CAM |
| CHI3L2    | -1.938837 | 0.0281398 | CD4\_T\_Naive\_CM\_CAM |
| MED15     | -2.192359 | 0.0109083 | CD4\_T\_Naive\_CM\_CAM |
| TRBV19    | -2.521269 | 0.0002293 | CD4\_T\_Naive\_CM\_CAM |
| ANXA5     | -2.717096 | 0.0011829 | CD4\_T\_Naive\_CM\_CAM |
| TRAV29DV5 | -2.820045 | 0.0391099 | CD4\_T\_Naive\_CM\_CAM |
| TRBV12-3  | -2.829849 | 0.0301214 | CD4\_T\_Naive\_CM\_CAM |
| RCBTB2    | -3.134366 | 0.0171996 | CD4\_T\_Naive\_CM\_CAM |
| TRBV6-5   | -3.137120 | 0.0020025 | CD4\_T\_Naive\_CM\_CAM |
| RBMX2     | -3.144391 | 0.0281398 | CD4\_T\_Naive\_CM\_CAM |
| TRBV12-4  | -3.372057 | 0.0103547 | CD4\_T\_Naive\_CM\_CAM |
| TRAV8-1   | -3.499255 | 0.0170781 | CD4\_T\_Naive\_CM\_CAM |
| TRBV7-9   | -3.961439 | 0.0002293 | CD4\_T\_Naive\_CM\_CAM |
| TRBV5-6   | -4.294513 | 0.0170781 | CD4\_T\_Naive\_CM\_CAM |
| TRBV6-1   | -4.974279 | 0.0023292 | CD4\_T\_Naive\_CM\_CAM |
| TRBV3-1   | -5.167653 | 0.0006423 | CD4\_T\_Naive\_CM\_CAM |
| TRAV8-4   | -5.249429 | 0.0002293 | CD4\_T\_Naive\_CM\_CAM |
| TRBV28    | -5.270779 | 0.0002293 | CD4\_T\_Naive\_CM\_CAM |
| GATD3A    | -6.101101 | 0.0170781 | CD4\_T\_Naive\_CM\_CAM |
| CRIP2     | -7.185925 | 0.0000004 | CD4\_T\_Naive\_CM\_CAM |
| TRBV30    | -7.199817 | 0.0170781 | CD4\_T\_Naive\_CM\_CAM |
| TRBV4-1   | -7.393158 | 0.0013838 | CD4\_T\_Naive\_CM\_CAM |

\[1\] “Results for celltype:” \[1\] “CD4\_Th\_CAM” \[1\] “Genes
upregulated in disease”

| genes |    logFC |      padj | cell\_type   |
|:------|---------:|----------:|:-------------|
| SOCS2 | 2.015560 | 0.0239304 | CD4\_Th\_CAM |
| HLA-C | 1.172837 | 0.0296902 | CD4\_Th\_CAM |

\[1\] “Genes downregulated in disease”

| genes    |     logFC |      padj | cell\_type   |
|:---------|----------:|----------:|:-------------|
| SELL     | -1.111503 | 0.0363707 | CD4\_Th\_CAM |
| TRAC     | -1.439802 | 0.0067672 | CD4\_Th\_CAM |
| TRBC1    | -1.449897 | 0.0168309 | CD4\_Th\_CAM |
| LGALS1   | -1.527621 | 0.0335454 | CD4\_Th\_CAM |
| HSPB1    | -1.706348 | 0.0008379 | CD4\_Th\_CAM |
| CD74     | -2.143536 | 0.0005794 | CD4\_Th\_CAM |
| IGF2R    | -2.966972 | 0.0167543 | CD4\_Th\_CAM |
| TRAV8-2  | -3.785172 | 0.0008379 | CD4\_Th\_CAM |
| TRAV12-2 | -3.918954 | 0.0067672 | CD4\_Th\_CAM |
| LGALS9   | -5.875214 | 0.0239304 | CD4\_Th\_CAM |
| TRAV8-4  | -6.549860 | 0.0168309 | CD4\_Th\_CAM |
| TRAV3    | -6.607086 | 0.0037002 | CD4\_Th\_CAM |
| TRBV6-2  | -6.888250 | 0.0020281 | CD4\_Th\_CAM |
| TRBV6-1  | -7.258807 | 0.0160778 | CD4\_Th\_CAM |
| TRBV6-6  | -7.306876 | 0.0008379 | CD4\_Th\_CAM |
| TRBV12-3 | -7.463529 | 0.0020281 | CD4\_Th\_CAM |
| TRBV7-9  | -7.601095 | 0.0011404 | CD4\_Th\_CAM |
| TRBV12-4 | -7.648138 | 0.0024301 | CD4\_Th\_CAM |

\[1\] “Results for celltype:” \[1\] “FoxP3-Treg\_CAM” \[1\] “Genes
upregulated in disease”

| genes  |    logFC |      padj | cell\_type      |
|:-------|---------:|----------:|:----------------|
| MTSS1  | 3.305912 | 0.0297197 | FoxP3-Treg\_CAM |
| ACTG2  | 3.030478 | 0.0015407 | FoxP3-Treg\_CAM |
| ATP1B1 | 2.315764 | 0.0029557 | FoxP3-Treg\_CAM |
| ELK3   | 1.673669 | 0.0067395 | FoxP3-Treg\_CAM |
| HDAC4  | 1.667549 | 0.0183836 | FoxP3-Treg\_CAM |

\[1\] “Genes downregulated in disease”

| genes      |     logFC |      padj | cell\_type      |
|:-----------|----------:|----------:|:----------------|
| RTKN2      | -1.325763 | 0.0336772 | FoxP3-Treg\_CAM |
| GSTP1      | -1.419004 | 0.0297197 | FoxP3-Treg\_CAM |
| CITED4     | -2.119445 | 0.0183836 | FoxP3-Treg\_CAM |
| SNHG25     | -2.789756 | 0.0183836 | FoxP3-Treg\_CAM |
| HLA-DRB5   | -2.809828 | 0.0045034 | FoxP3-Treg\_CAM |
| TRAV17     | -4.298352 | 0.0136605 | FoxP3-Treg\_CAM |
| AL627171.2 | -4.445038 | 0.0064998 | FoxP3-Treg\_CAM |
| TRAV26-1   | -4.823497 | 0.0029948 | FoxP3-Treg\_CAM |
| TRBV6-5    | -5.810043 | 0.0001424 | FoxP3-Treg\_CAM |
| TRAV13-1   | -6.509035 | 0.0259648 | FoxP3-Treg\_CAM |
| TRBV5-5    | -6.759540 | 0.0259648 | FoxP3-Treg\_CAM |
| TRAV27     | -6.791445 | 0.0180754 | FoxP3-Treg\_CAM |
| TRBV6-2    | -6.891983 | 0.0023641 | FoxP3-Treg\_CAM |
| CCL15      | -6.997224 | 0.0297197 | FoxP3-Treg\_CAM |
| TRBV11-1   | -7.050697 | 0.0277598 | FoxP3-Treg\_CAM |
| TRBV28     | -7.211573 | 0.0015407 | FoxP3-Treg\_CAM |
| TRAV14DV4  | -7.267412 | 0.0013595 | FoxP3-Treg\_CAM |
| TRBV12-3   | -7.279290 | 0.0045034 | FoxP3-Treg\_CAM |
| TRBV24-1   | -7.304115 | 0.0004719 | FoxP3-Treg\_CAM |
| TRBV3-1    | -7.432565 | 0.0064982 | FoxP3-Treg\_CAM |
| TRAV8-2    | -7.446478 | 0.0001817 | FoxP3-Treg\_CAM |
| TRAV29DV5  | -7.533647 | 0.0012739 | FoxP3-Treg\_CAM |
| TRAV12-2   | -7.669714 | 0.0002320 | FoxP3-Treg\_CAM |
| TRBV2      | -7.713291 | 0.0015407 | FoxP3-Treg\_CAM |
| TRBV11-2   | -7.763314 | 0.0001817 | FoxP3-Treg\_CAM |

\[1\] “Results for celltype:” \[1\] “CD4\_TEM\_CAM” \[1\] “Genes
upregulated in disease”

| genes      |     logFC |      padj | cell\_type    |
|:-----------|----------:|----------:|:--------------|
| FMN1       | 6.1288126 | 0.0348395 | CD4\_TEM\_CAM |
| DLGAP3     | 3.4499185 | 0.0146521 | CD4\_TEM\_CAM |
| DHTKD1     | 2.7184366 | 0.0348395 | CD4\_TEM\_CAM |
| NCAPG2     | 2.5839985 | 0.0177158 | CD4\_TEM\_CAM |
| AL136456.1 | 1.9957849 | 0.0318707 | CD4\_TEM\_CAM |
| ALAD       | 1.7077196 | 0.0318707 | CD4\_TEM\_CAM |
| DGKD       | 1.6633072 | 0.0172701 | CD4\_TEM\_CAM |
| CEBPD      | 1.5799187 | 0.0282062 | CD4\_TEM\_CAM |
| CRACR2A    | 1.5471664 | 0.0348395 | CD4\_TEM\_CAM |
| ATP8A1     | 1.5233528 | 0.0348395 | CD4\_TEM\_CAM |
| FAM193A    | 1.3705356 | 0.0453551 | CD4\_TEM\_CAM |
| PCAT1      | 1.3656413 | 0.0180735 | CD4\_TEM\_CAM |
| XYLT1      | 1.3148559 | 0.0146521 | CD4\_TEM\_CAM |
| PLCB1      | 1.1993139 | 0.0441492 | CD4\_TEM\_CAM |
| AKT3       | 1.1965737 | 0.0095423 | CD4\_TEM\_CAM |
| RUNX2      | 1.1949940 | 0.0282062 | CD4\_TEM\_CAM |
| PCF11      | 1.1766573 | 0.0318707 | CD4\_TEM\_CAM |
| PDE3B      | 1.1312824 | 0.0419255 | CD4\_TEM\_CAM |
| ELK3       | 1.0964341 | 0.0318707 | CD4\_TEM\_CAM |
| CERK       | 1.0536950 | 0.0348395 | CD4\_TEM\_CAM |
| HLA-C      | 1.0283064 | 0.0386693 | CD4\_TEM\_CAM |
| ZEB1       | 1.0196424 | 0.0458725 | CD4\_TEM\_CAM |
| FBXL17     | 1.0182777 | 0.0476252 | CD4\_TEM\_CAM |
| LARP1      | 1.0080886 | 0.0458500 | CD4\_TEM\_CAM |
| TBC1D22A   | 0.9979368 | 0.0254163 | CD4\_TEM\_CAM |
| RBM33      | 0.9747973 | 0.0377318 | CD4\_TEM\_CAM |
| VPS13B     | 0.9632457 | 0.0423606 | CD4\_TEM\_CAM |
| CDC14A     | 0.9130486 | 0.0414104 | CD4\_TEM\_CAM |
| RPS6KA3    | 0.9117998 | 0.0146768 | CD4\_TEM\_CAM |
| SYTL3      | 0.8854149 | 0.0451970 | CD4\_TEM\_CAM |
| JAK1       | 0.8838815 | 0.0095423 | CD4\_TEM\_CAM |
| TLK1       | 0.8824614 | 0.0348395 | CD4\_TEM\_CAM |
| ANKRD11    | 0.8731654 | 0.0348395 | CD4\_TEM\_CAM |
| USP3       | 0.8507536 | 0.0318707 | CD4\_TEM\_CAM |
| ACIN1      | 0.8507233 | 0.0494539 | CD4\_TEM\_CAM |
| MAML2      | 0.7786000 | 0.0386693 | CD4\_TEM\_CAM |
| FNBP4      | 0.7774340 | 0.0432892 | CD4\_TEM\_CAM |
| CBLB       | 0.7174881 | 0.0476252 | CD4\_TEM\_CAM |
| RIPOR2     | 0.7054428 | 0.0348395 | CD4\_TEM\_CAM |
| DDX17      | 0.6717504 | 0.0355916 | CD4\_TEM\_CAM |

\[1\] “Genes downregulated in disease”

| genes      |      logFC |      padj | cell\_type    |
|:-----------|-----------:|----------:|:--------------|
| KRTCAP2    | -0.7391641 | 0.0351443 | CD4\_TEM\_CAM |
| CIAO2B     | -0.7483198 | 0.0430764 | CD4\_TEM\_CAM |
| SERF2      | -0.7637491 | 0.0318707 | CD4\_TEM\_CAM |
| WDR83OS    | -0.7824104 | 0.0348395 | CD4\_TEM\_CAM |
| MYL6       | -0.7861345 | 0.0348395 | CD4\_TEM\_CAM |
| NDUFB8     | -0.7866780 | 0.0303892 | CD4\_TEM\_CAM |
| ELOB       | -0.7885713 | 0.0348395 | CD4\_TEM\_CAM |
| PSMA5      | -0.7944616 | 0.0348395 | CD4\_TEM\_CAM |
| MYL12A     | -0.8053803 | 0.0410689 | CD4\_TEM\_CAM |
| TXN        | -0.8088751 | 0.0359101 | CD4\_TEM\_CAM |
| TMSB10     | -0.8277985 | 0.0476252 | CD4\_TEM\_CAM |
| NOP10      | -0.8300975 | 0.0476252 | CD4\_TEM\_CAM |
| IFITM2     | -0.8436096 | 0.0372643 | CD4\_TEM\_CAM |
| UBL5       | -0.8457090 | 0.0318707 | CD4\_TEM\_CAM |
| ATP5MC1    | -0.8504148 | 0.0385862 | CD4\_TEM\_CAM |
| ATP6V1F    | -0.8540521 | 0.0348395 | CD4\_TEM\_CAM |
| MICOS10    | -0.8582937 | 0.0146768 | CD4\_TEM\_CAM |
| PFN1       | -0.8660401 | 0.0146521 | CD4\_TEM\_CAM |
| OSTC       | -0.8844435 | 0.0348395 | CD4\_TEM\_CAM |
| COX7A2     | -0.8913065 | 0.0146521 | CD4\_TEM\_CAM |
| TRMT112    | -0.8914877 | 0.0235913 | CD4\_TEM\_CAM |
| POLR2L     | -0.9025744 | 0.0348395 | CD4\_TEM\_CAM |
| SH3BGRL3   | -0.9025887 | 0.0182142 | CD4\_TEM\_CAM |
| RBCK1      | -0.9139863 | 0.0318707 | CD4\_TEM\_CAM |
| CKLF       | -0.9273174 | 0.0348395 | CD4\_TEM\_CAM |
| GSTO1      | -0.9322101 | 0.0265250 | CD4\_TEM\_CAM |
| ZNHIT1     | -0.9342107 | 0.0348395 | CD4\_TEM\_CAM |
| MRPL33     | -0.9387422 | 0.0476252 | CD4\_TEM\_CAM |
| SDF2L1     | -0.9479817 | 0.0318707 | CD4\_TEM\_CAM |
| MZT2A      | -0.9498367 | 0.0348395 | CD4\_TEM\_CAM |
| BANF1      | -0.9553266 | 0.0146521 | CD4\_TEM\_CAM |
| HSBP1      | -0.9588953 | 0.0375312 | CD4\_TEM\_CAM |
| THYN1      | -0.9655174 | 0.0424637 | CD4\_TEM\_CAM |
| COX8A      | -0.9659335 | 0.0282062 | CD4\_TEM\_CAM |
| HSPE1      | -0.9676423 | 0.0146768 | CD4\_TEM\_CAM |
| GZMA       | -0.9869870 | 0.0143689 | CD4\_TEM\_CAM |
| ATP5MPL    | -0.9971709 | 0.0221553 | CD4\_TEM\_CAM |
| ORMDL2     | -1.0007248 | 0.0318707 | CD4\_TEM\_CAM |
| MOSPD3     | -1.0026924 | 0.0483250 | CD4\_TEM\_CAM |
| POMP       | -1.0068219 | 0.0172173 | CD4\_TEM\_CAM |
| LTB        | -1.0099834 | 0.0254163 | CD4\_TEM\_CAM |
| CD48       | -1.0103916 | 0.0146521 | CD4\_TEM\_CAM |
| PLP2       | -1.0138949 | 0.0348395 | CD4\_TEM\_CAM |
| ALOX5AP    | -1.0158263 | 0.0318707 | CD4\_TEM\_CAM |
| ATP5F1E    | -1.0177106 | 0.0294757 | CD4\_TEM\_CAM |
| APOBEC3H   | -1.0205416 | 0.0346933 | CD4\_TEM\_CAM |
| BOLA3      | -1.0320776 | 0.0386693 | CD4\_TEM\_CAM |
| APEX1      | -1.0360240 | 0.0318707 | CD4\_TEM\_CAM |
| CISD3      | -1.0390279 | 0.0146521 | CD4\_TEM\_CAM |
| GSTP1      | -1.0531880 | 0.0235913 | CD4\_TEM\_CAM |
| PDLIM2     | -1.0536862 | 0.0362421 | CD4\_TEM\_CAM |
| CAPZA2     | -1.0703945 | 0.0235913 | CD4\_TEM\_CAM |
| SUMO3      | -1.1033925 | 0.0441492 | CD4\_TEM\_CAM |
| H2AFZ      | -1.1165978 | 0.0095423 | CD4\_TEM\_CAM |
| CITED4     | -1.1450842 | 0.0318707 | CD4\_TEM\_CAM |
| CTSB       | -1.1571361 | 0.0458500 | CD4\_TEM\_CAM |
| LTA        | -1.1611226 | 0.0348395 | CD4\_TEM\_CAM |
| PHPT1      | -1.1817845 | 0.0294757 | CD4\_TEM\_CAM |
| CD40LG     | -1.1913938 | 0.0000394 | CD4\_TEM\_CAM |
| NDUFA3     | -1.2109871 | 0.0095423 | CD4\_TEM\_CAM |
| HIGD1A     | -1.2120867 | 0.0348395 | CD4\_TEM\_CAM |
| MRPL54     | -1.2380159 | 0.0009543 | CD4\_TEM\_CAM |
| YIF1A      | -1.2518670 | 0.0318707 | CD4\_TEM\_CAM |
| SPAG7      | -1.2577260 | 0.0282062 | CD4\_TEM\_CAM |
| CCDC167    | -1.2692999 | 0.0348395 | CD4\_TEM\_CAM |
| LGALS1     | -1.2705557 | 0.0146521 | CD4\_TEM\_CAM |
| LINC00892  | -1.2990400 | 0.0065041 | CD4\_TEM\_CAM |
| AKR1B1     | -1.3583168 | 0.0287550 | CD4\_TEM\_CAM |
| ADA        | -1.3622405 | 0.0348395 | CD4\_TEM\_CAM |
| ACP5       | -1.3727830 | 0.0453551 | CD4\_TEM\_CAM |
| RNF181     | -1.3729306 | 0.0009614 | CD4\_TEM\_CAM |
| C17orf49   | -1.4178656 | 0.0017497 | CD4\_TEM\_CAM |
| HSPB1      | -1.4571921 | 0.0003937 | CD4\_TEM\_CAM |
| WARS       | -1.5049372 | 0.0424637 | CD4\_TEM\_CAM |
| CD74       | -1.5274595 | 0.0146521 | CD4\_TEM\_CAM |
| AC004585.1 | -1.5984507 | 0.0458500 | CD4\_TEM\_CAM |
| PTTG1      | -1.6030591 | 0.0124062 | CD4\_TEM\_CAM |
| CTSZ       | -1.6759380 | 0.0214582 | CD4\_TEM\_CAM |
| NUDT1      | -1.7195084 | 0.0059562 | CD4\_TEM\_CAM |
| HLA-DMA    | -1.8083663 | 0.0483250 | CD4\_TEM\_CAM |
| C1orf162   | -2.0154786 | 0.0069074 | CD4\_TEM\_CAM |
| XCL1       | -2.0749163 | 0.0000013 | CD4\_TEM\_CAM |
| GZMH       | -2.0799973 | 0.0180735 | CD4\_TEM\_CAM |
| IFNG       | -2.9307335 | 0.0019734 | CD4\_TEM\_CAM |
| TRAV1-2    | -3.4572171 | 0.0348395 | CD4\_TEM\_CAM |
| QSOX1      | -3.4797980 | 0.0282062 | CD4\_TEM\_CAM |
| TRAV19     | -4.3116006 | 0.0146521 | CD4\_TEM\_CAM |
| TRBV7-3    | -4.6367487 | 0.0214582 | CD4\_TEM\_CAM |
| FUCA2      | -6.7164476 | 0.0146521 | CD4\_TEM\_CAM |
| IGHV1-3    | -6.7409047 | 0.0172173 | CD4\_TEM\_CAM |
| HBEGF      | -6.8148022 | 0.0146521 | CD4\_TEM\_CAM |
| TRBV10-2   | -7.0874527 | 0.0318707 | CD4\_TEM\_CAM |
| TRAV14DV4  | -7.5689801 | 0.0221294 | CD4\_TEM\_CAM |
| HLA-DRB5   | -8.3131394 | 0.0000394 | CD4\_TEM\_CAM |

\[1\] “Results for celltype:” \[1\] “CD8\_TEM\_CAM” \[1\] “Genes
upregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Genes downregulated in disease”

| genes  |     logFC |      padj | cell\_type    |
|:-------|----------:|----------:|:--------------|
| PECAM1 | -2.414285 | 0.0070344 | CD8\_TEM\_CAM |
| STX7   | -2.668147 | 0.0070344 | CD8\_TEM\_CAM |

\[1\] “Results for celltype:” \[1\] “NKT\_CAM” \[1\] “Genes upregulated
in disease”

| genes    |     logFC |      padj | cell\_type |
|:---------|----------:|----------:|:-----------|
| NEGR1    | 3.8544096 | 0.0417053 | NKT\_CAM   |
| HNRNPUL1 | 0.9183659 | 0.0181287 | NKT\_CAM   |

\[1\] “Genes downregulated in disease”

| genes   |     logFC |      padj | cell\_type |
|:--------|----------:|----------:|:-----------|
| RREB1   | -2.502390 | 0.0417053 | NKT\_CAM   |
| CCL4L2  | -2.950748 | 0.0181287 | NKT\_CAM   |
| PPP1R3B | -4.246331 | 0.0181287 | NKT\_CAM   |
| MYL9    | -7.782021 | 0.0417053 | NKT\_CAM   |

\[1\] “Results for celltype:” \[1\] “CD14\_Monocyte\_CAM” \[1\] “Genes
upregulated in disease”

| genes      |    logFC |      padj | cell\_type          |
|:-----------|---------:|----------:|:--------------------|
| IL1R2      | 5.384123 | 0.0348342 | CD14\_Monocyte\_CAM |
| ADAMTS2    | 4.468815 | 0.0097794 | CD14\_Monocyte\_CAM |
| VSIG4      | 2.709228 | 0.0348342 | CD14\_Monocyte\_CAM |
| AC084871.1 | 2.077111 | 0.0071528 | CD14\_Monocyte\_CAM |

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type          |
|:------|----------:|----------:|:--------------------|
| C1QB  | -6.406533 | 0.0348342 | CD14\_Monocyte\_CAM |

\[1\] “Results for celltype:” \[1\] “CD16\_Monocyte\_CAM” \[1\] “Genes
upregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type          |
|:------|----------:|----------:|:--------------------|
| IL4I1 | -2.778376 | 0.0060756 | CD16\_Monocyte\_CAM |
| HAMP  | -4.476269 | 0.0060756 | CD16\_Monocyte\_CAM |

\[1\] “Results for celltype:” \[1\] “EVT\_CAM” \[1\] “Genes upregulated
in disease”

| genes      |    logFC |      padj | cell\_type |
|:-----------|---------:|----------:|:-----------|
| AC063944.1 | 3.670199 | 0.0471515 | EVT\_CAM   |

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type |
|:------|----------:|----------:|:-----------|
| FTL   | -7.858888 | 0.0072908 | EVT\_CAM   |
| MT2A  | -8.628985 | 0.0377687 | EVT\_CAM   |
| RPL39 | -9.498683 | 0.0471515 | EVT\_CAM   |

\[1\] “Results for celltype:” \[1\] “STB\_CAM” \[1\] “Genes upregulated
in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Genes downregulated in disease”

| genes  |     logFC |      padj | cell\_type |
|:-------|----------:|----------:|:-----------|
| ADGRL3 | -5.875056 | 0.0159121 | STB\_CAM   |
| PLCB1  | -6.458392 | 0.0023946 | STB\_CAM   |
| GPNMB  | -6.522743 | 0.0023946 | STB\_CAM   |

\[1\] “Results for celltype:” \[1\] “Smooth\_muscle\_CAM” \[1\] “Genes
upregulated in disease”

| genes      |    logFC |      padj | cell\_type          |
|:-----------|---------:|----------:|:--------------------|
| S100A4     | 7.383460 | 0.0020250 | Smooth\_muscle\_CAM |
| FP671120.4 | 5.622551 | 0.0026750 | Smooth\_muscle\_CAM |
| HFM1       | 5.342155 | 0.0000699 | Smooth\_muscle\_CAM |
| IGFBP6     | 5.139661 | 0.0002674 | Smooth\_muscle\_CAM |
| MTRNR2L12  | 4.782766 | 0.0000699 | Smooth\_muscle\_CAM |
| RPL39      | 4.256679 | 0.0086271 | Smooth\_muscle\_CAM |
| MT-ATP6    | 4.147909 | 0.0084012 | Smooth\_muscle\_CAM |
| FBN2       | 3.950063 | 0.0086271 | Smooth\_muscle\_CAM |
| LINC01091  | 3.793060 | 0.0014142 | Smooth\_muscle\_CAM |
| RPS27      | 3.780748 | 0.0073255 | Smooth\_muscle\_CAM |
| TAFA1      | 3.750532 | 0.0441189 | Smooth\_muscle\_CAM |
| RPL36      | 3.503978 | 0.0109426 | Smooth\_muscle\_CAM |
| TPT1       | 3.261437 | 0.0073255 | Smooth\_muscle\_CAM |
| S100A6     | 3.259460 | 0.0084012 | Smooth\_muscle\_CAM |
| TMSB10     | 3.250926 | 0.0382889 | Smooth\_muscle\_CAM |
| MT-ND6     | 3.148891 | 0.0074484 | Smooth\_muscle\_CAM |
| RPLP1      | 3.146351 | 0.0098495 | Smooth\_muscle\_CAM |
| RPL28      | 3.067523 | 0.0382889 | Smooth\_muscle\_CAM |
| ELN        | 3.057608 | 0.0157586 | Smooth\_muscle\_CAM |
| CTNNA3     | 2.775498 | 0.0375681 | Smooth\_muscle\_CAM |
| PLA2G5     | 2.717782 | 0.0370975 | Smooth\_muscle\_CAM |

\[1\] “Genes downregulated in disease”

| genes  |     logFC |      padj | cell\_type          |
|:-------|----------:|----------:|:--------------------|
| H19    | -3.235860 | 0.0317502 | Smooth\_muscle\_CAM |
| NABP1  | -3.590423 | 0.0443669 | Smooth\_muscle\_CAM |
| ABCC9  | -3.870639 | 0.0084012 | Smooth\_muscle\_CAM |
| PDE8B  | -4.145213 | 0.0084012 | Smooth\_muscle\_CAM |
| RGS5   | -4.156567 | 0.0071034 | Smooth\_muscle\_CAM |
| FRMD4B | -4.356230 | 0.0084012 | Smooth\_muscle\_CAM |
| STEAP4 | -4.720173 | 0.0022903 | Smooth\_muscle\_CAM |
| KCNIP4 | -6.293910 | 0.0252622 | Smooth\_muscle\_CAM |
| PI15   | -6.578137 | 0.0084012 | Smooth\_muscle\_CAM |
| FGF14  | -6.613583 | 0.0014142 | Smooth\_muscle\_CAM |
| REN    | -6.716534 | 0.0016191 | Smooth\_muscle\_CAM |
| DNAH11 | -6.892572 | 0.0008606 | Smooth\_muscle\_CAM |
| ADARB2 | -8.605282 | 0.0001383 | Smooth\_muscle\_CAM |

# Volcano plots in averagePE

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes      |    logFC |      padj | cell\_type    |
|:-----------|---------:|----------:|:--------------|
| IGHV5-10-1 | 4.717552 | 0.0063942 | Naive\_B\_CAM |
| IGHV6-1    | 4.364089 | 0.0307936 | Naive\_B\_CAM |
| MT-ND6     | 1.862188 | 0.0127668 | Naive\_B\_CAM |

\[1\] “Genes downregulated in disease”

| genes    |     logFC |      padj | cell\_type    |
|:---------|----------:|----------:|:--------------|
| IGLV1-40 | -2.690356 | 0.0127668 | Naive\_B\_CAM |
| IGHV2-5  | -3.870863 | 0.0030617 | Naive\_B\_CAM |
| IGKV1-8  | -4.286473 | 0.0025884 | Naive\_B\_CAM |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes  |    logFC |      padj | cell\_type     |
|:-------|---------:|----------:|:---------------|
| MT-ND6 | 2.242537 | 0.0272034 | Memory\_B\_CAM |

\[1\] “Genes downregulated in disease”

| genes    |     logFC |      padj | cell\_type     |
|:---------|----------:|----------:|:---------------|
| ID3      | -1.210575 | 0.0272034 | Memory\_B\_CAM |
| JPT1     | -1.312769 | 0.0272034 | Memory\_B\_CAM |
| ANXA4    | -2.340144 | 0.0180142 | Memory\_B\_CAM |
| IGKV1-39 | -4.883969 | 0.0272034 | Memory\_B\_CAM |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes  |     logFC |      padj | cell\_type                    |
|:-------|----------:|----------:|:------------------------------|
| ERAP2  | 2.2703007 | 0.0000819 | Fetal\_CD4\_T\_Naive\_CM\_CAM |
| FKBP5  | 1.9850984 | 0.0216858 | Fetal\_CD4\_T\_Naive\_CM\_CAM |
| MT-CYB | 1.2836113 | 0.0475563 | Fetal\_CD4\_T\_Naive\_CM\_CAM |
| BCL2   | 0.9811713 | 0.0432809 | Fetal\_CD4\_T\_Naive\_CM\_CAM |

\[1\] “Genes downregulated in disease”

| genes    |     logFC |      padj | cell\_type                    |
|:---------|----------:|----------:|:------------------------------|
| DNPH1    | -1.191209 | 0.0133051 | Fetal\_CD4\_T\_Naive\_CM\_CAM |
| RGS1     | -2.885664 | 0.0000819 | Fetal\_CD4\_T\_Naive\_CM\_CAM |
| TRBV10-3 | -3.367925 | 0.0000691 | Fetal\_CD4\_T\_Naive\_CM\_CAM |
| TRBV30   | -4.260392 | 0.0173448 | Fetal\_CD4\_T\_Naive\_CM\_CAM |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes   |     logFC |      padj | cell\_type             |
|:--------|----------:|----------:|:-----------------------|
| FMN1    | 6.3355004 | 0.0000030 | CD4\_T\_Naive\_CM\_CAM |
| TRBV6-2 | 3.6669146 | 0.0357875 | CD4\_T\_Naive\_CM\_CAM |
| MT-ND6  | 1.4755540 | 0.0268463 | CD4\_T\_Naive\_CM\_CAM |
| MT-CYB  | 1.4087794 | 0.0268463 | CD4\_T\_Naive\_CM\_CAM |
| MT-ND3  | 1.3459895 | 0.0268463 | CD4\_T\_Naive\_CM\_CAM |
| HLA-C   | 0.9752889 | 0.0460393 | CD4\_T\_Naive\_CM\_CAM |

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type             |
|:------|----------:|----------:|:-----------------------|
| CRIP2 | -3.092692 | 0.0357875 | CD4\_T\_Naive\_CM\_CAM |
| HES4  | -4.682519 | 0.0000916 | CD4\_T\_Naive\_CM\_CAM |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes     |     logFC |      padj | cell\_type   |
|:----------|----------:|----------:|:-------------|
| TRBV20-1  | 4.3226204 | 0.0128701 | CD4\_Th\_CAM |
| ANKRD27   | 1.7924659 | 0.0105059 | CD4\_Th\_CAM |
| ATP8A1    | 1.6338894 | 0.0148762 | CD4\_Th\_CAM |
| RBM28     | 1.5403128 | 0.0413684 | CD4\_Th\_CAM |
| ERAP2     | 1.5051277 | 0.0083821 | CD4\_Th\_CAM |
| SOCS2     | 1.5015222 | 0.0139306 | CD4\_Th\_CAM |
| MT-ND6    | 1.4771506 | 0.0179417 | CD4\_Th\_CAM |
| PIK3R5    | 1.3969550 | 0.0091710 | CD4\_Th\_CAM |
| MT-ND5    | 1.3667333 | 0.0316498 | CD4\_Th\_CAM |
| MT-ATP8   | 1.3357125 | 0.0186761 | CD4\_Th\_CAM |
| PATJ      | 1.2255987 | 0.0150281 | CD4\_Th\_CAM |
| IVNS1ABP  | 1.2179013 | 0.0283428 | CD4\_Th\_CAM |
| LINC01934 | 1.1745117 | 0.0495249 | CD4\_Th\_CAM |
| TRERF1    | 1.1359419 | 0.0458590 | CD4\_Th\_CAM |
| LINC01619 | 1.1212107 | 0.0209561 | CD4\_Th\_CAM |
| PDXDC1    | 1.1018241 | 0.0490511 | CD4\_Th\_CAM |
| VPS13B    | 1.0503448 | 0.0475805 | CD4\_Th\_CAM |
| UBQLN2    | 0.9767986 | 0.0432429 | CD4\_Th\_CAM |
| TAF15     | 0.8712943 | 0.0316498 | CD4\_Th\_CAM |

\[1\] “Genes downregulated in disease”

| genes      |      logFC |      padj | cell\_type   |
|:-----------|-----------:|----------:|:-------------|
| COX7B      | -0.7103103 | 0.0432429 | CD4\_Th\_CAM |
| COX7A2     | -0.7184952 | 0.0432429 | CD4\_Th\_CAM |
| MYL6       | -0.7867745 | 0.0444865 | CD4\_Th\_CAM |
| UPF2       | -0.8041596 | 0.0475805 | CD4\_Th\_CAM |
| UBL5       | -0.8252123 | 0.0134568 | CD4\_Th\_CAM |
| SEC61G     | -0.8344256 | 0.0440647 | CD4\_Th\_CAM |
| AC138123.1 | -0.8365479 | 0.0432429 | CD4\_Th\_CAM |
| SDF2L1     | -0.8640645 | 0.0471502 | CD4\_Th\_CAM |
| NDUFA3     | -0.8944363 | 0.0179417 | CD4\_Th\_CAM |
| SNRPD3     | -0.8999723 | 0.0177685 | CD4\_Th\_CAM |
| TRBC1      | -0.9007194 | 0.0269205 | CD4\_Th\_CAM |
| PSME2      | -0.9307750 | 0.0083821 | CD4\_Th\_CAM |
| HSPE1      | -0.9767643 | 0.0017836 | CD4\_Th\_CAM |
| CCDC167    | -1.0160397 | 0.0239844 | CD4\_Th\_CAM |
| MGST3      | -1.0463463 | 0.0017836 | CD4\_Th\_CAM |
| CKLF       | -1.1317170 | 0.0134568 | CD4\_Th\_CAM |
| DNPH1      | -1.1489260 | 0.0177685 | CD4\_Th\_CAM |
| HSPB1      | -1.1850838 | 0.0003198 | CD4\_Th\_CAM |
| FKBP2      | -1.2061004 | 0.0316498 | CD4\_Th\_CAM |
| CAPG       | -1.2711484 | 0.0300960 | CD4\_Th\_CAM |
| HSPA1A     | -1.2745015 | 0.0413684 | CD4\_Th\_CAM |
| CCL5       | -1.2894287 | 0.0413684 | CD4\_Th\_CAM |
| HOPX       | -1.3912608 | 0.0416851 | CD4\_Th\_CAM |
| LINC00892  | -1.4758142 | 0.0490511 | CD4\_Th\_CAM |
| CD74       | -1.5002865 | 0.0002320 | CD4\_Th\_CAM |
| CSRP1      | -1.5333125 | 0.0170578 | CD4\_Th\_CAM |
| LGALS9     | -3.1370427 | 0.0239628 | CD4\_Th\_CAM |
| TRBV6-2    | -3.1736034 | 0.0413684 | CD4\_Th\_CAM |
| TRAV3      | -4.4187818 | 0.0017836 | CD4\_Th\_CAM |
| HSPA6      | -6.0085924 | 0.0017836 | CD4\_Th\_CAM |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes       |     logFC |      padj | cell\_type      |
|:------------|----------:|----------:|:----------------|
| TRBV4-2     | 4.8433526 | 0.0183449 | FoxP3-Treg\_CAM |
| TRBV7-3     | 4.0788064 | 0.0183449 | FoxP3-Treg\_CAM |
| TRBV7-9     | 3.8867200 | 0.0412205 | FoxP3-Treg\_CAM |
| TRBV20-1    | 3.8255874 | 0.0167599 | FoxP3-Treg\_CAM |
| TRAV38-2DV8 | 3.6050297 | 0.0202452 | FoxP3-Treg\_CAM |
| TRBV7-6     | 3.5060317 | 0.0492866 | FoxP3-Treg\_CAM |
| CADM1       | 2.4540800 | 0.0149197 | FoxP3-Treg\_CAM |
| TRBV19      | 2.1901633 | 0.0312118 | FoxP3-Treg\_CAM |
| MT-CYB      | 2.0169812 | 0.0186468 | FoxP3-Treg\_CAM |
| MT-ND6      | 1.8535477 | 0.0073021 | FoxP3-Treg\_CAM |
| MT-ATP8     | 1.6197343 | 0.0371953 | FoxP3-Treg\_CAM |
| PAN3        | 1.4369515 | 0.0371953 | FoxP3-Treg\_CAM |
| GAS5        | 1.4235573 | 0.0088081 | FoxP3-Treg\_CAM |
| PIK3R5      | 1.3695758 | 0.0263618 | FoxP3-Treg\_CAM |
| BCAS3       | 1.3632868 | 0.0396313 | FoxP3-Treg\_CAM |
| MT-ND5      | 1.3283073 | 0.0492866 | FoxP3-Treg\_CAM |
| CXCR4       | 1.2113191 | 0.0088081 | FoxP3-Treg\_CAM |
| TXNIP       | 0.9371201 | 0.0422812 | FoxP3-Treg\_CAM |

\[1\] “Genes downregulated in disease”

| genes      |     logFC |      padj | cell\_type      |
|:-----------|----------:|----------:|:----------------|
| CTSB       | -1.109939 | 0.0183449 | FoxP3-Treg\_CAM |
| ATOX1      | -1.325961 | 0.0181960 | FoxP3-Treg\_CAM |
| AC093512.2 | -1.350904 | 0.0114938 | FoxP3-Treg\_CAM |
| ARG2       | -1.580555 | 0.0223537 | FoxP3-Treg\_CAM |
| HLA-DRB5   | -1.983633 | 0.0073021 | FoxP3-Treg\_CAM |
| CCL4L2     | -2.811218 | 0.0088081 | FoxP3-Treg\_CAM |
| TRAV26-1   | -3.236527 | 0.0043603 | FoxP3-Treg\_CAM |
| TRBV6-2    | -4.005152 | 0.0058008 | FoxP3-Treg\_CAM |
| TRBV28     | -4.132702 | 0.0073021 | FoxP3-Treg\_CAM |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes      |     logFC |      padj | cell\_type    |
|:-----------|----------:|----------:|:--------------|
| TRBV20-1   | 4.3898113 | 0.0241529 | CD4\_TEM\_CAM |
| TRBV7-9    | 4.2706804 | 0.0324167 | CD4\_TEM\_CAM |
| TRAV20     | 3.9138235 | 0.0465708 | CD4\_TEM\_CAM |
| TRAV12-2   | 3.8381431 | 0.0496661 | CD4\_TEM\_CAM |
| TRBV27     | 3.5937953 | 0.0114373 | CD4\_TEM\_CAM |
| FMN1       | 3.4887529 | 0.0465708 | CD4\_TEM\_CAM |
| TRBV6-5    | 3.4580335 | 0.0465708 | CD4\_TEM\_CAM |
| TRAV12-1   | 2.9518925 | 0.0251778 | CD4\_TEM\_CAM |
| DHTKD1     | 2.3075300 | 0.0040495 | CD4\_TEM\_CAM |
| MMS22L     | 2.2867033 | 0.0137792 | CD4\_TEM\_CAM |
| MT-CYB     | 1.9591102 | 0.0452955 | CD4\_TEM\_CAM |
| MT-ND6     | 1.9576008 | 0.0010670 | CD4\_TEM\_CAM |
| GSAP       | 1.8997775 | 0.0050565 | CD4\_TEM\_CAM |
| NCAPG2     | 1.7482100 | 0.0227896 | CD4\_TEM\_CAM |
| MT-ND5     | 1.6817450 | 0.0046829 | CD4\_TEM\_CAM |
| FGFR1      | 1.6490054 | 0.0496661 | CD4\_TEM\_CAM |
| MT-ATP8    | 1.6031477 | 0.0162793 | CD4\_TEM\_CAM |
| SPATA20    | 1.5223593 | 0.0236521 | CD4\_TEM\_CAM |
| AC016831.7 | 1.5003135 | 0.0137905 | CD4\_TEM\_CAM |
| ERAP2      | 1.4623938 | 0.0010670 | CD4\_TEM\_CAM |
| TRIM23     | 1.4502009 | 0.0137905 | CD4\_TEM\_CAM |
| AL136456.1 | 1.4426069 | 0.0245624 | CD4\_TEM\_CAM |
| PFKFB3     | 1.4171136 | 0.0208126 | CD4\_TEM\_CAM |
| BCL2L11    | 1.4074434 | 0.0332180 | CD4\_TEM\_CAM |
| ATP8A1     | 1.4004916 | 0.0074819 | CD4\_TEM\_CAM |
| MFN1       | 1.3782265 | 0.0200450 | CD4\_TEM\_CAM |
| C15orf39   | 1.3604199 | 0.0262517 | CD4\_TEM\_CAM |
| SIK2       | 1.3553360 | 0.0145922 | CD4\_TEM\_CAM |
| FBRSL1     | 1.2921573 | 0.0145922 | CD4\_TEM\_CAM |
| LIMD1      | 1.2917392 | 0.0416508 | CD4\_TEM\_CAM |
| TMEM168    | 1.2851341 | 0.0132415 | CD4\_TEM\_CAM |
| CRACR2A    | 1.2786888 | 0.0109532 | CD4\_TEM\_CAM |
| CFH        | 1.2759836 | 0.0263677 | CD4\_TEM\_CAM |
| APBB1      | 1.2570424 | 0.0026844 | CD4\_TEM\_CAM |
| OBSCN      | 1.2514768 | 0.0482317 | CD4\_TEM\_CAM |
| PCAT1      | 1.2491280 | 0.0010670 | CD4\_TEM\_CAM |
| MT-ND4L    | 1.2409073 | 0.0495693 | CD4\_TEM\_CAM |
| WDR36      | 1.2157664 | 0.0128509 | CD4\_TEM\_CAM |
| NR3C2      | 1.1688297 | 0.0123638 | CD4\_TEM\_CAM |
| PIGB       | 1.1614648 | 0.0241108 | CD4\_TEM\_CAM |
| GAS5       | 1.1491902 | 0.0208126 | CD4\_TEM\_CAM |
| SUSD6      | 1.1252527 | 0.0114373 | CD4\_TEM\_CAM |
| MT-ND2     | 1.0913617 | 0.0425940 | CD4\_TEM\_CAM |
| CHD7       | 1.0899259 | 0.0434934 | CD4\_TEM\_CAM |
| TET2       | 1.0769024 | 0.0224194 | CD4\_TEM\_CAM |
| DCAF17     | 1.0674479 | 0.0483819 | CD4\_TEM\_CAM |
| PTPRM      | 1.0577482 | 0.0144823 | CD4\_TEM\_CAM |
| MAP2K5     | 1.0505016 | 0.0478620 | CD4\_TEM\_CAM |
| PTCD3      | 1.0484599 | 0.0245624 | CD4\_TEM\_CAM |
| AC010609.1 | 1.0478116 | 0.0416508 | CD4\_TEM\_CAM |
| DIS3L2     | 1.0314170 | 0.0276235 | CD4\_TEM\_CAM |
| SLAMF6     | 1.0189053 | 0.0195615 | CD4\_TEM\_CAM |
| NUP107     | 1.0155272 | 0.0251778 | CD4\_TEM\_CAM |
| ZNF718     | 0.9921847 | 0.0469034 | CD4\_TEM\_CAM |
| DGKD       | 0.9854770 | 0.0396966 | CD4\_TEM\_CAM |
| KANSL3     | 0.9817526 | 0.0381249 | CD4\_TEM\_CAM |
| MADD       | 0.9786027 | 0.0251778 | CD4\_TEM\_CAM |
| FYCO1      | 0.9781519 | 0.0108039 | CD4\_TEM\_CAM |
| USP19      | 0.9630885 | 0.0417574 | CD4\_TEM\_CAM |
| KLF10      | 0.9610961 | 0.0422350 | CD4\_TEM\_CAM |
| KDM3B      | 0.9607772 | 0.0241067 | CD4\_TEM\_CAM |
| XYLT1      | 0.9589158 | 0.0055866 | CD4\_TEM\_CAM |
| TXNIP      | 0.9570941 | 0.0108039 | CD4\_TEM\_CAM |
| ZNF236     | 0.9535283 | 0.0387421 | CD4\_TEM\_CAM |
| TYK2       | 0.9457302 | 0.0367061 | CD4\_TEM\_CAM |
| RUNX2      | 0.9412711 | 0.0074819 | CD4\_TEM\_CAM |
| PDE3B      | 0.9214374 | 0.0206826 | CD4\_TEM\_CAM |
| AUTS2      | 0.9209237 | 0.0100385 | CD4\_TEM\_CAM |
| AKT3       | 0.9101886 | 0.0019693 | CD4\_TEM\_CAM |
| TMEM259    | 0.9032330 | 0.0033135 | CD4\_TEM\_CAM |
| TRAF3      | 0.8995765 | 0.0442094 | CD4\_TEM\_CAM |
| RFX3       | 0.8837793 | 0.0465708 | CD4\_TEM\_CAM |
| MAPK8      | 0.8768373 | 0.0478620 | CD4\_TEM\_CAM |
| IL7R       | 0.8756571 | 0.0230573 | CD4\_TEM\_CAM |
| SLMAP      | 0.8733478 | 0.0145922 | CD4\_TEM\_CAM |
| ABCC1      | 0.8719263 | 0.0284369 | CD4\_TEM\_CAM |
| MAML2      | 0.8509256 | 0.0004864 | CD4\_TEM\_CAM |
| FNDC3A     | 0.8324059 | 0.0417574 | CD4\_TEM\_CAM |
| GNAQ       | 0.8318905 | 0.0398969 | CD4\_TEM\_CAM |
| CBLB       | 0.8316066 | 0.0003164 | CD4\_TEM\_CAM |
| RASA1      | 0.8292982 | 0.0244083 | CD4\_TEM\_CAM |
| AC068587.4 | 0.8151523 | 0.0128509 | CD4\_TEM\_CAM |
| LMBR1      | 0.8147842 | 0.0244083 | CD4\_TEM\_CAM |
| PCF11      | 0.8082042 | 0.0208126 | CD4\_TEM\_CAM |
| PIK3R5     | 0.8079602 | 0.0200075 | CD4\_TEM\_CAM |
| VPS13B     | 0.8079348 | 0.0100385 | CD4\_TEM\_CAM |
| FBXL17     | 0.8036590 | 0.0137905 | CD4\_TEM\_CAM |
| VWA8       | 0.8018194 | 0.0245624 | CD4\_TEM\_CAM |
| ARL15      | 0.7616004 | 0.0416508 | CD4\_TEM\_CAM |
| CDC14A     | 0.7572101 | 0.0145922 | CD4\_TEM\_CAM |
| KIAA2013   | 0.7525087 | 0.0434934 | CD4\_TEM\_CAM |
| SMYD3      | 0.7475526 | 0.0046829 | CD4\_TEM\_CAM |
| RASGRP2    | 0.7333971 | 0.0442457 | CD4\_TEM\_CAM |
| SYTL3      | 0.7303800 | 0.0182198 | CD4\_TEM\_CAM |
| CXCR4      | 0.7277860 | 0.0473854 | CD4\_TEM\_CAM |
| DPYD       | 0.7259633 | 0.0040748 | CD4\_TEM\_CAM |
| RICTOR     | 0.7253902 | 0.0213518 | CD4\_TEM\_CAM |
| OGT        | 0.7083455 | 0.0059367 | CD4\_TEM\_CAM |
| VPS13D     | 0.7080964 | 0.0358836 | CD4\_TEM\_CAM |
| RALGAPA1   | 0.7043905 | 0.0271012 | CD4\_TEM\_CAM |
| HECTD1     | 0.7029922 | 0.0434934 | CD4\_TEM\_CAM |
| VTI1A      | 0.6955127 | 0.0241529 | CD4\_TEM\_CAM |
| CPSF6      | 0.6768833 | 0.0378251 | CD4\_TEM\_CAM |
| CAMK2D     | 0.6733162 | 0.0400122 | CD4\_TEM\_CAM |
| RIPOR2     | 0.6706396 | 0.0016971 | CD4\_TEM\_CAM |
| AP3B1      | 0.6618165 | 0.0202503 | CD4\_TEM\_CAM |
| LRBA       | 0.6607503 | 0.0208126 | CD4\_TEM\_CAM |
| XIST       | 0.6307739 | 0.0473829 | CD4\_TEM\_CAM |
| PRDM2      | 0.6266970 | 0.0370476 | CD4\_TEM\_CAM |
| CAPRIN1    | 0.6011134 | 0.0243436 | CD4\_TEM\_CAM |
| ASH1L      | 0.5864646 | 0.0365762 | CD4\_TEM\_CAM |
| TTC17      | 0.5513321 | 0.0326986 | CD4\_TEM\_CAM |
| KDM5A      | 0.5437709 | 0.0450467 | CD4\_TEM\_CAM |
| PPP3CC     | 0.5344085 | 0.0451563 | CD4\_TEM\_CAM |
| RPS6KA3    | 0.5194528 | 0.0385965 | CD4\_TEM\_CAM |
| FOXN3      | 0.5166050 | 0.0208126 | CD4\_TEM\_CAM |
| ZC3HAV1    | 0.4648978 | 0.0497097 | CD4\_TEM\_CAM |

\[1\] “Genes downregulated in disease”

| genes       |      logFC |      padj | cell\_type    |
|:------------|-----------:|----------:|:--------------|
| STK17A      | -0.4620643 | 0.0483819 | CD4\_TEM\_CAM |
| GUK1        | -0.4710418 | 0.0416508 | CD4\_TEM\_CAM |
| RNF167      | -0.4720628 | 0.0463524 | CD4\_TEM\_CAM |
| HCST        | -0.4721817 | 0.0483819 | CD4\_TEM\_CAM |
| NDUFB8      | -0.4802584 | 0.0426966 | CD4\_TEM\_CAM |
| GNG5        | -0.4902376 | 0.0497097 | CD4\_TEM\_CAM |
| CYBA        | -0.4949073 | 0.0416508 | CD4\_TEM\_CAM |
| SELENOH     | -0.4983602 | 0.0416508 | CD4\_TEM\_CAM |
| TSPO        | -0.5072921 | 0.0450467 | CD4\_TEM\_CAM |
| CALM3       | -0.5084453 | 0.0483889 | CD4\_TEM\_CAM |
| NDUFA4      | -0.5146353 | 0.0422350 | CD4\_TEM\_CAM |
| MAGOH       | -0.5160036 | 0.0442094 | CD4\_TEM\_CAM |
| TMEM50A     | -0.5204247 | 0.0309649 | CD4\_TEM\_CAM |
| CMTM3       | -0.5228712 | 0.0363234 | CD4\_TEM\_CAM |
| MRPS6       | -0.5251009 | 0.0434934 | CD4\_TEM\_CAM |
| WDR83OS     | -0.5281718 | 0.0320354 | CD4\_TEM\_CAM |
| PSME1       | -0.5323872 | 0.0484060 | CD4\_TEM\_CAM |
| PSMB10      | -0.5403477 | 0.0199761 | CD4\_TEM\_CAM |
| PGAM1       | -0.5404028 | 0.0365455 | CD4\_TEM\_CAM |
| NDUFAB1     | -0.5446976 | 0.0418472 | CD4\_TEM\_CAM |
| COX14       | -0.5463305 | 0.0395672 | CD4\_TEM\_CAM |
| TXNDC12     | -0.5465858 | 0.0416508 | CD4\_TEM\_CAM |
| SERF2       | -0.5478230 | 0.0241067 | CD4\_TEM\_CAM |
| NDUFA13     | -0.5502896 | 0.0434934 | CD4\_TEM\_CAM |
| NDUFB3      | -0.5511878 | 0.0268514 | CD4\_TEM\_CAM |
| CIB1        | -0.5512515 | 0.0147304 | CD4\_TEM\_CAM |
| COPZ1       | -0.5544050 | 0.0499450 | CD4\_TEM\_CAM |
| NDUFA6      | -0.5586589 | 0.0465708 | CD4\_TEM\_CAM |
| NDUFA1      | -0.5605261 | 0.0200450 | CD4\_TEM\_CAM |
| COX6A1      | -0.5606209 | 0.0245624 | CD4\_TEM\_CAM |
| ACTB        | -0.5627989 | 0.0367061 | CD4\_TEM\_CAM |
| IAH1        | -0.5628958 | 0.0483819 | CD4\_TEM\_CAM |
| PSMA5       | -0.5630311 | 0.0185578 | CD4\_TEM\_CAM |
| UBB         | -0.5670859 | 0.0200075 | CD4\_TEM\_CAM |
| SEM1        | -0.5696699 | 0.0230573 | CD4\_TEM\_CAM |
| NDUFA12     | -0.5714998 | 0.0422350 | CD4\_TEM\_CAM |
| OCIAD2      | -0.5717781 | 0.0123196 | CD4\_TEM\_CAM |
| CIAO2B      | -0.5741706 | 0.0145922 | CD4\_TEM\_CAM |
| ACTG1       | -0.5767667 | 0.0332180 | CD4\_TEM\_CAM |
| AURKAIP1    | -0.5788928 | 0.0162793 | CD4\_TEM\_CAM |
| STOML2      | -0.5811190 | 0.0463524 | CD4\_TEM\_CAM |
| MIF         | -0.5818341 | 0.0233146 | CD4\_TEM\_CAM |
| ATP6V1F     | -0.5830302 | 0.0232114 | CD4\_TEM\_CAM |
| JPT1        | -0.5832734 | 0.0244083 | CD4\_TEM\_CAM |
| NDUFC1      | -0.5865619 | 0.0221027 | CD4\_TEM\_CAM |
| SNRPG       | -0.5872970 | 0.0205445 | CD4\_TEM\_CAM |
| BRK1        | -0.5891746 | 0.0255555 | CD4\_TEM\_CAM |
| CLIC1       | -0.5922665 | 0.0195615 | CD4\_TEM\_CAM |
| BANF1       | -0.5928195 | 0.0147304 | CD4\_TEM\_CAM |
| SUMO1       | -0.5939966 | 0.0344623 | CD4\_TEM\_CAM |
| CD3D        | -0.5952928 | 0.0145922 | CD4\_TEM\_CAM |
| ATP5MC1     | -0.5982282 | 0.0202503 | CD4\_TEM\_CAM |
| ABRACL      | -0.5984874 | 0.0221748 | CD4\_TEM\_CAM |
| COMMD8      | -0.5991500 | 0.0152808 | CD4\_TEM\_CAM |
| DCTN3       | -0.5995706 | 0.0431931 | CD4\_TEM\_CAM |
| VDAC1       | -0.6021860 | 0.0208126 | CD4\_TEM\_CAM |
| APEX1       | -0.6023806 | 0.0382593 | CD4\_TEM\_CAM |
| MRPL54      | -0.6043154 | 0.0206826 | CD4\_TEM\_CAM |
| PSMB6       | -0.6044074 | 0.0173540 | CD4\_TEM\_CAM |
| DAD1        | -0.6058119 | 0.0230573 | CD4\_TEM\_CAM |
| RBCK1       | -0.6090630 | 0.0236521 | CD4\_TEM\_CAM |
| TIMM17B     | -0.6099755 | 0.0358489 | CD4\_TEM\_CAM |
| NDUFA2      | -0.6105173 | 0.0199761 | CD4\_TEM\_CAM |
| LAMTOR5     | -0.6107890 | 0.0252697 | CD4\_TEM\_CAM |
| MYL12A      | -0.6114652 | 0.0244083 | CD4\_TEM\_CAM |
| CD40LG      | -0.6119271 | 0.0075598 | CD4\_TEM\_CAM |
| UCHL3       | -0.6128560 | 0.0469034 | CD4\_TEM\_CAM |
| SELENOK     | -0.6139424 | 0.0319829 | CD4\_TEM\_CAM |
| PSMA2       | -0.6141906 | 0.0200450 | CD4\_TEM\_CAM |
| MRPL24      | -0.6150401 | 0.0441182 | CD4\_TEM\_CAM |
| ATP6V0B     | -0.6150585 | 0.0251778 | CD4\_TEM\_CAM |
| ANAPC11     | -0.6151024 | 0.0137905 | CD4\_TEM\_CAM |
| DSTN        | -0.6176721 | 0.0241529 | CD4\_TEM\_CAM |
| ELOC        | -0.6194318 | 0.0244083 | CD4\_TEM\_CAM |
| UBL5        | -0.6199812 | 0.0200450 | CD4\_TEM\_CAM |
| ELOB        | -0.6236044 | 0.0137792 | CD4\_TEM\_CAM |
| HSBP1       | -0.6242896 | 0.0252635 | CD4\_TEM\_CAM |
| NOP10       | -0.6258729 | 0.0232114 | CD4\_TEM\_CAM |
| COX8A       | -0.6285293 | 0.0350818 | CD4\_TEM\_CAM |
| APOBEC3H    | -0.6294251 | 0.0422350 | CD4\_TEM\_CAM |
| S100A11     | -0.6303534 | 0.0350818 | CD4\_TEM\_CAM |
| MZT2A       | -0.6304551 | 0.0416508 | CD4\_TEM\_CAM |
| COX17       | -0.6316992 | 0.0247947 | CD4\_TEM\_CAM |
| HSPE1       | -0.6334951 | 0.0147304 | CD4\_TEM\_CAM |
| CAPZA2      | -0.6345698 | 0.0245624 | CD4\_TEM\_CAM |
| POLR2C      | -0.6357825 | 0.0389018 | CD4\_TEM\_CAM |
| METTL5      | -0.6413016 | 0.0350633 | CD4\_TEM\_CAM |
| MIEN1       | -0.6445956 | 0.0244083 | CD4\_TEM\_CAM |
| IFITM2      | -0.6559776 | 0.0200450 | CD4\_TEM\_CAM |
| PIN4        | -0.6559969 | 0.0241067 | CD4\_TEM\_CAM |
| CD48        | -0.6569943 | 0.0195615 | CD4\_TEM\_CAM |
| PSME2       | -0.6587522 | 0.0029011 | CD4\_TEM\_CAM |
| TXN         | -0.6600279 | 0.0108039 | CD4\_TEM\_CAM |
| ARPC5L      | -0.6600372 | 0.0062424 | CD4\_TEM\_CAM |
| COTL1       | -0.6602595 | 0.0116279 | CD4\_TEM\_CAM |
| CNIH4       | -0.6647200 | 0.0416508 | CD4\_TEM\_CAM |
| H2AFZ       | -0.6668756 | 0.0166967 | CD4\_TEM\_CAM |
| SEC61G      | -0.6690062 | 0.0253721 | CD4\_TEM\_CAM |
| TRMT112     | -0.6698329 | 0.0100385 | CD4\_TEM\_CAM |
| GLIPR2      | -0.6709802 | 0.0108039 | CD4\_TEM\_CAM |
| ATP5MD      | -0.6711177 | 0.0289660 | CD4\_TEM\_CAM |
| JOSD2       | -0.6732083 | 0.0134447 | CD4\_TEM\_CAM |
| SDF2L1      | -0.6735728 | 0.0106021 | CD4\_TEM\_CAM |
| COX7A2      | -0.6736442 | 0.0055866 | CD4\_TEM\_CAM |
| FKBP2       | -0.6755594 | 0.0200450 | CD4\_TEM\_CAM |
| BTF3L4      | -0.6768334 | 0.0194815 | CD4\_TEM\_CAM |
| ZNHIT1      | -0.6771994 | 0.0137792 | CD4\_TEM\_CAM |
| S100A6      | -0.6790246 | 0.0294530 | CD4\_TEM\_CAM |
| LSM6        | -0.6820164 | 0.0055866 | CD4\_TEM\_CAM |
| TPI1        | -0.6825824 | 0.0078734 | CD4\_TEM\_CAM |
| TMEM134     | -0.6903394 | 0.0326986 | CD4\_TEM\_CAM |
| MRPL33      | -0.6903614 | 0.0241529 | CD4\_TEM\_CAM |
| HPCAL1      | -0.6961123 | 0.0115697 | CD4\_TEM\_CAM |
| OSTC        | -0.6986240 | 0.0074819 | CD4\_TEM\_CAM |
| CHCHD1      | -0.6990002 | 0.0344623 | CD4\_TEM\_CAM |
| ATP5MPL     | -0.7008426 | 0.0137905 | CD4\_TEM\_CAM |
| MITD1       | -0.7021353 | 0.0483819 | CD4\_TEM\_CAM |
| MICOS10     | -0.7023969 | 0.0015952 | CD4\_TEM\_CAM |
| DBI         | -0.7037288 | 0.0208126 | CD4\_TEM\_CAM |
| PDCL3       | -0.7043611 | 0.0208126 | CD4\_TEM\_CAM |
| RNASEK      | -0.7048035 | 0.0163680 | CD4\_TEM\_CAM |
| THYN1       | -0.7074797 | 0.0108039 | CD4\_TEM\_CAM |
| RGL4        | -0.7078076 | 0.0326986 | CD4\_TEM\_CAM |
| PLAAT4      | -0.7117288 | 0.0241067 | CD4\_TEM\_CAM |
| ATP5F1E     | -0.7127938 | 0.0245624 | CD4\_TEM\_CAM |
| POLR3K      | -0.7143602 | 0.0230573 | CD4\_TEM\_CAM |
| GAPDH       | -0.7186862 | 0.0055866 | CD4\_TEM\_CAM |
| PET100      | -0.7190826 | 0.0307774 | CD4\_TEM\_CAM |
| UBE2S       | -0.7192095 | 0.0081849 | CD4\_TEM\_CAM |
| MGST3       | -0.7289394 | 0.0108039 | CD4\_TEM\_CAM |
| PHPT1       | -0.7298188 | 0.0489090 | CD4\_TEM\_CAM |
| CISD2       | -0.7323032 | 0.0202503 | CD4\_TEM\_CAM |
| SELENOM     | -0.7326330 | 0.0483889 | CD4\_TEM\_CAM |
| GSTP1       | -0.7416449 | 0.0125821 | CD4\_TEM\_CAM |
| SH3BGRL3    | -0.7433313 | 0.0046191 | CD4\_TEM\_CAM |
| YIF1A       | -0.7453221 | 0.0243436 | CD4\_TEM\_CAM |
| MYL6        | -0.7462453 | 0.0026844 | CD4\_TEM\_CAM |
| SUMO3       | -0.7538630 | 0.0241529 | CD4\_TEM\_CAM |
| LINC00892   | -0.7637405 | 0.0145922 | CD4\_TEM\_CAM |
| POMP        | -0.7665437 | 0.0046829 | CD4\_TEM\_CAM |
| FIS1        | -0.7699725 | 0.0128509 | CD4\_TEM\_CAM |
| DNPH1       | -0.7714263 | 0.0066103 | CD4\_TEM\_CAM |
| C4orf48     | -0.7742553 | 0.0128509 | CD4\_TEM\_CAM |
| ADA         | -0.7781889 | 0.0483819 | CD4\_TEM\_CAM |
| ATP5ME      | -0.7861153 | 0.0336962 | CD4\_TEM\_CAM |
| LDHA        | -0.7891181 | 0.0059367 | CD4\_TEM\_CAM |
| CISD3       | -0.7903855 | 0.0020221 | CD4\_TEM\_CAM |
| ALOX5AP     | -0.7951292 | 0.0123638 | CD4\_TEM\_CAM |
| GSTO1       | -0.8097212 | 0.0004864 | CD4\_TEM\_CAM |
| MRPL52      | -0.8230279 | 0.0145922 | CD4\_TEM\_CAM |
| NDUFB1      | -0.8400938 | 0.0145922 | CD4\_TEM\_CAM |
| BOLA3       | -0.8475568 | 0.0055866 | CD4\_TEM\_CAM |
| ERG28       | -0.8493720 | 0.0151090 | CD4\_TEM\_CAM |
| SNRNP35     | -0.8562881 | 0.0255555 | CD4\_TEM\_CAM |
| TMEM191C    | -0.8768603 | 0.0434934 | CD4\_TEM\_CAM |
| HSPB1       | -0.8908568 | 0.0015952 | CD4\_TEM\_CAM |
| DPH3        | -0.8935576 | 0.0227896 | CD4\_TEM\_CAM |
| HIGD1A      | -0.9058351 | 0.0109532 | CD4\_TEM\_CAM |
| NDUFA3      | -0.9068489 | 0.0033135 | CD4\_TEM\_CAM |
| CKLF        | -0.9082492 | 0.0013902 | CD4\_TEM\_CAM |
| LTA         | -0.9162338 | 0.0046829 | CD4\_TEM\_CAM |
| NKG7        | -0.9270639 | 0.0132415 | CD4\_TEM\_CAM |
| TMEM205     | -0.9354914 | 0.0245624 | CD4\_TEM\_CAM |
| RNF181      | -0.9467911 | 0.0004864 | CD4\_TEM\_CAM |
| NUDT1       | -0.9603207 | 0.0108039 | CD4\_TEM\_CAM |
| CCDC167     | -0.9787274 | 0.0051253 | CD4\_TEM\_CAM |
| HIST2H2AA4  | -1.0094930 | 0.0241067 | CD4\_TEM\_CAM |
| CBWD1       | -1.0118468 | 0.0465708 | CD4\_TEM\_CAM |
| C17orf49    | -1.0513911 | 0.0004515 | CD4\_TEM\_CAM |
| LGALS1      | -1.0994522 | 0.0010670 | CD4\_TEM\_CAM |
| IGFBP4      | -1.1004672 | 0.0190410 | CD4\_TEM\_CAM |
| ISOC2       | -1.1154199 | 0.0074819 | CD4\_TEM\_CAM |
| PTTG1       | -1.1348833 | 0.0015952 | CD4\_TEM\_CAM |
| HLA-DPA1    | -1.1397841 | 0.0139494 | CD4\_TEM\_CAM |
| CITED4      | -1.2100881 | 0.0000018 | CD4\_TEM\_CAM |
| RGS16       | -1.2220482 | 0.0145922 | CD4\_TEM\_CAM |
| AC004585.1  | -1.2502695 | 0.0064762 | CD4\_TEM\_CAM |
| CD74        | -1.2524634 | 0.0021784 | CD4\_TEM\_CAM |
| BHLHE40-AS1 | -1.3004534 | 0.0312675 | CD4\_TEM\_CAM |
| C1orf162    | -1.3314270 | 0.0012131 | CD4\_TEM\_CAM |
| IFNG        | -1.4141282 | 0.0307774 | CD4\_TEM\_CAM |
| HLA-DMA     | -1.4484486 | 0.0027291 | CD4\_TEM\_CAM |
| GZMH        | -1.5862962 | 0.0051253 | CD4\_TEM\_CAM |
| XCL1        | -1.6148174 | 0.0000000 | CD4\_TEM\_CAM |
| MXRA7       | -1.6255841 | 0.0004515 | CD4\_TEM\_CAM |
| HLA-DRB1    | -1.6595028 | 0.0362143 | CD4\_TEM\_CAM |
| ZNF683      | -1.8649956 | 0.0465708 | CD4\_TEM\_CAM |
| QSOX1       | -1.9069279 | 0.0208126 | CD4\_TEM\_CAM |
| CCR4        | -1.9798058 | 0.0037699 | CD4\_TEM\_CAM |
| CCL4L2      | -2.6655841 | 0.0004864 | CD4\_TEM\_CAM |
| FUCA2       | -3.2484646 | 0.0245624 | CD4\_TEM\_CAM |
| IGHV1-3     | -4.1160541 | 0.0004515 | CD4\_TEM\_CAM |
| HLA-DRB5    | -5.2233252 | 0.0000000 | CD4\_TEM\_CAM |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes    |     logFC |      padj | cell\_type |
|:---------|----------:|----------:|:-----------|
| ZSWIM3   | 4.2268053 | 0.0141522 | NKT\_CAM   |
| MT-ND6   | 1.8010189 | 0.0141522 | NKT\_CAM   |
| L3MBTL4  | 1.2408644 | 0.0362879 | NKT\_CAM   |
| SLC1A5   | 1.0701875 | 0.0271170 | NKT\_CAM   |
| HNRNPUL1 | 0.6782167 | 0.0303383 | NKT\_CAM   |

\[1\] “Genes downregulated in disease”

| genes      |      logFC |      padj | cell\_type |
|:-----------|-----------:|----------:|:-----------|
| CKS1B      | -0.9772416 | 0.0271170 | NKT\_CAM   |
| TYMS       | -1.0227065 | 0.0141522 | NKT\_CAM   |
| CCL4L2     | -2.0942372 | 0.0141522 | NKT\_CAM   |
| AC034238.1 | -2.5142030 | 0.0419187 | NKT\_CAM   |
| MYL9       | -4.1184665 | 0.0419187 | NKT\_CAM   |
| HLA-DRB5   | -5.0492280 | 0.0141522 | NKT\_CAM   |
| DEFA3      | -8.6197282 | 0.0141522 | NKT\_CAM   |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes      |    logFC |      padj | cell\_type          |
|:-----------|---------:|----------:|:--------------------|
| IL1R2      | 5.898331 | 0.0031443 | CD14\_Monocyte\_CAM |
| AL139042.1 | 5.621903 | 0.0014369 | CD14\_Monocyte\_CAM |
| SPRY1      | 5.399208 | 0.0074563 | CD14\_Monocyte\_CAM |
| MFGE8      | 5.059921 | 0.0081399 | CD14\_Monocyte\_CAM |
| ADAMTS2    | 4.604609 | 0.0000078 | CD14\_Monocyte\_CAM |
| ALOX15B    | 4.192387 | 0.0154709 | CD14\_Monocyte\_CAM |
| AMPH       | 4.015612 | 0.0129254 | CD14\_Monocyte\_CAM |
| VSIG4      | 3.237043 | 0.0000016 | CD14\_Monocyte\_CAM |
| TPST1      | 3.075667 | 0.0154709 | CD14\_Monocyte\_CAM |
| FMN1       | 2.414963 | 0.0496641 | CD14\_Monocyte\_CAM |
| FLT3       | 2.019883 | 0.0431139 | CD14\_Monocyte\_CAM |
| MMP17      | 1.813243 | 0.0079562 | CD14\_Monocyte\_CAM |
| MT-ND6     | 1.773858 | 0.0014369 | CD14\_Monocyte\_CAM |
| RBMS2      | 1.396969 | 0.0246518 | CD14\_Monocyte\_CAM |
| MS4A4E     | 1.366648 | 0.0154709 | CD14\_Monocyte\_CAM |
| SNTB1      | 1.122528 | 0.0129254 | CD14\_Monocyte\_CAM |
| SLC16A6    | 1.100043 | 0.0332511 | CD14\_Monocyte\_CAM |
| DOCK10     | 1.094543 | 0.0246518 | CD14\_Monocyte\_CAM |
| NDRG1      | 1.052717 | 0.0129254 | CD14\_Monocyte\_CAM |

\[1\] “Genes downregulated in disease”

| genes  |      logFC |      padj | cell\_type          |
|:-------|-----------:|----------:|:--------------------|
| SNRPD3 | -0.8773568 | 0.0496641 | CD14\_Monocyte\_CAM |
| BATF   | -1.1199919 | 0.0300671 | CD14\_Monocyte\_CAM |
| GYPC   | -1.2398782 | 0.0188255 | CD14\_Monocyte\_CAM |
| FCGR1B | -1.2918467 | 0.0129254 | CD14\_Monocyte\_CAM |
| ANPEP  | -1.2971653 | 0.0246518 | CD14\_Monocyte\_CAM |
| GK     | -1.3079640 | 0.0285060 | CD14\_Monocyte\_CAM |
| CCR5AS | -1.7207819 | 0.0431139 | CD14\_Monocyte\_CAM |
| RASA4  | -1.8795184 | 0.0129254 | CD14\_Monocyte\_CAM |
| HSPB1  | -1.9675383 | 0.0000078 | CD14\_Monocyte\_CAM |
| G0S2   | -1.9968194 | 0.0054052 | CD14\_Monocyte\_CAM |
| CCL4L2 | -3.0127782 | 0.0249589 | CD14\_Monocyte\_CAM |
| CLU    | -4.5489010 | 0.0327595 | CD14\_Monocyte\_CAM |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes  |    logFC |      padj | cell\_type          |
|:-------|---------:|----------:|:--------------------|
| DIRC3  | 2.655244 | 0.0165255 | CD16\_Monocyte\_CAM |
| TPST1  | 2.599893 | 0.0003853 | CD16\_Monocyte\_CAM |
| MT-ND6 | 2.406923 | 0.0019576 | CD16\_Monocyte\_CAM |

\[1\] “Genes downregulated in disease”

| genes    |     logFC |      padj | cell\_type          |
|:---------|----------:|----------:|:--------------------|
| IL2RG    | -1.816887 | 0.0013826 | CD16\_Monocyte\_CAM |
| CCL3L1   | -1.887989 | 0.0159343 | CD16\_Monocyte\_CAM |
| SLC39A8  | -1.922702 | 0.0177458 | CD16\_Monocyte\_CAM |
| IL4I1    | -1.940714 | 0.0002147 | CD16\_Monocyte\_CAM |
| PPA1     | -1.958321 | 0.0003252 | CD16\_Monocyte\_CAM |
| CCL3     | -2.062079 | 0.0062713 | CD16\_Monocyte\_CAM |
| VMO1     | -2.392184 | 0.0465401 | CD16\_Monocyte\_CAM |
| HAMP     | -2.450716 | 0.0087030 | CD16\_Monocyte\_CAM |
| G0S2     | -3.408864 | 0.0062713 | CD16\_Monocyte\_CAM |
| CCL2     | -3.449248 | 0.0007592 | CD16\_Monocyte\_CAM |
| C15orf48 | -3.457075 | 0.0003321 | CD16\_Monocyte\_CAM |
| CCL4L2   | -3.815135 | 0.0003321 | CD16\_Monocyte\_CAM |
| CCL4     | -4.391263 | 0.0002053 | CD16\_Monocyte\_CAM |
| CXCL9    | -4.893074 | 0.0088486 | CD16\_Monocyte\_CAM |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes |    logFC |      padj | cell\_type      |
|:------|---------:|----------:|:----------------|
| SOX5  | 3.113059 | 0.0231115 | Macrophage\_CAM |
| PDE7B | 2.915812 | 0.0239968 | Macrophage\_CAM |

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type      |
|:------|----------:|----------:|:----------------|
| FLT1  | -1.849048 | 0.0239968 | Macrophage\_CAM |
| CCL5  | -1.959507 | 0.0019012 | Macrophage\_CAM |
| MMP12 | -4.163902 | 0.0019012 | Macrophage\_CAM |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes |    logFC |      padj | cell\_type |
|:------|---------:|----------:|:-----------|
| DERL3 | 5.791319 | 0.0353442 | STB\_CAM   |
| FSTL3 | 5.412260 | 0.0448367 | STB\_CAM   |
| LPL   | 4.794174 | 0.0353442 | STB\_CAM   |

\[1\] “Genes downregulated in disease”

| genes      |     logFC |      padj | cell\_type |
|:-----------|----------:|----------:|:-----------|
| RNF213     | -2.939173 | 0.0353442 | STB\_CAM   |
| AC011287.1 | -3.699988 | 0.0128849 | STB\_CAM   |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes     |    logFC |      padj | cell\_type          |
|:----------|---------:|----------:|:--------------------|
| S100A4    | 5.579660 | 0.0264541 | Smooth\_muscle\_CAM |
| CORIN     | 4.074064 | 0.0106007 | Smooth\_muscle\_CAM |
| MTRNR2L12 | 3.261961 | 0.0022351 | Smooth\_muscle\_CAM |
| GRIA1     | 2.982155 | 0.0154891 | Smooth\_muscle\_CAM |

\[1\] “Genes downregulated in disease”

| genes     |     logFC |      padj | cell\_type          |
|:----------|----------:|----------:|:--------------------|
| TNFAIP3   | -3.105069 | 0.0008647 | Smooth\_muscle\_CAM |
| LINC01705 | -3.147153 | 0.0008157 | Smooth\_muscle\_CAM |
| FGF14     | -3.215604 | 0.0322359 | Smooth\_muscle\_CAM |
| NAMPT     | -3.478537 | 0.0004219 | Smooth\_muscle\_CAM |
| ADARB2    | -3.999782 | 0.0129698 | Smooth\_muscle\_CAM |
| PABPC1    | -4.092390 | 0.0028759 | Smooth\_muscle\_CAM |
| REN       | -4.530367 | 0.0003448 | Smooth\_muscle\_CAM |
| DNAH11    | -4.556316 | 0.0003154 | Smooth\_muscle\_CAM |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes    |    logFC |      padj | cell\_type       |
|:---------|---------:|----------:|:-----------------|
| PLXNA4   | 5.426034 | 0.0440968 | Endothelial\_CAM |
| ERAP2    | 2.238797 | 0.0440968 | Endothelial\_CAM |
| PRICKLE1 | 2.164776 | 0.0440968 | Endothelial\_CAM |

\[1\] “Genes downregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

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
