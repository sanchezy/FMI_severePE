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

> Tissue: Myometrium

``` r
# This notebook is for tissue and annotation level:
tissue= "Myometrium" # "CAM", "PBMC", "Myometrium"
annotation_level <- "CellTypeManual.l3"

# List of cell types
cellgroups$CellTypeManual.l3 <- paste0(cellgroups$CellTypeManual.l3,"_",tissue)
```

# Volcano plots in Early disease

\[1\] “Results for celltype:” \[1\] “Naive\_B\_Myometrium” \[1\] “Genes
upregulated in disease”

| genes  |     logFC |      padj | cell\_type           |
|:-------|----------:|----------:|:---------------------|
| MT-CO1 | 12.475472 | 0.0006748 | Naive\_B\_Myometrium |
| MT-ND3 | 11.183738 | 0.0006748 | Naive\_B\_Myometrium |
| MT-ND1 |  9.999432 | 0.0062615 | Naive\_B\_Myometrium |
| TPST1  |  6.166617 | 0.0257361 | Naive\_B\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Results for celltype:” \[1\] “Memory\_B\_Myometrium” \[1\] “Genes
upregulated in disease”

| genes   |    logFC |      padj | cell\_type            |
|:--------|---------:|----------:|:----------------------|
| MT-CO1  | 12.79300 | 0.0001107 | Memory\_B\_Myometrium |
| MT-ND3  | 11.08903 | 0.0001856 | Memory\_B\_Myometrium |
| MT-ND4L | 10.92389 | 0.0006056 | Memory\_B\_Myometrium |
| MT-ND1  | 10.25693 | 0.0049583 | Memory\_B\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Results for celltype:” \[1\] “CD4\_T\_Naive\_CM\_Myometrium”
\[1\] “Genes upregulated in disease”

| genes   |     logFC |      padj | cell\_type                    |
|:--------|----------:|----------:|:------------------------------|
| MT-CO1  | 13.155091 | 0.0002657 | CD4\_T\_Naive\_CM\_Myometrium |
| MT-ND3  | 12.122263 | 0.0001927 | CD4\_T\_Naive\_CM\_Myometrium |
| MT-ND4L | 11.820322 | 0.0022662 | CD4\_T\_Naive\_CM\_Myometrium |
| MT-ND1  | 11.350126 | 0.0043101 | CD4\_T\_Naive\_CM\_Myometrium |
| ACTN1   |  8.127173 | 0.0007842 | CD4\_T\_Naive\_CM\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes  |     logFC |      padj | cell\_type                    |
|:-------|----------:|----------:|:------------------------------|
| HSPA1A | -2.308478 | 0.0307107 | CD4\_T\_Naive\_CM\_Myometrium |
| NABP1  | -3.422327 | 0.0179016 | CD4\_T\_Naive\_CM\_Myometrium |

\[1\] “Results for celltype:” \[1\] “CD4\_Th\_Myometrium” \[1\] “Genes
upregulated in disease”

| genes      |      logFC |      padj | cell\_type          |
|:-----------|-----------:|----------:|:--------------------|
| MT-ND3     | 12.0370162 | 0.0000476 | CD4\_Th\_Myometrium |
| MT-ND4L    | 11.8142573 | 0.0004221 | CD4\_Th\_Myometrium |
| MT-ND1     | 11.2624697 | 0.0004432 | CD4\_Th\_Myometrium |
| MT-CO1     |  9.6347120 | 0.0000788 | CD4\_Th\_Myometrium |
| TRBV29-1   |  5.3998350 | 0.0004432 | CD4\_Th\_Myometrium |
| EDA        |  3.7157720 | 0.0104883 | CD4\_Th\_Myometrium |
| ACTN1      |  3.2881726 | 0.0008839 | CD4\_Th\_Myometrium |
| KLF2       |  2.5587780 | 0.0000672 | CD4\_Th\_Myometrium |
| LINC02762  |  2.5155879 | 0.0297638 | CD4\_Th\_Myometrium |
| ERAP2      |  2.2482461 | 0.0001618 | CD4\_Th\_Myometrium |
| AL035530.2 |  1.9292515 | 0.0076587 | CD4\_Th\_Myometrium |
| PIM2       |  1.9043450 | 0.0494101 | CD4\_Th\_Myometrium |
| LPAR6      |  1.8017203 | 0.0084234 | CD4\_Th\_Myometrium |
| TXK        |  1.7875363 | 0.0148354 | CD4\_Th\_Myometrium |
| CDC25B     |  1.7589409 | 0.0032969 | CD4\_Th\_Myometrium |
| DUSP1      |  1.7277478 | 0.0084123 | CD4\_Th\_Myometrium |
| MYADM      |  1.5611055 | 0.0084723 | CD4\_Th\_Myometrium |
| RBL1       |  1.5515558 | 0.0032969 | CD4\_Th\_Myometrium |
| TIMP1      |  1.5229620 | 0.0475070 | CD4\_Th\_Myometrium |
| UBQLN2     |  1.4966429 | 0.0492860 | CD4\_Th\_Myometrium |
| TCP11L2    |  1.4876182 | 0.0056071 | CD4\_Th\_Myometrium |
| CA5B       |  1.4774674 | 0.0056071 | CD4\_Th\_Myometrium |
| JADE1      |  1.3675196 | 0.0001171 | CD4\_Th\_Myometrium |
| FOXJ3      |  1.3560858 | 0.0014051 | CD4\_Th\_Myometrium |
| CDC14A     |  1.3529131 | 0.0407415 | CD4\_Th\_Myometrium |
| OXR1       |  1.3446624 | 0.0110906 | CD4\_Th\_Myometrium |
| TSPAN18    |  1.3108180 | 0.0469192 | CD4\_Th\_Myometrium |
| MKRN1      |  1.2067045 | 0.0055410 | CD4\_Th\_Myometrium |
| LDLRAP1    |  1.2008286 | 0.0444360 | CD4\_Th\_Myometrium |
| COG3       |  1.1814531 | 0.0132230 | CD4\_Th\_Myometrium |
| MYC        |  1.1670663 | 0.0071580 | CD4\_Th\_Myometrium |
| TRAT1      |  1.1610692 | 0.0444360 | CD4\_Th\_Myometrium |
| STK38      |  1.1577712 | 0.0020285 | CD4\_Th\_Myometrium |
| C12orf75   |  1.1493346 | 0.0486350 | CD4\_Th\_Myometrium |
| FBXL14     |  1.1371667 | 0.0148354 | CD4\_Th\_Myometrium |
| EGLN1      |  1.1307163 | 0.0056071 | CD4\_Th\_Myometrium |
| STT3B      |  1.0814171 | 0.0444360 | CD4\_Th\_Myometrium |
| IPO7       |  1.0680473 | 0.0014691 | CD4\_Th\_Myometrium |
| ZBTB44     |  1.0275114 | 0.0452938 | CD4\_Th\_Myometrium |
| WDR36      |  1.0041124 | 0.0108167 | CD4\_Th\_Myometrium |
| ALKBH5     |  0.9888786 | 0.0059403 | CD4\_Th\_Myometrium |
| ADD3       |  0.9855459 | 0.0147758 | CD4\_Th\_Myometrium |
| SORBS3     |  0.9819598 | 0.0040190 | CD4\_Th\_Myometrium |
| MOCS2      |  0.9817281 | 0.0294089 | CD4\_Th\_Myometrium |
| LEF1       |  0.9501366 | 0.0386568 | CD4\_Th\_Myometrium |
| CAPRIN1    |  0.9352957 | 0.0014691 | CD4\_Th\_Myometrium |
| VIRMA      |  0.8952733 | 0.0262903 | CD4\_Th\_Myometrium |
| ZNF800     |  0.8499292 | 0.0274768 | CD4\_Th\_Myometrium |
| METTL16    |  0.8248594 | 0.0444360 | CD4\_Th\_Myometrium |
| LNPEP      |  0.8034579 | 0.0274768 | CD4\_Th\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes     |      logFC |      padj | cell\_type          |
|:----------|-----------:|----------:|:--------------------|
| IDH2      | -0.6948338 | 0.0092517 | CD4\_Th\_Myometrium |
| EIF6      | -0.7750767 | 0.0274768 | CD4\_Th\_Myometrium |
| AKR1A1    | -0.8666867 | 0.0096167 | CD4\_Th\_Myometrium |
| ARPC5L    | -0.8764658 | 0.0160185 | CD4\_Th\_Myometrium |
| POMP      | -0.8769291 | 0.0165832 | CD4\_Th\_Myometrium |
| CISD3     | -0.8908939 | 0.0069729 | CD4\_Th\_Myometrium |
| TRAF1     | -0.9403920 | 0.0492860 | CD4\_Th\_Myometrium |
| IDH3G     | -0.9670296 | 0.0463502 | CD4\_Th\_Myometrium |
| PDCL3     | -0.9758514 | 0.0123055 | CD4\_Th\_Myometrium |
| MMP24OS   | -0.9814572 | 0.0218528 | CD4\_Th\_Myometrium |
| TMEM141   | -1.0814298 | 0.0040806 | CD4\_Th\_Myometrium |
| EML4      | -1.1611254 | 0.0031893 | CD4\_Th\_Myometrium |
| GIPC1     | -1.1664801 | 0.0274768 | CD4\_Th\_Myometrium |
| CKLF      | -1.1689757 | 0.0148354 | CD4\_Th\_Myometrium |
| CREM      | -1.1987681 | 0.0000000 | CD4\_Th\_Myometrium |
| NME1      | -1.2294102 | 0.0050805 | CD4\_Th\_Myometrium |
| HIF1A     | -1.2584171 | 0.0160185 | CD4\_Th\_Myometrium |
| PLEKHF1   | -1.2919198 | 0.0023801 | CD4\_Th\_Myometrium |
| LINC01871 | -1.3332773 | 0.0007750 | CD4\_Th\_Myometrium |
| IL2RB     | -1.4356219 | 0.0492860 | CD4\_Th\_Myometrium |
| DNPH1     | -1.5336549 | 0.0160185 | CD4\_Th\_Myometrium |
| CAPG      | -1.6038985 | 0.0047272 | CD4\_Th\_Myometrium |
| LUZP1     | -1.6348214 | 0.0026612 | CD4\_Th\_Myometrium |
| ANKRD39   | -1.6841191 | 0.0196103 | CD4\_Th\_Myometrium |
| SH2D2A    | -1.7347205 | 0.0009937 | CD4\_Th\_Myometrium |
| ETV6      | -2.0910991 | 0.0000007 | CD4\_Th\_Myometrium |
| TRBV24-1  | -5.9898830 | 0.0218528 | CD4\_Th\_Myometrium |

\[1\] “Results for celltype:” \[1\] “CD4\_TEM\_Myometrium” \[1\] “Genes
upregulated in disease”

| genes   |     logFC |      padj | cell\_type           |
|:--------|----------:|----------:|:---------------------|
| MT-CO2  | 13.249060 | 0.0007097 | CD4\_TEM\_Myometrium |
| MT-ND3  | 12.179430 | 0.0007097 | CD4\_TEM\_Myometrium |
| MT-ND4L | 11.884163 | 0.0016224 | CD4\_TEM\_Myometrium |
| MT-ND1  | 11.604084 | 0.0050870 | CD4\_TEM\_Myometrium |
| TPP1    |  8.089229 | 0.0032257 | CD4\_TEM\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type           |
|:------|----------:|----------:|:---------------------|
| GZMB  | -4.424879 | 0.0060433 | CD4\_TEM\_Myometrium |

\[1\] “Results for celltype:” \[1\] “CD8\_T\_Naive\_CM\_Myometrium”
\[1\] “Genes upregulated in disease”

| genes   |     logFC |      padj | cell\_type                    |
|:--------|----------:|----------:|:------------------------------|
| MT-CO1  | 12.826156 | 0.0174205 | CD8\_T\_Naive\_CM\_Myometrium |
| MT-ND3  | 11.797657 | 0.0012579 | CD8\_T\_Naive\_CM\_Myometrium |
| MT-ND4L | 11.504966 | 0.0174205 | CD8\_T\_Naive\_CM\_Myometrium |
| MT-ND1  | 10.659671 | 0.0174205 | CD8\_T\_Naive\_CM\_Myometrium |
| AHNAK   |  7.689288 | 0.0174205 | CD8\_T\_Naive\_CM\_Myometrium |
| CNOT1   |  7.347821 | 0.0283096 | CD8\_T\_Naive\_CM\_Myometrium |
| ACTN1   |  4.220623 | 0.0174205 | CD8\_T\_Naive\_CM\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type                    |
|:------|----------:|----------:|:------------------------------|
| CD7   | -1.501792 | 0.0174205 | CD8\_T\_Naive\_CM\_Myometrium |

\[1\] “Results for celltype:” \[1\] “CD8\_TEM\_Myometrium” \[1\] “Genes
upregulated in disease”

| genes    |    logFC |      padj | cell\_type           |
|:---------|---------:|----------:|:---------------------|
| TRAV12-2 | 8.723446 | 0.0089005 | CD8\_TEM\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes  |     logFC |      padj | cell\_type           |
|:-------|----------:|----------:|:---------------------|
| CD300A | -2.031445 | 0.0072514 | CD8\_TEM\_Myometrium |

\[1\] “Results for celltype:” \[1\] “NKT\_Myometrium” \[1\] “Genes
upregulated in disease”

| genes   |     logFC |      padj | cell\_type      |
|:--------|----------:|----------:|:----------------|
| MT-CO1  | 15.777173 | 0.0000102 | NKT\_Myometrium |
| MT-CO2  | 15.059226 | 0.0000271 | NKT\_Myometrium |
| MT-ND3  | 13.512589 | 0.0002139 | NKT\_Myometrium |
| MT-ND4L | 13.259050 | 0.0004086 | NKT\_Myometrium |
| MT-ND1  | 12.918068 | 0.0061112 | NKT\_Myometrium |
| TMEM163 |  8.644997 | 0.0241875 | NKT\_Myometrium |
| AMPD3   |  8.475363 | 0.0000146 | NKT\_Myometrium |
| BCL2    |  1.905195 | 0.0339477 | NKT\_Myometrium |
| PRDM2   |  1.626177 | 0.0339477 | NKT\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes   |     logFC |      padj | cell\_type      |
|:--------|----------:|----------:|:----------------|
| H2AFV   | -1.039473 | 0.0061112 | NKT\_Myometrium |
| GSTO1   | -1.268993 | 0.0014281 | NKT\_Myometrium |
| PLAAT4  | -1.359089 | 0.0498714 | NKT\_Myometrium |
| TMEM107 | -1.558402 | 0.0187798 | NKT\_Myometrium |
| NENF    | -2.464513 | 0.0005573 | NKT\_Myometrium |

\[1\] “Results for celltype:” \[1\] “CD16\_NK\_Myometrium” \[1\] “Genes
upregulated in disease”

| genes   |    logFC |      padj | cell\_type           |
|:--------|---------:|----------:|:---------------------|
| MT-CO2  | 12.87848 | 0.0075725 | CD16\_NK\_Myometrium |
| MT-ND4L | 11.19620 | 0.0129097 | CD16\_NK\_Myometrium |
| MT-ND3  | 11.11822 | 0.0090863 | CD16\_NK\_Myometrium |
| MT-ND1  | 10.24793 | 0.0129097 | CD16\_NK\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Results for celltype:” \[1\] “CD14\_Monocyte\_Myometrium” \[1\]
“Genes upregulated in disease”

| genes   |    logFC |      padj | cell\_type                 |
|:--------|---------:|----------:|:---------------------------|
| FAR2    | 1.994471 | 0.0395067 | CD14\_Monocyte\_Myometrium |
| IL13RA1 | 1.548063 | 0.0079428 | CD14\_Monocyte\_Myometrium |
| MYADM   | 1.537120 | 0.0017390 | CD14\_Monocyte\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes    |     logFC |      padj | cell\_type                 |
|:---------|----------:|----------:|:---------------------------|
| SIGLEC10 | -1.679904 | 0.0017390 | CD14\_Monocyte\_Myometrium |
| CASP5    | -2.456328 | 0.0079428 | CD14\_Monocyte\_Myometrium |
| MCTP2    | -2.458933 | 0.0219469 | CD14\_Monocyte\_Myometrium |
| SNED1    | -3.651625 | 0.0219469 | CD14\_Monocyte\_Myometrium |

\[1\] “Results for celltype:” \[1\] “CD16\_Monocyte\_Myometrium” \[1\]
“Genes upregulated in disease”

| genes   |     logFC |      padj | cell\_type                 |
|:--------|----------:|----------:|:---------------------------|
| MT-ND3  | 12.562767 | 0.0000005 | CD16\_Monocyte\_Myometrium |
| MT-ND4L | 11.513410 | 0.0027828 | CD16\_Monocyte\_Myometrium |
| MT-CO2  |  9.074818 | 0.0002111 | CD16\_Monocyte\_Myometrium |
| MT-CO1  |  4.431828 | 0.0151151 | CD16\_Monocyte\_Myometrium |
| MS4A4E  |  2.407244 | 0.0033312 | CD16\_Monocyte\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type                 |
|:------|----------:|----------:|:---------------------------|
| TIMP1 | -2.886463 | 0.0033312 | CD16\_Monocyte\_Myometrium |

\[1\] “Results for celltype:” \[1\] “Nonclassical-monocyte\_Myometrium”
\[1\] “Genes upregulated in disease”

| genes   |     logFC |      padj | cell\_type                        |
|:--------|----------:|----------:|:----------------------------------|
| MT-ND3  | 11.464300 | 0.0000051 | Nonclassical-monocyte\_Myometrium |
| MT-ND4L | 10.352679 | 0.0000250 | Nonclassical-monocyte\_Myometrium |
| RNASE1  |  7.933523 | 0.0138774 | Nonclassical-monocyte\_Myometrium |
| MT-ND1  |  6.035167 | 0.0000250 | Nonclassical-monocyte\_Myometrium |
| ERAP2   |  2.703235 | 0.0050768 | Nonclassical-monocyte\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Results for celltype:” \[1\] “Macrophage\_Myometrium” \[1\]
“Genes upregulated in disease”

| genes      |    logFC |      padj | cell\_type             |
|:-----------|---------:|----------:|:-----------------------|
| MT-ND3     | 5.708658 | 0.0016832 | Macrophage\_Myometrium |
| MACROD2    | 4.322958 | 0.0404171 | Macrophage\_Myometrium |
| PNPLA7     | 3.452565 | 0.0463330 | Macrophage\_Myometrium |
| ERAP2      | 2.616627 | 0.0252194 | Macrophage\_Myometrium |
| AC025171.2 | 2.396697 | 0.0430607 | Macrophage\_Myometrium |
| MT-ND2     | 2.152351 | 0.0404171 | Macrophage\_Myometrium |
| PON2       | 1.791052 | 0.0430607 | Macrophage\_Myometrium |
| DPEP2      | 1.750447 | 0.0430607 | Macrophage\_Myometrium |
| HACD4      | 1.556192 | 0.0463330 | Macrophage\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes    |     logFC |      padj | cell\_type             |
|:---------|----------:|----------:|:-----------------------|
| CWF19L1  | -1.211997 | 0.0463330 | Macrophage\_Myometrium |
| PTPN7    | -1.561623 | 0.0478157 | Macrophage\_Myometrium |
| U62317.4 | -1.660282 | 0.0228084 | Macrophage\_Myometrium |
| DYNC2H1  | -1.660407 | 0.0463330 | Macrophage\_Myometrium |
| PROCR    | -2.502808 | 0.0478157 | Macrophage\_Myometrium |
| SLC39A8  | -2.513208 | 0.0255465 | Macrophage\_Myometrium |
| C1orf21  | -3.372247 | 0.0016832 | Macrophage\_Myometrium |
| DPP4     | -3.577007 | 0.0430607 | Macrophage\_Myometrium |
| PDPN     | -4.005679 | 0.0017010 | Macrophage\_Myometrium |
| CCL23    | -4.553187 | 0.0016832 | Macrophage\_Myometrium |
| RAI14    | -5.197667 | 0.0022203 | Macrophage\_Myometrium |

\[1\] “Results for celltype:” \[1\] “LED\_Myometrium” \[1\] “Genes
upregulated in disease”

| genes  |    logFC |      padj | cell\_type      |
|:-------|---------:|----------:|:----------------|
| ERAP2  | 4.583615 | 0.0000072 | LED\_Myometrium |
| LTC4S  | 4.229713 | 0.0047860 | LED\_Myometrium |
| NTN1   | 3.495635 | 0.0047860 | LED\_Myometrium |
| KLF4   | 2.843027 | 0.0035655 | LED\_Myometrium |
| GGT5   | 2.772426 | 0.0266574 | LED\_Myometrium |
| ADIRF  | 2.141052 | 0.0487177 | LED\_Myometrium |
| MITF   | 2.107503 | 0.0487177 | LED\_Myometrium |
| OAS3   | 1.996895 | 0.0310756 | LED\_Myometrium |
| EPHX1  | 1.761417 | 0.0442784 | LED\_Myometrium |
| FAM49A | 1.074763 | 0.0468790 | LED\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes      |     logFC |      padj | cell\_type      |
|:-----------|----------:|----------:|:----------------|
| GMDS       | -1.196876 | 0.0364533 | LED\_Myometrium |
| TMEM237    | -1.316927 | 0.0468790 | LED\_Myometrium |
| TSPAN15    | -1.472304 | 0.0185096 | LED\_Myometrium |
| CYBC1      | -1.509273 | 0.0335535 | LED\_Myometrium |
| CSF1       | -1.600185 | 0.0047860 | LED\_Myometrium |
| SOX4       | -1.662448 | 0.0335535 | LED\_Myometrium |
| LDLRAD3    | -1.719323 | 0.0134897 | LED\_Myometrium |
| KRBA2      | -1.783611 | 0.0052341 | LED\_Myometrium |
| DTNB       | -1.830204 | 0.0468790 | LED\_Myometrium |
| RDH10      | -2.165191 | 0.0389931 | LED\_Myometrium |
| ABHD17C    | -2.198118 | 0.0168882 | LED\_Myometrium |
| NEO1       | -3.325944 | 0.0487177 | LED\_Myometrium |
| TIMP1      | -3.749763 | 0.0035655 | LED\_Myometrium |
| AL691447.2 | -3.782056 | 0.0487177 | LED\_Myometrium |
| CCDC85A    | -3.859531 | 0.0457120 | LED\_Myometrium |
| CEMIP      | -5.214791 | 0.0035655 | LED\_Myometrium |
| CYP4Z1     | -6.022642 | 0.0487177 | LED\_Myometrium |
| PADI2      | -6.597656 | 0.0335535 | LED\_Myometrium |

\[1\] “Results for celltype:” \[1\] “Endothelial\_Myometrium” \[1\]
“Genes upregulated in disease”

| genes  |    logFC |     padj | cell\_type              |
|:-------|---------:|---------:|:------------------------|
| MT-ND3 | 3.510357 | 2.62e-05 | Endothelial\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes |      logFC |      padj | cell\_type              |
|:------|-----------:|----------:|:------------------------|
| RBP5  |  -3.197249 | 0.0062853 | Endothelial\_Myometrium |
| SMOC1 |  -6.912114 | 0.0449854 | Endothelial\_Myometrium |
| CGA   |  -9.340251 | 0.0164900 | Endothelial\_Myometrium |
| NPY   | -10.912157 | 0.0013329 | Endothelial\_Myometrium |

# Volcano plots in Late disease

\[1\] “Results for celltype:” \[1\] “CD4\_T\_Naive\_CM\_Myometrium”
\[1\] “Genes upregulated in disease”

| genes |    logFC |      padj | cell\_type                    |
|:------|---------:|----------:|:------------------------------|
| MYB   | 2.037973 | 0.0052471 | CD4\_T\_Naive\_CM\_Myometrium |
| SOCS2 | 1.392152 | 0.0006518 | CD4\_T\_Naive\_CM\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Results for celltype:” \[1\] “CD16\_Monocyte\_Myometrium” \[1\]
“Genes upregulated in disease”

| genes |    logFC |      padj | cell\_type                 |
|:------|---------:|----------:|:---------------------------|
| TCN2  | 1.490004 | 0.0387576 | CD16\_Monocyte\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes  |     logFC |      padj | cell\_type                 |
|:-------|----------:|----------:|:---------------------------|
| SLC2A3 | -2.295233 | 0.0181648 | CD16\_Monocyte\_Myometrium |
| OLR1   | -2.737356 | 0.0171015 | CD16\_Monocyte\_Myometrium |
| HAMP   | -4.590769 | 0.0387576 | CD16\_Monocyte\_Myometrium |

\[1\] “Results for celltype:” \[1\] “Endothelial\_Myometrium” \[1\]
“Genes upregulated in disease”

| genes   |    logFC |      padj | cell\_type              |
|:--------|---------:|----------:|:------------------------|
| PROM1   | 6.734898 | 0.0048122 | Endothelial\_Myometrium |
| ACSM3   | 4.384707 | 0.0002655 | Endothelial\_Myometrium |
| TRIM54  | 2.911032 | 0.0052782 | Endothelial\_Myometrium |
| ITGB3   | 2.878558 | 0.0406823 | Endothelial\_Myometrium |
| DNAJC5G | 1.977710 | 0.0021731 | Endothelial\_Myometrium |
| ZNF704  | 1.699443 | 0.0406823 | Endothelial\_Myometrium |
| PDE2A   | 1.379478 | 0.0300078 | Endothelial\_Myometrium |
| MYH10   | 1.316748 | 0.0307494 | Endothelial\_Myometrium |
| LRP5    | 1.156695 | 0.0352195 | Endothelial\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes   |     logFC |      padj | cell\_type              |
|:--------|----------:|----------:|:------------------------|
| STARD4  | -1.394332 | 0.0406823 | Endothelial\_Myometrium |
| UQCC2   | -1.609536 | 0.0002655 | Endothelial\_Myometrium |
| PLD1    | -1.816980 | 0.0300078 | Endothelial\_Myometrium |
| C2CD4B  | -2.264148 | 0.0300078 | Endothelial\_Myometrium |
| EGR1    | -2.541244 | 0.0473387 | Endothelial\_Myometrium |
| PLAUR   | -3.320991 | 0.0233383 | Endothelial\_Myometrium |
| GRAMD1B | -4.344889 | 0.0119218 | Endothelial\_Myometrium |

# Volcano plots in averagePE

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes  |    logFC |      padj | cell\_type           |
|:-------|---------:|----------:|:---------------------|
| MT-CO1 | 5.940716 | 0.0038178 | Naive\_B\_Myometrium |
| MT-ND3 | 5.366504 | 0.0038178 | Naive\_B\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes   |    logFC |      padj | cell\_type            |
|:--------|---------:|----------:|:----------------------|
| MT-CO1  | 6.190660 | 0.0006508 | Memory\_B\_Myometrium |
| MT-ND3  | 5.465695 | 0.0008386 | Memory\_B\_Myometrium |
| MT-ND4L | 5.285548 | 0.0033099 | Memory\_B\_Myometrium |
| MT-ND1  | 4.760249 | 0.0387563 | Memory\_B\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes   |    logFC |      padj | cell\_type                    |
|:--------|---------:|----------:|:------------------------------|
| MT-CO1  | 6.155837 | 0.0028924 | CD4\_T\_Naive\_CM\_Myometrium |
| MT-ND3  | 5.753414 | 0.0028924 | CD4\_T\_Naive\_CM\_Myometrium |
| MT-ND4L | 5.477872 | 0.0273328 | CD4\_T\_Naive\_CM\_Myometrium |
| MT-ND1  | 5.240773 | 0.0459380 | CD4\_T\_Naive\_CM\_Myometrium |
| ACTN1   | 4.046184 | 0.0028924 | CD4\_T\_Naive\_CM\_Myometrium |
| BBC3    | 3.953858 | 0.0415545 | CD4\_T\_Naive\_CM\_Myometrium |
| SOCS2   | 1.431917 | 0.0162327 | CD4\_T\_Naive\_CM\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes      |     logFC |      padj | cell\_type          |
|:-----------|----------:|----------:|:--------------------|
| MT-ND3     | 5.9185537 | 0.0001871 | CD4\_Th\_Myometrium |
| MT-ND4L    | 5.7435558 | 0.0017446 | CD4\_Th\_Myometrium |
| MT-ND1     | 5.3298967 | 0.0032640 | CD4\_Th\_Myometrium |
| MT-CO1     | 4.6525739 | 0.0004537 | CD4\_Th\_Myometrium |
| ACTN1      | 1.5210095 | 0.0308657 | CD4\_Th\_Myometrium |
| KLF2       | 1.3152991 | 0.0007353 | CD4\_Th\_Myometrium |
| TXK        | 1.2215254 | 0.0048447 | CD4\_Th\_Myometrium |
| AL035530.2 | 1.1436629 | 0.0459565 | CD4\_Th\_Myometrium |
| ERAP2      | 1.1201540 | 0.0048447 | CD4\_Th\_Myometrium |
| SORBS3     | 0.9721763 | 0.0001359 | CD4\_Th\_Myometrium |
| DUSP1      | 0.9390475 | 0.0497152 | CD4\_Th\_Myometrium |
| MYADM      | 0.9103872 | 0.0346044 | CD4\_Th\_Myometrium |
| CDC25B     | 0.8713727 | 0.0476196 | CD4\_Th\_Myometrium |
| JADE1      | 0.8381247 | 0.0048447 | CD4\_Th\_Myometrium |
| CA5B       | 0.8050985 | 0.0469529 | CD4\_Th\_Myometrium |
| EGLN1      | 0.8007432 | 0.0214027 | CD4\_Th\_Myometrium |
| MYC        | 0.6371841 | 0.0476196 | CD4\_Th\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes   |      logFC |      padj | cell\_type          |
|:--------|-----------:|----------:|:--------------------|
| AKR1A1  | -0.6112545 | 0.0441915 | CD4\_Th\_Myometrium |
| EML4    | -0.6452360 | 0.0469529 | CD4\_Th\_Myometrium |
| CREM    | -0.7189956 | 0.0003562 | CD4\_Th\_Myometrium |
| RHOC    | -0.8501934 | 0.0476196 | CD4\_Th\_Myometrium |
| POLR2J3 | -0.8502125 | 0.0476196 | CD4\_Th\_Myometrium |
| ETV6    | -1.2481392 | 0.0000116 | CD4\_Th\_Myometrium |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes   |    logFC |      padj | cell\_type           |
|:--------|---------:|----------:|:---------------------|
| MT-CO2  | 6.334635 | 0.0062819 | CD4\_TEM\_Myometrium |
| MT-ND3  | 5.854225 | 0.0062819 | CD4\_TEM\_Myometrium |
| MT-ND4L | 5.520691 | 0.0166791 | CD4\_TEM\_Myometrium |
| MT-ND1  | 5.394525 | 0.0439453 | CD4\_TEM\_Myometrium |
| TPP1    | 4.044220 | 0.0082101 | CD4\_TEM\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes    |     logFC |      padj | cell\_type           |
|:---------|----------:|----------:|:---------------------|
| GZMH     | -1.497816 | 0.0333816 | CD4\_TEM\_Myometrium |
| GZMB     | -2.795358 | 0.0082101 | CD4\_TEM\_Myometrium |
| HLA-DRB5 | -4.008987 | 0.0082101 | CD4\_TEM\_Myometrium |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes  |    logFC |      padj | cell\_type                    |
|:-------|---------:|----------:|:------------------------------|
| MT-ND3 | 5.583423 | 0.0166641 | CD8\_T\_Naive\_CM\_Myometrium |
| PDE4B  | 1.547370 | 0.0410794 | CD8\_T\_Naive\_CM\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes   |    logFC |      padj | cell\_type      |
|:--------|---------:|----------:|:----------------|
| MT-CO1  | 7.495685 | 0.0001480 | NKT\_Myometrium |
| MT-CO2  | 7.139903 | 0.0005187 | NKT\_Myometrium |
| MT-ND3  | 6.372875 | 0.0038586 | NKT\_Myometrium |
| MT-ND4L | 6.225793 | 0.0061756 | NKT\_Myometrium |
| AMPD3   | 3.865765 | 0.0308267 | NKT\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type      |
|:------|----------:|----------:|:----------------|
| GSTO1 | -0.755625 | 0.0102706 | NKT\_Myometrium |
| NENF  | -1.404692 | 0.0084062 | NKT\_Myometrium |
| TIGIT | -1.938208 | 0.0061756 | NKT\_Myometrium |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes  |    logFC |      padj | cell\_type           |
|:-------|---------:|----------:|:---------------------|
| MT-CO2 | 5.983853 | 0.0137511 | CD16\_NK\_Myometrium |
| MT-ND3 | 5.158626 | 0.0267501 | CD16\_NK\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes |      logFC |      padj | cell\_type           |
|:------|-----------:|----------:|:---------------------|
| MXD4  | -0.9616806 | 0.0156276 | CD16\_NK\_Myometrium |
| KLRC2 | -1.8764447 | 0.0137511 | CD16\_NK\_Myometrium |
| PTMS  | -2.5746688 | 0.0039671 | CD16\_NK\_Myometrium |
| LAG3  | -3.2438980 | 0.0137511 | CD16\_NK\_Myometrium |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes   |     logFC |      padj | cell\_type                 |
|:--------|----------:|----------:|:---------------------------|
| MFGE8   | 5.7656815 | 0.0357874 | CD14\_Monocyte\_Myometrium |
| TPST1   | 3.0345251 | 0.0437904 | CD14\_Monocyte\_Myometrium |
| VSIG4   | 2.8384551 | 0.0093440 | CD14\_Monocyte\_Myometrium |
| SLC40A1 | 1.8169681 | 0.0464402 | CD14\_Monocyte\_Myometrium |
| LHFPL2  | 1.3662576 | 0.0285114 | CD14\_Monocyte\_Myometrium |
| FAR2    | 1.1736852 | 0.0357874 | CD14\_Monocyte\_Myometrium |
| IL13RA1 | 0.9289398 | 0.0182842 | CD14\_Monocyte\_Myometrium |
| AHR     | 0.8804419 | 0.0285114 | CD14\_Monocyte\_Myometrium |
| MYADM   | 0.7383996 | 0.0357874 | CD14\_Monocyte\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes   |     logFC |      padj | cell\_type                 |
|:--------|----------:|----------:|:---------------------------|
| TNFRSF8 | -1.012851 | 0.0357874 | CD14\_Monocyte\_Myometrium |
| POLR2J3 | -1.201464 | 0.0357874 | CD14\_Monocyte\_Myometrium |
| CCR5AS  | -1.589785 | 0.0357874 | CD14\_Monocyte\_Myometrium |
| SNED1   | -2.549868 | 0.0182842 | CD14\_Monocyte\_Myometrium |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes   |    logFC |      padj | cell\_type                 |
|:--------|---------:|----------:|:---------------------------|
| MT-ND3  | 5.929815 | 0.0000087 | CD16\_Monocyte\_Myometrium |
| MT-ND4L | 5.247340 | 0.0096676 | CD16\_Monocyte\_Myometrium |
| MT-CO2  | 4.256478 | 0.0010789 | CD16\_Monocyte\_Myometrium |
| MS4A4E  | 1.883347 | 0.0000507 | CD16\_Monocyte\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes    |     logFC |      padj | cell\_type                 |
|:---------|----------:|----------:|:---------------------------|
| FABP5    | -1.508582 | 0.0333073 | CD16\_Monocyte\_Myometrium |
| MARCO    | -1.573589 | 0.0096676 | CD16\_Monocyte\_Myometrium |
| C15orf48 | -2.777806 | 0.0096676 | CD16\_Monocyte\_Myometrium |
| G0S2     | -3.695065 | 0.0223409 | CD16\_Monocyte\_Myometrium |
| HAMP     | -4.023321 | 0.0000478 | CD16\_Monocyte\_Myometrium |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes   |    logFC |      padj | cell\_type                        |
|:--------|---------:|----------:|:----------------------------------|
| MT-ND3  | 5.461946 | 0.0000815 | Nonclassical-monocyte\_Myometrium |
| MT-ND4L | 4.919404 | 0.0004758 | Nonclassical-monocyte\_Myometrium |
| MT-ND1  | 2.541841 | 0.0058931 | Nonclassical-monocyte\_Myometrium |
| ERAP2   | 1.615217 | 0.0070244 | Nonclassical-monocyte\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type                        |
|:------|----------:|----------:|:----------------------------------|
| LYPD2 | -3.054311 | 0.0016778 | Nonclassical-monocyte\_Myometrium |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes      |    logFC |      padj | cell\_type             |
|:-----------|---------:|----------:|:-----------------------|
| MT-ND3     | 2.563834 | 0.0110620 | Macrophage\_Myometrium |
| AC025171.2 | 1.426300 | 0.0094941 | Macrophage\_Myometrium |
| PON2       | 1.243275 | 0.0046393 | Macrophage\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes     |      logFC |      padj | cell\_type             |
|:----------|-----------:|----------:|:-----------------------|
| ICAM1     | -0.8948909 | 0.0046393 | Macrophage\_Myometrium |
| IL2RG     | -0.9941548 | 0.0127156 | Macrophage\_Myometrium |
| MT1F      | -1.1927346 | 0.0257411 | Macrophage\_Myometrium |
| PROCR     | -1.8861716 | 0.0094941 | Macrophage\_Myometrium |
| C1orf21   | -1.9168905 | 0.0083901 | Macrophage\_Myometrium |
| PDPN      | -2.0976816 | 0.0257411 | Macrophage\_Myometrium |
| SPARC     | -2.3515937 | 0.0046393 | Macrophage\_Myometrium |
| DPP4      | -2.6965385 | 0.0056484 | Macrophage\_Myometrium |
| MMP9      | -2.7898352 | 0.0094941 | Macrophage\_Myometrium |
| CXCL1     | -3.3647209 | 0.0094941 | Macrophage\_Myometrium |
| MMP12     | -4.6225825 | 0.0094941 | Macrophage\_Myometrium |
| LINC01933 | -5.7763271 | 0.0430935 | Macrophage\_Myometrium |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes |    logFC |      padj | cell\_type      |
|:------|---------:|----------:|:----------------|
| ERAP2 | 2.310666 | 0.0000825 | LED\_Myometrium |
| LTC4S | 2.262976 | 0.0163822 | LED\_Myometrium |
| GGT5  | 1.944141 | 0.0017825 | LED\_Myometrium |
| NTN1  | 1.749104 | 0.0376034 | LED\_Myometrium |
| KLF4  | 1.471912 | 0.0164104 | LED\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes |      logFC |      padj | cell\_type      |
|:------|-----------:|----------:|:----------------|
| KRBA2 | -0.9858624 | 0.0483724 | LED\_Myometrium |
| G0S2  | -5.7278670 | 0.0164104 | LED\_Myometrium |
| PADI2 | -6.2576862 | 0.0000626 | LED\_Myometrium |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes  |    logFC |      padj | cell\_type              |
|:-------|---------:|----------:|:------------------------|
| MT-ND3 | 1.663679 | 0.0009741 | Endothelial\_Myometrium |

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type              |
|:------|----------:|----------:|:------------------------|
| RBP5  | -1.912537 | 0.0121051 | Endothelial\_Myometrium |
| SMOC1 | -4.922118 | 0.0018082 | Endothelial\_Myometrium |
| CGA   | -7.758316 | 0.0000342 | Endothelial\_Myometrium |
| NPY   | -7.795989 | 0.0000342 | Endothelial\_Myometrium |

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
