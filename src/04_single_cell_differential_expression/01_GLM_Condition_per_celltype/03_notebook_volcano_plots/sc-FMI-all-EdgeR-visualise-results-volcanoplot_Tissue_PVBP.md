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

> Tissue: Placenta

``` r
# This notebook is for tissue and annotation level:
tissue= "Placenta" # "CAM", "PBMC", "Myometrium"
annotation_level <- "CellTypeManual.l3"

# List of cell types
cellgroups$CellTypeManual.l3 <- paste0(cellgroups$CellTypeManual.l3,"_",tissue)
```

# Volcano plots in Early disease

\[1\] “Results for celltype:” \[1\] “CD4\_T\_Naive\_CM\_Placenta” \[1\]
“Genes upregulated in disease”

| genes |    logFC |      padj | cell\_type                  |
|:------|---------:|----------:|:----------------------------|
| MPP7  | 1.436753 | 0.0175668 | CD4\_T\_Naive\_CM\_Placenta |

\[1\] “Genes downregulated in disease”

| genes    |     logFC |      padj | cell\_type                  |
|:---------|----------:|----------:|:----------------------------|
| TMEM80   | -1.706569 | 0.0175668 | CD4\_T\_Naive\_CM\_Placenta |
| TRBV21-1 | -3.381917 | 0.0175668 | CD4\_T\_Naive\_CM\_Placenta |

\[1\] “Results for celltype:” \[1\] “CD4\_TEM\_Placenta” \[1\] “Genes
upregulated in disease”

| genes      |    logFC |      padj | cell\_type         |
|:-----------|---------:|----------:|:-------------------|
| ERAP2      | 2.211053 | 0.0005653 | CD4\_TEM\_Placenta |
| AP001011.1 | 2.150855 | 0.0158773 | CD4\_TEM\_Placenta |
| OOEP       | 1.948529 | 0.0321968 | CD4\_TEM\_Placenta |
| DPEP2      | 1.938286 | 0.0414716 | CD4\_TEM\_Placenta |
| NELL2      | 1.740642 | 0.0011453 | CD4\_TEM\_Placenta |
| RIN3       | 1.700250 | 0.0288204 | CD4\_TEM\_Placenta |
| MT-CYB     | 1.595065 | 0.0414716 | CD4\_TEM\_Placenta |
| GAS5       | 1.588831 | 0.0098698 | CD4\_TEM\_Placenta |
| PIM2       | 1.520854 | 0.0018060 | CD4\_TEM\_Placenta |
| AC119396.1 | 1.421371 | 0.0196964 | CD4\_TEM\_Placenta |
| MXI1       | 1.419388 | 0.0425070 | CD4\_TEM\_Placenta |
| ZNF277     | 1.386618 | 0.0243057 | CD4\_TEM\_Placenta |
| ANKRD28    | 1.224478 | 0.0287724 | CD4\_TEM\_Placenta |
| SNHG6      | 1.175111 | 0.0121516 | CD4\_TEM\_Placenta |
| ITPKB      | 1.152903 | 0.0408931 | CD4\_TEM\_Placenta |

\[1\] “Genes downregulated in disease”

| genes       |      logFC |      padj | cell\_type         |
|:------------|-----------:|----------:|:-------------------|
| TXN         | -0.8892231 | 0.0439894 | CD4\_TEM\_Placenta |
| LDHA        | -1.1134448 | 0.0005152 | CD4\_TEM\_Placenta |
| CTSC        | -1.1419692 | 0.0005653 | CD4\_TEM\_Placenta |
| CKLF        | -1.1429099 | 0.0055809 | CD4\_TEM\_Placenta |
| LGALS1      | -1.5241032 | 0.0499841 | CD4\_TEM\_Placenta |
| MIR4435-2HG | -1.6618309 | 0.0045892 | CD4\_TEM\_Placenta |
| JAZF1       | -1.8541500 | 0.0121516 | CD4\_TEM\_Placenta |
| CAPG        | -2.0090399 | 0.0040172 | CD4\_TEM\_Placenta |
| XCL1        | -2.8610017 | 0.0408931 | CD4\_TEM\_Placenta |

\[1\] “Results for celltype:” \[1\] “CD8\_T\_Naive\_CM\_Placenta” \[1\]
“Genes upregulated in disease”

| genes    |    logFC |      padj | cell\_type                  |
|:---------|---------:|----------:|:----------------------------|
| TRAV12-1 | 6.680208 | 0.0402800 | CD8\_T\_Naive\_CM\_Placenta |
| TRBV7-2  | 5.117579 | 0.0402800 | CD8\_T\_Naive\_CM\_Placenta |
| TRBV2    | 3.781948 | 0.0264769 | CD8\_T\_Naive\_CM\_Placenta |
| ERAP2    | 1.812435 | 0.0264769 | CD8\_T\_Naive\_CM\_Placenta |

\[1\] “Genes downregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Results for celltype:” \[1\] “CD8\_TEM\_Placenta” \[1\] “Genes
upregulated in disease”

| genes      |    logFC |      padj | cell\_type         |
|:-----------|---------:|----------:|:-------------------|
| ERAP2      | 1.687349 | 0.0000109 | CD8\_TEM\_Placenta |
| AC022075.1 | 1.582868 | 0.0369246 | CD8\_TEM\_Placenta |

\[1\] “Genes downregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Results for celltype:” \[1\] “MAIT\_Placenta” \[1\] “Genes
upregulated in disease”

| genes    |    logFC |      padj | cell\_type     |
|:---------|---------:|----------:|:---------------|
| TRAV9-2  | 7.483492 | 0.0047142 | MAIT\_Placenta |
| TRBV6-5  | 3.743819 | 0.0481453 | MAIT\_Placenta |
| NR4A3    | 3.287802 | 0.0488013 | MAIT\_Placenta |
| RGS2     | 3.030092 | 0.0047142 | MAIT\_Placenta |
| ERAP2    | 2.628686 | 0.0000021 | MAIT\_Placenta |
| A2M      | 2.258525 | 0.0005882 | MAIT\_Placenta |
| PSAP     | 2.215604 | 0.0000000 | MAIT\_Placenta |
| RASGEF1B | 1.490579 | 0.0394317 | MAIT\_Placenta |
| MYADM    | 1.422058 | 0.0488013 | MAIT\_Placenta |
| PZP      | 1.398542 | 0.0268835 | MAIT\_Placenta |
| RGCC     | 1.239453 | 0.0268835 | MAIT\_Placenta |

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type     |
|:------|----------:|----------:|:---------------|
| KLRC4 | -4.056598 | 0.0447301 | MAIT\_Placenta |

\[1\] “Results for celltype:” \[1\] “CD56\_NK\_Placenta” \[1\] “Genes
upregulated in disease”

| genes      |    logFC |      padj | cell\_type         |
|:-----------|---------:|----------:|:-------------------|
| KRT72      | 6.458936 | 0.0086960 | CD56\_NK\_Placenta |
| USP9Y      | 6.074473 | 0.0118389 | CD56\_NK\_Placenta |
| ABHD15-AS1 | 3.859848 | 0.0000753 | CD56\_NK\_Placenta |
| CD3G       | 3.544173 | 0.0214521 | CD56\_NK\_Placenta |
| CD3D       | 2.829622 | 0.0219814 | CD56\_NK\_Placenta |
| TSPOAP1    | 2.735012 | 0.0217423 | CD56\_NK\_Placenta |
| IGF1R      | 2.503413 | 0.0069115 | CD56\_NK\_Placenta |
| MTSS1      | 2.469549 | 0.0069115 | CD56\_NK\_Placenta |
| SLCO4C1    | 2.349210 | 0.0449859 | CD56\_NK\_Placenta |
| DPEP2      | 2.282671 | 0.0306115 | CD56\_NK\_Placenta |
| RAP1GAP2   | 2.223510 | 0.0342386 | CD56\_NK\_Placenta |
| AC068587.4 | 2.068809 | 0.0331988 | CD56\_NK\_Placenta |
| ERAP2      | 2.013204 | 0.0069115 | CD56\_NK\_Placenta |
| ZNF831     | 1.925125 | 0.0086960 | CD56\_NK\_Placenta |
| RHBDF2     | 1.821519 | 0.0342386 | CD56\_NK\_Placenta |
| CRAMP1     | 1.816081 | 0.0495461 | CD56\_NK\_Placenta |
| SPON2      | 1.810841 | 0.0273801 | CD56\_NK\_Placenta |
| MYLIP      | 1.778306 | 0.0139004 | CD56\_NK\_Placenta |
| ARHGAP27   | 1.689661 | 0.0084387 | CD56\_NK\_Placenta |
| SYTL2      | 1.654370 | 0.0069115 | CD56\_NK\_Placenta |
| MGAT4A     | 1.581628 | 0.0292912 | CD56\_NK\_Placenta |
| MT-CYB     | 1.550612 | 0.0495461 | CD56\_NK\_Placenta |
| NLRP1      | 1.539276 | 0.0069115 | CD56\_NK\_Placenta |
| IQSEC1     | 1.501434 | 0.0449859 | CD56\_NK\_Placenta |
| C12orf75   | 1.433993 | 0.0296840 | CD56\_NK\_Placenta |
| SAMD4B     | 1.393628 | 0.0283484 | CD56\_NK\_Placenta |
| MT-CO3     | 1.366823 | 0.0283484 | CD56\_NK\_Placenta |

\[1\] “Genes downregulated in disease”

| genes   |     logFC |      padj | cell\_type         |
|:--------|----------:|----------:|:-------------------|
| EZR     | -1.264871 | 0.0273801 | CD56\_NK\_Placenta |
| S100A6  | -1.360180 | 0.0495607 | CD56\_NK\_Placenta |
| S100A10 | -1.409511 | 0.0086960 | CD56\_NK\_Placenta |
| S100A11 | -1.534349 | 0.0449859 | CD56\_NK\_Placenta |
| IGFLR1  | -1.543756 | 0.0099721 | CD56\_NK\_Placenta |
| NAA20   | -1.728926 | 0.0118389 | CD56\_NK\_Placenta |
| LGALS1  | -2.720529 | 0.0000007 | CD56\_NK\_Placenta |
| CD82    | -4.411512 | 0.0086960 | CD56\_NK\_Placenta |
| CSH2    | -7.149482 | 0.0086960 | CD56\_NK\_Placenta |

\[1\] “Results for celltype:” \[1\] “CD14\_Monocyte\_Placenta” \[1\]
“Genes upregulated in disease”

| genes  |     logFC |      padj | cell\_type               |
|:-------|----------:|----------:|:-------------------------|
| OCLN   | 4.2909720 | 0.0328817 | CD14\_Monocyte\_Placenta |
| ERAP2  | 1.8338619 | 0.0208118 | CD14\_Monocyte\_Placenta |
| HDAC9  | 1.2913047 | 0.0208118 | CD14\_Monocyte\_Placenta |
| SNTB1  | 1.0733209 | 0.0035727 | CD14\_Monocyte\_Placenta |
| GPD2   | 0.9320776 | 0.0224292 | CD14\_Monocyte\_Placenta |
| FOXN3  | 0.9242860 | 0.0360789 | CD14\_Monocyte\_Placenta |
| INPP5A | 0.8914372 | 0.0457403 | CD14\_Monocyte\_Placenta |
| NUP107 | 0.8790083 | 0.0457403 | CD14\_Monocyte\_Placenta |
| RICTOR | 0.7476156 | 0.0457403 | CD14\_Monocyte\_Placenta |

\[1\] “Genes downregulated in disease”

| genes      |      logFC |      padj | cell\_type               |
|:-----------|-----------:|----------:|:-------------------------|
| AC005280.2 | -0.8667925 | 0.0224292 | CD14\_Monocyte\_Placenta |
| PLA2G7     | -1.4541624 | 0.0457403 | CD14\_Monocyte\_Placenta |
| TNFRSF10C  | -1.5204251 | 0.0457403 | CD14\_Monocyte\_Placenta |
| SLC22A18AS | -1.6063603 | 0.0457403 | CD14\_Monocyte\_Placenta |

\[1\] “Results for celltype:” \[1\] “Fetal-Monocyte\_Placenta” \[1\]
“Genes upregulated in disease”

| genes  |    logFC |      padj | cell\_type               |
|:-------|---------:|----------:|:-------------------------|
| HBG2   | 4.863492 | 0.0031023 | Fetal-Monocyte\_Placenta |
| IFITM3 | 2.836201 | 0.0248347 | Fetal-Monocyte\_Placenta |
| SASH1  | 2.149919 | 0.0041138 | Fetal-Monocyte\_Placenta |

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type               |
|:------|----------:|----------:|:-------------------------|
| G0S2  | -2.045024 | 0.0248347 | Fetal-Monocyte\_Placenta |
| HP    | -2.311218 | 0.0138161 | Fetal-Monocyte\_Placenta |

\[1\] “Results for celltype:” \[1\] “CD16\_Monocyte\_Placenta” \[1\]
“Genes upregulated in disease”

| genes      |     logFC |      padj | cell\_type               |
|:-----------|----------:|----------:|:-------------------------|
| KIAA1958   | 2.6333833 | 0.0422489 | CD16\_Monocyte\_Placenta |
| RNF213-AS1 | 2.5795200 | 0.0474113 | CD16\_Monocyte\_Placenta |
| IFIT1      | 2.3952966 | 0.0474113 | CD16\_Monocyte\_Placenta |
| ERAP2      | 2.3584928 | 0.0000008 | CD16\_Monocyte\_Placenta |
| ZBP1       | 2.3433129 | 0.0422489 | CD16\_Monocyte\_Placenta |
| RSAD2      | 2.2248267 | 0.0381885 | CD16\_Monocyte\_Placenta |
| IFIT3      | 2.1068346 | 0.0230253 | CD16\_Monocyte\_Placenta |
| IFI44L     | 1.9551965 | 0.0228814 | CD16\_Monocyte\_Placenta |
| MX1        | 1.8724871 | 0.0136130 | CD16\_Monocyte\_Placenta |
| SAMD4A     | 1.8507138 | 0.0198376 | CD16\_Monocyte\_Placenta |
| DDX58      | 1.8307051 | 0.0198376 | CD16\_Monocyte\_Placenta |
| MX2        | 1.8133713 | 0.0308700 | CD16\_Monocyte\_Placenta |
| CMPK2      | 1.7382278 | 0.0153483 | CD16\_Monocyte\_Placenta |
| IFI44      | 1.6862796 | 0.0323721 | CD16\_Monocyte\_Placenta |
| OAS3       | 1.6519239 | 0.0474113 | CD16\_Monocyte\_Placenta |
| IFI6       | 1.5804633 | 0.0119025 | CD16\_Monocyte\_Placenta |
| DDX60      | 1.5791035 | 0.0136130 | CD16\_Monocyte\_Placenta |
| SAMD9L     | 1.5301945 | 0.0128695 | CD16\_Monocyte\_Placenta |
| IFIH1      | 1.4688285 | 0.0474113 | CD16\_Monocyte\_Placenta |
| DDX60L     | 1.4530458 | 0.0198376 | CD16\_Monocyte\_Placenta |
| EIF2AK2    | 1.4494568 | 0.0153483 | CD16\_Monocyte\_Placenta |
| TDRD7      | 1.4450317 | 0.0422489 | CD16\_Monocyte\_Placenta |
| RUFY4      | 1.4272005 | 0.0474113 | CD16\_Monocyte\_Placenta |
| SNTB1      | 1.3756370 | 0.0334419 | CD16\_Monocyte\_Placenta |
| PTK2B      | 1.3682113 | 0.0136130 | CD16\_Monocyte\_Placenta |
| C2         | 1.3521797 | 0.0334419 | CD16\_Monocyte\_Placenta |
| STAT2      | 1.3364781 | 0.0153483 | CD16\_Monocyte\_Placenta |
| SP140      | 1.3343702 | 0.0198376 | CD16\_Monocyte\_Placenta |
| PARP14     | 1.3304783 | 0.0474113 | CD16\_Monocyte\_Placenta |
| HELZ2      | 1.3081710 | 0.0491376 | CD16\_Monocyte\_Placenta |
| SRGAP2B    | 1.2599949 | 0.0474113 | CD16\_Monocyte\_Placenta |
| FMNL2      | 1.2592182 | 0.0440024 | CD16\_Monocyte\_Placenta |
| NLRC5      | 1.2563782 | 0.0307320 | CD16\_Monocyte\_Placenta |
| GPD2       | 1.2458584 | 0.0483175 | CD16\_Monocyte\_Placenta |
| XAF1       | 1.2443385 | 0.0378948 | CD16\_Monocyte\_Placenta |
| BTN3A3     | 1.1977522 | 0.0422489 | CD16\_Monocyte\_Placenta |
| ARHGEF3    | 1.1713545 | 0.0474113 | CD16\_Monocyte\_Placenta |
| UTRN       | 1.1629012 | 0.0483175 | CD16\_Monocyte\_Placenta |
| APOL6      | 1.1300997 | 0.0474113 | CD16\_Monocyte\_Placenta |
| XRN1       | 1.1266602 | 0.0277474 | CD16\_Monocyte\_Placenta |
| MCU        | 1.1261503 | 0.0474113 | CD16\_Monocyte\_Placenta |
| PML        | 1.1108008 | 0.0448984 | CD16\_Monocyte\_Placenta |
| SAMD9      | 1.0633838 | 0.0491920 | CD16\_Monocyte\_Placenta |
| PARP9      | 1.0421847 | 0.0307320 | CD16\_Monocyte\_Placenta |
| C1GALT1    | 0.9573305 | 0.0316786 | CD16\_Monocyte\_Placenta |
| PAG1       | 0.9454555 | 0.0370493 | CD16\_Monocyte\_Placenta |
| ESYT2      | 0.9122061 | 0.0422489 | CD16\_Monocyte\_Placenta |
| LY6E       | 0.8868426 | 0.0474113 | CD16\_Monocyte\_Placenta |
| COG5       | 0.8423529 | 0.0474113 | CD16\_Monocyte\_Placenta |
| DENND1A    | 0.8420733 | 0.0474113 | CD16\_Monocyte\_Placenta |
| NSMCE2     | 0.8035747 | 0.0474113 | CD16\_Monocyte\_Placenta |
| RANBP9     | 0.6810907 | 0.0474113 | CD16\_Monocyte\_Placenta |

\[1\] “Genes downregulated in disease”

| genes      |      logFC |      padj | cell\_type               |
|:-----------|-----------:|----------:|:-------------------------|
| UBAC1      | -0.7044889 | 0.0474113 | CD16\_Monocyte\_Placenta |
| PA2G4      | -0.7771543 | 0.0474113 | CD16\_Monocyte\_Placenta |
| PPA1       | -1.0399688 | 0.0422489 | CD16\_Monocyte\_Placenta |
| G0S2       | -3.2892885 | 0.0228814 | CD16\_Monocyte\_Placenta |
| AC004556.3 | -5.7965714 | 0.0422489 | CD16\_Monocyte\_Placenta |

\[1\] “Results for celltype:” \[1\] “Nonclassical-monocyte\_Placenta”
\[1\] “Genes upregulated in disease”

| genes   |    logFC |      padj | cell\_type                      |
|:--------|---------:|----------:|:--------------------------------|
| ERAP2   | 2.478030 | 0.0000104 | Nonclassical-monocyte\_Placenta |
| DOCK4   | 1.909138 | 0.0339971 | Nonclassical-monocyte\_Placenta |
| IFI6    | 1.508496 | 0.0230962 | Nonclassical-monocyte\_Placenta |
| TBC1D2  | 1.467106 | 0.0339971 | Nonclassical-monocyte\_Placenta |
| MAML2   | 1.357901 | 0.0230962 | Nonclassical-monocyte\_Placenta |
| RBM47   | 1.315087 | 0.0286528 | Nonclassical-monocyte\_Placenta |
| ARMH3   | 1.267504 | 0.0322591 | Nonclassical-monocyte\_Placenta |
| ARHGEF3 | 1.189951 | 0.0475429 | Nonclassical-monocyte\_Placenta |
| ZFAND3  | 1.065558 | 0.0458928 | Nonclassical-monocyte\_Placenta |

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type                      |
|:------|----------:|----------:|:--------------------------------|
| TESC  | -1.195234 | 0.0458928 | Nonclassical-monocyte\_Placenta |

\[1\] “Results for celltype:” \[1\] “Macrophage\_Placenta” \[1\] “Genes
upregulated in disease”

| genes   |    logFC |      padj | cell\_type           |
|:--------|---------:|----------:|:---------------------|
| GSDMA   | 6.817857 | 0.0000003 | Macrophage\_Placenta |
| ZNF385D | 4.719401 | 0.0001532 | Macrophage\_Placenta |
| ERAP2   | 2.656507 | 0.0000074 | Macrophage\_Placenta |
| HSD3B7  | 1.454662 | 0.0015560 | Macrophage\_Placenta |
| ELL2    | 1.391345 | 0.0264884 | Macrophage\_Placenta |
| DYSF    | 1.315364 | 0.0109291 | Macrophage\_Placenta |
| RMDN3   | 1.225602 | 0.0003062 | Macrophage\_Placenta |

\[1\] “Genes downregulated in disease”

| genes |     logFC |     padj | cell\_type           |
|:------|----------:|---------:|:---------------------|
| FOLR2 | -1.986148 | 2.77e-05 | Macrophage\_Placenta |

\[1\] “Results for celltype:” \[1\] “CTB\_Placenta” \[1\] “Genes
upregulated in disease”

| genes   |     logFC |      padj | cell\_type    |
|:--------|----------:|----------:|:--------------|
| KRT17   | 8.0591945 | 0.0074890 | CTB\_Placenta |
| CD63    | 3.0557487 | 0.0031943 | CTB\_Placenta |
| TPM1    | 2.7822088 | 0.0278002 | CTB\_Placenta |
| KRT19   | 2.4626154 | 0.0007084 | CTB\_Placenta |
| EFHD2   | 2.4159743 | 0.0016580 | CTB\_Placenta |
| IGF2    | 2.3913220 | 0.0074890 | CTB\_Placenta |
| CYP11A1 | 1.9426768 | 0.0111666 | CTB\_Placenta |
| GPX3    | 1.8811178 | 0.0184842 | CTB\_Placenta |
| KRT8    | 1.7773337 | 0.0267115 | CTB\_Placenta |
| NEK6    | 1.7691264 | 0.0007696 | CTB\_Placenta |
| S100P   | 1.7046765 | 0.0207815 | CTB\_Placenta |
| BASP1   | 1.5624219 | 0.0275635 | CTB\_Placenta |
| HMGB3   | 1.4895600 | 0.0207350 | CTB\_Placenta |
| TPT1    | 1.3857254 | 0.0287547 | CTB\_Placenta |
| RPL37A  | 1.3825471 | 0.0441543 | CTB\_Placenta |
| KRT18   | 1.3753116 | 0.0424593 | CTB\_Placenta |
| RPL36A  | 1.1705273 | 0.0111666 | CTB\_Placenta |
| RPL21   | 1.1256549 | 0.0253531 | CTB\_Placenta |
| RPS8    | 1.1184134 | 0.0142013 | CTB\_Placenta |
| RPL37   | 1.0812174 | 0.0176589 | CTB\_Placenta |
| RPS14   | 0.9699441 | 0.0423775 | CTB\_Placenta |
| RPL22   | 0.8723196 | 0.0267115 | CTB\_Placenta |
| RPS3A   | 0.8614950 | 0.0111666 | CTB\_Placenta |
| RPL13   | 0.8405242 | 0.0111666 | CTB\_Placenta |
| RPL15   | 0.8101662 | 0.0305816 | CTB\_Placenta |
| RPS6    | 0.7697740 | 0.0267115 | CTB\_Placenta |
| MIF     | 0.7171788 | 0.0424593 | CTB\_Placenta |
| RPL26   | 0.6839925 | 0.0441543 | CTB\_Placenta |
| RPS27A  | 0.6536248 | 0.0267115 | CTB\_Placenta |
| RPL24   | 0.6399412 | 0.0441543 | CTB\_Placenta |

\[1\] “Genes downregulated in disease”

| genes   |      logFC |      padj | cell\_type    |
|:--------|-----------:|----------:|:--------------|
| NDUFA6  | -0.8989871 | 0.0111666 | CTB\_Placenta |
| SMARCA2 | -1.1356401 | 0.0423775 | CTB\_Placenta |
| SEPTIN9 | -1.1755554 | 0.0163931 | CTB\_Placenta |
| PGD     | -1.4512298 | 0.0074890 | CTB\_Placenta |
| PHYH    | -1.4584708 | 0.0424593 | CTB\_Placenta |
| PIR     | -2.0282161 | 0.0016580 | CTB\_Placenta |
| MPP7    | -2.0712687 | 0.0036892 | CTB\_Placenta |

\[1\] “Results for celltype:” \[1\] “STB\_Placenta” \[1\] “Genes
upregulated in disease”

| genes      |    logFC |      padj | cell\_type    |
|:-----------|---------:|----------:|:--------------|
| IGFBP3     | 9.268516 | 0.0305549 | STB\_Placenta |
| XACT       | 7.943229 | 0.0219478 | STB\_Placenta |
| FSTL3      | 7.876156 | 0.0001295 | STB\_Placenta |
| ADAM2      | 7.654413 | 0.0102377 | STB\_Placenta |
| TMEM45A    | 7.092075 | 0.0099802 | STB\_Placenta |
| AC087857.1 | 6.174981 | 0.0307443 | STB\_Placenta |
| AC026167.1 | 6.163577 | 0.0307443 | STB\_Placenta |
| AC106798.1 | 5.952874 | 0.0305549 | STB\_Placenta |
| LEP        | 5.298699 | 0.0022922 | STB\_Placenta |
| MIR210HG   | 5.262350 | 0.0016703 | STB\_Placenta |
| ARMS2      | 4.933283 | 0.0022922 | STB\_Placenta |
| PNCK       | 4.874359 | 0.0147225 | STB\_Placenta |
| IFT27      | 4.799281 | 0.0001295 | STB\_Placenta |
| HTRA4      | 4.737491 | 0.0000352 | STB\_Placenta |
| FLT1       | 4.657075 | 0.0003673 | STB\_Placenta |
| BHLHE40    | 3.690681 | 0.0307443 | STB\_Placenta |
| EGLN3      | 3.608907 | 0.0267081 | STB\_Placenta |
| IGF2       | 3.552155 | 0.0286189 | STB\_Placenta |
| LPL        | 3.510080 | 0.0120950 | STB\_Placenta |
| AIF1L      | 3.470693 | 0.0012923 | STB\_Placenta |
| CRH        | 3.057630 | 0.0248603 | STB\_Placenta |
| KIAA0753   | 3.008544 | 0.0305549 | STB\_Placenta |
| HTRA1      | 2.986958 | 0.0107169 | STB\_Placenta |
| PPP1R1C    | 2.911112 | 0.0305549 | STB\_Placenta |
| LIMS2      | 2.867949 | 0.0332625 | STB\_Placenta |
| SH3BP5     | 2.788337 | 0.0244618 | STB\_Placenta |
| KRT86      | 2.695616 | 0.0120950 | STB\_Placenta |
| DERL3      | 2.552763 | 0.0307443 | STB\_Placenta |
| GALK1      | 2.515990 | 0.0317350 | STB\_Placenta |
| LY6D       | 2.439876 | 0.0062696 | STB\_Placenta |
| KRT81      | 2.419357 | 0.0332625 | STB\_Placenta |
| SPAG4      | 2.404235 | 0.0305549 | STB\_Placenta |
| SMIM3      | 2.304228 | 0.0332625 | STB\_Placenta |
| MIR193BHG  | 2.253136 | 0.0120950 | STB\_Placenta |
| AL731684.1 | 2.160598 | 0.0219279 | STB\_Placenta |
| CCDC183    | 2.150428 | 0.0307443 | STB\_Placenta |
| GDPD3      | 2.141044 | 0.0120950 | STB\_Placenta |
| TGFB1      | 2.088513 | 0.0180621 | STB\_Placenta |
| INHA       | 2.085179 | 0.0098246 | STB\_Placenta |
| ARMCX6     | 2.074872 | 0.0352423 | STB\_Placenta |
| FAAP20     | 1.997003 | 0.0028930 | STB\_Placenta |
| SH3PXD2A   | 1.680326 | 0.0420469 | STB\_Placenta |
| EBLN3P     | 1.639120 | 0.0305549 | STB\_Placenta |
| RPL5       | 1.450148 | 0.0440828 | STB\_Placenta |
| EIF2S3     | 1.342206 | 0.0305549 | STB\_Placenta |
| GDI2       | 1.268377 | 0.0307443 | STB\_Placenta |
| PLPP5      | 1.247648 | 0.0377406 | STB\_Placenta |

\[1\] “Genes downregulated in disease”

| genes      |     logFC |      padj | cell\_type    |
|:-----------|----------:|----------:|:--------------|
| CASP4      | -1.052883 | 0.0305549 | STB\_Placenta |
| DNMT3A     | -1.156918 | 0.0495340 | STB\_Placenta |
| RABGAP1L   | -1.218552 | 0.0041003 | STB\_Placenta |
| MORC4      | -1.343658 | 0.0305549 | STB\_Placenta |
| MTAP       | -1.366813 | 0.0248603 | STB\_Placenta |
| CDYL2      | -1.442092 | 0.0016178 | STB\_Placenta |
| TWSG1      | -1.567760 | 0.0303146 | STB\_Placenta |
| RASSF6     | -1.618461 | 0.0420469 | STB\_Placenta |
| MRAS       | -1.670037 | 0.0093339 | STB\_Placenta |
| ABCG2      | -1.714422 | 0.0167089 | STB\_Placenta |
| CDKN2B     | -1.724545 | 0.0120950 | STB\_Placenta |
| PYGL       | -1.764234 | 0.0420469 | STB\_Placenta |
| ALDH4A1    | -1.769686 | 0.0019607 | STB\_Placenta |
| PLCD3      | -1.803937 | 0.0219669 | STB\_Placenta |
| SYNPO2L    | -1.861298 | 0.0305549 | STB\_Placenta |
| NAV1       | -1.940413 | 0.0495340 | STB\_Placenta |
| AMPD3      | -1.992828 | 0.0102377 | STB\_Placenta |
| TBC1D9     | -1.999732 | 0.0307443 | STB\_Placenta |
| NCF4-AS1   | -2.113655 | 0.0003673 | STB\_Placenta |
| STON2      | -2.154496 | 0.0120950 | STB\_Placenta |
| GSTA3      | -2.177000 | 0.0226497 | STB\_Placenta |
| MYO1B      | -2.198687 | 0.0164600 | STB\_Placenta |
| SNTB1      | -2.436976 | 0.0307443 | STB\_Placenta |
| MCPH1-AS1  | -2.446668 | 0.0012923 | STB\_Placenta |
| NEDD4L     | -2.478179 | 0.0001295 | STB\_Placenta |
| KIAA1211   | -2.558052 | 0.0120950 | STB\_Placenta |
| RAB27B     | -2.609235 | 0.0078698 | STB\_Placenta |
| ASB2       | -2.647818 | 0.0028930 | STB\_Placenta |
| ACOT1      | -2.940030 | 0.0093339 | STB\_Placenta |
| P3H2       | -3.339030 | 0.0307443 | STB\_Placenta |
| NFIA-AS2   | -3.513977 | 0.0082286 | STB\_Placenta |
| ADAMTSL1   | -3.594766 | 0.0049143 | STB\_Placenta |
| CPQ        | -4.362170 | 0.0102377 | STB\_Placenta |
| LINC02109  | -8.250863 | 0.0018887 | STB\_Placenta |
| AC008825.1 | -8.393533 | 0.0062696 | STB\_Placenta |

# Volcano plots in Late disease

\[1\] “Results for celltype:” \[1\] “Naive\_B\_Placenta” \[1\] “Genes
upregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Genes downregulated in disease”

| genes    |     logFC |      padj | cell\_type         |
|:---------|----------:|----------:|:-------------------|
| IGLV1-44 | -2.524552 | 0.0007507 | Naive\_B\_Placenta |
| BACE2    | -2.946636 | 0.0007507 | Naive\_B\_Placenta |
| IGKV1-8  | -3.101062 | 0.0007173 | Naive\_B\_Placenta |
| TCL6     | -4.660098 | 0.0007173 | Naive\_B\_Placenta |
| HSPB1    | -4.922312 | 0.0099073 | Naive\_B\_Placenta |
| CLTB     | -6.494036 | 0.0007173 | Naive\_B\_Placenta |
| FDX1     | -6.641065 | 0.0036568 | Naive\_B\_Placenta |

\[1\] “Results for celltype:” \[1\] “CD8\_T\_Naive\_CM\_Placenta” \[1\]
“Genes upregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Genes downregulated in disease”

| genes      |    logFC |      padj | cell\_type                  |
|:-----------|---------:|----------:|:----------------------------|
| AC013264.1 | -1.45493 | 0.0351599 | CD8\_T\_Naive\_CM\_Placenta |
| HBG2       | -6.31538 | 0.0351599 | CD8\_T\_Naive\_CM\_Placenta |

\[1\] “Results for celltype:” \[1\] “MAIT\_Placenta” \[1\] “Genes
upregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Genes downregulated in disease”

| genes      |     logFC |      padj | cell\_type     |
|:-----------|----------:|----------:|:---------------|
| NCR3       | -1.043648 | 0.0115725 | MAIT\_Placenta |
| AC010967.1 | -4.863378 | 0.0088866 | MAIT\_Placenta |

\[1\] “Results for celltype:” \[1\] “NKT\_Placenta” \[1\] “Genes
upregulated in disease”

| genes  |    logFC |      padj | cell\_type    |
|:-------|---------:|----------:|:--------------|
| INPP4B | 1.561871 | 0.0109038 | NKT\_Placenta |

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type    |
|:------|----------:|----------:|:--------------|
| STMN1 | -2.790383 | 0.0245263 | NKT\_Placenta |

\[1\] “Results for celltype:” \[1\] “CD16\_Monocyte\_Placenta” \[1\]
“Genes upregulated in disease”

| genes |     logFC |      padj | cell\_type               |
|:------|----------:|----------:|:-------------------------|
| FLNB  | 1.7144214 | 0.0263878 | CD16\_Monocyte\_Placenta |
| HECA  | 0.7971445 | 0.0405733 | CD16\_Monocyte\_Placenta |

\[1\] “Genes downregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Results for celltype:” \[1\] “Nonclassical-monocyte\_Placenta”
\[1\] “Genes upregulated in disease”

| genes  |    logFC |      padj | cell\_type                      |
|:-------|---------:|----------:|:--------------------------------|
| CD163  | 5.035168 | 0.0074705 | Nonclassical-monocyte\_Placenta |
| RNASE2 | 3.969727 | 0.0302151 | Nonclassical-monocyte\_Placenta |

\[1\] “Genes downregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Results for celltype:” \[1\] “Intermediate-macrophage\_Placenta”
\[1\] “Genes upregulated in disease”

| genes    |    logFC |      padj | cell\_type                        |
|:---------|---------:|----------:|:----------------------------------|
| TNFRSF21 | 2.228944 | 0.0407418 | Intermediate-macrophage\_Placenta |

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type                        |
|:------|----------:|----------:|:----------------------------------|
| DUSP2 | -3.057637 | 0.0019314 | Intermediate-macrophage\_Placenta |

\[1\] “Results for celltype:” \[1\] “Macrophage\_Placenta” \[1\] “Genes
upregulated in disease”

| genes |     logFC |      padj | cell\_type           |
|:------|----------:|----------:|:---------------------|
| CSH1  | 5.9571447 | 0.0293425 | Macrophage\_Placenta |
| RMDN3 | 0.8400423 | 0.0293425 | Macrophage\_Placenta |

\[1\] “Genes downregulated in disease”

| genes   |     logFC |      padj | cell\_type           |
|:--------|----------:|----------:|:---------------------|
| PLA2G2D | -4.300843 | 0.0293425 | Macrophage\_Placenta |

\[1\] “Results for celltype:” \[1\] “STB\_Placenta” \[1\] “Genes
upregulated in disease”

| genes |    logFC |      padj | cell\_type    |
|:------|---------:|----------:|:--------------|
| LEP   | 4.127213 | 0.0007515 | STB\_Placenta |
| HTRA4 | 3.112228 | 0.0000189 | STB\_Placenta |

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type    |
|:------|----------:|----------:|:--------------|
| TRIM5 | -1.558006 | 0.0321418 | STB\_Placenta |

# Volcano plots in averagePE

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes |    logFC |      padj | cell\_type         |
|:------|---------:|----------:|:-------------------|
| THRB  | 1.618126 | 0.0078699 | Naive\_B\_Placenta |

\[1\] “Genes downregulated in disease”

| genes     |     logFC |      padj | cell\_type         |
|:----------|----------:|----------:|:-------------------|
| IGLV1-44  | -1.707892 | 0.0016365 | Naive\_B\_Placenta |
| IGKV1-8   | -1.992690 | 0.0016365 | Naive\_B\_Placenta |
| HSPB1     | -2.742736 | 0.0421591 | Naive\_B\_Placenta |
| TCL6      | -2.766002 | 0.0009479 | Naive\_B\_Placenta |
| IGKV1-16  | -3.370036 | 0.0441048 | Naive\_B\_Placenta |
| CLTB      | -3.609773 | 0.0016365 | Naive\_B\_Placenta |
| FDX1      | -3.988542 | 0.0024688 | Naive\_B\_Placenta |
| IGHV7-4-1 | -6.931084 | 0.0009479 | Naive\_B\_Placenta |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes     |    logFC |     padj | cell\_type          |
|:----------|---------:|---------:|:--------------------|
| IGKV1D-16 | 6.365848 | 0.019093 | Memory\_B\_Placenta |

\[1\] “Genes downregulated in disease”

| genes    |     logFC |      padj | cell\_type          |
|:---------|----------:|----------:|:--------------------|
| HOPX     | -2.742520 | 0.0481495 | Memory\_B\_Placenta |
| IGHV3-72 | -5.066492 | 0.0008680 | Memory\_B\_Placenta |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes |     logFC |      padj | cell\_type                  |
|:------|----------:|----------:|:----------------------------|
| MPP7  | 0.9945184 | 0.0009368 | CD4\_T\_Naive\_CM\_Placenta |

\[1\] “Genes downregulated in disease”

| genes  |     logFC |      padj | cell\_type                  |
|:-------|----------:|----------:|:----------------------------|
| TMEM80 | -1.386272 | 0.0004162 | CD4\_T\_Naive\_CM\_Placenta |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes |     logFC |     padj | cell\_type           |
|:------|----------:|---------:|:---------------------|
| CASP8 | 0.9827241 | 0.021518 | FoxP3-Treg\_Placenta |

\[1\] “Genes downregulated in disease”

| genes |     logFC |     padj | cell\_type           |
|:------|----------:|---------:|:---------------------|
| CSH1  | -4.128652 | 0.021518 | FoxP3-Treg\_Placenta |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes      |     logFC |      padj | cell\_type         |
|:-----------|----------:|----------:|:-------------------|
| PRKN       | 1.8030679 | 0.0288538 | CD4\_TEM\_Placenta |
| ERAP2      | 1.1687885 | 0.0121909 | CD4\_TEM\_Placenta |
| FP671120.4 | 1.0982372 | 0.0442610 | CD4\_TEM\_Placenta |
| GAS5       | 1.0494673 | 0.0116275 | CD4\_TEM\_Placenta |
| AC119396.1 | 0.9358105 | 0.0257464 | CD4\_TEM\_Placenta |
| HECA       | 0.8941041 | 0.0032246 | CD4\_TEM\_Placenta |
| PIM2       | 0.7599829 | 0.0442610 | CD4\_TEM\_Placenta |
| RASA3      | 0.7280707 | 0.0442610 | CD4\_TEM\_Placenta |
| INTS6      | 0.7215686 | 0.0442610 | CD4\_TEM\_Placenta |
| RICTOR     | 0.6945980 | 0.0442610 | CD4\_TEM\_Placenta |
| ANK3       | 0.6138591 | 0.0442610 | CD4\_TEM\_Placenta |

\[1\] “Genes downregulated in disease”

| genes       |      logFC |      padj | cell\_type         |
|:------------|-----------:|----------:|:-------------------|
| CLIC1       | -0.5692723 | 0.0442610 | CD4\_TEM\_Placenta |
| NCF1        | -0.7492597 | 0.0442610 | CD4\_TEM\_Placenta |
| CKLF        | -0.8148055 | 0.0051300 | CD4\_TEM\_Placenta |
| LDHA        | -0.8303667 | 0.0000326 | CD4\_TEM\_Placenta |
| MT2A        | -0.8339731 | 0.0390427 | CD4\_TEM\_Placenta |
| RGS1        | -1.0209854 | 0.0442610 | CD4\_TEM\_Placenta |
| LGALS1      | -1.0465160 | 0.0442610 | CD4\_TEM\_Placenta |
| MIR4435-2HG | -1.0642982 | 0.0288538 | CD4\_TEM\_Placenta |
| AGPAT2      | -1.2062907 | 0.0433396 | CD4\_TEM\_Placenta |
| CAPG        | -1.2658453 | 0.0151541 | CD4\_TEM\_Placenta |
| NDFIP2      | -1.9905382 | 0.0121909 | CD4\_TEM\_Placenta |
| TRBV5-6     | -2.2517119 | 0.0442610 | CD4\_TEM\_Placenta |
| XCL1        | -2.9821562 | 0.0001958 | CD4\_TEM\_Placenta |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes |    logFC |      padj | cell\_type                  |
|:------|---------:|----------:|:----------------------------|
| ERAP2 | 1.152982 | 0.0073921 | CD8\_T\_Naive\_CM\_Placenta |

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type                  |
|:------|----------:|----------:|:----------------------------|
| GLRX  | -2.028561 | 0.0266611 | CD8\_T\_Naive\_CM\_Placenta |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes |     logFC |      padj | cell\_type             |
|:------|----------:|----------:|:-----------------------|
| EDA   | 2.1778092 | 0.0170352 | GZMK\_CD8\_T\_Placenta |
| TXK   | 0.9203278 | 0.0170352 | GZMK\_CD8\_T\_Placenta |

\[1\] “Genes downregulated in disease”

| genes |      logFC |     padj | cell\_type             |
|:------|-----------:|---------:|:-----------------------|
| NKG7  | -0.7586009 | 0.017065 | GZMK\_CD8\_T\_Placenta |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes      |     logFC |      padj | cell\_type         |
|:-----------|----------:|----------:|:-------------------|
| AC022075.1 | 0.9881724 | 0.0366754 | CD8\_TEM\_Placenta |
| ERAP2      | 0.8949720 | 0.0022397 | CD8\_TEM\_Placenta |

\[1\] “Genes downregulated in disease”

| genes  |     logFC |      padj | cell\_type         |
|:-------|----------:|----------:|:-------------------|
| TRAV27 | -2.749687 | 0.0069387 | CD8\_TEM\_Placenta |
| TRBV27 | -2.847242 | 0.0105434 | CD8\_TEM\_Placenta |
| TRBV28 | -3.403608 | 0.0013980 | CD8\_TEM\_Placenta |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes |    logFC |      padj | cell\_type     |
|:------|---------:|----------:|:---------------|
| ERAP2 | 1.415390 | 0.0003423 | MAIT\_Placenta |
| A2M   | 1.227352 | 0.0141320 | MAIT\_Placenta |
| PSAP  | 1.059560 | 0.0001092 | MAIT\_Placenta |

\[1\] “Genes downregulated in disease”

| genes      |     logFC |      padj | cell\_type     |
|:-----------|----------:|----------:|:---------------|
| AC010967.1 | -2.972952 | 0.0047448 | MAIT\_Placenta |
| LINC02446  | -3.078081 | 0.0258222 | MAIT\_Placenta |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes |    logFC |      padj | cell\_type    |
|:------|---------:|----------:|:--------------|
| EXOC4 | 1.026525 | 0.0370525 | NKT\_Placenta |

\[1\] “Genes downregulated in disease”

| genes |     logFC |      padj | cell\_type    |
|:------|----------:|----------:|:--------------|
| ITM2A | -1.067548 | 0.0445781 | NKT\_Placenta |
| MT2A  | -1.072850 | 0.0370525 | NKT\_Placenta |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes      |     logFC |      padj | cell\_type         |
|:-----------|----------:|----------:|:-------------------|
| USP9Y      | 4.5508841 | 0.0000509 | CD56\_NK\_Placenta |
| KRT72      | 3.5623192 | 0.0110470 | CD56\_NK\_Placenta |
| UTY        | 2.8505479 | 0.0180430 | CD56\_NK\_Placenta |
| PRKY       | 2.6148101 | 0.0262408 | CD56\_NK\_Placenta |
| ABHD15-AS1 | 2.0991918 | 0.0003530 | CD56\_NK\_Placenta |
| CHSY3      | 1.9562616 | 0.0180385 | CD56\_NK\_Placenta |
| CD3G       | 1.9265124 | 0.0207034 | CD56\_NK\_Placenta |
| MTSS1      | 1.9115293 | 0.0000383 | CD56\_NK\_Placenta |
| CHD7       | 1.8707024 | 0.0251939 | CD56\_NK\_Placenta |
| IGF1R      | 1.7641617 | 0.0003530 | CD56\_NK\_Placenta |
| CDHR3      | 1.7355620 | 0.0387297 | CD56\_NK\_Placenta |
| DNHD1      | 1.5962531 | 0.0394588 | CD56\_NK\_Placenta |
| LINGO2     | 1.5742145 | 0.0233358 | CD56\_NK\_Placenta |
| SLCO4C1    | 1.5277956 | 0.0139174 | CD56\_NK\_Placenta |
| RAP1GAP2   | 1.4830665 | 0.0116126 | CD56\_NK\_Placenta |
| CD3D       | 1.4792381 | 0.0341229 | CD56\_NK\_Placenta |
| MIR646HG   | 1.4743517 | 0.0268397 | CD56\_NK\_Placenta |
| ATXN7L1    | 1.4633127 | 0.0241612 | CD56\_NK\_Placenta |
| AC119396.1 | 1.4211309 | 0.0447150 | CD56\_NK\_Placenta |
| SBK1       | 1.3744716 | 0.0403385 | CD56\_NK\_Placenta |
| PPM1L      | 1.3099922 | 0.0131153 | CD56\_NK\_Placenta |
| MYLIP      | 1.3047188 | 0.0011129 | CD56\_NK\_Placenta |
| RUFY2      | 1.2823407 | 0.0262408 | CD56\_NK\_Placenta |
| DPEP2      | 1.2746019 | 0.0301900 | CD56\_NK\_Placenta |
| BNC2       | 1.1997038 | 0.0455731 | CD56\_NK\_Placenta |
| TLE1       | 1.1960890 | 0.0288238 | CD56\_NK\_Placenta |
| CAMK2D     | 1.1739021 | 0.0116126 | CD56\_NK\_Placenta |
| AMMECR1    | 1.1392208 | 0.0475205 | CD56\_NK\_Placenta |
| CAMK1D     | 1.1386508 | 0.0185106 | CD56\_NK\_Placenta |
| LRIG2      | 1.1194775 | 0.0185106 | CD56\_NK\_Placenta |
| WWOX       | 1.1159841 | 0.0380111 | CD56\_NK\_Placenta |
| CRAMP1     | 1.1132939 | 0.0224375 | CD56\_NK\_Placenta |
| ZNF831     | 1.1078026 | 0.0116126 | CD56\_NK\_Placenta |
| LAIR1      | 1.0900993 | 0.0352425 | CD56\_NK\_Placenta |
| ERAP2      | 1.0837823 | 0.0131153 | CD56\_NK\_Placenta |
| NOTCH1     | 1.0767936 | 0.0077167 | CD56\_NK\_Placenta |
| AC068587.4 | 1.0580556 | 0.0421648 | CD56\_NK\_Placenta |
| LCORL      | 1.0531838 | 0.0263441 | CD56\_NK\_Placenta |
| CDR2       | 1.0448405 | 0.0160010 | CD56\_NK\_Placenta |
| SNX29      | 1.0385011 | 0.0116126 | CD56\_NK\_Placenta |
| TMEM181    | 1.0384938 | 0.0150905 | CD56\_NK\_Placenta |
| PDE4D      | 1.0380824 | 0.0425430 | CD56\_NK\_Placenta |
| ARHGAP12   | 1.0303640 | 0.0357740 | CD56\_NK\_Placenta |
| OSBPL5     | 1.0215867 | 0.0483532 | CD56\_NK\_Placenta |
| TGFBR3     | 1.0199721 | 0.0405207 | CD56\_NK\_Placenta |
| FAM193B    | 1.0069652 | 0.0387297 | CD56\_NK\_Placenta |
| UST        | 1.0010398 | 0.0387297 | CD56\_NK\_Placenta |
| EXT1       | 0.9983924 | 0.0483532 | CD56\_NK\_Placenta |
| ARL15      | 0.9949508 | 0.0214200 | CD56\_NK\_Placenta |
| VWA8       | 0.9943326 | 0.0166677 | CD56\_NK\_Placenta |
| PXN        | 0.9851729 | 0.0272101 | CD56\_NK\_Placenta |
| MT-ATP8    | 0.9802692 | 0.0290256 | CD56\_NK\_Placenta |
| IQSEC1     | 0.9792868 | 0.0132755 | CD56\_NK\_Placenta |
| NCOA6      | 0.9785537 | 0.0185106 | CD56\_NK\_Placenta |
| NLRP1      | 0.9712759 | 0.0029665 | CD56\_NK\_Placenta |
| MGAT4A     | 0.9642271 | 0.0177753 | CD56\_NK\_Placenta |
| HDAC7      | 0.9624116 | 0.0132755 | CD56\_NK\_Placenta |
| AL645568.1 | 0.9582343 | 0.0204932 | CD56\_NK\_Placenta |
| MT-CYB     | 0.9471342 | 0.0258986 | CD56\_NK\_Placenta |
| SAMD4B     | 0.9427491 | 0.0111238 | CD56\_NK\_Placenta |
| PCCA       | 0.9407113 | 0.0411679 | CD56\_NK\_Placenta |
| ARHGAP27   | 0.9385312 | 0.0131153 | CD56\_NK\_Placenta |
| SLC12A6    | 0.9344685 | 0.0355530 | CD56\_NK\_Placenta |
| CLASP1     | 0.9312348 | 0.0215711 | CD56\_NK\_Placenta |
| SPON2      | 0.9296332 | 0.0407041 | CD56\_NK\_Placenta |
| KANSL1     | 0.9288893 | 0.0077167 | CD56\_NK\_Placenta |
| MAML2      | 0.9284749 | 0.0298018 | CD56\_NK\_Placenta |
| CWC27      | 0.9266298 | 0.0263441 | CD56\_NK\_Placenta |
| NCALD      | 0.9264896 | 0.0421648 | CD56\_NK\_Placenta |
| EXOC4      | 0.9246305 | 0.0158571 | CD56\_NK\_Placenta |
| RHBDF2     | 0.9196564 | 0.0462922 | CD56\_NK\_Placenta |
| ABR        | 0.9163329 | 0.0203983 | CD56\_NK\_Placenta |
| COG5       | 0.9133216 | 0.0218707 | CD56\_NK\_Placenta |
| PAM        | 0.9106896 | 0.0251939 | CD56\_NK\_Placenta |
| KDM7A      | 0.8987172 | 0.0483532 | CD56\_NK\_Placenta |
| GBF1       | 0.8920351 | 0.0387297 | CD56\_NK\_Placenta |
| RMDN3      | 0.8859705 | 0.0483532 | CD56\_NK\_Placenta |
| TNKS       | 0.8832460 | 0.0357740 | CD56\_NK\_Placenta |
| FAM172A    | 0.8662580 | 0.0340836 | CD56\_NK\_Placenta |
| ST6GAL1    | 0.8662176 | 0.0483532 | CD56\_NK\_Placenta |
| RFFL       | 0.8485369 | 0.0407041 | CD56\_NK\_Placenta |
| VAV1       | 0.8464553 | 0.0372747 | CD56\_NK\_Placenta |
| GAB3       | 0.8420960 | 0.0342620 | CD56\_NK\_Placenta |
| RIPOR2     | 0.8236972 | 0.0483532 | CD56\_NK\_Placenta |
| PSME4      | 0.8203740 | 0.0387297 | CD56\_NK\_Placenta |
| PRR5       | 0.8202325 | 0.0270892 | CD56\_NK\_Placenta |
| RORA-AS1   | 0.8193173 | 0.0378188 | CD56\_NK\_Placenta |
| NCBP1      | 0.8182979 | 0.0463114 | CD56\_NK\_Placenta |
| ZBTB16     | 0.8145112 | 0.0378188 | CD56\_NK\_Placenta |
| EPS15L1    | 0.8144035 | 0.0387297 | CD56\_NK\_Placenta |
| SUMF1      | 0.8120751 | 0.0475205 | CD56\_NK\_Placenta |
| SIDT1      | 0.8085921 | 0.0279084 | CD56\_NK\_Placenta |
| CREBBP     | 0.7999577 | 0.0185106 | CD56\_NK\_Placenta |
| USP33      | 0.7909158 | 0.0131153 | CD56\_NK\_Placenta |
| RAB27B     | 0.7882585 | 0.0270892 | CD56\_NK\_Placenta |
| SSH2       | 0.7853514 | 0.0185106 | CD56\_NK\_Placenta |
| RABGAP1    | 0.7850972 | 0.0357740 | CD56\_NK\_Placenta |
| ATM        | 0.7742634 | 0.0216501 | CD56\_NK\_Placenta |
| SYTL2      | 0.7732484 | 0.0257740 | CD56\_NK\_Placenta |
| ARMC8      | 0.7696084 | 0.0368635 | CD56\_NK\_Placenta |
| BCAS3      | 0.7673819 | 0.0342620 | CD56\_NK\_Placenta |
| IKZF1      | 0.7660476 | 0.0256194 | CD56\_NK\_Placenta |
| MT-CO3     | 0.7640618 | 0.0298018 | CD56\_NK\_Placenta |
| MCTP2      | 0.7523491 | 0.0407041 | CD56\_NK\_Placenta |
| SHLD2      | 0.7468036 | 0.0494793 | CD56\_NK\_Placenta |
| PIK3R5     | 0.7357563 | 0.0483532 | CD56\_NK\_Placenta |
| CELF2      | 0.7327456 | 0.0130221 | CD56\_NK\_Placenta |
| ZFPM1      | 0.7307325 | 0.0378188 | CD56\_NK\_Placenta |
| RABGAP1L   | 0.7253259 | 0.0298018 | CD56\_NK\_Placenta |
| FNIP1      | 0.7148024 | 0.0496699 | CD56\_NK\_Placenta |
| SRPK2      | 0.7013360 | 0.0475205 | CD56\_NK\_Placenta |
| CARD11     | 0.6948321 | 0.0233690 | CD56\_NK\_Placenta |
| TBC1D5     | 0.6824771 | 0.0369224 | CD56\_NK\_Placenta |
| MKLN1      | 0.6752199 | 0.0453948 | CD56\_NK\_Placenta |
| CELF1      | 0.6732094 | 0.0368635 | CD56\_NK\_Placenta |
| LRBA       | 0.6719425 | 0.0483532 | CD56\_NK\_Placenta |
| ARID1A     | 0.6656806 | 0.0338293 | CD56\_NK\_Placenta |
| DPYD       | 0.6548339 | 0.0218969 | CD56\_NK\_Placenta |
| RICTOR     | 0.6496985 | 0.0342151 | CD56\_NK\_Placenta |
| FOXN3      | 0.6373144 | 0.0475205 | CD56\_NK\_Placenta |
| ARRB1      | 0.6366378 | 0.0462922 | CD56\_NK\_Placenta |
| SEC31A     | 0.6351438 | 0.0462922 | CD56\_NK\_Placenta |
| DCAF5      | 0.6316167 | 0.0483532 | CD56\_NK\_Placenta |
| TAF15      | 0.6239331 | 0.0475205 | CD56\_NK\_Placenta |
| JMJD1C     | 0.6191515 | 0.0483532 | CD56\_NK\_Placenta |
| CCDC88C    | 0.5989916 | 0.0416299 | CD56\_NK\_Placenta |
| MYO1F      | 0.5870588 | 0.0483532 | CD56\_NK\_Placenta |
| AKNA       | 0.5576196 | 0.0455945 | CD56\_NK\_Placenta |
| MYO9B      | 0.5325005 | 0.0483532 | CD56\_NK\_Placenta |

\[1\] “Genes downregulated in disease”

| genes    |      logFC |      padj | cell\_type         |
|:---------|-----------:|----------:|:-------------------|
| GNAS     | -0.5198083 | 0.0462922 | CD56\_NK\_Placenta |
| CD37     | -0.5324259 | 0.0483532 | CD56\_NK\_Placenta |
| EEF1A1   | -0.5449312 | 0.0296405 | CD56\_NK\_Placenta |
| TRIR     | -0.5497888 | 0.0483532 | CD56\_NK\_Placenta |
| PSMA6    | -0.5504083 | 0.0483532 | CD56\_NK\_Placenta |
| SNRPB    | -0.5625637 | 0.0463114 | CD56\_NK\_Placenta |
| MYL6     | -0.5758872 | 0.0387297 | CD56\_NK\_Placenta |
| ITGB1    | -0.5800629 | 0.0371944 | CD56\_NK\_Placenta |
| LMAN2    | -0.5837294 | 0.0483532 | CD56\_NK\_Placenta |
| CLIC1    | -0.5882670 | 0.0416299 | CD56\_NK\_Placenta |
| PARK7    | -0.5910805 | 0.0378188 | CD56\_NK\_Placenta |
| RPL5     | -0.5924893 | 0.0311649 | CD56\_NK\_Placenta |
| RPL12    | -0.6020269 | 0.0478986 | CD56\_NK\_Placenta |
| ATP5F1E  | -0.6155460 | 0.0483532 | CD56\_NK\_Placenta |
| RPL10    | -0.6177159 | 0.0462922 | CD56\_NK\_Placenta |
| CD160    | -0.6242479 | 0.0483532 | CD56\_NK\_Placenta |
| RPL8     | -0.6249310 | 0.0462922 | CD56\_NK\_Placenta |
| SH3BGRL3 | -0.6250433 | 0.0116126 | CD56\_NK\_Placenta |
| EEF1D    | -0.6271600 | 0.0483532 | CD56\_NK\_Placenta |
| NDUFS5   | -0.6332419 | 0.0463114 | CD56\_NK\_Placenta |
| RPL29    | -0.6344081 | 0.0474892 | CD56\_NK\_Placenta |
| COMMD7   | -0.6374317 | 0.0483532 | CD56\_NK\_Placenta |
| NDUFA1   | -0.6386361 | 0.0407041 | CD56\_NK\_Placenta |
| NDUFA11  | -0.6418488 | 0.0475205 | CD56\_NK\_Placenta |
| ERH      | -0.6448058 | 0.0386036 | CD56\_NK\_Placenta |
| RABAC1   | -0.6487527 | 0.0483532 | CD56\_NK\_Placenta |
| MAGOH    | -0.6491255 | 0.0477874 | CD56\_NK\_Placenta |
| HMGN1    | -0.6522716 | 0.0380248 | CD56\_NK\_Placenta |
| HMGN2    | -0.6567978 | 0.0335840 | CD56\_NK\_Placenta |
| POLR2J   | -0.6588807 | 0.0483532 | CD56\_NK\_Placenta |
| FAU      | -0.6624176 | 0.0381725 | CD56\_NK\_Placenta |
| RPL24    | -0.6625069 | 0.0342620 | CD56\_NK\_Placenta |
| CCDC107  | -0.6626018 | 0.0405207 | CD56\_NK\_Placenta |
| RPL14    | -0.6628097 | 0.0335840 | CD56\_NK\_Placenta |
| ARPC3    | -0.6663952 | 0.0150905 | CD56\_NK\_Placenta |
| TBCB     | -0.6723298 | 0.0378188 | CD56\_NK\_Placenta |
| TRBC2    | -0.6775190 | 0.0462922 | CD56\_NK\_Placenta |
| SEC11A   | -0.6778227 | 0.0311561 | CD56\_NK\_Placenta |
| ENO1     | -0.6779739 | 0.0251939 | CD56\_NK\_Placenta |
| RPL19    | -0.6833138 | 0.0349368 | CD56\_NK\_Placenta |
| LDHA     | -0.6850704 | 0.0462922 | CD56\_NK\_Placenta |
| RPL7A    | -0.6883326 | 0.0352425 | CD56\_NK\_Placenta |
| RPS19    | -0.6897826 | 0.0475205 | CD56\_NK\_Placenta |
| NDUFB9   | -0.6937692 | 0.0483532 | CD56\_NK\_Placenta |
| FASLG    | -0.6993224 | 0.0263441 | CD56\_NK\_Placenta |
| MRPL41   | -0.7002195 | 0.0421648 | CD56\_NK\_Placenta |
| PSMD8    | -0.7027241 | 0.0357740 | CD56\_NK\_Placenta |
| FIS1     | -0.7076602 | 0.0462922 | CD56\_NK\_Placenta |
| RSL24D1  | -0.7087280 | 0.0288238 | CD56\_NK\_Placenta |
| RPS4X    | -0.7146079 | 0.0387297 | CD56\_NK\_Placenta |
| AGTRAP   | -0.7159078 | 0.0475205 | CD56\_NK\_Placenta |
| PSME2    | -0.7169989 | 0.0258941 | CD56\_NK\_Placenta |
| SLC25A3  | -0.7203773 | 0.0077686 | CD56\_NK\_Placenta |
| PCNP     | -0.7211235 | 0.0262408 | CD56\_NK\_Placenta |
| ITGB1BP1 | -0.7211415 | 0.0378188 | CD56\_NK\_Placenta |
| MIF      | -0.7230362 | 0.0387297 | CD56\_NK\_Placenta |
| RECQL    | -0.7249362 | 0.0475205 | CD56\_NK\_Placenta |
| RNF187   | -0.7267744 | 0.0462922 | CD56\_NK\_Placenta |
| LSM4     | -0.7312888 | 0.0353756 | CD56\_NK\_Placenta |
| RPL41    | -0.7370592 | 0.0373743 | CD56\_NK\_Placenta |
| NACA     | -0.7375955 | 0.0131153 | CD56\_NK\_Placenta |
| EZR      | -0.7434479 | 0.0239274 | CD56\_NK\_Placenta |
| HIGD2A   | -0.7435993 | 0.0366326 | CD56\_NK\_Placenta |
| ANAPC11  | -0.7436933 | 0.0246106 | CD56\_NK\_Placenta |
| UBE2L6   | -0.7497226 | 0.0407041 | CD56\_NK\_Placenta |
| SAP18    | -0.7500313 | 0.0150905 | CD56\_NK\_Placenta |
| FXYD5    | -0.7526944 | 0.0258986 | CD56\_NK\_Placenta |
| GNG5     | -0.7602651 | 0.0144172 | CD56\_NK\_Placenta |
| NUTF2    | -0.7744847 | 0.0270892 | CD56\_NK\_Placenta |
| ROMO1    | -0.7770060 | 0.0263441 | CD56\_NK\_Placenta |
| NME2     | -0.7913239 | 0.0311561 | CD56\_NK\_Placenta |
| SELENOW  | -0.8011473 | 0.0052464 | CD56\_NK\_Placenta |
| SEPTIN11 | -0.8152139 | 0.0376041 | CD56\_NK\_Placenta |
| GAPDH    | -0.8174843 | 0.0214639 | CD56\_NK\_Placenta |
| COMMD6   | -0.8200231 | 0.0185106 | CD56\_NK\_Placenta |
| EMC7     | -0.8230460 | 0.0357740 | CD56\_NK\_Placenta |
| VEGFB    | -0.8321416 | 0.0462922 | CD56\_NK\_Placenta |
| PLAAT3   | -0.8495389 | 0.0298018 | CD56\_NK\_Placenta |
| MZT2A    | -0.8605676 | 0.0462922 | CD56\_NK\_Placenta |
| PET100   | -0.8632184 | 0.0198046 | CD56\_NK\_Placenta |
| PRDX3    | -0.8713686 | 0.0465982 | CD56\_NK\_Placenta |
| TMEM208  | -0.8725069 | 0.0427229 | CD56\_NK\_Placenta |
| SNHG29   | -0.8763500 | 0.0116126 | CD56\_NK\_Placenta |
| RPL22L1  | -0.8815504 | 0.0112465 | CD56\_NK\_Placenta |
| POLR3K   | -0.8893062 | 0.0421648 | CD56\_NK\_Placenta |
| CCL5     | -0.8914240 | 0.0270471 | CD56\_NK\_Placenta |
| LAMTOR2  | -0.9108871 | 0.0355495 | CD56\_NK\_Placenta |
| CRIP1    | -0.9124258 | 0.0134871 | CD56\_NK\_Placenta |
| AP3S1    | -0.9240345 | 0.0447150 | CD56\_NK\_Placenta |
| CCDC47   | -0.9441162 | 0.0462922 | CD56\_NK\_Placenta |
| BECN1    | -0.9649495 | 0.0213196 | CD56\_NK\_Placenta |
| IFI27L2  | -0.9801653 | 0.0166677 | CD56\_NK\_Placenta |
| ANXA2    | -0.9819654 | 0.0180430 | CD56\_NK\_Placenta |
| NAA20    | -0.9841182 | 0.0150905 | CD56\_NK\_Placenta |
| NDUFA8   | -0.9843185 | 0.0342620 | CD56\_NK\_Placenta |
| TSPO     | -0.9930395 | 0.0369224 | CD56\_NK\_Placenta |
| S100A10  | -0.9988604 | 0.0020917 | CD56\_NK\_Placenta |
| TRDC     | -1.0143246 | 0.0148897 | CD56\_NK\_Placenta |
| SPRY2    | -1.0330518 | 0.0427229 | CD56\_NK\_Placenta |
| NUDT21   | -1.0348677 | 0.0110470 | CD56\_NK\_Placenta |
| VIM      | -1.0490377 | 0.0158571 | CD56\_NK\_Placenta |
| MT1X     | -1.0800345 | 0.0483532 | CD56\_NK\_Placenta |
| S100A6   | -1.1104845 | 0.0047189 | CD56\_NK\_Placenta |
| IGFLR1   | -1.1135350 | 0.0020917 | CD56\_NK\_Placenta |
| RANGRF   | -1.1188357 | 0.0295883 | CD56\_NK\_Placenta |
| RGL4     | -1.1590536 | 0.0112465 | CD56\_NK\_Placenta |
| MED8     | -1.1915941 | 0.0132755 | CD56\_NK\_Placenta |
| S100A11  | -1.2159882 | 0.0052464 | CD56\_NK\_Placenta |
| ITM2A    | -1.2198562 | 0.0222309 | CD56\_NK\_Placenta |
| OASL     | -1.2510441 | 0.0342620 | CD56\_NK\_Placenta |
| ANXA5    | -1.2612294 | 0.0416239 | CD56\_NK\_Placenta |
| CHCHD10  | -1.2878153 | 0.0416299 | CD56\_NK\_Placenta |
| FLT3LG   | -1.2907355 | 0.0131153 | CD56\_NK\_Placenta |
| APOBEC3H | -1.2992842 | 0.0298018 | CD56\_NK\_Placenta |
| PLP2     | -1.3138482 | 0.0357740 | CD56\_NK\_Placenta |
| LTB      | -1.3245067 | 0.0140768 | CD56\_NK\_Placenta |
| LGALS3   | -1.4730994 | 0.0475205 | CD56\_NK\_Placenta |
| DNPH1    | -1.5932825 | 0.0154311 | CD56\_NK\_Placenta |
| CDKN2A   | -1.6965427 | 0.0190263 | CD56\_NK\_Placenta |
| CAPG     | -1.7831186 | 0.0185106 | CD56\_NK\_Placenta |
| CSRP1    | -1.8423446 | 0.0203983 | CD56\_NK\_Placenta |
| MT1F     | -1.8836746 | 0.0185942 | CD56\_NK\_Placenta |
| CD82     | -1.9942423 | 0.0387297 | CD56\_NK\_Placenta |
| LGALS1   | -2.0269092 | 0.0000000 | CD56\_NK\_Placenta |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes      |     logFC |      padj | cell\_type               |
|:-----------|----------:|----------:|:-------------------------|
| IL1R2      | 4.8680464 | 0.0158397 | CD14\_Monocyte\_Placenta |
| ADAMTS2    | 4.0793770 | 0.0008257 | CD14\_Monocyte\_Placenta |
| OCLN       | 3.5314120 | 0.0005874 | CD14\_Monocyte\_Placenta |
| VSIG4      | 2.9577605 | 0.0030059 | CD14\_Monocyte\_Placenta |
| TPST1      | 2.8665694 | 0.0091764 | CD14\_Monocyte\_Placenta |
| FMN1       | 2.4453802 | 0.0096860 | CD14\_Monocyte\_Placenta |
| SIGLEC1    | 1.2496750 | 0.0471105 | CD14\_Monocyte\_Placenta |
| AC004551.1 | 1.2436275 | 0.0293784 | CD14\_Monocyte\_Placenta |
| LHFPL2     | 1.1134394 | 0.0100604 | CD14\_Monocyte\_Placenta |
| GRB10      | 1.0566190 | 0.0100604 | CD14\_Monocyte\_Placenta |
| CALCRL     | 1.0265603 | 0.0293784 | CD14\_Monocyte\_Placenta |
| HDAC9      | 0.9431794 | 0.0007341 | CD14\_Monocyte\_Placenta |
| FAR2       | 0.8276572 | 0.0155072 | CD14\_Monocyte\_Placenta |
| SNTB1      | 0.7491801 | 0.0005874 | CD14\_Monocyte\_Placenta |
| EIF2AK2    | 0.7179928 | 0.0293784 | CD14\_Monocyte\_Placenta |
| CSF2RA     | 0.6359028 | 0.0471105 | CD14\_Monocyte\_Placenta |
| FOXN3      | 0.6237381 | 0.0100604 | CD14\_Monocyte\_Placenta |
| GPD2       | 0.5708229 | 0.0293784 | CD14\_Monocyte\_Placenta |
| ZC3HAV1    | 0.5216000 | 0.0293784 | CD14\_Monocyte\_Placenta |

\[1\] “Genes downregulated in disease”

| genes      |      logFC |      padj | cell\_type               |
|:-----------|-----------:|----------:|:-------------------------|
| ASAH1      | -0.5184657 | 0.0489031 | CD14\_Monocyte\_Placenta |
| RAB32      | -0.5573265 | 0.0268922 | CD14\_Monocyte\_Placenta |
| RAB4A      | -0.5635786 | 0.0327692 | CD14\_Monocyte\_Placenta |
| AC005280.2 | -0.6184572 | 0.0100604 | CD14\_Monocyte\_Placenta |
| AP1S2      | -0.6251650 | 0.0222559 | CD14\_Monocyte\_Placenta |
| FAM13A     | -0.6348051 | 0.0293784 | CD14\_Monocyte\_Placenta |
| SLC46A2    | -0.6585324 | 0.0222559 | CD14\_Monocyte\_Placenta |
| LTB        | -0.7996356 | 0.0293784 | CD14\_Monocyte\_Placenta |
| VNN2       | -0.9462968 | 0.0293784 | CD14\_Monocyte\_Placenta |
| TNFRSF10C  | -0.9764314 | 0.0489031 | CD14\_Monocyte\_Placenta |
| CCR5AS     | -0.9897805 | 0.0471105 | CD14\_Monocyte\_Placenta |
| CLEC6A     | -1.0991582 | 0.0461639 | CD14\_Monocyte\_Placenta |
| PLA2G7     | -1.1704736 | 0.0036368 | CD14\_Monocyte\_Placenta |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes      |    logFC |      padj | cell\_type               |
|:-----------|---------:|----------:|:-------------------------|
| HBG2       | 3.597640 | 0.0063701 | Fetal-Monocyte\_Placenta |
| CXCL8      | 2.397114 | 0.0479689 | Fetal-Monocyte\_Placenta |
| IFITM3     | 2.071215 | 0.0022035 | Fetal-Monocyte\_Placenta |
| ABCA1      | 2.066536 | 0.0087679 | Fetal-Monocyte\_Placenta |
| FP671120.4 | 1.783184 | 0.0119542 | Fetal-Monocyte\_Placenta |
| HMOX1      | 1.719291 | 0.0174258 | Fetal-Monocyte\_Placenta |
| SASH1      | 1.372600 | 0.0022035 | Fetal-Monocyte\_Placenta |

\[1\] “Genes downregulated in disease”

| genes | logFC | padj | cell\_type |
|:------|------:|-----:|:-----------|

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes   |     logFC |      padj | cell\_type               |
|:--------|----------:|----------:|:-------------------------|
| ADAMTS2 | 2.6124728 | 0.0258235 | CD16\_Monocyte\_Placenta |
| TPST1   | 2.5481374 | 0.0051888 | CD16\_Monocyte\_Placenta |
| ZBP1    | 2.2921182 | 0.0006567 | CD16\_Monocyte\_Placenta |
| FMN1    | 1.8814912 | 0.0253495 | CD16\_Monocyte\_Placenta |
| IFIT1   | 1.8136451 | 0.0076019 | CD16\_Monocyte\_Placenta |
| USP18   | 1.7557727 | 0.0224437 | CD16\_Monocyte\_Placenta |
| RSAD2   | 1.7348765 | 0.0034194 | CD16\_Monocyte\_Placenta |
| VSIG4   | 1.6130441 | 0.0133539 | CD16\_Monocyte\_Placenta |
| IFIT3   | 1.5132419 | 0.0034194 | CD16\_Monocyte\_Placenta |
| MX1     | 1.3618676 | 0.0012988 | CD16\_Monocyte\_Placenta |
| CMPK2   | 1.3529124 | 0.0013350 | CD16\_Monocyte\_Placenta |
| IFI44L  | 1.2683862 | 0.0102127 | CD16\_Monocyte\_Placenta |
| ERAP2   | 1.2653031 | 0.0004119 | CD16\_Monocyte\_Placenta |
| MS4A4E  | 1.2079315 | 0.0437586 | CD16\_Monocyte\_Placenta |
| PRKCE   | 1.1870795 | 0.0181316 | CD16\_Monocyte\_Placenta |
| OAS3    | 1.1720105 | 0.0164518 | CD16\_Monocyte\_Placenta |
| TCN2    | 1.1708464 | 0.0110542 | CD16\_Monocyte\_Placenta |
| IFI6    | 1.1685575 | 0.0004119 | CD16\_Monocyte\_Placenta |
| MX2     | 1.1664998 | 0.0171523 | CD16\_Monocyte\_Placenta |
| OASL    | 1.1392773 | 0.0102887 | CD16\_Monocyte\_Placenta |
| IFI44   | 1.1283917 | 0.0137599 | CD16\_Monocyte\_Placenta |
| C2      | 1.1037486 | 0.0018777 | CD16\_Monocyte\_Placenta |
| FMNL2   | 1.0889307 | 0.0012988 | CD16\_Monocyte\_Placenta |
| DDX60   | 1.0861329 | 0.0034194 | CD16\_Monocyte\_Placenta |
| DDX58   | 1.0732494 | 0.0296869 | CD16\_Monocyte\_Placenta |
| DDX60L  | 1.0685076 | 0.0027912 | CD16\_Monocyte\_Placenta |
| ISG20   | 1.0641921 | 0.0391551 | CD16\_Monocyte\_Placenta |
| HELZ2   | 1.0640746 | 0.0077261 | CD16\_Monocyte\_Placenta |
| SNTB1   | 1.0493280 | 0.0040132 | CD16\_Monocyte\_Placenta |
| SAMD4A  | 1.0490521 | 0.0428538 | CD16\_Monocyte\_Placenta |
| SAMD9L  | 1.0434666 | 0.0022017 | CD16\_Monocyte\_Placenta |
| CPEB3   | 1.0355434 | 0.0052931 | CD16\_Monocyte\_Placenta |
| OAS2    | 1.0035250 | 0.0053126 | CD16\_Monocyte\_Placenta |
| CCND3   | 0.9966044 | 0.0317702 | CD16\_Monocyte\_Placenta |
| ISG15   | 0.9949836 | 0.0284895 | CD16\_Monocyte\_Placenta |
| ARHGEF3 | 0.9928244 | 0.0027912 | CD16\_Monocyte\_Placenta |
| PPM1L   | 0.9339283 | 0.0430940 | CD16\_Monocyte\_Placenta |
| ZC3HAV1 | 0.9318726 | 0.0429442 | CD16\_Monocyte\_Placenta |
| EIF2AK2 | 0.9292534 | 0.0102127 | CD16\_Monocyte\_Placenta |
| STAT2   | 0.8878776 | 0.0052931 | CD16\_Monocyte\_Placenta |
| XAF1    | 0.8686827 | 0.0122661 | CD16\_Monocyte\_Placenta |
| UTRN    | 0.8556445 | 0.0133539 | CD16\_Monocyte\_Placenta |
| SP140   | 0.8465110 | 0.0171523 | CD16\_Monocyte\_Placenta |
| GPD2    | 0.8285410 | 0.0429442 | CD16\_Monocyte\_Placenta |
| DIP2B   | 0.8265429 | 0.0409925 | CD16\_Monocyte\_Placenta |
| PTK2B   | 0.7816738 | 0.0188289 | CD16\_Monocyte\_Placenta |
| NLRC5   | 0.7740852 | 0.0328437 | CD16\_Monocyte\_Placenta |
| APOL6   | 0.7711525 | 0.0273339 | CD16\_Monocyte\_Placenta |
| KIF13A  | 0.7677355 | 0.0272017 | CD16\_Monocyte\_Placenta |
| ITGAM   | 0.7546734 | 0.0184810 | CD16\_Monocyte\_Placenta |
| PARP9   | 0.7530230 | 0.0066166 | CD16\_Monocyte\_Placenta |
| LY6E    | 0.7378671 | 0.0034194 | CD16\_Monocyte\_Placenta |
| PAG1    | 0.7183676 | 0.0051888 | CD16\_Monocyte\_Placenta |
| PML     | 0.7130592 | 0.0448390 | CD16\_Monocyte\_Placenta |
| SPIDR   | 0.7022549 | 0.0391551 | CD16\_Monocyte\_Placenta |
| XRN1    | 0.6879290 | 0.0317702 | CD16\_Monocyte\_Placenta |
| ESYT2   | 0.6651955 | 0.0158755 | CD16\_Monocyte\_Placenta |
| DTX3L   | 0.6439028 | 0.0428538 | CD16\_Monocyte\_Placenta |
| GNGT2   | 0.6378296 | 0.0296869 | CD16\_Monocyte\_Placenta |
| WDFY2   | 0.6347974 | 0.0476533 | CD16\_Monocyte\_Placenta |
| RHBDD1  | 0.6101505 | 0.0448390 | CD16\_Monocyte\_Placenta |
| LHFPL2  | 0.6084938 | 0.0206097 | CD16\_Monocyte\_Placenta |
| ADAM17  | 0.6013916 | 0.0348571 | CD16\_Monocyte\_Placenta |
| HECA    | 0.5624862 | 0.0051888 | CD16\_Monocyte\_Placenta |
| SLC38A6 | 0.5414111 | 0.0373482 | CD16\_Monocyte\_Placenta |
| NSMCE2  | 0.5291229 | 0.0429442 | CD16\_Monocyte\_Placenta |

\[1\] “Genes downregulated in disease”

| genes     |      logFC |      padj | cell\_type               |
|:----------|-----------:|----------:|:-------------------------|
| CREBL2    | -0.5005128 | 0.0430940 | CD16\_Monocyte\_Placenta |
| TNFAIP8L1 | -0.7599193 | 0.0325494 | CD16\_Monocyte\_Placenta |
| PPA1      | -0.8480139 | 0.0051888 | CD16\_Monocyte\_Placenta |
| PLA2G7    | -0.9002636 | 0.0157169 | CD16\_Monocyte\_Placenta |
| TTC39C    | -0.9546111 | 0.0253495 | CD16\_Monocyte\_Placenta |
| C15orf48  | -1.3090721 | 0.0171523 | CD16\_Monocyte\_Placenta |
| MT1X      | -1.3551961 | 0.0102127 | CD16\_Monocyte\_Placenta |
| TIFAB     | -1.6193322 | 0.0284895 | CD16\_Monocyte\_Placenta |
| RASA4B    | -1.6593631 | 0.0315242 | CD16\_Monocyte\_Placenta |
| EBI3      | -3.0262467 | 0.0039163 | CD16\_Monocyte\_Placenta |
| MT1G      | -3.9274379 | 0.0071234 | CD16\_Monocyte\_Placenta |
| FABP4     | -4.5047346 | 0.0253495 | CD16\_Monocyte\_Placenta |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes   |     logFC |      padj | cell\_type                      |
|:--------|----------:|----------:|:--------------------------------|
| CD163   | 3.8045005 | 0.0060407 | Nonclassical-monocyte\_Placenta |
| VSIG4   | 3.5302316 | 0.0189940 | Nonclassical-monocyte\_Placenta |
| FMN1    | 3.3194366 | 0.0189940 | Nonclassical-monocyte\_Placenta |
| DOCK4   | 1.2793849 | 0.0189940 | Nonclassical-monocyte\_Placenta |
| IFI6    | 1.1943084 | 0.0006966 | Nonclassical-monocyte\_Placenta |
| ERAP2   | 1.1679744 | 0.0189940 | Nonclassical-monocyte\_Placenta |
| TBC1D2  | 1.0585530 | 0.0189940 | Nonclassical-monocyte\_Placenta |
| RBM47   | 0.8461219 | 0.0189940 | Nonclassical-monocyte\_Placenta |
| ARHGEF3 | 0.8141359 | 0.0202278 | Nonclassical-monocyte\_Placenta |

\[1\] “Genes downregulated in disease”

| genes  |      logFC |      padj | cell\_type                      |
|:-------|-----------:|----------:|:--------------------------------|
| TESC   | -0.8575627 | 0.0270781 | Nonclassical-monocyte\_Placenta |
| NEURL1 | -1.4176469 | 0.0192220 | Nonclassical-monocyte\_Placenta |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes  |    logFC |      padj | cell\_type                        |
|:-------|---------:|----------:|:----------------------------------|
| PDGFC  | 1.848506 | 0.0480196 | Intermediate-macrophage\_Placenta |
| TCN2   | 1.761327 | 0.0108433 | Intermediate-macrophage\_Placenta |
| OAS3   | 1.494528 | 0.0480196 | Intermediate-macrophage\_Placenta |
| CMPK2  | 1.458390 | 0.0480196 | Intermediate-macrophage\_Placenta |
| IFIT3  | 1.436581 | 0.0480196 | Intermediate-macrophage\_Placenta |
| IFI44L | 1.371236 | 0.0480196 | Intermediate-macrophage\_Placenta |
| SPRED1 | 1.259161 | 0.0480196 | Intermediate-macrophage\_Placenta |
| SNTB1  | 1.176942 | 0.0480196 | Intermediate-macrophage\_Placenta |

\[1\] “Genes downregulated in disease”

| genes     |     logFC |      padj | cell\_type                        |
|:----------|----------:|----------:|:----------------------------------|
| CST6      | -1.164354 | 0.0480196 | Intermediate-macrophage\_Placenta |
| LINC00937 | -1.216497 | 0.0480196 | Intermediate-macrophage\_Placenta |
| DUSP2     | -2.072237 | 0.0017863 | Intermediate-macrophage\_Placenta |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes   |     logFC |      padj | cell\_type           |
|:--------|----------:|----------:|:---------------------|
| GSDMA   | 3.4711845 | 0.0000202 | Macrophage\_Placenta |
| ZNF385D | 2.4860007 | 0.0010785 | Macrophage\_Placenta |
| TPST1   | 2.2848901 | 0.0072908 | Macrophage\_Placenta |
| DPYSL3  | 1.9622885 | 0.0361572 | Macrophage\_Placenta |
| TBC1D16 | 1.7157658 | 0.0338316 | Macrophage\_Placenta |
| P2RY1   | 1.5018669 | 0.0042004 | Macrophage\_Placenta |
| SNTB1   | 1.4421326 | 0.0003287 | Macrophage\_Placenta |
| MCEMP1  | 1.2101676 | 0.0349125 | Macrophage\_Placenta |
| ERAP2   | 1.1487414 | 0.0081256 | Macrophage\_Placenta |
| SLC11A1 | 1.1015205 | 0.0479365 | Macrophage\_Placenta |
| PPARG   | 1.0736162 | 0.0081256 | Macrophage\_Placenta |
| RMDN3   | 1.0328220 | 0.0000001 | Macrophage\_Placenta |
| ELL2    | 0.9821808 | 0.0023178 | Macrophage\_Placenta |
| PLIN2   | 0.9629528 | 0.0322435 | Macrophage\_Placenta |
| FAM20C  | 0.9597376 | 0.0043354 | Macrophage\_Placenta |
| IFI6    | 0.9445222 | 0.0107840 | Macrophage\_Placenta |
| TDP2    | 0.8631547 | 0.0396070 | Macrophage\_Placenta |
| STK39   | 0.8627038 | 0.0442220 | Macrophage\_Placenta |
| ZCCHC2  | 0.8326033 | 0.0222021 | Macrophage\_Placenta |
| CBLB    | 0.8032932 | 0.0321228 | Macrophage\_Placenta |
| FMNL2   | 0.8003592 | 0.0453378 | Macrophage\_Placenta |
| IRAK1   | 0.7977131 | 0.0297363 | Macrophage\_Placenta |
| HSD3B7  | 0.7615371 | 0.0214890 | Macrophage\_Placenta |
| APOBR   | 0.7612176 | 0.0072908 | Macrophage\_Placenta |
| CLIP4   | 0.7577123 | 0.0062074 | Macrophage\_Placenta |
| EIF2AK2 | 0.7399629 | 0.0297363 | Macrophage\_Placenta |
| KCTD13  | 0.7183159 | 0.0482362 | Macrophage\_Placenta |
| TXNDC11 | 0.6627755 | 0.0222021 | Macrophage\_Placenta |
| OLR1    | 0.6212141 | 0.0222021 | Macrophage\_Placenta |
| HECA    | 0.5662908 | 0.0450916 | Macrophage\_Placenta |
| UBE3A   | 0.5382939 | 0.0396070 | Macrophage\_Placenta |
| DOP1B   | 0.5104190 | 0.0249599 | Macrophage\_Placenta |

\[1\] “Genes downregulated in disease”

| genes      |      logFC |      padj | cell\_type           |
|:-----------|-----------:|----------:|:---------------------|
| GCA        | -0.6523135 | 0.0222021 | Macrophage\_Placenta |
| CREBL2     | -0.7075917 | 0.0297363 | Macrophage\_Placenta |
| CD14       | -0.7377568 | 0.0110279 | Macrophage\_Placenta |
| SRGAP3     | -0.7846689 | 0.0179689 | Macrophage\_Placenta |
| CD48       | -0.8161427 | 0.0361572 | Macrophage\_Placenta |
| PNRC1      | -0.8163011 | 0.0297363 | Macrophage\_Placenta |
| GIMAP7     | -0.8199462 | 0.0396070 | Macrophage\_Placenta |
| CLEC4A     | -0.8789758 | 0.0475418 | Macrophage\_Placenta |
| RNF113A    | -0.8875025 | 0.0110279 | Macrophage\_Placenta |
| JAML       | -0.9471128 | 0.0004534 | Macrophage\_Placenta |
| LYZ        | -1.0499147 | 0.0222021 | Macrophage\_Placenta |
| AC005280.2 | -1.0996702 | 0.0044392 | Macrophage\_Placenta |
| FOLR2      | -1.1635456 | 0.0004180 | Macrophage\_Placenta |
| RASA4      | -1.2368351 | 0.0253833 | Macrophage\_Placenta |
| PLD4       | -1.2908554 | 0.0123913 | Macrophage\_Placenta |
| PPA1       | -1.3062162 | 0.0322435 | Macrophage\_Placenta |
| ITGB2-AS1  | -1.3506827 | 0.0240881 | Macrophage\_Placenta |
| METTL7B    | -1.4764937 | 0.0297363 | Macrophage\_Placenta |
| TSPAN33    | -1.6193214 | 0.0109951 | Macrophage\_Placenta |
| DUSP5      | -1.6320470 | 0.0453378 | Macrophage\_Placenta |
| RARRES1    | -1.9015132 | 0.0114220 | Macrophage\_Placenta |
| MT1E       | -1.9590500 | 0.0218309 | Macrophage\_Placenta |
| CASP5      | -1.9800657 | 0.0453378 | Macrophage\_Placenta |
| CLEC4D     | -2.0663308 | 0.0451818 | Macrophage\_Placenta |
| TMEM176B   | -2.1215385 | 0.0042976 | Macrophage\_Placenta |
| ENPP2      | -2.2211863 | 0.0245782 | Macrophage\_Placenta |
| IGF1       | -2.3856036 | 0.0042976 | Macrophage\_Placenta |
| NOTCH2NLB  | -2.7874445 | 0.0222021 | Macrophage\_Placenta |
| CXCL12     | -3.3060058 | 0.0000392 | Macrophage\_Placenta |
| HSD11B1    | -3.7595723 | 0.0361572 | Macrophage\_Placenta |
| TNFRSF4    | -3.8485926 | 0.0240074 | Macrophage\_Placenta |
| PLA2G2D    | -3.8968861 | 0.0003287 | Macrophage\_Placenta |
| UBD        | -5.5377777 | 0.0088289 | Macrophage\_Placenta |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes  |     logFC |      padj | cell\_type    |
|:-------|----------:|----------:|:--------------|
| EFHD2  | 2.1018278 | 0.0000424 | CTB\_Placenta |
| IGF2   | 1.8561850 | 0.0002200 | CTB\_Placenta |
| S100P  | 1.5328196 | 0.0002200 | CTB\_Placenta |
| KRT19  | 1.3739805 | 0.0020960 | CTB\_Placenta |
| GPX3   | 1.3528219 | 0.0099302 | CTB\_Placenta |
| SH3GL1 | 1.2434043 | 0.0003671 | CTB\_Placenta |
| NEK6   | 0.9995955 | 0.0107673 | CTB\_Placenta |
| HMGB3  | 0.9404534 | 0.0478951 | CTB\_Placenta |
| TCEAL4 | 0.9171384 | 0.0063767 | CTB\_Placenta |
| RPS8   | 0.6758377 | 0.0382546 | CTB\_Placenta |
| RPS6   | 0.5130453 | 0.0382546 | CTB\_Placenta |

\[1\] “Genes downregulated in disease”

| genes  |      logFC |      padj | cell\_type    |
|:-------|-----------:|----------:|:--------------|
| PGD    | -0.9795341 | 0.0099302 | CTB\_Placenta |
| NFE2L3 | -1.3110479 | 0.0284361 | CTB\_Placenta |

\[1\] “Results for celltype:” \[1\] “Genes upregulated in disease”

| genes      |     logFC |      padj | cell\_type    |
|:-----------|----------:|----------:|:--------------|
| IGFBP3     | 6.1490467 | 0.0033461 | STB\_Placenta |
| LTF        | 6.0120211 | 0.0072238 | STB\_Placenta |
| FSTL3      | 5.6689714 | 0.0000002 | STB\_Placenta |
| XACT       | 5.6028197 | 0.0008932 | STB\_Placenta |
| SERPINA3   | 5.3947439 | 0.0305344 | STB\_Placenta |
| ADAM2      | 5.1024483 | 0.0004246 | STB\_Placenta |
| LEP        | 4.7129557 | 0.0000002 | STB\_Placenta |
| HTRA4      | 3.9248596 | 0.0000000 | STB\_Placenta |
| TMEM45A    | 3.7129572 | 0.0133145 | STB\_Placenta |
| MUC1       | 3.6661006 | 0.0176728 | STB\_Placenta |
| ARMS2      | 3.6571989 | 0.0000209 | STB\_Placenta |
| AC087857.1 | 3.6005977 | 0.0162097 | STB\_Placenta |
| AC026167.1 | 3.5706205 | 0.0173924 | STB\_Placenta |
| MIR210HG   | 3.5703566 | 0.0000090 | STB\_Placenta |
| FLT1       | 3.3560253 | 0.0000017 | STB\_Placenta |
| AC004704.1 | 3.0357807 | 0.0062270 | STB\_Placenta |
| QPCT       | 2.9984822 | 0.0173924 | STB\_Placenta |
| PLOD2      | 2.8769227 | 0.0160767 | STB\_Placenta |
| IGF2       | 2.7762790 | 0.0009686 | STB\_Placenta |
| AL513164.1 | 2.7553913 | 0.0172993 | STB\_Placenta |
| PNCK       | 2.7334979 | 0.0172179 | STB\_Placenta |
| PAPSS2     | 2.6758205 | 0.0436316 | STB\_Placenta |
| FN1        | 2.6091173 | 0.0499650 | STB\_Placenta |
| CAMK2N1    | 2.3767860 | 0.0446402 | STB\_Placenta |
| IFT27      | 2.3275786 | 0.0009805 | STB\_Placenta |
| LPL        | 2.2707431 | 0.0031327 | STB\_Placenta |
| CRH        | 2.2574833 | 0.0018825 | STB\_Placenta |
| H19        | 2.2417164 | 0.0119237 | STB\_Placenta |
| TMEM91     | 2.1457179 | 0.0156129 | STB\_Placenta |
| EGLN3      | 2.1154698 | 0.0178619 | STB\_Placenta |
| P4HA1      | 2.0998810 | 0.0264222 | STB\_Placenta |
| AIF1L      | 2.0832377 | 0.0005384 | STB\_Placenta |
| LIMS2      | 1.9936215 | 0.0062270 | STB\_Placenta |
| CST6       | 1.9161430 | 0.0004246 | STB\_Placenta |
| MAOB       | 1.9078398 | 0.0207673 | STB\_Placenta |
| HTRA1      | 1.8872400 | 0.0049624 | STB\_Placenta |
| TREM1      | 1.8415146 | 0.0027634 | STB\_Placenta |
| SH3BP5     | 1.8152533 | 0.0091563 | STB\_Placenta |
| BIN2       | 1.7893715 | 0.0131028 | STB\_Placenta |
| PHYHIP     | 1.7849821 | 0.0162097 | STB\_Placenta |
| DERL3      | 1.7355201 | 0.0099387 | STB\_Placenta |
| PPP1R1C    | 1.7275751 | 0.0240954 | STB\_Placenta |
| GUCA2A     | 1.6965686 | 0.0426190 | STB\_Placenta |
| INHA       | 1.6849473 | 0.0000616 | STB\_Placenta |
| SPAG4      | 1.6694511 | 0.0062270 | STB\_Placenta |
| PPIA       | 1.6675856 | 0.0128942 | STB\_Placenta |
| KIAA0753   | 1.6472109 | 0.0393055 | STB\_Placenta |
| KRT86      | 1.6345817 | 0.0097080 | STB\_Placenta |
| GALK1      | 1.6283281 | 0.0156129 | STB\_Placenta |
| LY6D       | 1.6044800 | 0.0012642 | STB\_Placenta |
| GDPD3      | 1.5939173 | 0.0006160 | STB\_Placenta |
| LMCD1      | 1.5641022 | 0.0436316 | STB\_Placenta |
| SLC6A8     | 1.5450017 | 0.0492948 | STB\_Placenta |
| AP003307.1 | 1.5057586 | 0.0365330 | STB\_Placenta |
| PPL        | 1.5049631 | 0.0212323 | STB\_Placenta |
| SMIM3      | 1.4818045 | 0.0158382 | STB\_Placenta |
| PAPPA2     | 1.4718337 | 0.0267598 | STB\_Placenta |
| FAAP20     | 1.4560437 | 0.0000474 | STB\_Placenta |
| CCDC183    | 1.3707946 | 0.0158382 | STB\_Placenta |
| FLNB       | 1.3664608 | 0.0413214 | STB\_Placenta |
| GCHFR      | 1.3630111 | 0.0021624 | STB\_Placenta |
| KRT81      | 1.3555951 | 0.0480741 | STB\_Placenta |
| AL731684.1 | 1.3364742 | 0.0122217 | STB\_Placenta |
| MIR193BHG  | 1.3302199 | 0.0140855 | STB\_Placenta |
| AC116424.1 | 1.3143786 | 0.0212323 | STB\_Placenta |
| RPL13      | 1.2896098 | 0.0221116 | STB\_Placenta |
| CDCA4      | 1.2623685 | 0.0119237 | STB\_Placenta |
| BCL3       | 1.1954395 | 0.0377889 | STB\_Placenta |
| AC108690.1 | 1.1909145 | 0.0358092 | STB\_Placenta |
| HMGN3      | 1.1813909 | 0.0131028 | STB\_Placenta |
| TGFB1      | 1.1778869 | 0.0264222 | STB\_Placenta |
| RPS3A      | 1.1454463 | 0.0172179 | STB\_Placenta |
| ARMCX6     | 1.1331807 | 0.0446402 | STB\_Placenta |
| RAB34      | 1.1148141 | 0.0174911 | STB\_Placenta |
| SH3PXD2A   | 1.1141128 | 0.0172179 | STB\_Placenta |
| RPL36A     | 1.0953180 | 0.0273170 | STB\_Placenta |
| TPT1       | 1.0776486 | 0.0358092 | STB\_Placenta |
| RPL37      | 1.0574104 | 0.0492948 | STB\_Placenta |
| ZFAS1      | 1.0052215 | 0.0128942 | STB\_Placenta |
| RPL26      | 1.0017275 | 0.0273014 | STB\_Placenta |
| AP3S1      | 0.9803850 | 0.0069761 | STB\_Placenta |
| PTMS       | 0.9618632 | 0.0473729 | STB\_Placenta |
| RPL5       | 0.9523012 | 0.0212323 | STB\_Placenta |
| RPS27L     | 0.9495825 | 0.0470829 | STB\_Placenta |
| OST4       | 0.9483017 | 0.0399602 | STB\_Placenta |
| BSCL2      | 0.9305392 | 0.0172179 | STB\_Placenta |
| RPS25      | 0.9283933 | 0.0413116 | STB\_Placenta |
| FKBP8      | 0.8732989 | 0.0366201 | STB\_Placenta |
| TUSC2      | 0.8659268 | 0.0182357 | STB\_Placenta |
| CCDC115    | 0.8534488 | 0.0349281 | STB\_Placenta |
| SMAGP      | 0.8531216 | 0.0497983 | STB\_Placenta |
| PLPP5      | 0.8438830 | 0.0129103 | STB\_Placenta |
| FIBP       | 0.8418058 | 0.0281595 | STB\_Placenta |
| C11orf24   | 0.7958018 | 0.0488765 | STB\_Placenta |
| TOMM20     | 0.7827701 | 0.0497983 | STB\_Placenta |
| TOMM6      | 0.7431493 | 0.0447576 | STB\_Placenta |
| SNHG16     | 0.7412182 | 0.0373535 | STB\_Placenta |
| SORBS3     | 0.6518493 | 0.0485169 | STB\_Placenta |

\[1\] “Genes downregulated in disease”

| genes      |      logFC |      padj | cell\_type    |
|:-----------|-----------:|----------:|:--------------|
| SNX27      | -0.6502717 | 0.0061811 | STB\_Placenta |
| RABGAP1L   | -0.6887417 | 0.0155816 | STB\_Placenta |
| GCLC       | -0.7074628 | 0.0427827 | STB\_Placenta |
| CASP4      | -0.7491741 | 0.0099387 | STB\_Placenta |
| DYRK2      | -0.7713425 | 0.0365330 | STB\_Placenta |
| KCTD3      | -0.7747330 | 0.0158382 | STB\_Placenta |
| TRIM33     | -0.7978145 | 0.0338116 | STB\_Placenta |
| ABCD3      | -0.8100797 | 0.0189633 | STB\_Placenta |
| CPEB4      | -0.8484995 | 0.0410946 | STB\_Placenta |
| MORC4      | -0.8512844 | 0.0273170 | STB\_Placenta |
| CDYL2      | -0.8621808 | 0.0033461 | STB\_Placenta |
| PXK        | -0.8739235 | 0.0212323 | STB\_Placenta |
| PARVA      | -0.9190432 | 0.0276850 | STB\_Placenta |
| WWC1       | -0.9288144 | 0.0160767 | STB\_Placenta |
| EVL        | -0.9362932 | 0.0446402 | STB\_Placenta |
| ABAT       | -0.9428266 | 0.0011759 | STB\_Placenta |
| ALDH4A1    | -0.9459933 | 0.0151917 | STB\_Placenta |
| MYO18A     | -0.9483422 | 0.0240954 | STB\_Placenta |
| CDKN2B     | -0.9576318 | 0.0212323 | STB\_Placenta |
| ELMO1      | -0.9931417 | 0.0410946 | STB\_Placenta |
| PPP2R3A    | -1.0040760 | 0.0264222 | STB\_Placenta |
| PIK3C2B    | -1.0224506 | 0.0473729 | STB\_Placenta |
| GAB1       | -1.0320690 | 0.0475512 | STB\_Placenta |
| NCF4-AS1   | -1.0773481 | 0.0106455 | STB\_Placenta |
| AMPD3      | -1.1046380 | 0.0264222 | STB\_Placenta |
| ABCG2      | -1.1887930 | 0.0062270 | STB\_Placenta |
| STON2      | -1.2009845 | 0.0358092 | STB\_Placenta |
| NAV1       | -1.2265740 | 0.0410946 | STB\_Placenta |
| SYNPO2L    | -1.2850862 | 0.0074533 | STB\_Placenta |
| MCPH1-AS1  | -1.3009974 | 0.0128942 | STB\_Placenta |
| STON1      | -1.3141410 | 0.0042348 | STB\_Placenta |
| TBC1D9     | -1.3260890 | 0.0156129 | STB\_Placenta |
| TRIM5      | -1.3588146 | 0.0008932 | STB\_Placenta |
| RAB27B     | -1.3645712 | 0.0238518 | STB\_Placenta |
| PLCXD2     | -1.3729524 | 0.0358092 | STB\_Placenta |
| IMMP2L     | -1.3801183 | 0.0358092 | STB\_Placenta |
| ADCY7      | -1.4156148 | 0.0252270 | STB\_Placenta |
| ASB2       | -1.4288526 | 0.0172179 | STB\_Placenta |
| AGPAT5     | -1.4333054 | 0.0436316 | STB\_Placenta |
| MYO1B      | -1.4376650 | 0.0113436 | STB\_Placenta |
| NEDD4L     | -1.4759799 | 0.0003631 | STB\_Placenta |
| WNT7A      | -1.5359869 | 0.0495028 | STB\_Placenta |
| TCHHL1     | -1.6188871 | 0.0162220 | STB\_Placenta |
| PIP5K1C    | -1.6544057 | 0.0421086 | STB\_Placenta |
| RADX       | -1.6795549 | 0.0193401 | STB\_Placenta |
| LRIG3      | -1.6839292 | 0.0308936 | STB\_Placenta |
| HEG1       | -1.8008218 | 0.0320462 | STB\_Placenta |
| LINC01484  | -1.8036813 | 0.0160767 | STB\_Placenta |
| KIAA1211   | -1.8685877 | 0.0012642 | STB\_Placenta |
| P3H2       | -1.9708199 | 0.0413214 | STB\_Placenta |
| NFIA-AS2   | -1.9993554 | 0.0172993 | STB\_Placenta |
| AR         | -2.1729805 | 0.0251433 | STB\_Placenta |
| ADAMTSL1   | -2.2325845 | 0.0047012 | STB\_Placenta |
| AC099792.1 | -2.4781760 | 0.0172179 | STB\_Placenta |
| CPQ        | -2.4789859 | 0.0176728 | STB\_Placenta |
| AC104248.1 | -2.7261044 | 0.0160767 | STB\_Placenta |
| AC079760.1 | -2.9415404 | 0.0339449 | STB\_Placenta |
| APOC3      | -3.3814734 | 0.0497983 | STB\_Placenta |

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
