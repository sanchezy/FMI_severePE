Spatial transcriptomics: SC-CTB and EVT are spatially patterned
================

Author: Yara E. Sanchez Corrales

# Introduction

SC-CTB is an trophoblast population in the CAM. This population was
recently described in [Marsh et al
2022](https://elifesciences.org/articles/78829). Here I calculate the
pearson coeffient per donor and make a plot. The plots of two cell types
are done in python cell2location. <br>

Note: The deconvolution for donor “LR” is not great.

> SC-CTB and EVT are spatially anti-correlated

<br>

# Load libraries and object

``` r
dirsave <- "~/Projects/FMI-all-Spatial-20240312/plots/20250110_Tropho_decidua_parietalis/"
```

``` r
# read the object
spatial.1 <- readRDS(file = "~/Projects/FMI-all-Spatial-20240312/objects/FMI-CAM-Spatial-20241120.rds")
# create a new metadata column
spatial.1@meta.data$donor2 <- paste0("F", str_sub(spatial.1@meta.data$donor, start= -2))
# create a composite column with the condition
spatial.1@meta.data$donor.Condition2 <- paste0(spatial.1@meta.data$Condition,"_",spatial.1@meta.data$donor2)

spatial.1
```

    ## An object of class Seurat 
    ## 36601 features across 47771 samples within 1 assay 
    ## Active assay: Spatial (36601 features, 0 variable features)
    ##  19 images present: FVQ.CAM, FLR.CAM, FVB.CAM, FCM.CAM, FAM.CAM, FGS.CAM, FEP.CAM, FVS.CAM, FLJ.CAM, FJD.CAM, FVT.CAM, FYI.CAM, FNS.CAM, FRK.CAM, FFM.CAM, FRH.CAM, FSH.CAM, FND.CAM, FMM.CAM

``` r
# subset by the chorion region
spatial.1 <- subset(spatial.1, subset = regions %in% c("chorion"))
```

``` r
# Subset by donor, filter the SC.CTB and EVT abundance by 0.5 

df <- data.frame()
df_temp <- data.frame()

# loop through donors

for (i in 1:length(donors_to_plot)) {

# Subset the data for the current donor
temp <- subset(spatial.1, subset = donor2 %in% donors_to_plot[i])
# filter by expression 
temp <- subset(temp, subset = c2loc_SC.CTB > 0.5)
# temp <- subset(temp, subset = c2loc_EVT > 0.5)

# Check if temp has rows; skip if not
  if (nrow(temp) > 0) {
    # Extract metadata and correlation coefficient
    donor <- unique(temp@meta.data$donor2)
    GA_Condition <- unique(temp@meta.data$GA_Condition)
    pearson <- cor(temp$c2loc_SC.CTB, temp$c2loc_EVT, method = "pearson")
    test <- cor.test(temp$c2loc_SC.CTB, temp$c2loc_EVT, method = "pearson") 
    
    GA_Cat <- unique(temp@meta.data$GA_Cat)
      
    # Create a temporary data frame for this donor
    df_temp <- data.frame(
      donor = donor,
      GA_Condition = GA_Condition,
      pearson = pearson,
      GA_Cat = GA_Cat
    )
    
    # Append to the final data frame
    df <- rbind(df, df_temp)
  } else {
    # If no data for this donor, optionally print a warning
    print(paste("No data for donor:", donors_to_plot[i]))
  }
  
print(donors_to_plot[i])
print(test)

  # Remove temporary variable
  rm(temp)
}
```

    ## [1] "FVQ"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = -7.1619, df = 677, p-value = 2.081e-12
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.3339572 -0.1940182
    ## sample estimates:
    ##        cor 
    ## -0.2653848 
    ## 
    ## [1] "FLR"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = 0.47893, df = 3, p-value = 0.6647
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.8050508  0.9300830
    ## sample estimates:
    ##       cor 
    ## 0.2665082 
    ## 
    ## [1] "FVB"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = -3.6374, df = 502, p-value = 0.0003039
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.24417006 -0.07393913
    ## sample estimates:
    ##        cor 
    ## -0.1602459 
    ## 
    ## [1] "FCM"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = -6.1769, df = 324, p-value = 1.958e-09
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.4184536 -0.2238497
    ## sample estimates:
    ##        cor 
    ## -0.3245823 
    ## 
    ## [1] "FAM"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = -12.506, df = 603, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.5149099 -0.3881533
    ## sample estimates:
    ##        cor 
    ## -0.4538244 
    ## 
    ## [1] "FGS"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = -2.493, df = 828, p-value = 0.01286
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.15346213 -0.01837284
    ## sample estimates:
    ##         cor 
    ## -0.08631422 
    ## 
    ## [1] "FEP"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = -20.925, df = 414, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.7606464 -0.6667451
    ## sample estimates:
    ##       cor 
    ## -0.716932 
    ## 
    ## [1] "FVS"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = -6.8367, df = 166, p-value = 1.476e-10
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.5790429 -0.3415584
    ## sample estimates:
    ##        cor 
    ## -0.4687278 
    ## 
    ## [1] "FLJ"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = -2.298, df = 440, p-value = 0.02203
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.2001400 -0.0157891
    ## sample estimates:
    ##        cor 
    ## -0.1089008 
    ## 
    ## [1] "FJD"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = -13.462, df = 458, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.5949108 -0.4635837
    ## sample estimates:
    ##        cor 
    ## -0.5324438 
    ## 
    ## [1] "FVT"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = -5.2776, df = 578, p-value = 1.855e-07
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.2907513 -0.1353620
    ## sample estimates:
    ##        cor 
    ## -0.2144129 
    ## 
    ## [1] "FYI"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = -29.099, df = 666, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.7797624 -0.7127522
    ## sample estimates:
    ##        cor 
    ## -0.7481588 
    ## 
    ## [1] "FRK"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = 1.4348, df = 438, p-value = 0.1521
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.0252478  0.1608532
    ## sample estimates:
    ##        cor 
    ## 0.06839767 
    ## 
    ## [1] "FNS"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = 3.2868, df = 407, p-value = 0.001101
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.06484649 0.25380944
    ## sample estimates:
    ##       cor 
    ## 0.1608011 
    ## 
    ## [1] "FFM"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = -15.553, df = 645, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.5761345 -0.4638380
    ## sample estimates:
    ##        cor 
    ## -0.5222465 
    ## 
    ## [1] "FRH"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = -5.2539, df = 366, p-value = 2.533e-07
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.3573782 -0.1671163
    ## sample estimates:
    ##        cor 
    ## -0.2648227 
    ## 
    ## [1] "FND"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = -4.205, df = 338, p-value = 3.347e-05
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.3216987 -0.1194358
    ## sample estimates:
    ##        cor 
    ## -0.2229656 
    ## 
    ## [1] "FSH"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = -16.754, df = 806, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.5576662 -0.4552320
    ## sample estimates:
    ##        cor 
    ## -0.5082445 
    ## 
    ## [1] "FMM"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = 5.7666, df = 619, p-value = 1.278e-08
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.1497756 0.2991570
    ## sample estimates:
    ##       cor 
    ## 0.2257932

``` r
# View the final results
# print(df)
```

``` r
p1
```

![](/home/ssd/ysanchez/Projects/FMI-all-Spatial-20240312/scripts/github_repository/09_Spatial_patterns_extraplacental_trophoblasts/04_notebook_cell2loc_CAM_trophoblasts_anticorrelation-2025-05-28_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
# ggsave(paste0(dirsave,"Pearson_coeffient_SC.CTB_EVT_20250112.png"), p1, width=3,height=4, bg = "white",  units = 'in', dpi = 300)
```

# Stats

To test for a statistical difference in correlation coefficients between
two groups of donors, we can use Fisher’s Z-transformation. This
transformation converts Pearson correlation coefficients into normally
distributed values, allowing for comparison between groups.

Zobserved = (z1 – z2) / (square root of \[ (1 / N1 – 3) + (1 / N2 – 3)
\]

``` r
# Fisher's Z-transformation function
fisher_z <- function(r) {
  0.5 * log((1 + r) / (1 - r))
}

# Apply transformation to each group
df$z_cor <- fisher_z(df$pearson)

# Compute means and sample sizes per group
group_means <- tapply(df$z_cor, df$GA_Condition, mean, na.rm = TRUE)
group_ns <- tapply(df$z_cor, df$GA_Condition, function(x) sum(!is.na(x)))
```

``` r
# Extract values
Z1 <- group_means[1]
Z2 <- group_means[2]
n1 <- group_ns[1]
n2 <- group_ns[2]

# Compute Z-score for difference
Z_diff <- (Z1 - Z2) / sqrt((1 / (n1 - 3)) + (1 / (n2 - 3)))

# Compute p-value (two-tailed test)
p_value <- 2 * (1 - pnorm(abs(Z_diff)))

# Print results
cat("Early disease: Z-score:", Z_diff, "\nP-value:", p_value, "\n")
```

    ## Early disease: Z-score: 0 
    ## P-value: 1

``` r
# Extract values
Z1 <- group_means[3]
Z2 <- group_means[4]
n1 <- group_ns[3]
n2 <- group_ns[4]

# Compute Z-score for difference
Z_diff <- (Z1 - Z2) / sqrt((1 / (n1 - 3)) + (1 / (n2 - 3)))

# Compute p-value (two-tailed test)
p_value <- 2 * (1 - pnorm(abs(Z_diff)))

# Print results
cat("Late disease: Z-score:", Z_diff, "\nP-value:", p_value, "\n")
```

    ## Late disease: Z-score: 0.005658 
    ## P-value: 0.9954856

# Correlation test per condition

``` r
# Subset by condition, filter the SC.CTB and EVT abundance by 0.5 
# Test for association between paired samples

df <- data.frame()
df_temp <- data.frame()

# loop through donors
conditions <- unique(spatial.1@meta.data$GA_Condition)

for (i in 1:length(conditions)) {

# Subset the data for the current condition
temp <- subset(spatial.1, subset = GA_Condition %in% conditions[i])

# filter by expression 
temp <- subset(temp, subset = c2loc_SC.CTB > 0.5)


# Check if temp has rows; skip if not
  if (nrow(temp) > 0) {
    # Extract metadata and correlation coefficient
    # donor <- unique(temp@meta.data$donor2)
    # GA_Condition <- unique(temp@meta.data$GA_Condition)
    pearson <- cor(temp$c2loc_SC.CTB, temp$c2loc_EVT, method = "pearson")
    test <- cor.test(temp$c2loc_SC.CTB, temp$c2loc_EVT, method = "pearson") 
    
    GA_Cat <- unique(temp@meta.data$GA_Cat)
      
    # Create a temporary data frame for this donor
    df_temp <- data.frame(
      GA_Condition = GA_Condition,
      pearson = pearson,
      GA_Cat = GA_Cat
    )
    
  
  
print(conditions[i])
print(test)

  # Remove temporary variable
  rm(temp)
} }
```

    ## [1] "Early_Control"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = -8.4898, df = 1186, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.2922517 -0.1849990
    ## sample estimates:
    ##        cor 
    ## -0.2393554 
    ## 
    ## [1] "Early_PE"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = -17.419, df = 2343, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.3739529 -0.3022654
    ## sample estimates:
    ##        cor 
    ## -0.3386004 
    ## 
    ## [1] "Late_Control"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = -13.209, df = 2997, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.2680909 -0.2004393
    ## sample estimates:
    ##        cor 
    ## -0.2345491 
    ## 
    ## [1] "Late_PE"
    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$c2loc_SC.CTB and temp$c2loc_EVT
    ## t = -15.561, df = 2782, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.3167841 -0.2484275
    ## sample estimates:
    ##        cor 
    ## -0.2829651

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
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggExtra_0.10.0     lmerTest_3.1-3     lme4_1.1-27.1      Matrix_1.5-4.1     ggsignif_0.6.3     ggpubr_0.4.0       cowplot_1.1.1      paletteer_1.6.0   
    ##  [9] magrittr_2.0.1     forcats_0.5.1      purrr_1.0.2        readr_2.1.1        tidyr_1.1.4        tibble_3.1.6       tidyverse_1.3.1    stringr_1.4.0     
    ## [17] ggplot2_3.3.5      knitr_1.45         patchwork_1.1.1    SeuratObject_4.1.3 Seurat_4.2.1       dplyr_1.0.7       
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] readxl_1.3.1           backports_1.4.1        plyr_1.8.6             igraph_1.2.11          lazyeval_0.2.2         sp_1.5-1               splines_4.1.1         
    ##   [8] listenv_0.8.0          scattermore_0.7        digest_0.6.29          htmltools_0.5.8.1      fansi_1.0.2            tensor_1.5             cluster_2.1.2         
    ##  [15] ROCR_1.0-11            tzdb_0.2.0             globals_0.14.0         modelr_0.1.8           matrixStats_0.61.0     spatstat.sparse_3.0-0  colorspace_2.0-2      
    ##  [22] rvest_1.0.2            ggrepel_0.9.1          haven_2.4.3            xfun_0.41              prismatic_1.1.2        crayon_1.4.2           jsonlite_1.7.3        
    ##  [29] progressr_0.10.0       spatstat.data_3.0-0    survival_3.2-13        zoo_1.8-9              glue_1.6.1             polyclip_1.10-0        gtable_0.3.0          
    ##  [36] leiden_0.3.9           car_3.0-12             future.apply_1.8.1     abind_1.4-5            scales_1.1.1           DBI_1.1.2              rstatix_0.7.0         
    ##  [43] spatstat.random_3.0-1  miniUI_0.1.1.1         Rcpp_1.0.8             viridisLite_0.4.0      xtable_1.8-4           reticulate_1.24        htmlwidgets_1.5.4     
    ##  [50] httr_1.4.2             RColorBrewer_1.1-2     ellipsis_0.3.2         ica_1.0-2              farver_2.1.0           pkgconfig_2.0.3        uwot_0.1.14           
    ##  [57] dbplyr_2.1.1           deldir_1.0-6           utf8_1.2.2             labeling_0.4.2         tidyselect_1.1.1       rlang_1.1.1            reshape2_1.4.4        
    ##  [64] later_1.3.0            munsell_0.5.0          cellranger_1.1.0       tools_4.1.1            cli_3.6.1              generics_0.1.1         broom_0.7.11          
    ##  [71] ggridges_0.5.3         evaluate_0.23          fastmap_1.1.1          yaml_2.2.2             goftest_1.2-3          rematch2_2.1.2         fs_1.5.2              
    ##  [78] fitdistrplus_1.1-6     RANN_2.6.1             pbapply_1.5-0          future_1.23.0          nlme_3.1-155           mime_0.12              xml2_1.3.3            
    ##  [85] compiler_4.1.1         rstudioapi_0.13        plotly_4.10.0          png_0.1-7              spatstat.utils_3.0-1   reprex_2.0.1           stringi_1.7.6         
    ##  [92] highr_0.9              lattice_0.20-45        nloptr_1.2.2.3         vctrs_0.6.5            pillar_1.6.5           lifecycle_1.0.4        spatstat.geom_3.0-3   
    ##  [99] lmtest_0.9-39          RcppAnnoy_0.0.19       data.table_1.14.2      irlba_2.3.5            httpuv_1.6.5           R6_2.5.1               promises_1.2.0.1      
    ## [106] KernSmooth_2.23-20     gridExtra_2.3          parallelly_1.30.0      codetools_0.2-18       boot_1.3-28            MASS_7.3-55            assertthat_0.2.1      
    ## [113] withr_2.5.0            sctransform_0.3.5      parallel_4.1.1         hms_1.1.1              grid_4.1.1             minqa_1.2.4            rmarkdown_2.25        
    ## [120] carData_3.0-5          Rtsne_0.15             spatstat.explore_3.0-5 numDeriv_2016.8-1.1    shiny_1.7.1            lubridate_1.8.0

``` r
# rmarkdown::render(
#   input = "/home/ssd/ysanchez/Projects/FMI-all-Spatial-20240312/scripts/github_repository/09_Spatial_patterns_extraplacental_trophoblasts/04_notebook_cell2loc_CAM_trophoblasts_anticorrelation.Rmd",
#   output_format = "github_document",
# output_file = paste0("/home/ssd/ysanchez/Projects/FMI-all-Spatial-20240312/scripts/github_repository/09_Spatial_patterns_extraplacental_trophoblasts/04_notebook_cell2loc_CAM_trophoblasts_anticorrelation-", Sys.Date(), ".md"),
# )
```
