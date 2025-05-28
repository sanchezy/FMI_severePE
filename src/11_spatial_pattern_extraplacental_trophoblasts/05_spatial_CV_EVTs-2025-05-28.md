Spatial transcriptomics: Coefficient of variation
================

Author: Yara E. Sanchez Corrales

# Introduction

Here quantify the CV per capture area in EVT. Plot it per condition and
do some stats.

# Load libraries

``` r
dirsave <- "/Users/ysanchez/Documents/Projects-analysis/FMI-Spatial-transcriptomics-20240316/plots/20250318_EVT_spatial_variation/"
```

# Read the object

    ## Processing slice: FVQ.PVBP1 
    ## [1] 0.4317103
    ## Processing slice: FVQ.PVBP2 
    ## [1] 0.3772355
    ## Processing slice: FLR.PVBP 
    ## [1] 0.2912858
    ## Processing slice: FVB.PVBP 
    ## [1] 0.2357514
    ## Processing slice: FCM.PVBP1 
    ## [1] 0.8156569
    ## Processing slice: FCM.PVBP2 
    ## [1] 0.8830684
    ## Processing slice: FAM.PVBP 
    ## [1] 0.4825008
    ## Processing slice: FEP.PVBP 
    ## [1] 0.5029176
    ## Processing slice: FLJ.PVBP1 
    ## [1] 0.3258394
    ## Processing slice: FLJ.PVBP2 
    ## [1] 0.6529948
    ## Processing slice: FJD.PVBP 
    ## [1] 0.5489463
    ## Processing slice: FVT.PVBP 
    ## [1] 0.5019533
    ## Processing slice: FYI.PVBP 
    ## [1] 0.5474003
    ## Processing slice: FNS.PVBP 
    ## [1] 1.289241
    ## Processing slice: FRK.PVBP1 
    ## [1] 0.5358485
    ## Processing slice: FFM.PVBP 
    ## [1] 0.4169061
    ## Processing slice: FRH.PVBP1 
    ## [1] 0.4687007
    ## Processing slice: FSH.PVBP 
    ## [1] 0.3013543
    ## Processing slice: FND.PVBP 
    ## [1] 0.3063303
    ## Processing slice: FMM.PVBP 
    ## [1] 0.5218167

    ## Processing slice: FVQ.MYO 
    ## [1] 0.6075197
    ## Processing slice: FVB.MYO 
    ## [1] 0.4654599
    ## Processing slice: FCM.MYO 
    ## [1] 1.321458
    ## Processing slice: FAM.MYO 
    ## [1] 1.360072
    ## Processing slice: FGS.MYO 
    ## [1] 0.6729663
    ## Processing slice: FEP.MYO 
    ## [1] 2.135942
    ## Processing slice: FVS.MYO 
    ## [1] 1.414201
    ## Processing slice: FLJ.MYO 
    ## [1] 1.02069
    ## Processing slice: FJD.MYO 
    ## [1] 1.495754
    ## Processing slice: FVT.MYO 
    ## [1] 1.174067
    ## Processing slice: FRK.MYO 
    ## [1] 0.4244746
    ## Processing slice: FFM.MYO 
    ## [1] 0.7153717
    ## Processing slice: FRH.MYO 
    ## [1] 0.7050134
    ## Processing slice: FSH.MYO 
    ## [1] 0.3644783
    ## Processing slice: FND.MYO 
    ## [1] 1.074453
    ## Processing slice: FMM.MYO 
    ## [1] 0.3871452

    ## Processing slice: FVQ.MYO 
    ## [1] 2.545864
    ## Processing slice: FLR.MYO1 
    ## [1] 0.7517317
    ## Processing slice: FLR.MYO2 
    ## [1] 0.6315411
    ## Processing slice: FVB.MYO 
    ## [1] 2.401406
    ## Processing slice: FCM.MYO 
    ## [1] 2.573039
    ## Processing slice: FAM.MYO 
    ## [1] 2.469685
    ## Processing slice: FGS.MYO 
    ## [1] 2.155911
    ## Processing slice: FEP.MYO 
    ## [1] 0.8130542
    ## Processing slice: FVS.MYO 
    ## [1] 1.171086
    ## Processing slice: FLJ.MYO 
    ## [1] 1.39653
    ## Processing slice: FJD.MYO 
    ## [1] 2.358099
    ## Processing slice: FVT.MYO 
    ## [1] 1.504043
    ## Processing slice: FYI.MYO 
    ## [1] 0.645825
    ## Processing slice: FNS.MYO 
    ## [1] 0.6576659
    ## Processing slice: FRK.MYO 
    ## [1] 0.8981047
    ## Processing slice: FFM.MYO 
    ## [1] 0.481908
    ## Processing slice: FRH.MYO 
    ## [1] 0.6320986
    ## Processing slice: FSH.MYO 
    ## [1] 1.824247
    ## Processing slice: FND.MYO 
    ## [1] 0.9696552
    ## Processing slice: FMM.MYO 
    ## [1] 0.2725343

### CV

``` r
barplot <-  ggplot(cv_placenta, aes(x=GA_Condition, y=CV, fill = GA_Condition)) +
  geom_boxplot(color="black") +  
 scale_fill_manual(values = colors_conditions2)  + 
 theme(
  panel.background = element_blank(),
  axis.line = element_line(colour = "grey"),
   strip.background =element_rect(fill="white"),
    axis.text.x = element_text(angle = 0, size = 12),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 12, color = "black"),
    legend.position = "left",
    legend.text = element_text(size = 18),
    legend.title= element_blank()) +
    labs(x= "", y = "Coefficent of variation")  + NoLegend() +
    scale_x_discrete(labels=c("Early_Control" = "Control", "Early_PE" = "PE",
                              "Late_Control" = "Control", "Late_PE" = "PE")) +
     facet_wrap(~ GA_Cat, scales = "free_x") + 
# Add significance bracket with manual asterisks
  geom_signif(comparisons = list(c("Early_Control", "Early_PE")), 
              annotations = "**",  # Change this to "**" or "*" depending on your p-value
              y_position = max(cv_placenta$CV) * 1.6,  # Adjust the position
              tip_length = 0.02,  # Length of the bracket tips
              textsize = 6) + # Size of the asterisks +
  
    scale_x_discrete(labels=c("Early_Control" = "Control", "Early_PE" = "PE",
                              "Late_Control" = "Control", "Late_PE" = "PE")) # Facet by Condition 
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

``` r
 barplot <-barplot +  ylim(0.2,2.2)
 
# ggsave(paste0(dirsave,"PVBP_decidua_CV_c2loc_EVT_20250318.png"), barplot, width=3,height=4, bg = "white",  units = 'in', dpi = 300) 
# barplot
```

``` r
# Convert GA_Condition to a factor with explicit levels
cv_placenta$GA_Condition <- factor(cv_placenta$GA_Condition, 
                                   levels = c("Early_Control", "Early_PE", "Late_Control", "Late_PE"))


# Updated barplot code
barplot <- ggplot(cv_placenta, aes(x = GA_Condition, y = CV, fill = GA_Condition)) +
  geom_boxplot(color = "black") +  
  scale_fill_manual(values = colors_conditions2) + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 0, size = 12),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 12, color = "black"),
    legend.position = "left",
    legend.text = element_text(size = 18),
    legend.title = element_blank()
  ) +
  labs(x = "", y = "Coefficient of variation") +
  NoLegend() +
  scale_x_discrete(labels = c("Early_Control" = "Control", "Early_PE" = "PE",
                              "Late_Control" = "Control", "Late_PE" = "PE")) +
  facet_wrap(~ GA_Cat, scales = "free_x") +

  # Add multiple significance brackets with geom_signif
  geom_signif(
    comparisons = list(c("Early_Control", "Early_PE"), c("Late_Control", "Late_PE")), 
    annotations = c("*", "*"),  # Customize significance level
    map_signif_level = TRUE,  # Automatically map significance level to symbol
    y_position = c(max(cv_placenta$CV) * 0.1, max(cv_placenta$CV) * 0.1),  # Adjust bracket heights
    tip_length = 0.02,  # Length of bracket tips
    textsize = 6  # Size of the annotation text
  )

 barplot <- barplot +  ylim(0.2,2.2)
# Display the plot
# barplot
```

``` r
# Same but with the same scale
# p1 <- boxplot_gene(cv_placenta, CV,colors_conditions2) + ggtitle("PVBP decidua") + ylim(0.1,2.5)
# p2 <- boxplot_gene(cv_decidua, CV,colors_conditions2) + ggtitle("Myometrium decidua") + ylim(0.1,2.5)
# p3 <- boxplot_gene(cv_muscle, CV,colors_conditions2) + ggtitle("Myometrium muscle") + ylim(0,2.5)
# 
# p1 + p2 + p3
```

# Stats

### CV

``` r
# decidua of the PVBP
stats_cv(cv_placenta)
```

    ## 
    ##  Wilcoxon rank sum exact test
    ## 
    ## data:  data_early_p$CV by data_early_p$GA_Condition
    ## W = 0, p-value = 0.02857
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ##  Wilcoxon rank sum exact test
    ## 
    ## data:  data_late_p$CV by data_late_p$GA_Condition
    ## W = 31, p-value = 0.0303
    ## alternative hypothesis: true location shift is not equal to 0

``` r
# decidua of the Myometrium
stats_cv(cv_decidua)
```

    ## 
    ##  Wilcoxon rank sum exact test
    ## 
    ## data:  data_early_p$CV by data_early_p$GA_Condition
    ## W = 0, p-value = 0.09524
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ##  Wilcoxon rank sum exact test
    ## 
    ## data:  data_late_p$CV by data_late_p$GA_Condition
    ## W = 16, p-value = 0.1905
    ## alternative hypothesis: true location shift is not equal to 0

``` r
# muscle of myometrium
stats_cv(cv_muscle)
```

    ## 
    ##  Wilcoxon rank sum exact test
    ## 
    ## data:  data_early_p$CV by data_early_p$GA_Condition
    ## W = 7, p-value = 0.5556
    ## alternative hypothesis: true location shift is not equal to 0
    ## 
    ## 
    ##  Wilcoxon rank sum exact test
    ## 
    ## data:  data_late_p$CV by data_late_p$GA_Condition
    ## W = 22, p-value = 0.2468
    ## alternative hypothesis: true location shift is not equal to 0

``` r
sessionInfo()
```

    ## R version 4.3.2 (2023-10-31)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS Sonoma 14.3
    ## 
    ## Matrix products: default
    ## BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: Europe/London
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] spdep_1.3-10       sf_1.0-19          spData_2.3.4       ggsignif_0.6.4     paletteer_1.6.0    magrittr_2.0.3     lubridate_1.9.3    forcats_1.0.0     
    ##  [9] purrr_1.0.4        readr_2.1.5        tidyr_1.3.1        tibble_3.2.1       tidyverse_2.0.0    stringr_1.5.1      ggplot2_3.5.2      knitr_1.45        
    ## [17] patchwork_1.2.0    Seurat_5.0.1       SeuratObject_5.0.1 sp_2.1-3           dplyr_1.1.4       
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3     wk_0.9.4               rstudioapi_0.17.1      jsonlite_1.8.8         spatstat.utils_3.0-4   rmarkdown_2.25         vctrs_0.6.5           
    ##   [8] ROCR_1.0-11            spatstat.explore_3.2-6 htmltools_0.5.7        s2_1.1.7               sctransform_0.4.1      parallelly_1.36.0      KernSmooth_2.23-22    
    ##  [15] htmlwidgets_1.6.4      ica_1.0-3              plyr_1.8.9             plotly_4.10.4          zoo_1.8-12             igraph_2.0.1.1         mime_0.12             
    ##  [22] lifecycle_1.0.4        pkgconfig_2.0.3        Matrix_1.6-5           R6_2.5.1               fastmap_1.1.1          fitdistrplus_1.1-11    future_1.33.1         
    ##  [29] shiny_1.8.0            digest_0.6.34          colorspace_2.1-0       rematch2_2.1.2         tensor_1.5             RSpectra_0.16-1        irlba_2.3.5.1         
    ##  [36] progressr_0.14.0       fansi_1.0.6            spatstat.sparse_3.0-3  timechange_0.3.0       httr_1.4.7             polyclip_1.10-6        abind_1.4-5           
    ##  [43] compiler_4.3.2         proxy_0.4-27           withr_3.0.1            DBI_1.2.1              fastDummies_1.7.3      MASS_7.3-60.0.1        classInt_0.4-11       
    ##  [50] units_0.8-7            tools_4.3.2            lmtest_0.9-40          httpuv_1.6.14          future.apply_1.11.1    goftest_1.2-3          glue_1.7.0            
    ##  [57] nlme_3.1-164           promises_1.2.1         grid_4.3.2             Rtsne_0.17             cluster_2.1.6          reshape2_1.4.4         generics_0.1.3        
    ##  [64] gtable_0.3.4           spatstat.data_3.0-4    tzdb_0.4.0             class_7.3-22           data.table_1.16.0      hms_1.1.3              utf8_1.2.4            
    ##  [71] spatstat.geom_3.2-8    RcppAnnoy_0.0.22       ggrepel_0.9.5          RANN_2.6.1             pillar_1.9.0           spam_2.10-0            RcppHNSW_0.6.0        
    ##  [78] later_1.3.2            splines_4.3.2          lattice_0.22-5         survival_3.5-7         deldir_2.0-2           tidyselect_1.2.1       miniUI_0.1.1.1        
    ##  [85] pbapply_1.7-2          gridExtra_2.3          scattermore_1.2        xfun_0.42              matrixStats_1.2.0      stringi_1.8.4          boot_1.3-28.1         
    ##  [92] lazyeval_0.2.2         yaml_2.3.8             evaluate_0.23          codetools_0.2-19       cli_3.6.5              uwot_0.1.16            xtable_1.8-4          
    ##  [99] reticulate_1.35.0      munsell_0.5.0          Rcpp_1.0.13            globals_0.16.2         spatstat.random_3.2-2  png_0.1-8              parallel_4.3.2        
    ## [106] ellipsis_0.3.2         dotCall64_1.1-1        listenv_0.9.1          viridisLite_0.4.2      e1071_1.7-14           scales_1.3.0           ggridges_0.5.6        
    ## [113] leiden_0.4.3.1         rlang_1.1.6            cowplot_1.1.3

