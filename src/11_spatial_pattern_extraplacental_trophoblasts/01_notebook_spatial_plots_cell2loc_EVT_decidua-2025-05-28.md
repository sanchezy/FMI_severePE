Spatial transcriptomics: Plot EVT deconvoluted at PVBP and Myometrium
================

Author: Yara E. Sanchez Corrales

# Load libraries

# Introduction

Here I plot the spatial deconvoluted plots of EVTs. In this notebook,
the images are commented.

``` r
# order the donor: Early_Control, Early_PE, Late_Control, Late_PE,

dirsave <- "~/Projects/FMI-all-Spatial-20240312/plots/20250124_Spatial_maps_EVTs_decidua_basalis/"

donors_to_plot <- c("FVQ", 
                 "FLR",
                 "FVB",

                 "FCM",
                 "FAM",
                 "FGS",
                 "FEP",
                 "FVS",
                 
                 "FLJ",
                 "FJD",
                 "FVT",
                 "FYI",
                 "FRK",
                 "FNS",
                 
                 "FFM",
                 "FRH",
                 "FND",
                 "FSH",
                 "FMM")
```

**Spatial domains to plot** <br>

    ## [1] "decidua-myo" "muscle"      "decidua-cam"

# Spatial maps

``` r
# transparency of the HE image
alpha <- 0.8

# point-size-factor
spotsize <- 1.3

feature_to_plot <- c("c2loc_EVT")

# width and lentght for plots

w = 7
h = 7 
  
# limits
l <- c(-2,4)
```

## PVBP

``` r
# make and empty list of plots
tissue <- "PVBP"
# spotsize <- 1.0

plots1 = list()
plots5 = list()
plots = list()

donors_to_plot <- c("FVB","FAM") 

metadata <- spatial.placenta@meta.data


# loop through donors
for (i in 1:length(donors_to_plot)) {

  metadatatemp <- metadata %>% filter(donor2 %in% donors_to_plot[i])

  images_to_plot <- unique(metadatatemp$Image_name)
  condition <- unique(metadatatemp$GA_Condition)
    
plots1[[i]]  <- SpatialFeaturePlot(spatial.placenta, features = feature_to_plot[1], pt.size.factor = spotsize, ncol = 1, crop = FALSE, images = images_to_plot,image.alpha = 0.3)  &  paletteer::scale_fill_paletteer_c("viridis::plasma", limits = c(0, 3)) &  plot_annotation(title = paste0(condition," - ",donors_to_plot[i], " - ", feature_to_plot[1]))

plots5[[i]]  <- SpatialDimPlot(spatial.placenta, group.by = 'regions', ncol = 1, crop = FALSE, pt.size.factor = spotsize, images = images_to_plot, cols = region_colors)

}
```

``` r
for (i in 1:length(donors_to_plot)) {
  # print(plots1[[i]])
  # 
  # ggsave(paste0(dirsave, tissue, "_", donors_to_plot[i],"_region_decidua_c2loc_EVT_score_20250124.png"), plots1[[i]], width=7,height=7, bg = "white",  units = 'in', dpi = 300)
}
```

``` r
for (i in 1:length(donors_to_plot)) {
  
  # print(plots5[[i]])
  # ggsave(paste0(dirsave, tissue, "_", donors_to_plot[i],"_regions_20250124.png"), plots5[[i]], width=7,height=7, bg = "white",  units = 'in', dpi = 300)
}
```

``` r
# legend for regions
p1 <- plots5[[1]] + theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) + guides(fill = guide_legend(override.aes = list(size=22), guide_legend(title="Spatial domain")))

# Extract the legend. Returns a gtable
leg <- get_legend(p1)

# Convert to a ggplot and print - ggpubr 
legend <- as_ggplot(leg)
# legend
# 
# ggsave(paste0(dirsave,tissue,"_legend_spatial_domain_20250124.png"), legend, width=2,height=3, bg = "white",  units = 'in', dpi = 300)
```

``` r
# legend for ifn1 score
p1 <- plots1[[1]] + theme(legend.key.width=unit(3,"cm"),legend.key.height = unit(1.5,"cm"),legend.title=element_text(size=20), legend.text=element_text(size=20) )  #+ theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) + guides(guide_legend(title="IFN-I score"))

# Extract the legend. Returns a gtable
leg <- get_legend(p1)

# Convert to a ggplot and print - ggpubr 
legend <- as_ggplot(leg)
# legend
# 
# ggsave(paste0(dirsave,tissue,"_legend_c2loc_EVT_20250124.png"), legend, width=8,height=2, bg = "white",  units = 'in', dpi = 300)
```

## Myometrium

``` r
for (i in 1:length(donors_to_plot)) {
  # print(plots1[[i]])
  # 
  # ggsave(paste0(dirsave, tissue, "_", donors_to_plot[i],"_region_decidua_c2loc_EVT_score_20250124.png"), plots1[[i]], width=7,height=7, bg = "white",  units = 'in', dpi = 300)
  # print(plots5[[i]])
}
```

``` r
for (i in 1:length(donors_to_plot)) {
  # # print(plots1[[i]])
  # print(plots5[[i]])
  # ggsave(paste0(dirsave, tissue, "_", donors_to_plot[i],"_regions_20250124.png"), plots5[[i]], width=7,height=7, bg = "white",  units = 'in', dpi = 300)
}
```

``` r
# legend for regions
p1 <- plots5[[1]] + theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) + guides(fill = guide_legend(override.aes = list(size=22), guide_legend(title="Spatial domain")))

# Extract the legend. Returns a gtable
leg <- get_legend(p1)

# Convert to a ggplot and print - ggpubr 
legend <- as_ggplot(leg)
# legend
# 
# ggsave(paste0(dirsave,tissue,"_legend_spatial_domain_20250124.png"), legend, width=2,height=3, bg = "white",  units = 'in', dpi = 300)
```

``` r
# legend for ifn1 score
p1 <- plots1[[1]] + theme(legend.key.width=unit(3,"cm"),legend.key.height = unit(1.5,"cm"),legend.title=element_text(size=20), legend.text=element_text(size=20) )  #+ theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) + guides(guide_legend(title="IFN-I score"))

# Extract the legend. Returns a gtable
leg <- get_legend(p1)

# Convert to a ggplot and print - ggpubr 
legend <- as_ggplot(leg)
# legend
# 
# ggsave(paste0(dirsave,tissue,"_legend_c2loc_EVT_20250124.png"), legend, width=8,height=2, bg = "white",  units = 'in', dpi = 300)
```

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
    ##  [1] ggpubr_0.4.0       cowplot_1.1.1      paletteer_1.6.0    magrittr_2.0.1     forcats_0.5.1      purrr_1.0.2        readr_2.1.1        tidyr_1.1.4       
    ##  [9] tibble_3.1.6       tidyverse_1.3.1    stringr_1.4.0      ggplot2_3.3.5      knitr_1.45         patchwork_1.1.1    SeuratObject_4.1.3 Seurat_4.2.1      
    ## [17] dplyr_1.0.7       
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
    ##  [85] compiler_4.1.1         rstudioapi_0.13        plotly_4.10.0          png_0.1-7              ggsignif_0.6.3         spatstat.utils_3.0-1   reprex_2.0.1          
    ##  [92] stringi_1.7.6          highr_0.9              lattice_0.20-45        Matrix_1.5-4.1         vctrs_0.6.5            pillar_1.6.5           lifecycle_1.0.4       
    ##  [99] spatstat.geom_3.0-3    lmtest_0.9-39          RcppAnnoy_0.0.19       data.table_1.14.2      irlba_2.3.5            httpuv_1.6.5           R6_2.5.1              
    ## [106] promises_1.2.0.1       KernSmooth_2.23-20     gridExtra_2.3          parallelly_1.30.0      codetools_0.2-18       MASS_7.3-55            assertthat_0.2.1      
    ## [113] withr_2.5.0            sctransform_0.3.5      parallel_4.1.1         hms_1.1.1              grid_4.1.1             rmarkdown_2.25         carData_3.0-5         
    ## [120] Rtsne_0.15             spatstat.explore_3.0-5 shiny_1.7.1            lubridate_1.8.0


