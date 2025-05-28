Spatial transcriptomics: Calculate the mean abundance of EVTs
deconvoluted per spatial region
================

Author: Yara E. Sanchez-Corrales

# Introduction

Here I want to plot the abundance of EVT per spot per region and per
donor. The relevant regions are: <br> 1.- Decidua (Placenta) <br> 2.-
Decidua (Myometrium) <br>

We detected very few cells in the myometrium - almost indistinguisable.

The hypothesis is that EVTs are less invasive in PE, so we expect less
EVTs in disease.

*Note*: Update in November: new deconvolution done using a new
single-cell reference

For the stats, we use Linear Mixed Model (LMM) using lme4. It allows for
random intercepts per donor, accounting for repeated measures.

# Load libraries

``` r
library(dplyr)
library(Seurat)
library(patchwork)
library(knitr)
library(ggplot2)
library(stringr)
library(tidyverse)
library(magrittr)
library(paletteer)
library(ggsignif)
```

``` r
dirsave <- "~/Projects/FMI-all-Spatial-20240312/plots/20241204_EVT_decidua_basalis/"
```

# Read the object

``` r
# read the object
spatial.myo <- readRDS(file = "~/Projects/FMI-all-Spatial-20240312/objects/FMI-Myometrium-Spatial-20241120.rds")

spatial.placenta <- readRDS(file = "~/Projects/FMI-all-Spatial-20240312/objects/FMI-Placenta-Spatial-20241120.rds")

# create a new metadata column
spatial.myo@meta.data$donor2 <- paste0("F", str_sub(spatial.myo@meta.data$donor, start= -2))
spatial.placenta@meta.data$donor2 <- paste0("F", str_sub(spatial.placenta@meta.data$donor, start= -2))

# create a composite column with the condition
spatial.placenta@meta.data$donor.Condition2 <- paste0(spatial.placenta@meta.data$Condition,"_",spatial.placenta@meta.data$donor2)
spatial.myo@meta.data$donor.Condition2 <- paste0(spatial.myo@meta.data$Condition,"_",spatial.myo@meta.data$donor2)
```

``` r
# Remove the background 
spatial.myo <- subset(spatial.myo, subset = regions %in% c("background", "villi"), invert = TRUE)
spatial.placenta <- subset(spatial.placenta, subset = regions %in% c("background"), invert = TRUE)
unique(spatial.myo@meta.data$regions)
```

    ## [1] "muscle"  "decidua"

``` r
unique(spatial.placenta@meta.data$regions)
```

    ## [1] "villi"                          "vasculature_intermediate_villi" "decidua"

    ## # A tibble: 6 × 3
    ##   donor   regions total_spots
    ##   <chr>   <chr>         <int>
    ## 1 F1668RK decidua         769
    ## 2 F1668RK muscle         1481
    ## 3 F1676VQ decidua         117
    ## 4 F1676VQ muscle         2607
    ## 5 F1678CM decidua         289
    ## 6 F1678CM muscle         1183

``` r
# dataframe with the information about number of spots per region per donor
number_of_spots_per_region_placenta <- spatial.placenta@meta.data %>% group_by(donor,regions) %>%
summarise(total_spots = n(), .groups = 'drop')
# number_of_spots_per_region_placenta
```

# Placenta: spatial region: decidua-placenta

``` r
# loop through donors and append the df.

# List of donors
donors <- unique(spatial.placenta@meta.data$donor)

results_c2loc_decidua_p <- data.frame()
results_c2loc_decidua_p <- df_c2loc_region_donor(spatial.placenta,number_of_spots_per_region_placenta, donors[1], "decidua")
```

    ## Joining, by = "donor"
    ## Joining, by = "CellTypeManual.l2"
    ## Joining, by = "donor"

``` r
for (d in 2:length(donors)) {
  temp <- data.frame()
  temp <- df_c2loc_region_donor(spatial.placenta, number_of_spots_per_region_placenta, donors[d], "decidua")
  results_c2loc_decidua_p <- rbind(temp,results_c2loc_decidua_p)
}
```

    ## Joining, by = "donor"
    ## Joining, by = "CellTypeManual.l2"
    ## Joining, by = "donor"
    ## Joining, by = "donor"
    ## Joining, by = "CellTypeManual.l2"
    ## Joining, by = "donor"
    ## Joining, by = "donor"
    ## Joining, by = "CellTypeManual.l2"
    ## Joining, by = "donor"
    ## Joining, by = "donor"
    ## Joining, by = "CellTypeManual.l2"
    ## Joining, by = "donor"
    ## Joining, by = "donor"
    ## Joining, by = "CellTypeManual.l2"
    ## Joining, by = "donor"
    ## Joining, by = "donor"
    ## Joining, by = "CellTypeManual.l2"
    ## Joining, by = "donor"
    ## Joining, by = "donor"
    ## Joining, by = "CellTypeManual.l2"
    ## Joining, by = "donor"
    ## Joining, by = "donor"
    ## Joining, by = "CellTypeManual.l2"
    ## Joining, by = "donor"
    ## Joining, by = "donor"
    ## Joining, by = "CellTypeManual.l2"
    ## Joining, by = "donor"
    ## Joining, by = "donor"
    ## Joining, by = "CellTypeManual.l2"
    ## Joining, by = "donor"
    ## Joining, by = "donor"
    ## Joining, by = "CellTypeManual.l2"
    ## Joining, by = "donor"
    ## Joining, by = "donor"
    ## Joining, by = "CellTypeManual.l2"
    ## Joining, by = "donor"
    ## Joining, by = "donor"
    ## Joining, by = "CellTypeManual.l2"
    ## Joining, by = "donor"
    ## Joining, by = "donor"
    ## Joining, by = "CellTypeManual.l2"
    ## Joining, by = "donor"
    ## Joining, by = "donor"
    ## Joining, by = "CellTypeManual.l2"
    ## Joining, by = "donor"
    ## Joining, by = "donor"
    ## Joining, by = "CellTypeManual.l2"
    ## Joining, by = "donor"
    ## Joining, by = "donor"
    ## Joining, by = "CellTypeManual.l2"
    ## Joining, by = "donor"
    ## Joining, by = "donor"
    ## Joining, by = "CellTypeManual.l2"
    ## Joining, by = "donor"

``` r
# # # append the information about GA_Condition
results_c2loc_decidua_p <- results_c2loc_decidua_p %>% inner_join(donor.order.condition)
```

    ## Joining, by = c("donor", "donor2", "condition", "GA_Condition")

``` r
head(results_c2loc_decidua_p)
```

    ## # A tibble: 6 × 12
    ##   CellTypeManual.l2 predicted_abund… donor regions total_spots abundance_per_s… proportion_regi… Cumulative_norm… cellgroup
    ##   <chr>                        <dbl> <chr> <chr>         <int>            <dbl>            <dbl>            <dbl> <chr>    
    ## 1 EVT                           844. F202… decidua         721            1.17            0.155             0.155 Trophobl…
    ## 2 Decidual                      428. F202… decidua         721            0.593           0.0786            0.234 Stromal  
    ## 3 CTB                           383. F202… decidua         721            0.531           0.0704            0.304 Trophobl…
    ## 4 STB                           294. F202… decidua         721            0.408           0.0541            0.428 Trophobl…
    ## 5 Dendritic                     285. F202… decidua         721            0.395           0.0524            0.481 Myeloid  
    ## 6 NKT                           281. F202… decidua         721            0.390           0.0517            0.532 T/NK     
    ## # … with 3 more variables: donor2 <chr>, condition <chr>, GA_Condition <chr>

# Plot all the values (per spot) EVT in PVBP and Myometrium

``` r
# metadata, column, region, colors_conditions
# boxplot_gene_donor(spatial.myo@meta.data,!!sym(c2loc_EVT),"muscle",colors_conditions2) 
##### PVBP 
region <- "decidua"
column <- "c2loc_EVT"
metadata <- spatial.placenta@meta.data
colors_conditions <- colors_conditions2
# 
  
 filtered_data <- metadata %>%
    filter(regions %in% region)
    
     # Filter rows where the selected 'gene' column is > 0
filtered_data <- filtered_data %>%
select(donor.Condition2,"c2loc_EVT", donor2, donor,GA_Condition,GA_Cat) 
  
  # Filter rows where the selected 'gene' column is > 0
  filtered_data <- filtered_data %>%
    filter(c2loc_EVT > 0)
  
  
  # reorder by donor gestational age
  filtered_data$donor.Condition2 <- factor(filtered_data$donor.Condition2, levels =  donor.order.condition$donor2, ordered = TRUE)
  
  # Create the bar plot using ggplot2
  
plot_decidua_p_all <-  ggplot(filtered_data, aes(x=GA_Condition, y=c2loc_EVT, fill = GA_Condition)) +
  geom_boxplot(color="black",outlier.size = 0.3) +  
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), dotsize = 0.0) +
 scale_fill_manual(values = colors_conditions)  + 
 theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    strip.background =element_rect(fill="white"),
    axis.text.x = element_text(angle = 0, size = 12),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 12, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 18),
    legend.title= element_blank()) +
labs(y = "Mean abundance at decidua", x = "")  +   # Add significance bracket with manual asterisks
  geom_signif(comparisons = list(c("Early_Control", "Early_PE")), 
              annotations = "*",  # Change this to "**" or "*" depending on your p-value
              y_position = max(filtered_data$c2loc_EVT) * 1.05,  # Adjust the position
              tip_length = 0.02,  # Length of the bracket tips
              textsize = 6) + NoLegend() +
    scale_x_discrete(labels=c("Early_Control" = "Control", "Early_PE" = "PE",
                              "Late_Control" = "Control", "Late_PE" = "PE"))  +
     facet_grid(.~ GA_Cat, scales = "free_x")  # Facet by Condition 

filtered_data_p <- filtered_data

# plot_decidua_p_all
```

``` r
# metadata, column, region, colors_conditions
# boxplot_gene_donor(spatial.myo@meta.data,!!sym(c2loc_EVT),"muscle",colors_conditions2) 

##### Myometrium
region <- "decidua"
column <- "c2loc_EVT"
metadata <- spatial.myo@meta.data
colors_conditions <- colors_conditions2
# 
  
 filtered_data <- metadata %>%
    filter(regions %in% region)
    
     # Filter rows where the selected 'gene' column is > 0
filtered_data <- filtered_data %>%
select(donor.Condition2,"c2loc_EVT", donor2, donor,GA_Condition,GA_Cat) 
  
  # Filter rows where the selected 'gene' column is > 0
  filtered_data <- filtered_data %>%
    filter(c2loc_EVT > 0)
  
  # Add the info about 
    # filtered_data <- filtered_data %>% inner_join(donor.order.condition)
  
  # reorder by donor gestational age
  filtered_data$donor.Condition2 <- factor(filtered_data$donor.Condition2, levels =  donor.order.condition$donor2, ordered = TRUE)
  
  # Create the bar plot using ggplot2
  
plot_decidua_myo_all <-  ggplot(filtered_data, aes(x=GA_Condition, y=c2loc_EVT, fill = GA_Condition)) +
  geom_boxplot(color="black",outlier.size = 0.3) +  
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), dotsize = 0.0) +
 scale_fill_manual(values = colors_conditions)  + 
 theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    strip.background =element_rect(fill="white"),
    axis.text.x = element_text(angle = 0, size = 12),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 12, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 18),
    legend.title= element_blank()) +
labs(y = "Mean abundance at decidua", x = "")  + NoLegend() +
  

 # Add significance bracket with manual asterisks
  geom_signif(comparisons = list(c("Early_Control", "Early_PE")), 
              annotations = "**",  # Change this to "**" or "*" depending on your p-value
              y_position = max(filtered_data$c2loc_EVT) * 1.25,  # Adjust the position
              tip_length = 0.02,  # Length of the bracket tips
              textsize = 6) + # Size of the asterisks +
  
    scale_x_discrete(labels=c("Early_Control" = "Control", "Early_PE" = "PE",
                              "Late_Control" = "Control", "Late_PE" = "PE")) +
     facet_grid(.~ GA_Cat, scales = "free_x")  # Facet by Condition 

filtered_data_m <- filtered_data

# plot_decidua_myo_all
```

``` r
plot_decidua_p_all <- plot_decidua_p_all + ylim(0, 7) + NoLegend()
plot_decidua_myo_all <- plot_decidua_myo_all + ylim(0, 7) + NoLegend()

plot_decidua_p_all + plot_decidua_myo_all
```

    ## Bin width defaults to 1/30 of the range of the data. Pick better value with `binwidth`.

    ## Warning: Computation failed in `stat_signif()`:
    ## missing value where TRUE/FALSE needed

    ## Bin width defaults to 1/30 of the range of the data. Pick better value with `binwidth`.

    ## Warning: Computation failed in `stat_signif()`:
    ## missing value where TRUE/FALSE needed

![](/home/ssd/ysanchez/Projects/FMI-all-Spatial-20240312/scripts/github_repository/09_Spatial_patterns_extraplacental_trophoblasts/02_notebook_spatial_mean_abundance_EVT_cell2loc_results-2025-05-28_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

# Stats

## Stats - EVT decidua at PVBP

``` r
# histogram - is it normally distributed? No.
data_p <- filtered_data_p

hist <- ggplot(data_p, aes(x=c2loc_EVT, color=GA_Condition)) +
  geom_histogram(fill="white") 
# hist
```

### Early - PVBP

``` r
# Assuming 'data' is a data frame with columns 'abundance_per_spot'  and 'condition'
data_early_p <- data_p %>% filter(GA_Cat == "Early")
data_late_p <- data_p %>% filter(GA_Cat == "Late")
```

``` r
library(lme4)
```

    ## Loading required package: Matrix

    ## Warning: package 'Matrix' was built under R version 4.2.0

    ## 
    ## Attaching package: 'Matrix'

    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

``` r
library(lmerTest) 
```

    ## 
    ## Attaching package: 'lmerTest'

    ## The following object is masked from 'package:lme4':
    ## 
    ##     lmer

    ## The following object is masked from 'package:stats':
    ## 
    ##     step

``` r
#  # Enhances lmer() with p-values
model <- lmer((c2loc_EVT) ~ GA_Condition + (1 | donor2), data = data_early_p)
summary(model)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
    ## Formula: (c2loc_EVT) ~ GA_Condition + (1 | donor2)
    ##    Data: data_early_p
    ## 
    ## REML criterion at convergence: 4928.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8421 -0.7113 -0.0410  0.6631  8.0895 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  donor2   (Intercept) 0.0642   0.2534  
    ##  Residual             0.3866   0.6217  
    ## Number of obs: 2597, groups:  donor2, 6
    ## 
    ## Fixed effects:
    ##                      Estimate Std. Error      df t value Pr(>|t|)    
    ## (Intercept)            1.7381     0.1481  4.1143  11.738 0.000257 ***
    ## GA_ConditionEarly_PE  -0.6319     0.2089  4.0736  -3.025 0.038051 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## GA_CndtE_PE -0.709

``` r
# qqnorm(resid(model))
# qqline(resid(model))
```

``` r
# plot(fitted(model), resid(model))
```

### Late - PVBP

``` r
model <- lmer((c2loc_EVT) ~ GA_Condition + (1 | donor2), data = data_late_p)
summary(model) # includes p-values
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
    ## Formula: (c2loc_EVT) ~ GA_Condition + (1 | donor2)
    ##    Data: data_late_p
    ## 
    ## REML criterion at convergence: 9537.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.0919 -0.7726  0.0837  0.6469  4.7844 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  donor2   (Intercept) 0.1420   0.3768  
    ##  Residual             0.4517   0.6721  
    ## Number of obs: 4642, groups:  donor2, 11
    ## 
    ## Fixed effects:
    ##                      Estimate Std. Error        df t value Pr(>|t|)    
    ## (Intercept)          1.398604   0.154773  8.914041   9.037 8.81e-06 ***
    ## GA_ConditionLate_PE -0.001704   0.229632  8.923986  -0.007    0.994    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## GA_CndtL_PE -0.674

``` r
# qqnorm(resid(model))
# qqline(resid(model))
```

``` r
# plot(fitted(model), resid(model))
```

## Stats - EVT decidua at Myometrium

``` r
# histogram - is it normally distributed? No.
data_m <- filtered_data_m

hist <- ggplot(data_m, aes(x=c2loc_EVT, color=GA_Condition)) +
  geom_histogram(fill="white") 
# hist
```

``` r
# Assuming 'data' is a data frame with columns 'abundance_per_spot'  and 'condition'
data_early_m <- data_m %>% filter(GA_Cat == "Early")
data_late_m <- data_m %>% filter(GA_Cat == "Late")
```

### Early - Myometrium

``` r
model <- lmer((c2loc_EVT) ~ GA_Condition + (1 | donor2), data = data_early_m)
summary(model)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
    ## Formula: (c2loc_EVT) ~ GA_Condition + (1 | donor2)
    ##    Data: data_early_m
    ## 
    ## REML criterion at convergence: 4537.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.1886 -0.5399 -0.0529  0.0813  4.8251 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  donor2   (Intercept) 0.2398   0.4897  
    ##  Residual             0.5495   0.7413  
    ## Number of obs: 2012, groups:  donor2, 7
    ## 
    ## Fixed effects:
    ##                      Estimate Std. Error      df t value Pr(>|t|)   
    ## (Intercept)            2.2948     0.3490  5.0570   6.575  0.00117 **
    ## GA_ConditionEarly_PE  -1.7287     0.4125  5.0341  -4.191  0.00844 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## GA_CndtE_PE -0.846

``` r
# qqnorm(resid(model))
# qqline(resid(model))
```

``` r
# plot(fitted(model), resid(model))
```

### Late - Myometrium

``` r
model <- lmer((c2loc_EVT) ~ GA_Condition + (1 | donor2), data = data_late_m)
summary(model)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
    ## Formula: (c2loc_EVT) ~ GA_Condition + (1 | donor2)
    ##    Data: data_late_m
    ## 
    ## REML criterion at convergence: 4873.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.8298 -0.3312 -0.0253  0.0776  6.9216 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  donor2   (Intercept) 0.5529   0.7436  
    ##  Residual             0.2421   0.4921  
    ## Number of obs: 3393, groups:  donor2, 9
    ## 
    ## Fixed effects:
    ##                     Estimate Std. Error      df t value Pr(>|t|)  
    ## (Intercept)           0.7200     0.3721  6.9889   1.935   0.0943 .
    ## GA_ConditionLate_PE  -0.2880     0.4993  6.9897  -0.577   0.5821  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## GA_CndtL_PE -0.745

``` r
# qqnorm(resid(model))
# qqline(resid(model))
```

``` r
# plot(fitted(model), resid(model))
```

``` r
# save the plots
# ggsave(paste0(dirsave,"PVBP_decidua_EVT_per_spot_20250130.png"), plot_decidua_p_all, width=3,height=4, bg = "white",  units = 'in', dpi = 300)
# ggsave(paste0(dirsave,"Myometrium_decidua_EVT_per_spot_20250130.png"), plot_decidua_myo_all, width=3,height=4, bg = "white",  units = 'in', dpi = 300)
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
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] lmerTest_3.1-3     lme4_1.1-27.1      Matrix_1.5-4.1     ggsignif_0.6.3     ggpubr_0.4.0       cowplot_1.1.1     
    ##  [7] paletteer_1.6.0    magrittr_2.0.1     forcats_0.5.1      purrr_1.0.2        readr_2.1.1        tidyr_1.1.4       
    ## [13] tibble_3.1.6       tidyverse_1.3.1    stringr_1.4.0      ggplot2_3.3.5      knitr_1.45         patchwork_1.1.1   
    ## [19] SeuratObject_4.1.3 Seurat_4.2.1       dplyr_1.0.7       
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] readxl_1.3.1           backports_1.4.1        plyr_1.8.6             igraph_1.2.11          lazyeval_0.2.2        
    ##   [6] sp_1.5-1               splines_4.1.1          listenv_0.8.0          scattermore_0.7        digest_0.6.29         
    ##  [11] htmltools_0.5.8.1      fansi_1.0.2            tensor_1.5             cluster_2.1.2          ROCR_1.0-11           
    ##  [16] tzdb_0.2.0             globals_0.14.0         modelr_0.1.8           matrixStats_0.61.0     spatstat.sparse_3.0-0 
    ##  [21] colorspace_2.0-2       rvest_1.0.2            ggrepel_0.9.1          haven_2.4.3            xfun_0.41             
    ##  [26] prismatic_1.1.2        crayon_1.4.2           jsonlite_1.7.3         progressr_0.10.0       spatstat.data_3.0-0   
    ##  [31] survival_3.2-13        zoo_1.8-9              glue_1.6.1             polyclip_1.10-0        gtable_0.3.0          
    ##  [36] leiden_0.3.9           car_3.0-12             future.apply_1.8.1     abind_1.4-5            scales_1.1.1          
    ##  [41] DBI_1.1.2              rstatix_0.7.0          spatstat.random_3.0-1  miniUI_0.1.1.1         Rcpp_1.0.8            
    ##  [46] viridisLite_0.4.0      xtable_1.8-4           reticulate_1.24        htmlwidgets_1.5.4      httr_1.4.2            
    ##  [51] RColorBrewer_1.1-2     ellipsis_0.3.2         ica_1.0-2              farver_2.1.0           pkgconfig_2.0.3       
    ##  [56] uwot_0.1.14            dbplyr_2.1.1           deldir_1.0-6           utf8_1.2.2             labeling_0.4.2        
    ##  [61] tidyselect_1.1.1       rlang_1.1.1            reshape2_1.4.4         later_1.3.0            munsell_0.5.0         
    ##  [66] cellranger_1.1.0       tools_4.1.1            cli_3.6.1              generics_0.1.1         broom_0.7.11          
    ##  [71] ggridges_0.5.3         evaluate_0.23          fastmap_1.1.1          yaml_2.2.2             goftest_1.2-3         
    ##  [76] rematch2_2.1.2         fs_1.5.2               fitdistrplus_1.1-6     RANN_2.6.1             pbapply_1.5-0         
    ##  [81] future_1.23.0          nlme_3.1-155           mime_0.12              xml2_1.3.3             compiler_4.1.1        
    ##  [86] rstudioapi_0.13        plotly_4.10.0          png_0.1-7              spatstat.utils_3.0-1   reprex_2.0.1          
    ##  [91] stringi_1.7.6          highr_0.9              lattice_0.20-45        nloptr_1.2.2.3         vctrs_0.6.5           
    ##  [96] pillar_1.6.5           lifecycle_1.0.4        spatstat.geom_3.0-3    lmtest_0.9-39          RcppAnnoy_0.0.19      
    ## [101] data.table_1.14.2      irlba_2.3.5            httpuv_1.6.5           R6_2.5.1               promises_1.2.0.1      
    ## [106] KernSmooth_2.23-20     gridExtra_2.3          parallelly_1.30.0      codetools_0.2-18       boot_1.3-28           
    ## [111] MASS_7.3-55            assertthat_0.2.1       withr_2.5.0            sctransform_0.3.5      parallel_4.1.1        
    ## [116] hms_1.1.1              grid_4.1.1             minqa_1.2.4            rmarkdown_2.25         carData_3.0-5         
    ## [121] Rtsne_0.15             spatstat.explore_3.0-5 numDeriv_2016.8-1.1    shiny_1.7.1            lubridate_1.8.0


