# NDimClustering
 N-Dimensional MR Clustering

## Summary
This package is designed to cluster data on n-dimensional axis. The variable names correspond to clustering phenotype snp data but the clustering approach could be easily adapted to other problems. 

# Setup
Required packages:
* ggplot2
* tidyr
* tidyverse

### Package and Version details:
During development the following was used.
R version 4.3.1 (2023-06-16)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.3 LTS

Matrix products: default
BLAS:   libblas.so.3 
LAPACK: libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

attached base packages:
stats, graphics, grDevices, utils, datasets, methods, base     

other attached packages:
lubridate_1.9.2, forcats_1.0.0, stringr_1.5.0, dplyr_1.1.2, purrr_1.0.1    
readr_2.1.4, tibble_3.2.1, tidyverse_2.0.0, tidyr_1.3.0,ggplot2_3.4.2  

loaded via a namespace (and not attached):
gtable_0.3.3, compiler_4.3.1, tidyselect_1.2.0, systemfonts_1.0.4, scales_1.2.1     
textshaping_0.3.6, R6_2.5.1, labeling_0.4.2, generics_0.1.3, munsell_0.5.0    
pillar_1.9.0, tzdb_0.4.0, rlang_1.1.1, utf8_1.2.3, stringi_1.7.12, timechange_0.2.0
cli_3.6.1, withr_2.5.0, magrittr_2.0.3, grid_4.3.1, rstudioapi_0.14, hms_1.1.3         lifecycle_1.0.3, vctrs_0.6.3, glue_1.6.2, data.table_1.14.8, farver_2.1.1, ragg_1.2.5        fansi_1.0.4, colorspace_2.1-0, tools_4.3.1, pkgconfig_2.0.3  

# Tutorial
Tutorial on using some of the basic functions in `GettingStarted.qmd`

# Types of inputs
There are optional variations in which you can cluster with this package. This includes:
* `bp_on`: Weighting cluster contribution with a p-value assoicated with the initial scores.
* `clust_prob_on`: Weighting cluster contribution with the points distance to the cluster centre.
* `norm_typ`: The type of norm used in distance calculations.
* `clust_typ`: Types of clustering, either k-means or k-means minimising AIC. Also includes whether to transform the data prior to clustering to angles. Angle based clustering will find clusters that fit a line as apposed to balls.


