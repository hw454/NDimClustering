---
title: Testing Clustering methods.
author: Hayley Wragg
date: "21st Nov 2023"
output:
  html_document:
    keep_md: yes
  pdf_document: default
editor_options: 
  markdown: 
    wrap: sentence
---

# Install package

Load the "ndimclusteringR" package


```r
devtools::install("../ndimclusteringR/R")
```

```
## 
## ── R CMD build ─────────────────────────────────────────────────────────────────
##   
   checking for file ‘/home/hayley/Code/ICEP/NDimClustering/ndimclusteringR/DESCRIPTION’ ...
  
✔  checking for file ‘/home/hayley/Code/ICEP/NDimClustering/ndimclusteringR/DESCRIPTION’
## 
  
─  preparing ‘ndimclusteringR’:
##    checking DESCRIPTION meta-information ...
  
✔  checking DESCRIPTION meta-information
## 
  
─  checking for LF line-endings in source and make files and shell scripts
## 
  
─  checking for empty or unneeded directories
## 
  
─  building ‘ndimclusteringR_0.0.1.tar.gz’
## 
  
   
## 
Running /usr/lib/R/bin/R CMD INSTALL \
##   /tmp/RtmplYkigx/ndimclusteringR_0.0.1.tar.gz --install-tests 
## * installing to library ‘/home/hayley/R/x86_64-pc-linux-gnu-library/4.3’
## * installing *source* package ‘ndimclusteringR’ ...
## ** using staged installation
## ** R
## ** data
## *** moving datasets to lazyload DB
## ** byte-compile and prepare package for lazy loading
## ** help
## *** installing help indices
## ** building package indices
## ** testing if installed package can be loaded from temporary location
## ** testing if installed package can be loaded from final location
## ** testing if installed package keeps a record of temporary installation path
## * DONE (ndimclusteringR)
```

```r
library("ndimclusteringR")
```

# Source test

Load the functions for testing


```r
source("test_clust_kmeans.R")
```

# Step by step

If all the steps are included then they are performed in the following order.

1\.
Crop data to complete cases and lowest se.

2\.
Convert data to angles.

3\.
Perform PCA on angle data.

4\.
Cluster on the PCA

-   Set cluster centres
-   Converge cluster centres
-   if min then rerun for different number of clusters and choose set with min-aic.

### Test 1

Run the tests with:

-   No PCA
-   basic k-means clustering. k is set to the number of pathways.
-   no angles calculated before clustering.
-   number of pathways 2.


```r
pc_type <- "No_pca"
np <- 1
d <- 30
clust_typ <- "basic"
space_typ <- "regular"
test_clust_kmeans_function(d,
                           num_path = np,
                           pc_type = pc_type,
                           clust_typ = clust_typ,
                           space_typ = space_typ)
```

```
## [1] "Test No_pca with 2 pathways"
##   bp_on clust_prob_on         clust_typ ndim_typ how_cents pc_type num_paths
## 1  TRUE          TRUE test_basicregular test_all      rand  No_pca         1
##           res_dir
## 1 PC_TestResults/
## [1] "Clusters converged 4"
```

![](test_clust_kmeans_files/figure-html/unnamed-chunk-3-1.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-3-2.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-3-3.png)<!-- -->

### Test 2

Run the tests with:

-   No PCA
-   basic k-means clustering. k is set to the number of pathways.
-   no angles calculated before clustering.
-   number of pathways 4.


```r
pc_type <- "No_pca"
np <- 3
d <- 30
clust_typ <- "basic"
space_typ <- "regular"
test_clust_kmeans_function(d,
                           num_path = np,
                           pc_type = pc_type,
                           clust_typ = clust_typ,
                           space_typ = space_typ)
```

```
## [1] "Test No_pca with 4 pathways"
##   bp_on clust_prob_on         clust_typ ndim_typ how_cents pc_type num_paths
## 1  TRUE          TRUE test_basicregular test_all      rand  No_pca         3
##           res_dir
## 1 PC_TestResults/
## [1] "Clusters converged 6"
```

![](test_clust_kmeans_files/figure-html/unnamed-chunk-4-1.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-4-2.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-4-3.png)<!-- -->

### Test 3

Run the tests with:

-   No PCA
-   basic k-means clustering. k is set to the number of pathways.
-   angles calculated before pca and before clustering.
-   number of pathways 2.


```r
pc_type <- "No_pca"
np <- 1
d <- 30
clust_typ <- "basic"
space_typ <- "angle"
test_clust_kmeans_function(d,
                           num_path = np,
                           pc_type = pc_type,
                           clust_typ = clust_typ,
                           space_typ = space_typ)
```

```
## [1] "Test No_pca with 2 pathways"
##   bp_on clust_prob_on       clust_typ ndim_typ how_cents pc_type num_paths
## 1  TRUE          TRUE test_basicangle test_all      rand  No_pca         1
##           res_dir
## 1 PC_TestResults/
## [1] "Clusters converged 4"
```

![](test_clust_kmeans_files/figure-html/unnamed-chunk-5-1.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-5-2.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-5-3.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-5-4.png)<!-- -->

### Test 4

Run the tests with:

-   No PCA
-   basic k-means clustering. k is set to the number of pathways.
-   angles calculated before clustering.
-   number of pathways 4.


```r
pc_type <- "No_pca"
np <- 3
d <- 30
clust_typ <- "basic"
space_typ <- "angle"
test_clust_kmeans_function(d,
                           num_path = np,
                           pc_type = pc_type,
                           clust_typ = clust_typ,
                           space_typ = space_typ)
```

```
## [1] "Test No_pca with 4 pathways"
##   bp_on clust_prob_on       clust_typ ndim_typ how_cents pc_type num_paths
## 1  TRUE          TRUE test_basicangle test_all      rand  No_pca         3
##           res_dir
## 1 PC_TestResults/
## [1] "Clusters converged 3"
```

![](test_clust_kmeans_files/figure-html/unnamed-chunk-6-1.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-6-2.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-6-3.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-6-4.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-6-5.png)<!-- -->

### Test 5

Run the tests with:

-   No PCA
-   min-aic k-means clustering. Maximum k is set to the number of pathways.
-   no angles calculated before clustering.
-   number of pathways 2.


```r
pc_type <- "No_pca"
np <- 1
d <- 30
clust_typ <- "min"
space_typ <- "regular"
test_clust_kmeans_function(d,
                           num_path = np,
                           pc_type = pc_type,
                           clust_typ = clust_typ,
                           space_typ = space_typ)
```

```
## [1] "Test No_pca with 2 pathways"
##   bp_on clust_prob_on       clust_typ ndim_typ how_cents pc_type num_paths
## 1  TRUE          TRUE test_minregular test_all      rand  No_pca         1
##           res_dir
## 1 PC_TestResults/
## [1] "Clusters converged 2"
## [1] "Clusters converged 2"
## [1] "Clusters converged 3"
```

![](test_clust_kmeans_files/figure-html/unnamed-chunk-7-1.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-7-2.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-7-3.png)<!-- -->

### Test 6

Run the tests with:

-   No PCA
-   min-aic k-means clustering. Maximum k is set to the number of pathways.
-   no angles calculated before clustering.
-   number of pathways 4.


```r
pc_type <- "No_pca"
np <- 3
d <- 30
clust_typ <- "min"
space_typ <- "regular"
test_clust_kmeans_function(d,
                           num_path = np,
                           pc_type = pc_type,
                           clust_typ = clust_typ,
                           space_typ = space_typ)
```

```
## [1] "Test No_pca with 4 pathways"
##   bp_on clust_prob_on       clust_typ ndim_typ how_cents pc_type num_paths
## 1  TRUE          TRUE test_minregular test_all      rand  No_pca         3
##           res_dir
## 1 PC_TestResults/
## [1] "Clusters converged 2"
## [1] "Clusters converged 4"
## [1] "Clusters converged 3"
## [1] "Clusters converged 6"
## [1] "Clusters converged 5"
```

![](test_clust_kmeans_files/figure-html/unnamed-chunk-8-1.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-8-2.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-8-3.png)<!-- -->

### Test 7

Run the tests with:

-   No PCA
-   min k-means clustering. Maximum k is set to the number of pathways.
-   angles calculated before clustering.
-   number of pathways 2.


```r
pc_type <- "No_pca"
np <- 1
d <- 30
clust_typ <- "min"
space_typ <- "angle"
test_clust_kmeans_function(d,
                           num_path = np,
                           pc_type = pc_type,
                           clust_typ = clust_typ,
                           space_typ = space_typ)
```

```
## [1] "Test No_pca with 2 pathways"
##   bp_on clust_prob_on     clust_typ ndim_typ how_cents pc_type num_paths
## 1  TRUE          TRUE test_minangle test_all      rand  No_pca         1
##           res_dir
## 1 PC_TestResults/
## [1] "Clusters converged 2"
## [1] "Clusters converged 4"
## [1] "Clusters converged 3"
```

![](test_clust_kmeans_files/figure-html/unnamed-chunk-9-1.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-9-2.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-9-3.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-9-4.png)<!-- -->

### Test 8

Run the tests with:

-   No PCA
-   min-aic k-means clustering. Maximum k is set to the number of pathways.
-   angles calculated before clustering.
-   number of pathways 4.


```r
pc_type <- "No_pca"
np <- 3
d <- 30
clust_typ <- "min"
space_typ <- "angle"
test_clust_kmeans_function(d,
                           num_path = np,
                           pc_type = pc_type,
                           clust_typ = clust_typ,
                           space_typ = space_typ)
```

```
## [1] "Test No_pca with 4 pathways"
##   bp_on clust_prob_on     clust_typ ndim_typ how_cents pc_type num_paths
## 1  TRUE          TRUE test_minangle test_all      rand  No_pca         3
##           res_dir
## 1 PC_TestResults/
## [1] "Clusters converged 2"
## [1] "Clusters converged 4"
## [1] "Clusters converged 3"
## [1] "Clusters converged 6"
## [1] "Clusters converged 5"
```

![](test_clust_kmeans_files/figure-html/unnamed-chunk-10-1.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-10-2.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-10-3.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-10-4.png)<!-- -->

## Consider the location of the centroids.

The previous tests used centroids assigned randomly within the ranges on each axis.
Instead we will initialise the centroids using randomly selected points from out dataset.

### Test 9

Run the tests with:

-   No PCA
-   basic k-means clustering. k is set to the number of pathways.
-   no angles calculated before clustering.
-   number of pathways 4.
-   Centroids assigned using points


```r
pc_type <- "No_pca"
np <- 3
d <- 30
clust_typ <- "basic"
space_typ <- "regular"
how_cents <- "points"
test_clust_kmeans_function(d,
                           num_path = np,
                           pc_type = pc_type,
                           clust_typ = clust_typ,
                           space_typ = space_typ,
                           how_cents = how_cents)
```

```
## [1] "Test No_pca with 4 pathways"
##   bp_on clust_prob_on         clust_typ ndim_typ how_cents pc_type num_paths
## 1  TRUE          TRUE test_basicregular test_all    points  No_pca         3
##           res_dir
## 1 PC_TestResults/
## [1] "Clusters converged 6"
```

![](test_clust_kmeans_files/figure-html/unnamed-chunk-11-1.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-11-2.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-11-3.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-11-4.png)<!-- -->

### Test 10

Run the tests with: - No PCA - basic k-means clustering.
k is set to the number of pathways.
- angles calculated before clustering.
- number of pathways 4.
- Centroids assigned using points


```r
pc_type <- "No_pca"
np <- 3
d <- 30
clust_typ <- "basic"
space_typ <- "angle"
how_cents <- "points"
test_clust_kmeans_function(d,
                           num_path = np,
                           pc_type = pc_type,
                           clust_typ = clust_typ,
                           space_typ = space_typ,
                           how_cents = how_cents)
```

```
## [1] "Test No_pca with 4 pathways"
##   bp_on clust_prob_on       clust_typ ndim_typ how_cents pc_type num_paths
## 1  TRUE          TRUE test_basicangle test_all    points  No_pca         3
##           res_dir
## 1 PC_TestResults/
## [1] "Clusters converged 2"
```

![](test_clust_kmeans_files/figure-html/unnamed-chunk-12-1.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-12-2.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-12-3.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-12-4.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-12-5.png)<!-- -->

### Test 11

Run the tests with:

-   PCA with prcomp
-   basic k-means clustering. k is set to the number of pathways.
-   no angles calculated before clustering.
-   number of pathways 4.
-   Centroids assigned using points


```r
pc_type <- "prcomp"
np <- 3
d <- 30
clust_typ <- "basic"
space_typ <- "regular"
how_cents <- "points"
test_clust_kmeans_function(d,
                           num_path = np,
                           pc_type = pc_type,
                           clust_typ = clust_typ,
                           space_typ = space_typ,
                           how_cents = how_cents)
```

```
## [1] "Test prcomp with 4 pathways"
##   bp_on clust_prob_on         clust_typ ndim_typ how_cents pc_type num_paths
## 1  TRUE          TRUE test_basicregular test_all    points  prcomp         3
##           res_dir
## 1 PC_TestResults/
## [1] "Clusters converged 6"
```

![](test_clust_kmeans_files/figure-html/unnamed-chunk-13-1.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-13-2.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-13-3.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-13-4.png)<!-- -->

### Test 12

Run the tests with:

-   PCA with prcomp
-   basic k-means clustering. k is set to the number of pathways.
-   angles calculated before PCA and clustering.
-   number of pathways 4.
-   Centroids assigned using points


```r
pc_type <- "prcomp"
np <- 3
d <- 30
clust_typ <- "basic"
space_typ <- "angle"
how_cents <- "points"
test_clust_kmeans_function(d,
                           num_path = np,
                           pc_type = pc_type,
                           clust_typ = clust_typ,
                           space_typ = space_typ,
                           how_cents = how_cents)
```

```
## [1] "Test prcomp with 4 pathways"
##   bp_on clust_prob_on       clust_typ ndim_typ how_cents pc_type num_paths
## 1  TRUE          TRUE test_basicangle test_all    points  prcomp         3
##           res_dir
## 1 PC_TestResults/
## [1] "Clusters converged 5"
```

![](test_clust_kmeans_files/figure-html/unnamed-chunk-14-1.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-14-2.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-14-3.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-14-4.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-14-5.png)<!-- -->

### Test 13

Run the tests with:

-   No PCA
-   min-aic k-means clustering. k is set to the number of pathways.
-   no angles calculated before clustering.
-   number of pathways 4.
-   Centroids assigned using points


```r
pc_type <- "No_pca"
np <- 3
d <- 30
clust_typ <- "min"
space_typ <- "regular"
how_cents <- "points"
test_clust_kmeans_function(d,
                           num_path = np,
                           pc_type = pc_type,
                           clust_typ = clust_typ,
                           space_typ = space_typ,
                           how_cents = how_cents)
```

```
## [1] "Test No_pca with 4 pathways"
##   bp_on clust_prob_on       clust_typ ndim_typ how_cents pc_type num_paths
## 1  TRUE          TRUE test_minregular test_all    points  No_pca         3
##           res_dir
## 1 PC_TestResults/
## [1] "Clusters converged 2"
## [1] "Clusters converged 2"
## [1] "Clusters converged 4"
## [1] "Clusters converged 3"
## [1] "Clusters converged 4"
```

![](test_clust_kmeans_files/figure-html/unnamed-chunk-15-1.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-15-2.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-15-3.png)<!-- -->

### Test 14

Run the tests with:

-   No PCA
-   min-aic k-means clustering. k is set to the number of pathways.
-   angles calculated before clustering.
-   number of pathways 4.
-   Centroids assigned using points


```r
pc_type <- "No_pca"
np <- 3
d <- 30
clust_typ <- "min"
space_typ <- "angle"
how_cents <- "points"
test_clust_kmeans_function(d,
                           num_path = np,
                           pc_type = pc_type,
                           clust_typ = clust_typ,
                           space_typ = space_typ,
                           how_cents = how_cents)
```

```
## [1] "Test No_pca with 4 pathways"
##   bp_on clust_prob_on     clust_typ ndim_typ how_cents pc_type num_paths
## 1  TRUE          TRUE test_minangle test_all    points  No_pca         3
##           res_dir
## 1 PC_TestResults/
## [1] "Clusters converged 2"
## [1] "Clusters converged 7"
## [1] "Clusters converged 3"
## [1] "Clusters converged 3"
## [1] "Clusters converged 5"
```

![](test_clust_kmeans_files/figure-html/unnamed-chunk-16-1.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-16-2.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-16-3.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-16-4.png)<!-- -->

### Test 15

Run the tests with:

-   PCA with prcomp
-   min-aic k-means clustering. k is set to the number of pathways.
-   no angles calculated before clustering.
-   number of pathways 4.
-   Centroids assigned using points


```r
pc_type <- "prcomp"
np <- 3
d <- 30
clust_typ <- "min"
space_typ <- "regular"
how_cents <- "points"
test_clust_kmeans_function(d,
                           num_path = np,
                           pc_type = pc_type,
                           clust_typ = clust_typ,
                           space_typ = space_typ,
                           how_cents = how_cents)
```

```
## [1] "Test prcomp with 4 pathways"
##   bp_on clust_prob_on       clust_typ ndim_typ how_cents pc_type num_paths
## 1  TRUE          TRUE test_minregular test_all    points  prcomp         3
##           res_dir
## 1 PC_TestResults/
## [1] "Clusters converged 2"
## [1] "Clusters converged 6"
## [1] "Clusters converged 4"
## [1] "Clusters converged 5"
## [1] "Clusters converged 3"
```

![](test_clust_kmeans_files/figure-html/unnamed-chunk-17-1.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-17-2.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-17-3.png)<!-- -->

### Test 16

Run the tests with:

-   PCA with prcomp
-   min-aic k-means clustering. k is set to the number of pathways.
-   angles calculated before PCA and clustering.
-   number of pathways 4.
-   Centroids assigned using points


```r
pc_type <- "prcomp"
np <- 3
d <- 30
clust_typ <- "min"
space_typ <- "angle"
how_cents <- "points"
test_clust_kmeans_function(d,
                           num_path = np,
                           pc_type = pc_type,
                           clust_typ = clust_typ,
                           space_typ = space_typ,
                           how_cents = how_cents)
```

```
## [1] "Test prcomp with 4 pathways"
##   bp_on clust_prob_on     clust_typ ndim_typ how_cents pc_type num_paths
## 1  TRUE          TRUE test_minangle test_all    points  prcomp         3
##           res_dir
## 1 PC_TestResults/
## [1] "Clusters converged 2"
## [1] "Clusters converged 6"
## [1] "Clusters converged 4"
## [1] "Clusters converged 5"
## [1] "Clusters converged 3"
```

![](test_clust_kmeans_files/figure-html/unnamed-chunk-18-1.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-18-2.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-18-3.png)<!-- -->![](test_clust_kmeans_files/figure-html/unnamed-chunk-18-4.png)<!-- -->
