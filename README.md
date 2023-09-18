# NDimClustering
 N-Dimensional MR Clustering

## Summary
This package is designed to cluster data on n-dimensional axis. The variable names correspond to clustering phenotype snp data but the clustering approach could be easily adapted to other problems. 

# Setup
Required packages:
* ggplot2
* tidyr
* tidyverse

# Tutorial
Tutorial on using some of the basic functions in `GettingStarted.qmd`

# Types of inputs
There are optional variations in which you can cluster with this package. This includes:
* `bp_on`: Weighting cluster contribution with a p-value assoicated with the initial scores.
* `clust_prob_on`: Weighting cluster contribution with the points distance to the cluster centre.
* `norm_typ`: The type of norm used in distance calculations.
* `clust_typ`: Types of clustering, either k-means or k-means minimising AIC. Also includes whether to transform the data prior to clustering to angles. Angle based clustering will find clusters that fit a line as apposed to balls.


