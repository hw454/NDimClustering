devtools::install("ndimclusteringR")
library("ndimclusteringR")

# Script for running the cluster program setting the program inputs
data_dir <- "./TestData/paths"
for (i in 0:4){
  iter_traits <- list(
    dname = paste0(data_dir, i, "/"),
    clust_type = "basic",
    pca_type = "prcomp",
    n_clust = 3,
    n_pcs = 3,
    bin_angles = 1,
    bin_p_clust = 1,
    bin_p_score = 1,
    bin_d_score = 1,
    nan_rm = 1
  )
  clustering_program(iter_traits, test = 1)
}