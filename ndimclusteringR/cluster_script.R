devtools::install("ndimclusteringR")
library("ndimclusteringR")

# Script for running the cluster program setting the program inputs
data_dir <- "./TestData/paths"
res_dir <- "./TestResults"
for (i in 0:3){
  if (i > 1) {
    nc <- i
  }
  iter_traits <- list(
    dname = paste0(data_dir, i, "/"),
    clust_type = "basic",
    pca_type = "prcomp",
    nclust = nc,
    n_pcs = 3,
    bin_angles = 1,
    bin_p_clust = 1,
    bin_p_score = 1,
    bin_d_score = 1,
    nan_rm = 1,
    how_cents = "points"
  )
  res_dir <- math_path_label_str(iter_traits)
  iter_traits$res_dir <- res_dir
  clustering_program(iter_traits, test = 1)
}