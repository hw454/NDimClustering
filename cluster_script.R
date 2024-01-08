devtools::install("ndimclusteringR")
library("ndimclusteringR")

# Script for running the cluster program setting the program inputs
data_dir <- "./TestData/paths"
res_dir0 <- "./TestResults/"
if (!file.exists(res_dir0)){
  dir.create(res_dir0)
}
res_dir0 <- paste0(res_dir0, "paths")
for (i in 0:3){
  if (i > 1) {
    nc <- i
  } else { nc <- 1}
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
    point_eps = 0.4,
    how_cents = "points"
  )
  res_dir1 <- paste0(res_dir0, "paths", i)
  if (!file.exists(res_dir1)) {
    dir.create(res_dir1)
  }
  res_dir <- paste0(res_dir1, "/", make_path_label_str(iter_traits))
  iter_traits$res_dir <- res_dir
  clustering_program(iter_traits, test = 0)
}