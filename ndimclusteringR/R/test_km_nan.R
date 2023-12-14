#' Test the function for finding k clusters.
#'
#' @details Checks:
#' * there are k clusters found
#' * The number of points with clusters assigned is the number of points input.
#'
#' @export
#' @family tests
test_km_nan <- function() {
  make_clust_col <- function(i) {
    return(paste0("clust_", i))
  }
  dummy_traits <- c("T1", "T2", "T3", "T4", "T5", "T6")
  dummy_snps <- c("rs35662", "rs301884", "rs69696",
                  "rs4096", "rs646464", "rs1234")
  nclust <- 5
  num_axis <- length(dummy_traits)
  nsnps <- length(dummy_snps)
  print("Create mat")
  b_mat <- matrix(stats::runif(nsnps * num_axis, 0, 10), nrow = nsnps)
  colnames(b_mat) <- dummy_traits
  rownames(b_mat) <- dummy_snps
  cat_list <- rep("Exposure", ncol(b_mat))
  cat_list[1] <- "Outcome"
  print("km_nan")
  clust_out <- km_nan(b_mat, nclust = nclust)
  print(clust_out)
  expec_list <- c("clusters", "clust_dist", "centres")
  testit::assert("km_nan doesn't contain the right terms",
    all(expec_list %in% names(clust_out))
  )
  expec_ncents <- nclust
  expec_ncols <- ncol(b_mat)
  expec_nrows <- nrow(b_mat)
  expec_clust_cols <- c(make_clust_col(1),
    make_clust_col(2),
    make_clust_col(3),
    make_clust_col(4),
    make_clust_col(5)
  )
  expec_cent_rows <- seq_len(5)
  expec_clusterdf_cols <- c("clust_dist", "clust_num", "clust_prob", "ncents")
  testit::assert("Centres are wrong dimension",
    ncol(clust_out$centres) == expec_ncols
  )
  testit::assert("Wrong number of centres",
    nrow(clust_out$centres) == expec_ncents
  )
  testit::assert("Centres have wrong axis",
    all(colnames(clust_out$centres) %in% colnames(b_mat))
  )
  testit::assert("Centres has wrong rows",
    all(rownames(clust_out$centres) %in% expec_cent_rows)
  )
  testit::assert("Centres are wrong dimension",
    ncol(clust_out$clust_dist) == nclust
  )
  testit::assert("Clust_dist has wrong number of rows",
    nrow(clust_out$clust_dist) == expec_nrows
  )
  testit::assert("Clust_dist has wrong axis",
    all(colnames(clust_out$clust_dist) %in% expec_clust_cols)
  )
  testit::assert("Clust_dist has wrong rows",
    all(rownames(clust_out$clust_dist) %in% rownames(b_mat))
  )
  testit::assert("Cluster_df has wrong number of rows",
    nrow(clust_out$clusters) == expec_nrows
  )
  testit::assert("Cluster_df has wrong axis",
    all(colnames(clust_out$clusters) %in% expec_clusterdf_cols)
  )
  testit::assert("Clust_dist has wrong rows",
    all(rownames(clust_out$clusters) %in% rownames(b_mat))
  )
}