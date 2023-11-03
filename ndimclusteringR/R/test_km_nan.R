#' Test kmeans clustering ignoring NaN terms
#'
#' @description
#' Use common axes for distance calculations to igore NaNs and
#' retain data.
#'
#' @details
#' The cluster centres are randomly assigned using [make_rand_cent].
#' The cluster centres are check and reassigned using [check_clust_cent],
#' When the clusters are converged or the maximum number of iterations are
#' used then the clusters are return in the dataframe "cluster_df".
#' "clust_out" is the list containing the "clusters_df" dataframe
#' labelled "cluster" of the cluster membership with columns\:
#'   * "clust_num" the number of clusters
#'   * "clust_dist" the distance from the snp to the cluster centre
#'   (or distance between angles)
#'   * "clust_prob" probability the snp is in the cluster. Calculated
#'   using [calc_clust_prob].
#' "clust_out" also contains the "centroids_df" dataframe labelled "centres"
#' whose columns are the traits and rows are the cluster numbers. Additional
#' column "thresh_check" contains Bool indicating whether that cluster
#' converged.
#'
#' @return
#'
#' @family tests
#'
#' @export
test_km_nan <- function() {
  dummy_traits <- c("T1", "T2", "T3")
  dummy_snps <- c("rs35662", "rs301884", "rs69696",
                     "rs4096", "rs646464", "rs1234")
  num_axis <- length(dummy_traits)
  nsnps <- length(dummy_snps)
  b_mat <- matrix(runif(nsnps * num_axis, 0, 10), nrow = nsnps)
  colnames(b_mat) <- dummy_traits
  rownames(b_mat) <- dummy_snps
  nclust <- 5
  iter_max <- 300
  clust_threshold <- 1e-5
  norm_typ <- "F"
  na_rm <- TRUE
  cluster_df <- km_nan(b_mat,
                       nclust = nclust,
                       iter_max = iter_max,
                       clust_threshold = clust_threshold,
                       norm_typ = norm_typ,
                       na_rm = na_rm)
  print("... test with full matrix")
  expect_cols <- c("clusters", "centres", "clust_dist")
  testit::assert("... not dataframe clusters",
    is.data.frame(cluster_df$clusters))
  testit::assert("... not dataframe centres",
    is.data.frame(cluster_df$centres))
  testit::assert("... not dataframe clust_dist",
    is.data.frame(cluster_df$clust_dist))
  testit::assert("...km_nan has found the wrong list names",
    names(cluster_df) == expect_cols)
  testit::assert("...km_nan has found the wrong clust numbers",
    unique(cluster_df$clusters$clust_num) %in% rownames(cluster_df$centres))
  testit::assert("...km_nan has found the wrong centres column names",
    colnames(cluster_df$centres) == dummy_traits)
  testit::assert("...km_nan has found the wrong clusters row names",
    rownames(cluster_df$clusters) == dummy_snps)
  testit::assert("...km_nan has found the wrong clust_dist row names",
    rownames(cluster_df$clust_dist) == dummy_snps)
  nan_rows <- runif(4, 1, nsnps)
  nan_cols <- runif(4, 1, num_axis)
  for (i in 1:4){
    b_mat[nan_rows[i], nan_cols[i]] <- NaN
  }
  cluster_df <- km_nan(b_mat,
                       nclust = nclust,
                       iter_max = iter_max,
                       clust_threshold = clust_threshold,
                       norm_typ = norm_typ,
                       na_rm = na_rm)
  print("... test with NaN terms")
  testit::assert("...km_nan has found the wrong list names",
    names(cluster_df) == expect_cols)
  testit::assert("...km_nan has found the wrong clust numbers",
    unique(cluster_df$clusters$clust_num) %in% rownames(cluster_df$centres))
  testit::assert("...km_nan has found the wrong centres column names",
    colnames(cluster_df$centres) == dummy_traits)
  testit::assert("...km_nan has found the wrong clusters row names",
    rownames(cluster_df$clusters) == dummy_snps)
  testit::assert("...km_nan has found the wrong clust_dist row names",
    rownames(cluster_df$clust_dist) == dummy_snps)
  b_mat["T2"] <- NaN
  cluster_df <- km_nan(b_mat,
                       nclust = nclust,
                       iter_max = iter_max,
                       clust_threshold = clust_threshold,
                       norm_typ = norm_typ,
                       na_rm = na_rm)
  print("... test with NaN column")
  testit::assert("...km_nan has found the wrong list names",
    names(cluster_df) == expect_cols)
  testit::assert("...km_nan has found the wrong clust numbers",
    unique(cluster_df$clusters$clust_num) %in% rownames(cluster_df$centres))
  testit::assert("...km_nan has found the wrong centres column names",
    colnames(cluster_df$centres) == dummy_traits)
  testit::assert("...km_nan has found the wrong clusters row names",
    rownames(cluster_df$clusters) == dummy_snps)
  testit::assert("...km_nan has found the wrong clust_dist row names",
    rownames(cluster_df$clust_dist) == dummy_snps)
  return()
}