#' Find the cluster whose centroid is closest to snp_id
#'
#' @description
#'   Test the function which finds the distance between the snp and
#'   each cluster centroid. `find_closest_clust_snp`
#'
#'
#' @details
#'   The dataframe "snp_clust_df"  has a row labelled by the "snp_id".
#'   The columns are\:
#'   * "clust_num" The number of the closest cluster
#'   * "clust_dist" The distance from the snp to the centroid of that cluster.
#'   * "clust_prob" probability the snp is in the cluster it's been assigned.
#'   The dataframe "snp_clust_dist_df"  has a row labelled by the "snp_id".
#'   The columns are the clust numbers, each value is the distance from the
#'   snp to that cluster.
#'
#'   out_list = ("clusters" = snp_clust_df, "clust_dist" = snp_clust_dist_df)
#'
#' @return out_list
#'
#' @export
test_find_closest_clust_snp <- function() {
  # DUMMY INPUTS
  nclust <- 5
  dummy_traits <- c("T1", "T2", "T3")
  dummy_snps <- c("rs35662", "rs301884", "rs69696",
                     "rs4096", "rs646464", "rs1234")
  num_axis <- length(dummy_traits)
  nsnps <- length(dummy_snps)
  b_mat <- matrix(runif(nsnps * num_axis, 0, 10), nrow = nsnps)
  colnames(b_mat) <- dummy_traits
  rownames(b_mat) <- dummy_snps
  c_dist_rnd <- runif(nsnps, 0, 10)
  c_prob_rnd <- runif(nsnps, 0, 1)
  c_num_rnd <- sample(nsnps, 1, nclust)
  cluster_df <- data.frame(
    row.names = dummy_snps,
    clust_num = c_num_rnd,
    clust_dist = c_dist_rnd,
    clust_prob = c_prob_rnd
  )
  centres_df <- data.frame(row.names = 1:nclust)
  for (i in 1:num_axis){
    centres_df[dummy_traits[i]] <- runif(nclust, 0, 10)
  }
  snp_clust_list <- lapply(dummy_snps, find_closest_clust_snp,
                            b_mat = b_mat,
                            cluster_df = cluster_df,
                            centroids_df = centres_df,
                            )
  # Combine the list of dataframes into one dataframe.
  # Override Cluster_df with the new assignment
  cluster_df_list <- lapply(snp_clust_list,
                          df_cols,
                          col = "clusters")
  cluster_df <- Reduce(rbind, cluster_df_list)
  expect_cols_clust <- c("clust_num", "clust_dist", "clust_prob")
  clust_dist_df_list <- lapply(snp_clust_list,
                          df_cols,
                          col = "clust_dist")
  clust_dist_df <- Reduce(rbind, clust_dist_df_list)
  c_1 <- make_clust_col_name(1)
  c_2 <- make_clust_col_name(2)
  c_3 <- make_clust_col_name(3)
  c_4 <- make_clust_col_name(4)
  c_5 <- make_clust_col_name(5)
  expect_cols_dist <- c(c_1, c_2, c_3, c_4, c_5)
  testit::assert("find_closest_clust_snp wrong labels cluster_df",
    all(names(cluster_df) == expect_cols_clust))
  testit::assert("find_closest_clust_snp wrong labels clust_dist_df",
    all(names(clust_dist_df) == expect_cols_dist))
  testit::assert("find_closest_clust_snp wrong type",
    is.data.frame(cluster_df))
}
df_cols <- function(df, col) {
  testit::assert("not dataframe in function",
  is.data.frame(df[[col]]))
  return(df[[col]])
}