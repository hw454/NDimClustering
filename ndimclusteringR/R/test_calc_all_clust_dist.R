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
test_calc_all_clust_dist <- function() {
  make_clust_col_name <- function(i) {
    return(paste0("clust_", i))
  }
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
  centres_df <- data.frame(row.names = 1:nclust)
  for (i in 1:num_axis){
    centres_df[dummy_traits[i]] <- runif(nclust, 0, 10)
  }
  # Calculate the cluster distances
  clust_dist_df <- calc_all_clust_dist(b_mat, centres_df)
  c_1 <- make_clust_col_name(1)
  c_2 <- make_clust_col_name(2)
  c_3 <- make_clust_col_name(3)
  c_4 <- make_clust_col_name(4)
  c_5 <- make_clust_col_name(5)
  expect_cols_dist <- c(c_1, c_2, c_3, c_4, c_5)
  expec_rows <- rownames(b_mat)
  testit::assert("clust_dist_df wrong type",
                 is.data.frame(clust_dist_df))
  testit::assert("wrong number of columns clust_dist_df",
                 ncol(clust_dist_df) == length(expect_cols_dist))
  testit::assert("wrong number of rows clust_dist_df",
                 nrow(clust_dist_df) == nrow(b_mat))
  testit::assert("wrong labels clust_dist_df",
                 all(names(clust_dist_df) == expect_cols_dist))
  testit::assert("wrong rows clust_dist_df",
                 all(rownames(clust_dist_df) == expec_rows))
}