#' Test `calc_member_dist_cent` the function which finds the distance
#' from a point to the centre of a cluster numbered c_num
#'
#' @description
#'  * Check output is dataframe
#'  * Check the number of rows is same data_mat
#'  * Check there is a column for each cluster
#'  * Check the column names matches the cluster number labels
#'
#' @family tests
#'
#' @export
test_calc_member_dist_cent <- function() {
  nclust <- 3
  snp_id <- "rs35662"
  dummy_traits <- c("T1", "T2", "T3")
  dummy_snps <- c(snp_id, "rs301884", "rs69696",
                  "rs4096", "rs646464", "rs1234")
  expec_cols <- c("clust_dist")
  num_axis <- length(dummy_traits)
  nsnps <- length(dummy_snps)
  b_mat <- matrix(runif(nsnps * num_axis, 0, 10), nrow = nsnps)
  colnames(b_mat) <- dummy_traits
  rownames(b_mat) <- dummy_snps
  centres_df <- data.frame(row.names = 1:nclust)
  for (i in 1:num_axis){
    centres_df[dummy_traits[i]] <- runif(nclust, 0, 10)
  }
  c_dist_rnd <- runif(nsnps, 0, 10)
  c_prob_rnd <- runif(nsnps, 0, 1)
  c_num_rnd <- sample(nsnps, 1, nclust)
  cluster_df <- data.frame(
    row.names = dummy_snps,
    clust_num = c_num_rnd,
    clust_dist = c_dist_rnd,
    clust_prob = c_prob_rnd
  )
  dist_df <- calc_member_dist_cent(snp_id, b_mat, cluster_df, centres_df)
  testit::assert("... Output is not a dataframe",
                 is.data.frame(dist_df))
  testit::assert("... Wrong number of rows",
                 nrow(dist_df) == 1)
  testit::assert("... Wrong number of columns",
                 ncol(dist_df) == 1)
  testit::assert("... Wrong column label",
                 colnames(dist_df) == expec_cols)
}