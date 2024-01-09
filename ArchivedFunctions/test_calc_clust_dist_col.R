#' Test `calc_clust_dist_col` the function which finds the distance
#' from a point to the centre of a cluster numbered c_num
#'
#' @description
#'  * Check output is dataframe
#'  * Check the number of rows is same data_mat
#'  * Check there is one column
#'  * Check the column name matches the cluster number label
#'
#' @family tests
#'
#' @export
test_calc_clust_dist_col <- function() {
  nclust <- 3
  c_num <- 2
  dummy_traits <- c("T1", "T2", "T3")
  dummy_snps <- c("rs35662", "rs301884", "rs69696",
                  "rs4096", "rs646464", "rs1234")
  expec_cols <- c(make_clust_col_name(c_num))
  num_axis <- length(dummy_traits)
  nsnps <- length(dummy_snps)
  b_mat <- matrix(runif(nsnps * num_axis, 0, 10), nrow = nsnps)
  colnames(b_mat) <- dummy_traits
  rownames(b_mat) <- dummy_snps
  centres_df <- data.frame(row.names = 1:nclust)
  for (i in 1:num_axis){
    centres_df[dummy_traits[i]] <- runif(nclust, 0, 10)
  }
  dist_col <- calc_clust_dist_col(c_num, b_mat, centres_df)
  testit::assert("... Output is not a dataframe",
                 is.data.frame(dist_col))
  testit::assert("... Wrong number of rows",
                 nrow(dist_col) == nrow(b_mat))
  testit::assert("... Wrong number of columns",
                 ncol(dist_col) == 1)
  testit::assert("... Wrong column label",
                 colnames(dist_col) == expec_cols)
}