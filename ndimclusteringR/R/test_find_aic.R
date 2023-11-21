#' Test the function`find_all_ic`
#'
#' @description Test the column names after using `find_all_ic`
#' on a dummy dataframe.
#'
#' @family tests
#'
#' @export
test_find_aic <- function() {
  print("... Testing aic columns")
  # The row names will get overridden in rbind so assign rownames to column
  # then reassign when ncents has been chosen
  nclust <- 5
  dummy_traits <- c("T1", "T2", "T3")
  dummy_snps <- c("rs35662", "rs301884", "rs69696",
                  "rs4096", "rs646464", "rs1234")
  num_axis <- length(dummy_traits)
  nsnps <- length(dummy_snps)
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
  clust_dist_df <- data.frame(row.names = dummy_snps)
  clust_re <- list("clusters" = cluster_df,
                   "centres" = centres_df,
                   "clust_dist" = clust_dist_df)
  #clust_re, num_axis
  expect_cols <- c("snp_id", colnames(cluster_df), "aic", "bic")
  clust_re <- find_all_ic(clust_re, num_axis)
  testit::assert("find_all_aic has found the wrong row names",
                 colnames(clust_re) == expect_cols)
}