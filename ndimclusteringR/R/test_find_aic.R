#' Test the function`find_all_ic`
#'
#' @family tests
#' @export
#' #' @description Using [calc_clust_ic]
#'   for each cluster find the aic and bic and add this to the
#'   cluster results dataframe.
#'
#' @param clust_re cluster results dataframe found using [km_nan]
#' @param num_axis the number of axis
#'
#' @details
#'   The aic and bic are found for each cluster using [calc_clust_ic]
#'   this is added to clust_re\$clusters and clust_re\$clusters is output.
#'   Ncents is also added as an axis to clust_re\$clusters.
#'
#' @return clust_re\$clusters
#'
#' @family ic_functions
#'
#' @export
test_find_aic <- function() {
  print("... Testing aic columns")
  set.seed(240) # setting seed
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
  testit::assert("find_all_aic has found the wrong row names when empty",
    colnames(clust_re) == expect_cols)
  return()
}