#' Test the function for finding k clusters.
#'
#' @details Checks:
#' * there are k clusters found
#' * The number of points with clusters assigned is the number of points input.
#'
#' @export
#' @family tests
test_cluster_kmeans <- function() {
  dummy_traits <- c("T1", "T2", "T3", "T4", "T5", "T6")
  dummy_snps <- c("rs35662", "rs301884", "rs69696",
                  "rs4096", "rs646464", "rs1234")
  num_axis <- length(dummy_traits)
  nsnps <- length(dummy_snps)
  b_mat <- matrix(stats::runif(nsnps * num_axis, 0, 10), nrow = nsnps)
  se_mat <- matrix(stats::runif(nsnps * num_axis, 0, 10), nrow = nsnps)
  p_mat <- matrix(stats::runif(nsnps * num_axis, 0, 1), nrow = nsnps)
  colnames(b_mat) <- dummy_traits
  rownames(b_mat) <- dummy_snps
  colnames(se_mat) <- dummy_traits
  rownames(se_mat) <- dummy_snps
  cat_list <- rep("Exposure", ncol(dummy_beta))
  cat_list[1] <- "Outcome"
  trait_info <- list("pheno_category" = cat_list,
                     "phenotype" = trait_list)
  nang <- num_axis - 1
  mat_list <- list("beta" = b_mat,
                   "se" = se_mat,
                   "pval" = p_mat,
                   "trait_info" = trait_info,
                   "beta_pc" = b_mat[1:nang],
                   "se_pc" = se_mat[1:nang],
                   "pval_pc" = p_mat[1:nang],
                   "tranform" = diag(nang))
  clust_out <- cluster_kmeans(mat_list)
  expec_list <- c("clusters", "clust_dist", "centres")
}