#' Test the function for finding the aic of a set of clusters.
#'
#' @details Checks:
#' * there result is a cluster dataframe with an aic column
#' * the aic is numeric
#' * the aic is the same on all columns
#'
#' @export
#' @family tests
test_get_aic <- function() {
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
  colnames(p_mat) <- dummy_traits
  rownames(p_mat) <- dummy_snps
  cat_list <- rep("Exposure", ncol(b_mat))
  cat_list[1] <- "Outcome"
  trait_info <- list("pheno_category" = cat_list,
                     "phenotype" = dummy_traits)
  nang <- num_axis - 1
  mat_list <- list("beta" = b_mat,
                   "se" = se_mat,
                   "pval" = p_mat,
                   "trait_info" = trait_info,
                   "beta_pc" = b_mat[, 1:nang],
                   "se_pc" = se_mat[, 1:nang],
                   "pval_pc" = p_mat[, 1:nang],
                   "tranform" = diag(nang))
  nclust <- 3
  clust_out <- cluster_kmeans(mat_list, nclust = nclust)
  aic_df <- get_aic(clust_out)
  expec_nrows <- nrow(b_mat)
  expec_clusterdf_cols <- c("clust_dist", "clust_num",
                            "clust_prob", "ncents",
                            "num_axis", "aic")
  expec_ncols <- length(expec_clusterdf_cols)
  testit::assert("aic_df has wrong number of rows",
    nrow(aic_df) == expec_nrows
  )
  testit::assert("aic_df has wrong number of columns",
    ncol(aic_df) == expec_ncols
  )
  testit::assert("aic_df has wrong axis",
    all(colnames(aic_df) %in% expec_clusterdf_cols)
  )
  testit::assert("aic_df has wrong rows",
    all(rownames(aic_df) %in% rownames(b_mat))
  )
}