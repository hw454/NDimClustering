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
  colnames(p_mat) <- dummy_traits
  rownames(p_mat) <- dummy_snps
  cat_list <- rep("Exposure", ncol(b_mat))
  cat_list[1] <- "Outcome"
  trait_info <- list("pheno_category" = cat_list,
                     "phenotype" = dummy_traits)
  nang <- num_axis - 1
  b_mat_pc <- b_mat[, 1:nang]
  mat_list <- list("beta" = b_mat,
                   "se" = se_mat,
                   "pval" = p_mat,
                   "trait_info" = trait_info,
                   "beta_pc" = b_mat_pc,
                   "se_pc" = se_mat[, 1:nang],
                   "pval_pc" = p_mat[, 1:nang],
                   "tranform" = diag(nang))
  nclust <- 3
  clust_out <- cluster_kmeans(mat_list, nclust = nclust)
  expec_list <- c("clusters", "clust_dist", "centres")
  testit::assert("km_nan doesn't contain the right terms",
    all(expec_list %in% names(clust_out))
  )
  expec_ncents <- nclust
  expec_nrows <- nrow(b_mat)
  expec_clust_dist_cols <- c("snp_id",
    make_clust_col_name(1),
    make_clust_col_name(2),
    make_clust_col_name(3),
    "num_axis",
    "ncents"
  )
  expec_cent_rows <- seq_len(nclust)
  expec_clusterdf_cols <- c("snp_id",
                            "clust_dist",
                            "clust_num",
                            "clust_prob",
                            "ncents",
                            "num_axis")
  expec_ncols <- length(expec_clusterdf_cols)
  expec_cent_ncols <- nang + 2
  extra_label_cols <- 3
  expec_clust_dist_ncols <- nclust + extra_label_cols
  testit::assert("Centres are wrong dimension",
    ncol(clust_out$centres) == expec_cent_ncols
  )
  testit::assert("Wrong number of centres",
    nrow(clust_out$centres) == expec_ncents
  )
  testit::assert("Centres have wrong axis",
    all(colnames(b_mat_pc) %in% colnames(clust_out$centres))
  )
  testit::assert("Centres has wrong rows",
    all(rownames(clust_out$centres) %in% expec_cent_rows)
  )
  testit::assert("Clust_dist are wrong dimension",
    ncol(clust_out$clust_dist) == expec_clust_dist_ncols
  )
  testit::assert("Clust_dist has wrong number of rows",
    nrow(clust_out$clust_dist) == expec_nrows
  )
  testit::assert("Clust_dist and centres have different number of clusters",
    {ncol(clust_out$clust_dist) == nrow(clust_out$centres) + extra_label_cols
    }
  )
  testit::assert("Clust_dist has wrong axis",
    all(colnames(clust_out$clust_dist) %in% expec_clust_dist_cols)
  )
  testit::assert("Clust_dist has wrong rows",
    all(clust_out$clust_dist$snp_id %in% rownames(b_mat))
  )
  testit::assert("Cluster_df has wrong number of rows",
    nrow(clust_out$clusters) == expec_nrows
  )
  testit::assert("Cluster_df has wrong number of columns",
    ncol(clust_out$clusters) == expec_ncols
  )
  testit::assert("Cluster_df has wrong axis",
    all(colnames(clust_out$clusters) %in% expec_clusterdf_cols)
  )
  testit::assert("Cluster_df has wrong rows",
    all(clust_out$clusters$snp_id %in% rownames(b_mat))
  )
}