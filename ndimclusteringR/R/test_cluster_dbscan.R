#' Test the function for finding k clusters.
#'
#' @details Checks:
#' * there are k clusters found
#' * The number of points with clusters assigned is the number of points input.
#'
#' @export
#' @family tests
test_cluster_dbscan <- function() {
  dummy_traits <- c("T1", "T2", "T3", "T4", "T5", "T6")
  snp_id <- "rs35662"
  dummy_snps1 <- c(snp_id, "rs301884", "rs69696",
                   "rs4096", "rs646464", "rs1234",
                   "rs203052", "rs646416", "rs98765",
                   "rs124532", "rs161616", "rs001122",
                   "rs111144", "rs222224", "rs01234")
  dummy_snps2 <- c("rs409664", "rs323223", "rs998877",
                   "rs309664", "rs453223", "rs668877",
                   "rs209664", "rs673223", "rs558877",
                   "rs109664", "rs893223", "rs448877")
  num_axis <- length(dummy_traits)
  nsnps1 <- length(dummy_snps1)
  b_mat1 <- matrix(stats::runif(nsnps1 * num_axis, 5, 10), nrow = nsnps1)
  se_mat1 <- matrix(stats::runif(nsnps1 * num_axis, 5, 10), nrow = nsnps1)
  p_mat1 <- matrix(stats::runif(nsnps1 * num_axis, 0, 1), nrow = nsnps1)
  nsnps2 <- length(dummy_snps2)
  b_mat2 <- matrix(stats::runif(nsnps2 * num_axis, 0, 0.5), nrow = nsnps2)
  se_mat2 <- matrix(stats::runif(nsnps2 * num_axis, 0, 1), nrow = nsnps2)
  p_mat2 <- matrix(stats::runif(nsnps2 * num_axis, 0, 1), nrow = nsnps2)
  nsnps <- nsnps1 + nsnps2
  b_mat <- rbind(b_mat1, b_mat2)
  se_mat <- rbind(se_mat1, se_mat2)
  p_mat <- rbind(p_mat1, p_mat2)
  dummy_snps <- c(dummy_snps1, dummy_snps2)
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
  clust_out <- cluster_dbscan(mat_list)
  expec_list <- c("clusters", "clust_dist", "centres")
  testit::assert("cluster_dbscan doesn't contain the right terms",
    all(expec_list %in% names(clust_out))
  )
  expec_ncents <- length(unique(clust_out$clusters$clust_num))
  expec_nrows <- nrow(b_mat)
  extra_label_cols <- 3
  expec_cent_rows <- unique(clust_out$clusters$clust_num)
  clust_cols <- lapply(expec_cent_rows,
    make_clust_col_name
  )
  expec_clust_dist_cols <- c(list("snp_id",
      "num_axis",
      "ncents"
    ),
    clust_cols
  )
  expec_clusterdf_cols <- c("snp_id",
                            "clust_dist",
                            "clust_num",
                            "clust_prob",
                            "num_axis",
                            "ncents")
  expec_ncols <- length(expec_clusterdf_cols)
  expec_clust_dist_ncols <- expec_ncents + extra_label_cols
  expec_cents_ncols <- nang + 2
  print("...Test with complete data")
  testit::assert("Centres are wrong dimension",
    ncol(clust_out$centres) == expec_cents_ncols
  )
  testit::assert("Wrong number of centres",
    nrow(clust_out$centres) <= expec_ncents
  )
  testit::assert("Centres have wrong axis",
    all(colnames(b_mat_pc) %in% colnames(clust_out$centres))
  )
  testit::assert("Centres has wrong rows",
    all(rownames(clust_out$centres) %in% expec_cent_rows)
  )
  testit::assert("Clust_dist are wrong dimension",
    ncol(clust_out$clust_dist) <= expec_clust_dist_ncols
  )
  testit::assert("Clust_dist has wrong number of rows",
    nrow(clust_out$clust_dist) == expec_nrows
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
  testit::assert("Clust_dist has wrong rows",
    all(clust_out$clusters$snp_id %in% rownames(b_mat))
  )
  print("...Test with incomplete data")
  snp_id <- "rs35662"
  snp_id2 <- "rs646416"
  dummy_snps <- c(snp_id, "rs301884", "rs69696",
                  "rs4096", "rs646464", "rs1234",
                  "rs203052", snp_id2, "rs98765",
                  "rs124532", "rs161616", "rs001122",
                  "rs111144", "rs222224", "rs01234",
                  "rs409664", "rs323223", "rs998877",
                  "rs309664", "rs453223", "rs668877",
                  "rs209664", "rs673223", "rs558877",
                  "rs109664", "rs893223", "rs448877")
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
  # Create missing data
  b_mat[snp_id, "T2"] <- NA
  se_mat[snp_id, "T2"] <- NA
  p_mat[snp_id, "T2"] <- NA
  b_mat[snp_id, "T2"] <- NA
  se_mat[snp_id, "T2"] <- NA
  p_mat[snp_id, "T2"] <- NA
  b_mat[snp_id2, "T3"] <- NA
  se_mat[snp_id2, "T3"] <- NA
  p_mat[snp_id2, "T3"] <- NA
  b_mat[snp_id2, "T5"] <- NA
  se_mat[snp_id2, "T5"] <- NA
  p_mat[snp_id2, "T5"] <- NA
  b_mat_pc <- b_mat[, 1:nang]
  mat_list <- list("beta" = b_mat,
                   "se" = se_mat,
                   "pval" = p_mat,
                   "trait_info" = trait_info,
                   "beta_pc" = b_mat_pc,
                   "se_pc" = se_mat[, 1:nang],
                   "pval_pc" = p_mat[, 1:nang],
                   "tranform" = diag(nang))
  clust_out <- cluster_dbscan(mat_list)
  testit::assert("cluster_dbscan doesn't contain the right terms",
    all(expec_list %in% names(clust_out))
  )
  expec_ncents <- length(unique(clust_out$clusters$clust_num))
  expec_cent_rows <- unique(clust_out$clusters$clust_num)
  expec_clust_dist_ncols <- expec_ncents + extra_label_cols
  testit::assert("Centres are wrong dimension",
    ncol(clust_out$centres) == expec_cents_ncols
  )
  testit::assert("Wrong number of centres",
    nrow(clust_out$centres) <= expec_ncents
  )
  testit::assert("Centres have wrong axis",
    all(colnames(b_mat_pc) %in% colnames(clust_out$centres))
  )
  testit::assert("Centres has wrong rows",
    all(rownames(clust_out$centres) %in% expec_cent_rows)
  )
  testit::assert("Clust_dist are wrong dimension",
    ncol(clust_out$clust_dist) <= expec_clust_dist_ncols
  )
  testit::assert("Clust_dist has wrong number of rows",
    nrow(clust_out$clust_dist) == expec_nrows
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
  testit::assert("Clust_dist has wrong rows",
    all(clust_out$clusters$snp_id %in% rownames(b_mat))
  )
}