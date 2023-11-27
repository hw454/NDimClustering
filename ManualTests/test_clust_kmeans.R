number_row_col_names <- function(mat) {
  colnames(mat) <- seq_len(ncol(mat))
  rownames(mat) <- seq_len(nrow(mat))
  return(mat)
}
form_beta_corr <- function(n_path, d, rand_shift, a_list, b_list) {
  a <- a_list[n_path]
  b <- b_list[n_path]
  x_mat <- matrix(rep(1:d, times = d), nrow = d)
  beta <- (a * x_mat + b) + rand_shift
  beta[, 1] <- 1:d
  return(beta)
}
create_pc <- function(i) {
  return(paste0("PC", i))
}
test_clust_kmeans_function <- function(d = 10,
                                       num_path = 0,
                                       pc_type = "prcomp",
                                       n_pc = 2,
                                       space_typ = "angle",
                                       clust_typ = "min",
                                       how_cents = "rand",
                                       ndim_typ = "all") {
  np <- num_path + 1
  print(paste("Test", pc_type, "with", np, "pathways"))
  iter_traits <- data.frame(
    "bp_on" = TRUE,
    "clust_prob_on" = TRUE,
    "clust_typ" = paste0("test_", clust_typ, space_typ),
    "ndim_typ" = paste0("test_", ndim_typ),
    "how_cents" = how_cents,
    "pc_type" = pc_type,
    "num_paths" = num_path,
    "res_dir" = paste0("PC_TestResults/")
  )
  print(iter_traits)
  rand_mat <- matrix(runif(d * d, -1, 1), nrow = d)
  if (num_path > 0) {
    a_list <- runif(num_path + 1, 0, 3)
    b_list <- runif(num_path + 1, 0, 0)
    mat_list <- lapply(1:(num_path + 1),
      form_beta_corr,
      d = d,
      rand_shift = rand_mat,
      a_list = a_list,
      b_list = b_list
    )
    dummy_beta <- Reduce(rbind, mat_list)
  } else {
    dummy_beta <- rand_mat
  }
  dummy_se <- matrix(runif((num_path + 1) * d * d, 0, 1),
    nrow = (num_path + 1) * d
  )
  dummy_p <- matrix(runif((num_path + 1) * d * d),
    nrow = (num_path + 1) * d
  )
  dummy_beta <- number_row_col_names(dummy_beta)
  dummy_se <- number_row_col_names(dummy_se)
  dummy_p <- number_row_col_names(dummy_p)
  # Compute the Sample SE with na_rm
  num_axis <- ncol(dummy_beta)
  if (pc_type == "prcomp") {
    pca_beta <- stats::prcomp(dummy_beta,
      center = TRUE,
      scale = TRUE,
      retx = TRUE,
      rank = n_pc
    )
    t_mat <- pca_beta$rotation
    b_pc_mat <- pca_beta$x
    p_pc_mat <- ndimclusteringR::transform_coords_in_mat(dummy_p, t_mat)
    se_pc_mat <- ndimclusteringR::transform_coords_in_mat(dummy_se, t_mat)
    # Rescale se so it is within (0,1)
    se_pc_mat <- ndimclusteringR::calc_scale_mat(se_pc_mat)
    pc_cols <- lapply(1:n_pc, create_pc)
    colnames(se_pc_mat) <- pc_cols
    out_list <- list(
      "beta" = b_pc_mat,
      "pval" = p_pc_mat,
      "se" = se_pc_mat,
      "transform" = t_mat
    )
  } else {
    n <- ncol(dummy_beta)
    t_mat <- diag(n)
    out_list <- list(
      "beta" = dummy_beta,
      "pval" = dummy_p,
      "se" = dummy_se,
      "transform" = t_mat
    )
  }
  cluster_out <- ndimclusteringR::cluster_kmeans(out_list,
    iter_traits = iter_traits,
    nclust = num_path + 1,
    max_dist = 10.0,
    space_typ = space_typ,
    clust_typ = clust_typ,
    how_cents = how_cents
  )
  cluster_df <- cluster_out$clusters
  cluster_df["num_axis"] <- num_axis
  cluster_df <- tibble::rownames_to_column(cluster_df, "snp_id")
  c1 <- colnames(dummy_beta)[1]
  c2 <- colnames(dummy_beta)[2]
  nclust <- length(unique(cluster_df$clust_num))
  if (nclust >= 3) {
    ndimclusteringR::plot_clust_scatter_rgb_test(cluster_out$clust_dist,
      out_list$beta,
      out_list$se,
      iter_traits,
      num_axis = num_axis
    )
  }
  ndimclusteringR::plot_clust_scatter_test(cluster_df,
    out_list$beta,
    out_list$se,
    iter_traits,
    num_axis = num_axis
  )
  ndimclusteringR::plot_clust_expout_scatter_test(cluster_df,
    dummy_beta,
    dummy_se,
    iter_traits,
    exp_pheno = c1,
    out_pheno = c2,
    num_axis = num_axis
  )
}
