#' Single iteration for iterating through trait list.
#'
#' @param df_list the results from previous iterations to append to.
#' @param data_matrices - A list containing
#'   * "beta" main data matrix
#'   * "pval" probabilities associated with "beta" values
#'   * "se"  standard errors associated with "beta" values.
#'   * "trait_info" data frame of all the traits.
#' @param num_axis Number of trait axes in use.
#' @param thresholds - List of threshold related variables.
#'   * "threshmul" is multiplied by the data variance for cluster
#'   difference threshold.
#'   * "diff" is the the cluster difference threshold once found.
#'   * "clust" is the required distance between cluster centres for a
#'   cluster to be considered converged.
#' @param na_handling - List of the terms related to handling NaNs
#'   * "narm" - TRUE (default) if NaNs are to be removed in calculations.
#'   * "percent" - percentage of a column which must be not NaN. default\:0.95
#' @param iter_traits - list of terms indicating the type of program.
#'   * "iter" - integer, (default 0)
#'   * "bp_on" - TRUE (default) if probability of scores is to be used.
#'   * "clust_prob_on" - TRUE\:default switch for using prob of being in cluster
#'   * "clust_typ" - default\:"basic", the clustering method to use.
#' @param norm_typs - List of norm types
#'   * "clust" - default\:"F"
#'   * "thresh" - default\:"F"
#' @param nums - List of important numbers.
#'   * "max_dist" - Maximum distance between points, default=1
#'   * "np" - Number of principal components
#'   * "nc" - Number of clusters.
#'
#' @details
#'   df_list is the list of of the clusters, the principal component matrices,
#'   cluster scores, and the maximum difference between clusters.
#'   Each iteration generated using [clust_pca_compare_single].
#'
#' @family cluster_wrappers
#'
#' @return df_list
#'
#' @export
clust_pca_compare_single <- function(df_list, iter_traits,
                                            num_axis,
                                            data_matrices,
                                            na_handling,
                                            thresholds,
                                            norm_typs,
                                            nums) {
  # Extract the trait_df dataframe from df_list
  trait_df <- df_list$trait
  # If the trait is not all NaN then run clustering.
  print(paste("PCA on", num_axis, "axes"))
  # Get the data upto this axis
  b_iter_mat <- data_matrices$beta[, trait_df$label]
  se_iter_mat <- data_matrices$se[, trait_df$label]
  pval_iter_mat <- data_matrices$pval[, trait_df$label]
  # Cluster the data on these axes
  pca_list <- find_principal_components(b_iter_mat, pval_iter_mat, se_iter_mat,
                          nums$np, na_handling$narm)
  b_pc_mat    <- pca_list$beta
  p_pc_mat <- pca_list$pval
  t_mat       <- pca_list$transform
  # Store the matrices for result output
  df_list$e_mat <- t_mat
  df_list$b_pc <- b_pc_mat
  df_list$se_pc <- pca_list$se
  # Get column names for PCs
  pc_cols <- colnames(b_pc_mat)
  # Cluster the data on these axes
  if (grepl("angle", iter_traits$clust_typ, fixed = TRUE)) {
    st <- "angle"
  } else {
    st <- "regular"
  }
  if (grepl("min", iter_traits$clust_typ)) {
    cluster_out <- cluster_kmeans_min(pca_list,
                                      nums$nr,
                                      space_typ = st,
                                      clust_prob_on = iter_traits$clust_prob_on, # nolint
                                      norm_typ = norm_typs$clust,
                                      threshold = thresholds$clust,
                                      narm = na_handling$narm)
  } else if (grepl("basic", iter_traits$clust_typ)) {
    cluster_out <- cluster_kmeans_basic(pca_list,
                                        nums$nr,
                                        space_typ = st,
                                        clust_prob_on = iter_traits$clust_prob_on, # nolint
                                        threshold = thresholds$clust,
                                        norm_typ = norm_typs$clust,
                                        narm = na_handling$narm)
  }
  cluster_df <- cluster_out$clusters
  centroids_df <- cluster_out$centres
  # Calculate the distance to all the cluster centres
  clust_dist_df <- calc_clust_dist(df_list$b_pc, centroids_df)
  # Find the set of cluster numbers
  c_nums <- unique(cluster_df$clust_num)
  # Score the clustered data based on affiliations with axes.
  # Find the score for each PC
  c_score0 <- score_all_clusters(cluster_df,
                                beta_mat = b_pc_mat,
                                pval_mat = p_pc_mat,
                                bp_on = iter_traits$bp_on,
                                clust_prob_on = iter_traits$clust_prob_on,
                                num_axis = num_axis)
  # Iterate through each cluster and compare across the others to find if
  # any pair have a distinct difference.
  diff_score_list <- lapply(c_nums, compare_oneclust_tolist,
          c_nums = c_nums,
          c_score0 = c_score0,
          axis = pc_cols,
          clust_norm = norm_typs$clust
        )
  diff_score_list <- diff_score_list[!sapply(diff_score_list, is.null)]
  diff_scores <- Reduce(rbind, diff_score_list)
  # Find the pair with maximum cluster score diff and store in max_df0
  row <- which.max(diff_scores$diff)
  max_df0 <- diff_scores[row, ]
  max_df0["num_axis"] <- num_axis
  c_score0["num_axis"] <- num_axis
  clust_dist_df["num_axis"] <- num_axis
  cluster_df["num_axis"] <- num_axis
  clust_dist_df <- tibble::rownames_to_column(clust_dist_df, "snp_id")
  cluster_df <- tibble::rownames_to_column(cluster_df, "snp_id")
  df_list$clust_items <- rbind(df_list$clust_items, cluster_df)
  df_list$max_diff <- rbind(df_list$max_diff, max_df0)
  df_list$clust_scores <- rbind(df_list$clust_scores, c_score0)
  df_list$clust_membership <- dplyr::bind_rows(df_list$clust_membership,
                                              clust_dist_df)
  return(df_list)
}