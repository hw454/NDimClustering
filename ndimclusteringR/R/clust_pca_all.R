#' Find the principal components (PCs) and cluster the data.
#' Scores the clusters on the PCs
#'
#' @param data_matrices - A list containing
#'   * "beta" matrix corresponding to the data
#'   * "pval" probabilities associated with "beta" values.
#'   * "se" standard errors associated with "beta" values.
#'   * "trait_info" dataframe of all the traits.
#' @param out_pheno - Label for the outcome phenotype
#' @param thresholds - List for the  threshold related variables.
#'   * "threshmul" is multiplied by the data variance for cluster
#'   difference threshold.
#'   * "diff" is the the cluster difference threshold once found.
#'   * "clust" is the required distance between cluster centres for a
#'   cluster to be considered converged.
#' @param na_handling - List of the terms related to handling NaNs
#'   * "narm" - TRUE (default) if NaNs are to be removed in calculations.
#'   * "percent" - percentage of a column which must be not NaN. default\:0.95
#' @param iter_traits - list of terms indicating the type of program.
#'   * "iter" - default\:0
#'   * "bp_on" - default\:TRUE switch if pval scores is to be used.
#'   * "clust_prob_on" - default\:TRUE switch to use prob of being in cluster
#'   * "clust_typ" - default="basic", the clustering method to use.
#' @param norm_typs - List of norm types
#'   * "clust" - default="F"
#'   * "thresh" - default="F"
#' @param nums - List of important numbers.
#'   * "max_dist" - Maximum distance between points, default=1
#'   * "np" - Number of principal components
#'   * "nc" - Number of clusters.
#'
#' @details
#'   df_list is the list of of the clusters, the principal component matrices,
#'   cluster scores, and the maximum difference between clusters.
#'   Generated using [clust_pca_compare_single].
#'
#' @family cluster_wrappers
#'
#' @return df_list
#'
#' @export
clust_pca_all <- function(data_matrices,
                         out_pheno,
                         thresholds = list("threshmul" = 5,
                                            "diff" = 1e-5,
                                            "clust" = 1e-5),
                         iter_traits = list(iter = 0,
                                            "bp_on" = TRUE,
                                            "clust_prob_on" = TRUE,
                                            "clust_typ" = "basic"),
                         na_handling = list("narm" = TRUE, "percent" = 0.95),
                         norm_typs = list("clust" = "F", "thresh" = "F"),
                         nums = list("max_dist" = 1, "np" = 3, "nr" = 5)
                        ) {
  # Data frame for recording the cluster scores on the pcs
  c_scores_pc <- data.frame(
    num_axis = integer(),
    clust_num = integer()
  )
  # Add np columns for each PC
  c_scores_pc <- add_np_cols(c_scores_pc, nums$np)
  # Data frame for recording the cluster scores on the traits
  c_scores_tr <- data.frame(
    num_axis = integer(),
    clust_num = integer()
  )
  # Add np columns for each PC
  for (tr in data_matrices$trait_info$phenotype){
    c_scores_tr[tr] <- numeric()
  }
  # Initialise with outcome
  trait_df <- data.frame(
                       num_axis = 1,
                       label = out_pheno,
                       axes_ind = which(data_matrices$trait_info$phenotype == out_pheno)[1] #nolint
                       )
  max_df <- data.frame(num_axis = integer(),
                       cn1 = integer(),
                       cn2 = integer(),
                       diff = integer())
  cluster_df <- data.frame("snp_id" = character(),
                          "clust_num" = integer(),
                          "clust_prob" = numeric(),
                          "clust_dist" = numeric(),
                          "num_axis" = integer())
  cluster_df <- add_np_cols(cluster_df, nums$np)
  clust_dist_df <- data.frame("snp_id" = character())
  clust_dist_df <- add_nclust_cols(clust_dist_df, nums$nr)
  df_list <- list("clust_pc_scores" = c_scores_pc,
                  "clust_trait_scores" = c_scores_tr,
                  "max_diff" = max_df,
                  "e_list" = list(),
                  "b_pc_list" = list(),
                  "se_pc_list" = list(),
                  "trait" = trait_df,
                  "clust_items" = cluster_df,
                  "clust_membership" = clust_dist_df)
  # Fill trait_df with valid traits before running main program.
  trait_df <- make_trait_df(pheno_list = data_matrices$trait_info$phenotype,
                            data_mat = data_matrices$b_df,
                            na_percent = na_handling$percent
                            )
  num_axis <- nrow(trait_df)
  trait_df["num_axis"] <- num_axis
  df_list$trait <- trait_df
  out_list <- clust_pca_compare_single(df_list = df_list,
                                      iter_traits = iter_traits,
                                      num_axis = num_axis,
                                      data_matrices = data_matrices,
                                      na_handling = na_handling,
                                      thresholds = thresholds,
                                      norm_typs = norm_typs,
                                      nums = nums)
  return(out_list)
}