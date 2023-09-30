#' Cluster data on the principal components and score the clusters on
#' association. Compare the cluster scores and iteratively add a
#' trait each time.
#'
#' @param data_matrices - A list containing the matrices
#'   * "beta" matrix corresponding to the data.
#'   * "pval" the matrix of probabilities associated with each score.
#'   * "se" matrix of standard errors associated with each score.
#'   * "trait_info" which is the dataframe of all the traits.
#' @param out_pheno - Label for the outcome phenotype
#' @param thresholds - List for the three threshold related variables.
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
#'   * "clust_prob_on" - TRUE (default) if prob of being in cluster to be used.
#'   * "clust_typ" - default="basic", the clustering method to use.
#' @param norm_typs - List of norm types
#'   * "clust" - default="F"
#'   * "thresh" - default="F"
#' @param nums - List of important numbers.
#'   * "max_dist" - Maximum distance between points, default=1
#'   * "np" - Number of principal components
#'   * "nc" - Number of clusters.
#'
#' @family cluster_wrappers
#'
#' @description
#'   df_list is the list of of the clusters, the principal component matrices,
#'   cluster scores, and the maximum difference between clusters.
#'   Each iteration generated using [clust_pca_compare_single].
#'
#' @return df_list
#'
#' @export
clust_pca_compare_iterative <- function(data_matrices,
                         out_pheno,
                         thresholds = list("threshmul" = 5,
                                            "diff" = 1e-5,
                                            "clust" = 1e-5),
                         na_handling = list("narm" = TRUE, "percent" = 0.95),
                         iter_traits = list(iter = 0,
                                            "bp_on" = TRUE,
                                            "clust_prob_on" = TRUE,
                                            "clust_typ" = "basic"),
                         norm_typs = list("clust" = "F", "thresh" = "F"),
                         nums = list("max_dist" = 1, "np" = 3, "nr" = 5)
                        ) {
  # Iterate through the traits in trait_info contained in the
  # data_matrices list. Find the principal components.
  # Transform the data onto these.
  # Cluster the transformed data at each iteration.
  # If there is a distinct difference between two clusters exit.
  trait_df <- data.frame(
    label = character(),
    axes_ind = integer()
  )
  # Data frame for recording the cluster scores.
  c_scores <- data.frame(
    num_axis = integer(),
    clust_num = integer()
  )
  # Add np columns for each PC
  c_scores <- add_np_cols(c_scores, nums$np)
  # Initialise with outcome
  trait_df <- data.frame(label = out_pheno,
                       axes_ind = which(data_matrices$trait_info$phenotype == out_pheno)[1] #nolint
                       )
  max_df <- data.frame(num_axis = integer(),
                       cn1 = integer(),
                       cn2 = integer(),
                       max_diff = integer())
  cluster_df <- data.frame("snp_id" = character(),
                          "clust_num" = integer(),
                          "clust_prob" = numeric(),
                          "clust_dist" = numeric(),
                          "num_axis" = integer(),
                          "ncents" = integer())
  cluster_df <- add_np_cols(cluster_df, nums$np)
  df_list <- list("clust_scores" = c_scores,
                  "max_diff" = max_df,
                  "e_list" = list(),
                  "trait" = list(),
                  "cluster_items" = cluster_df)
  pheno_list <- data_matrices$trait_info$phenotype
  na <- length(pheno_list)
  for (ai in 1:na){
    # Update the traits forming the axis
    # Add the trait to the trait dataframe
    a <- pheno_list[ai]
    print(paste("New trait on axis", a))
    covered <- (a %in% trait_df$label)
    allna <- check_col_na(data_matrices$beta[, a], na_handling$percent)
    if (!covered && !allna) {
      # Add trait to trait dataframe
      a_ind <- which(pheno_list == a)[1]
      trait_row <- data.frame(label = a,
                              axes_ind = a_ind)
      trait_df <- rbind(trait_df, trait_row)
      df_list$trait <- rbind(df_list$trait, trait_df)
      df_list <- clust_pca_compare_single(df_list,
                                    iter_traits = iter_traits,
                                    num_axis = ai,
                                    data_matrices = data_matrices,
                                    na_handling = na_handling,
                                    thresholds = thresholds,
                                    norm_typs = norm_typs,
                                    nums = nums)
    }
  }
  return(df_list)
}