#' Function to run clustering, scoring and plot the results.
#'
#' @param data_matrices - A list containing
#'   * "beta" matrix corresponding to the data
#'   * "pval" probabilities associated with "beta" values.
#'   * "se" standard errors associated with "beta" values.
#'   * "trait_info" dataframe of all the traits.
#' @param exp_pheno - Label for the exposure phenotype
#' @param out_pheno - Label for the outcome phenotype
#' @param res_dir0 - The location of the base results directory.
#'  default \:"", the current directory.
#' @param thresholds - List for the  threshold related variables.
#'   * "threshmul" is multiplied by the data variance for cluster
#'   difference threshold.
#'   * "diff" is the the cluster difference threshold once found.
#'   * "clust" is the required distance between cluster centres for a
#'   cluster to be considered converged.
#' @param na_handling - List of the terms related to handling NaNs
#'   * "narm" - TRUE (default) if NaNs are to be removed in calculations.
#'   * "percent" - percentage of a column which must be not NaN. default\:0.95
#' @param iter The iteratio number
#' @param iter_df - list of terms indicating the type of program.
#'   * "iter" - default\:0
#'   * "bp_on" - default\:TRUE switch if pval scores is to be used.
#'   * "clust_prob_on" - default\:TRUE switch to use prob of being in cluster
#'   * "clust_typ" - default="basic", the clustering method to use.
#'   * "ndim_typ" - Whether to use all the traits or iterate through.
#' @param norm_typs - List of norm types
#'   * "clust" - default="F"
#'   * "thresh" - default="F"
#' @param nums - List of important numbers.
#'   * "max_dist" - Maximum distance between points, default=1
#'   * "np" - Number of principal components
#'   * "nc" - Number of clusters.
#'
#' @details
#'   based on iter_traits\:ndim_typ call the comparison and cluster function.
#'   If iter_trais\:ndim_typ=="all"{
#'     * run [clust_pca_all]
#'     * Plot the trait heatmap using [plot_trait_heatmap]
#'     * Plot the scatter plot using [plot_clust_scatter]
#'     * Plot the transform heatmap using [plot_transform_heatmap]
#'     } else{
#'     * run [clust_pca_compare_iterative]
#'     * Plot the heatmap for each iteration using [plot_trait_heatmap]
#'     * Plot the max difference through the iterations using
#'       [plot_max_diff]
#'     }
#'   df_list contains the results from all the iterations.
#'
#' @return df_list
#'
#' @family cluster_wrappers
#'
#' @export
full_cluster_and_plot <- function(data_matrices, exp_pheno, out_pheno,
  res_dir0 = "",
  iter = 1,
  iter_df = data.frame("iter" = 1,
                       "bp_on" = FALSE,
                       "clust_prob_on" = FALSE,
                       "clust_typ" = "basic",
                       "ndim_typ" = "all"),
  thresholds = list("threshmul" = 5,
                    "diff" = 1e-5,
                    "clust" = 1e-5),
  na_handling = list("narm" = TRUE, "percent" = 0.95),
  norm_typs = list("clust" = "F", "thresh" = "F"),
  nums = list("max_dist" = 1, "np" = 3, "nr" = 5)
) {
  iter_traits <- iter_df[iter, ]
  res_dir <- set_dir(res_dir0, iter_traits)
  iter_traits["res_dir"] <- res_dir
  # Find the distances between all points to initialise the threshold
  # for cluster difference.
  # FIXME
  # This needs to be set after PCA since axis change
  print("Begining algorithm for inputs")
  print(iter_traits)
  if (iter_traits$ndim_typ == "all") {
    out <- clust_pca_all(data_matrices = data_matrices,
                         out_pheno = out_pheno,
                         na_handling = na_handling,
                         iter_traits = iter_traits,
                         norm_typs = norm_typs,
                         nums = nums)
    print("Clust done")
    print("Plot pc scores")
    c_scores_pc <- out$clust_pc_scores
    plot_trait_heatmap(c_scores_pc, iter_traits, prefix = "pca")
    print("Plot trait scores")
    c_scores_tr <- out$clust_trait_scores
    plot_trait_heatmap(c_scores_tr, iter_traits, prefix = "pheno")
    print("Heatmap plot done")
    num_axis <- length(data_matrices$trait_info$phenotype)
    col1 <- colnames(out$b_pc)[1]
    col2 <- colnames(out$b_pc)[2]
    plot_clust_scatter(out$clust_items, out$b_pc, out$se_pc, iter_traits,
                       c1 = col1, c2 = col2, num_axis = num_axis)
    print("scatter plot done")
    plot_clust_exp_out_scatter(out$clust_items,
                               data_matrices$beta,
                               data_matrices$se,
                               iter_traits,
                               exp_pheno,
                               out_pheno,
                               num_axis = num_axis)
    print("exposure outcome scatter plot done")
    # Only plot max diff when iterating through the axis
    plot_clust_scatter_rgb(out$clust_membership,
                           out$b_pc,
                           out$se_pc,
                           iter_traits,
                           c1 = col1,
                           c2 = col2,
                           num_axis = num_axis)
    print("rgb pca scatter plot done")
    plot_clust_scatter_rgb(out$clust_membership,
                           data_matrices$beta,
                           data_matrices$se,
                           iter_traits,
                           c1 = out_pheno,
                           c2 = exp_pheno,
                           num_axis = num_axis)
    print("rgb exp out scatter plot done")
    # Plot the transform heatmap.
    plot_transform_heatmap(out$e_mat, iter_traits)
    out <- list("iter_df" = iter_df,
                "c_scores_pc" = c_scores_pc,
                "c_scores_tr" = c_scores_tr,
                "max_df" = out$max_diff)
  } else if (iter_traits$ndim_typ == "iterative") {
    out <- clust_pca_compare_iterative(data_matrices = data_matrices,
      out_pheno = out_pheno,
      na_handling = na_handling,
      iter_traits = iter_traits,
      norm_typs = norm_typs,
      nums = nums
    )
    print("Clust done")
    max_diff_df <- out$max_diff
    c_scores_tr <- out$clust_trait_scores
    plot_trait_heatmap(c_scores_tr, iter_traits)
    print("Heatmap plot done")
    plot_mxdf(max_diff_df, iter_traits)
    print("Diff plot done")
  }
  return(out)
}