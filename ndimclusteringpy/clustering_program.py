#' Program for preparing data then clustering it.
#'
#' @param iter_traits dataframe of the indicators for the types of
#'   calculations that will be done.
#'   Contains\:
#'   * `dname`\: data directory name
#'   * `clust_type`\: name of the clustering method
#'   * `pc_type`\: method for finding principal components,
#'                options are "manual", "prcomp", "none"
#'   * `n_clust`\: number of clusters
#'   * `n_pcs`\: number of PCs
#'   * `bin_angles`\: 1 or 0 use angles
#'   * `bin_p_clust`\: 1 or 0 use p-vals in clustering
#'   * `bin_p_score`\: 1 or 0 use p-vals in score
#'   * `bin_d_score`\: 1 or 0 use cluster distance for cluster scoring
#'   * `point_eps` \: the density of points for density based methods. default \:0.4
#'   * `nan_rm`\: how to handle NaN values. If 1 remove and use pairwise.
#' @param test Indicator for whether to run as a test or full program.
#'  default \: 1 (test case)
#'
#' @export
def clustering_program(iter_traits, test = 0):
  # Load data
  data_matrices <- setup_matrices(iter_traits.dname, test = test)

  # Reformulate the data
  # - If angles Y then calculate angles
  # - If PCs Y then calculate PCs
  # data_matrices = reformat_data(data_matrices,
  #   bin_angles = iter_traits.bin_angles,
  #   pca_type = iter_traits.pca_type,
  #   np = iter_traits.n_pcs
  # )
  # Cluster the data
  # - `basic`: kmeans clustering
  # - `min`: minimise aic
  # if (grepl("basic", iter_traits$clust_type,
  #           ignore.case = TRUE)) {
  #   clust_out <- cluster_kmeans(data_matrices,
  #     nclust = iter_traits$nclust,
  #     how_cents = iter_traits$how_cents,
  #     bin_p_clust = iter_traits$bin_p_clust
  #   )
  # } else if (grepl("min", iter_traits$clust_type,
  #                  ignore.case = TRUE)) {
  #   clust_out <- cluster_kmeans_min(data_matrices,
  #     nclust = iter_traits$nclust,
  #     how_cents = iter_traits$how_cents,
  #     bin_p_clust = iter_traits$bin_p_clust
  #   )
  # } else if (grepl("dbscan", iter_traits$clust_type,
  #                  ignore.case = TRUE)) {
  #   clust_out <- cluster_dbscan(data_matrices,
  #     eps = iter_traits$point_eps,
  #     bin_p_clust = FALSE
  #   )
  # }
  # num_axis <- ncol(data_matrices$beta_pc)
  # # Plot any final outputs
  # # -Plot the original data
  # c1 <- colnames(data_matrices$beta)[1]
  # c2 <- colnames(data_matrices$beta)[2]
  # plot_clust_scatter_test(clust_out$clusters,
  #   data_matrices$beta,
  #   iter_traits,
  #   c1 = c1,
  #   c2 = c2,
  #   num_axis = num_axis
  # )
  # # -Plot the reformatted data
  # c1 <- colnames(data_matrices$beta_pc)[1]
  # c2 <- colnames(data_matrices$beta_pc)[2]
  # plot_clust_scatter_test(clust_out$clusters,
  #   data_matrices$beta_pc,
  #   iter_traits,
  #   c1 = c1,
  #   c2 = c2,
  #   num_axis = num_axis,
  #   save_suffix = "_pc"
  # )
  # # If angles then plot the angles
  # if (iter_traits$bin_angles) {
  #   c1 <- colnames(data_matrices$beta_ang)[1]
  #   c2 <- colnames(data_matrices$beta_ang)[2]
  #   plot_clust_scatter_test(clust_out$clusters,
  #     data_matrices$beta_ang,
  #     iter_traits,
  #     c1 = c1,
  #     c2 = c2,
  #     num_axis = num_axis,
  #     save_suffix = "_ang"
  #   )