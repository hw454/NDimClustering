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
#'   * `nan_rm`\: how to handle NaN values. If 1 remove and use pairwise.
#' @param test Indicator for whether to run as a test or full program.
#'  default \: 1 (test case)
#'
#' @export
clustering_program <- function(iter_traits, test = 1) {
  # Load data
  data_matrices <- setup_matrices(iter_traits$dname, test = test)

  # Reformulate the data
  # - If angles Y then calculate angles
  # - If PCs Y then calculate PCs
  data_matrices <- reformat_data(data_matrices,
    bin_angles = iter_traits$bin_angles,
    pca_type = iter_traits$pca_type,
    np = iter_traits$n_pcs
  )
  print("after reformating")
  print(data_matrices)
  # Cluster the data
  # - `basic`: kmeans clustering
  # - `min`: minimise aic

  # Plot any final outputs
}