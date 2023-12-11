#' Program for preparing data then clustering it.
#'
#' @param iter_traits dataframe of the indicators for the types of
#'   calculations that will be done. Contains:
#'   * `dname`: data directory name
#'   * `clust_type`: name of the clustering method
#'   * `n_clust`: number of clusters
#'   * `n_pcs`: number of PCs
#'   * `bin_angles`: Y or N use angles
#'   * `bin_p_clust`: Y or N use p-vals in clustering
#'   * `bin_p_score`: Y or N use p-vals in score
#'   * `bin_d_score`: Y or N use cluster distance for cluster scoring
#'   * `nan_rm`: how to handle NaN values. If Y remove and use pairwise.
clustering_program <- function(iter_traits) {
  # Load data

  # Reformulate the data
  # - If angles Y then calculate angles
  # - If PCs Y then calculate PCs

  # Cluster the data
  # - `basic`: kmeans clustering
  # - `min`: minimise aic
  
  # Plot any final outputs
}