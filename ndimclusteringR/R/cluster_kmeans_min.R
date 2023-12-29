#' Min-aic clustering
#'
#' @description Run iterations of k-means clustering using [cluster_kmeans]
#' then take the clusters which minimise the aic.
#' AIC calculated using [get_aic].
#'
#' @param data_matrices List of the matrices of data.
#'  Conatins the\:
#'    * "beta" initial matrix of assocaition data
#'    * "se" initial matrix of standard error data
#'    * "pval" initial matrix of p-value data
#'    * "trait_info" list of informations about the traits
#'    * "beta_pc" reformatted matrix of assocaition data
#'    * "se" reformatted matrix of standard error data
#'    * "pval" reformatted matrix of p-value data
#'    * "transform" matrix used for the pca transformation
#'    * Rows correspond to snps, columns to traits
#' @param nclust The maximum number of clusters to consider. default\:10
#' @param max_dist The maximum distance between any two points.
#' @param how_cents How the centroids will be initialised. If "rand (default)"
#'   then the coordinates are created using a uniform distribution on the range
#'   of each axis. If "points" then the centroid coordinates are assign to
#'   nclust random points from the dataspace.
#' @param bin_p_clust Bool switch. If TRUE (default) then cluster
#'   probability if used to weight the cluster scores.
#' @param threshold The threshold for distance between cluster centres
#'   for clusters to be considered converged.
#'
#' @details
#'   clust_out = [clusters, clust_dist, centres].
#'   clust_out$clusters is cluster_df.
#'   The clust_dist dataframe has columns corresponding to the clusters,
#'     the rows are the datapoints and each value is the distance from that
#'     point to the cluster centre.
#'   The centres dataframe has columns corresponding to the axis and the rows
#'     corresponding to the clusters. Each row is the co-ordinate of the cluster
#'     centroid.
#'
#' @export
#' @family cluster_functions
#' @family k_means
cluster_kmeans_min <- function(data_matrices,
  nclust = 10,
  max_dist = 10.0,
  how_cents = "rand",
  bin_p_clust = TRUE,
  threshold = 1e-5
) {
  # Run kmeans clustering for clusters 1 to nclust
  clust_out_list <- lapply(1:nclust,
    cluster_kmeans,
    data_matrices = data_matrices,
    max_dist = max_dist,
    how_cents = how_cents,
    bin_p_clust = bin_p_clust,
    threshold = threshold
  )
  # Get the AIC for each cluster set
  aicdf_list <- lapply(clust_out_list,
    get_aic
  )
  # Combine all aic results into one dataframe
  ic_df <- Reduce(rbind, aicdf_list)
  # Find which cluster set minimises the AIC
  min_cents <- ic_df$ncents[which.min(ic_df$aic)]
  # Get the dataframes for the minimising cluster set.
  cluster_df <- clust_out_list[[min_cents]]$clusters
  centroids_df <- clust_out_list[[min_cents]]$centres
  clust_dist_df <- clust_out_list[[min_cents]]$clust_dist
  # Set the resulting dataframes to the output.
  clust_out <- list("clusters" = cluster_df,
                    "clust_dist" = clust_dist_df,
                    "centres" = centroids_df)
  return(clust_out)
}