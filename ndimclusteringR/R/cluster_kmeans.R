#' Cluster the data using kmeans then minimising aic.
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
#' @description
#'   If space_typ == "angle" then data is converted to angles.
#'   for i=1:(nr+1)
#'      The points are clustered using [km_nan] and i clusters
#'      find the aic for each cluster using [find_ic].
#'   The set of clusters which minimises the aic is chosen.
#'   The "clusters_df" dataframe of the cluster membership with columns\:
#'     * "clust_num" the number of clusters
#'     * "clust_dist" the distance from the snp to the cluster centre
#'       (or distance between angles)
#'     * "clust_prob" probability the snp is in the cluster. Calculated
#'       using [calc_clust_prob].
#'
#' @return clusters_df
#'
#' @family cluster_functions
#'
#' @export
cluster_kmeans <- function(data_matrices,
  nclust = 10,
  max_dist = 10.0,
  how_cents = "rand",
  bin_p_clust = TRUE,
  threshold = 1e-5
) {
  # Using the association scores for each SNP accross traits cluster the traits
  # using kmeans. Return the cluster setup which minimises AIC.
  b_mat <- data_matrices$beta
  se_mat <- data_matrices$se
  p_mat <- data_matrices$pval
  # Crop the data to the complete cases
  crop_se <- se_mat[stats::complete.cases(se_mat), ]
  crop_snp_list <- rownames(crop_se)

  # Crop the data to focus on the snps with the lowest standard error
  b_mat_crop <- b_mat[crop_snp_list, ]
  p_mat_crop <- p_mat[crop_snp_list, ]

  # Cluster the data
  clust_out <- km_nan(b_mat_crop,
    p_mat = p_mat_crop,
    nclust = nclust,
    clust_threshold = threshold,
    bin_p_clust = bin_p_clust,
    how_cents = how_cents
  )

  # Cluster number identification for each observation
  nan_snp_list <- lapply(setdiff(rownames(b_mat), crop_snp_list),
                         find_closest_clust_snp,
                         b_mat = b_mat,
                         cluster_df = clust_out$clusters,
                         centroids_df = clust_out$centres)
  f1 <- function(x) {
    x$clusters
  }
  f2 <- function(x) {
    x$clust_dist
  }
  snp_cluster_list <- lapply(nan_snp_list, f1)
  snp_c_dist_list <- lapply(nan_snp_list, f2)
  nan_cluster_df <- Reduce(rbind, snp_cluster_list)
  nan_clust_dist_df <- Reduce(rbind, snp_c_dist_list)
  if (bin_p_clust) {
    nan_cluster_df$clust_prob <- calc_clust_prob(nan_cluster_df$clust_dist)
  }
  clust_out$clusters <- rbind(clust_out$clusters, nan_cluster_df)
  clust_out$clust_dist <- rbind(clust_out$clust_dist, nan_clust_dist_df)
  # Add number of axis to the dataframe to combine outputs from iterative case.
  num_axis <- ncol(b_mat_crop)
  clust_out$clusters["num_axis"] <- num_axis
  clust_out$clust_dist["num_axis"] <- num_axis
  clust_out$centres["num_axis"] <- num_axis
  # ADDFEATURE - Assign junk clusters.
  return(clust_out)
}