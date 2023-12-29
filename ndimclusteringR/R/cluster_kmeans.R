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
#' 1. Data is cropped to the complete cases.
#' 2. The points are clustered using a k-means clustering [km_nan].
#' 3. The data removed in the crop is assigned to the closest cluster.
#' 4. The "clusters_df" dataframe of the cluster membership with columns\:
#'     * "clust_num" the number of clusters
#'     * "clust_dist" the distance from the snp to the cluster centre
#'       (or distance between angles)
#'     * "clust_prob" probability the snp is in the cluster. Calculated
#'       using [calc_clust_prob].
#'   clust_out = [clusters, clust_dist, centres].
#'   clust_out$clusters is cluster_df.
#'   The clust_dist dataframe has columns corresponding to the clusters,
#'     the rows are the datapoints and each value is the distance from that
#'     point to the cluster centre.
#'   The centres dataframe has columns corresponding to the axis and the rows
#'     corresponding to the clusters. Each row is the co-ordinate of the cluster
#'     centroid.
#'
#' @return clusters_df
#'
#' @family cluster_functions
#' @family k_means
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
  b_mat <- data_matrices$beta_pc
  se_mat <- data_matrices$se_pc
  p_mat <- data_matrices$pval_pc
  # Crop the data to the complete cases
  crop_se <- se_mat[stats::complete.cases(se_mat), ]
  crop_snp_list <- rownames(crop_se)
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
  # Get the cluster assignments
  snp_cluster_list <- lapply(nan_snp_list, f1)
  # Get the cluster distances
  snp_c_dist_list <- lapply(nan_snp_list, f2)
  # Convert the lists into dataframes
  nan_cluster_df <- Reduce(rbind, snp_cluster_list)
  nan_clust_dist_df <- Reduce(rbind, snp_c_dist_list)
  # Calculate the probability of being in each cluster.
  nan_cluster_df$clust_prob <- calc_clust_prob(nan_cluster_df$clust_dist)
  # Combine the kmeans cluster assignment with the additional clusters
  clust_out$clusters <- rbind(clust_out$clusters, nan_cluster_df)
  clust_out$clust_dist <- rbind(clust_out$clust_dist, nan_clust_dist_df)
  # Add number of axis to the dataframe to combine outputs from iterative case.
  num_axis <- ncol(b_mat_crop)
  clust_out$clusters["num_axis"] <- num_axis
  clust_out$clust_dist["num_axis"] <- num_axis
  clust_out$centres["num_axis"] <- num_axis
  clust_out$clusters["ncents"] <- nclust
  clust_out$clust_dist["ncents"] <- nclust
  clust_out$centres["ncents"] <- nclust
  # Convert row names to columns for stacking later results
  clust_out$clusters <- tibble::rownames_to_column(clust_out$clusters,
    var = "snp_id"
  )
  clust_out$clust_dist <- tibble::rownames_to_column(clust_out$clust_dist,
    var = "snp_id"
  )
  # ADDFEATURE - Assign junk clusters.
  return(clust_out)
}