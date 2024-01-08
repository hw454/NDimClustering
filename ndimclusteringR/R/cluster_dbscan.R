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
#' @param eps The density deining a cluster. default \: 0.4
#' @param bin_p_clust Should the cluster centroid be weighted by p-values. default \: FALSE
#'
#' @description
#'   If space_typ == "angle" then data is converted to angles.
#'   for i=1:(nr+1)
#'      The points are clustered using [km_nan] and i clusters
#'      find the aic for each cluster using [find_ic].
#'   The set of clusters which minimises the aic is chosen.
#'   Outputs: clustout = [clusters, clust_dist, centres]
#'   The "clusters_df" dataframe of the cluster membership with columns\:
#'     * "clust_num" the number of clusters
#'     * "clust_dist" the distance from the snp to the cluster centre
#'       (or distance between angles)
#'     * "clust_prob" probability the snp is in the cluster. Calculated
#'       using [calc_clust_prob].
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
#' @family dbscan_functions
#'
#' @export
cluster_dbscan <- function(data_matrices,
  eps = 0.4,
  bin_p_clust = FALSE
) {
  # Using the association scores for each SNP accross traits cluster the traits
  # using dbscan.
  b_mat <- data_matrices$beta_pc
  se_mat <- data_matrices$se_pc
  p_mat <- data_matrices$pval_pc
  # Crop the data to the complete cases
  crop_se <- se_mat[stats::complete.cases(se_mat), ]
  crop_snp_list <- rownames(crop_se)
  b_mat_crop <- b_mat[crop_snp_list, ]
  p_mat_crop <- p_mat[crop_snp_list, ]

  # Cluster the data
  clust_dbscan <- dbscan::dbscan(b_mat_crop, eps = eps)
  cluster_df <- data.frame(row.names = crop_snp_list,
    clust_num = clust_dbscan$cluster
  )
  # Find the number of clusters identified
  nclust <- length(unique(cluster_df$clust_num))
  # Calculate cluster centres
  centroids_df <- calc_clust_cent(cluster_df,
    b_mat = b_mat_crop,
    p_mat = p_mat_crop,
    bin_p_clust = bin_p_clust
  )
  # Calculate cluster distances
  clust_dist_df <- calc_all_clust_dist(centroids_df,
    b_mat = b_mat_crop
  )
  # Add member distance to clust dataframe
  mem_dist_list <- lapply(crop_snp_list,
    calc_member_dist_cent,
    b_mat = b_mat_crop,
    cluster_df = cluster_df,
    centroids_df = centroids_df
  )
  mem_dist_df <- Reduce(rbind, mem_dist_list)
  cluster_df <- cbind(
    rn = rownames(cluster_df),
    cluster_df,
    mem_dist_df,
    row.names = NULL
  )
  cluster_df <- tibble::column_to_rownames(cluster_df, var = "rn")
  cluster_df$clust_prob <- calc_clust_prob(cluster_df)
  # Assign all cluster info to clust_out
  clust_out <- list("clusters" = cluster_df,
    "centres" = centroids_df,
    "clust_dist" = clust_dist_df
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
  nan_cluster_df$clust_prob <- calc_clust_prob(nan_cluster_df)
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
  # Convert row names to columns to stacking later results
  clust_out$clusters <- tibble::rownames_to_column(clust_out$clusters,
    var = "snp_id"
  )
  clust_out$clust_dist <- tibble::rownames_to_column(clust_out$clust_dist,
    var = "snp_id"
  )
  return(clust_out)
}