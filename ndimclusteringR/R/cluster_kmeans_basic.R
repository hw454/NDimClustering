#' Cluster data with a standard kmeans approach
#'
#' @param data_list List of the matrices of data.
#'  Conatins the\:
#'    * "beta" matrix of assocaition data
#'    * "se" matrix of standard error data
#'    * "pval" matrix of p-value data
#'    * Rows correspond to snps, columns to traits
#' @param nclust Number of clusters to allocate
#' @param space_typ Whether to uses angles of spatial co-ordinates.
#'   If "angle" use angles found using [convert_point_mat_to_angle]
#' @param clust_prob_on Bool switch. If TRUE then calculate the snps
#'   probability of being in the cluster using [calc_clust_prob]
#' @param norm_typ The type of norm to use in distance calculations.
#'   The default is the Froebenius norm "F".
#' @param threshold The threshold with which clusters centres must differ
#'   by for clusters to be considered converged.
#' @param narm Boolean swithc. If TRUE (default) ignore NaNs in calculations.
#'
#' @description
#'   If space_typ == "angle" then data is converted to angles.
#'   The points are clustered using [km_nan] and "nclust" clusters.
#'   The "cluster_df" dataframe labelled
#'     "clusters" of the cluster membership with columns\:
#'     * "clust_num" the number of clusters
#'     * "clust_dist" the distance from the snp to the cluster centre
#'       (or distance between angles)
#'     * "clust_prob" probability the snp is in the cluster. Calculated
#'       using [calc_clust_prob].
#'
#' @return clusters_df
#'
#' @family clustering_functions
#'
#' @export
cluster_kmeans_basic <- function(data_list,
                                nclust = 10,
                                space_typ = "regular",
                                clust_prob_on = TRUE,
                                norm_typ = "F",
                                threshold = 1e-5,
                                narm = TRUE) {
  # Using the association scores for each SNP accross traits cluster the traits
  # using kmeans. Return the cluster setup which minimises AIC.

  b_mat <- data_list$beta
  se_mat <- data_list$se
  # Crop the data to those with the lowest standard error.
  # Add remaining terms to closest cluster once centres have been found.
  crop_se <- docore::lim(se_mat, 0, 2)
  crop_se <- crop_se[stats::complete.cases(crop_se), ]
  crop_snp_list <- rownames(crop_se)

  # Filter NaNs before clustering
  if (space_typ == "angle") {
    # For each point in b_df_comp convert the score to the angle
    # between the vectors to the origin and the unit vectors on the axis.
    b_mat_clust <- mat_to_angle_mat(b_mat)
  } else {
    b_mat_clust <- b_mat
  }
  # Crop the data to focus on the snps with the lowest standard error
  b_mat_crop <- b_mat_clust[crop_snp_list, ]

  # Initial cluster dataframe
  clust_out <- km_nan(b_mat_crop,
                    nclust = nclust,
                    clust_threshold = threshold,
                    norm_typ = norm_typ,
                    prob_on = clust_prob_on,
                    na_rm = narm)
  # cluster number identification for the snps with higher standard error.
  nan_snp_list <- lapply(setdiff(rownames(b_mat_clust), crop_snp_list),
                            find_closest_clust_snp,
                            b_mat = b_mat_clust,
                            cluster_df = clust_out$clusters,
                            centroids_df = clust_out$centres,
                            norm_typ = norm_typ)
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
  if (clust_prob_on) {
    nan_cluster_df$clust_prob <- calc_clust_prob(nan_cluster_df$clust_dist)
  }
  # ADDFEATURE - Assign junk clusters.
  clust_out$clusters <- rbind(clust_out$clusters, nan_cluster_df)
  clust_out$clust_dist <- rbind(clust_out$clust_dist, nan_clust_dist_df)
  return(clust_out)
}