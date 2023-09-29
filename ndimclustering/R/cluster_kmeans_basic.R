#' Cluster data with a standard kmeans approach
#'
#' @param b_mat Matrix of data. Rows correspond to snps.
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
cluster_kmeans_basic <- function(b_mat,
                                nclust = 10,
                                space_typ = "regular",
                                clust_prob_on = TRUE,
                                norm_typ = "F",
                                threshold = 1e-5,
                                narm = TRUE) {
  # Using the association scores for each SNP accross traits cluster the traits
  # using kmeans. Return the cluster setup which minimises AIC.

  # Filter NaNs before clustering
  if (space_typ == "angle") {
    # For each point in b_df_comp convert the score to the angle
    # between the vectors to the origin and the unit vectors on the axis.
    b_mat_clust <- mat_to_angle_mat(b_mat)
  } else {
    b_mat_clust <- b_mat
  }
  # max_dist <- max_dist_calc(b_mat_clust,
  #                          norm_typ = norm_typ,
  #                          na_rm = narm)
  # b_mat_comp <- b_mat_clust[complete.cases(b_mat_clust), ]

  # Initial cluster dataframe
  clust_out <- km_nan(b_mat_clust,
                    nclust = nclust,
                    clust_threshold = threshold,
                    norm_typ = norm_typ,
                    prob_on = clust_prob_on,
                    na_rm = narm)
  # cluster number identification for each observation
  # snp_cluster_list <- lapply(setdiff(names(b_mat_clust), names(b_mat_comp)),
  #                          find_closest_clust_snp,
  #                          b_mat = b_mat_clust,
  #                          cluster_df = clust_out$clusters,
  #                          centroids_df = clust_out$centres,
  #                          norm_typ = norm_typ,
  #                          max_dist = max_dist)
  # nan_cluster_df <- Reduce(rbind, snp_cluster_list)
  # if (clust_prob_on) {
  #  nan_cluster_df$clust_prob <- calc_clust_prob(nan_cluster_df$clust_dist)
  # }
  # clust_out <- rbind(clust_out, nan_cluster_df)
  return(clust_out)
}