#' Cluster the data using kmeans then minimising aic.
#' @param b_mat The data matrix
#' @param nr The maximum number of clusters to consider. default\:10
#' @param max_dist The maximum distance between any two points.
#' @param space_typ The spatial structure of the data when clustering.
#'   If "angle" then the data is converted to angles. If "regular" the
#'   data is left the same.
#' @param clust_prob_on Bool switch. If TRUE (default) then cluster
#'   probability if used to weight the cluster scores.
#' @param norm_typ The type of norm to be used for distance calculations.
#'   The default is the Froebenius norm "F".
#' @param threshold The threshold for distance between cluster centres
#'   for clusters to be considered converged.
#' @param narm Bool to indicate how to handle NaNs. If TRUE (default)
#'   NaNs are ignored in calculations.
#' @description
#' If space_typ == "angle" then data is converted to angles.
#' for i=1:(nr+1)
#'      The points are clustered using \link{km_nan} and i clusters
#'      find the aic for each cluster using \link{find_ic}
#' The set of clusters which minimises the aic is chosen.
#' The "clusters_df" dataframe of the cluster membership with columns\:
#'   * "clust_num" the number of clusters
#'   * "clust_dist" the distance from the snp to the cluster centre
#'   (or distance between angles)
#'   * "clust_prob" probability the snp is in the cluster. Calculated
#'   using \link{calc_clust_prob}.
#' @return clusters_df
#' @export
cluster_kmeans_min <- function(b_mat,
                              nr = 10,
                              max_dist = 10.0,
                              space_typ = "regular",
                              clust_prob_on = TRUE,
                              norm_typ = "F",
                              threshold = 1e-5,
                              narm = TRUE) {
  #' Using the association scores for each SNP accross traits cluster the traits
  #' using kmeans. Return the cluster setup which minimises AIC.

  if (space_typ == "angle") {
    # For each point in b_df_comp convert the score to the angle
    # between the vectors to the origin and the unit vectors on the axis.
    b_mat_clust <- mat_to_angle_mat(b_mat)
  } else {
    b_mat_clust <- b_mat
  }
  #max_dist <- max_dist_calc(b_mat_clust,
  #                          norm_typ = norm_typ,
  #                          na_rm = narm)
  # Filter complete cases
  #b_mat_comp <- b_mat_clust[complete.cases(b_mat_clust), ]

  # Initial cluster dataframe
  clust_re_list <- lapply(1:(nclust+1),km_nan,
                    data_mat = b_mat,
                    nclust = cents,
                    clust_threshold = clust_thres,
                    norm_typ = clust_norm,
                    na_rm = na_rm,
                    prob_on = clust_prob_on)
  ic_list <- lapply(clust_re_list, find_all_ic,
                    num_axis = ncol(b_mat))
  ic_df <- Reduce(rbind, ic_list)
  # Find the number of centres that minimizes the AIC
  min_cents <- ic_df$ncents[which.min(ic_df$aic)]
  # Use the number of centres to locate corresponding clusters since there
  # maybe variations due to machine precision in the aic values.
  clust_re_min_aic <- ic_df[which(ic_df$ncents == min_cents), ]
  if (clust_prob_on) {
    clust_out$clust_prob <- calc_clust_prob(clust_out$clust_dist)
  }
  return(clust_out)
}