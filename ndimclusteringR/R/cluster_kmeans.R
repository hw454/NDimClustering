#' Cluster the data using kmeans then minimising aic.
#'
#' @param data_list List of the matrices of data.
#'  Conatins the\:
#'    * "beta" matrix of assocaition data
#'    * "se" matrix of standard error data
#'    * "pval" matrix of p-value data
#'    * Rows correspond to snps, columns to traits
#' @param iter_traits - list of terms indicating the type of program.
#'   * "iter" - integer, (default 0)
#'   * "bp_on" - TRUE (default) if probability of scores is to be used.
#'   * "clust_prob_on" - TRUE\:default switch for using prob of being in cluster
#'   * "clust_typ" - default\:"basic", the clustering method to use.
#' @param nclust The maximum number of clusters to consider. default\:10
#' @param max_dist The maximum distance between any two points.
#' @param space_typ The spatial structure of the data when clustering.
#'   If "angle" then the data is converted to angles. If "regular" the
#'   data is left the same.
#' @param clust_typ String describing the clsutering method.
#'  * "basic" then find fixed number of clusters
#'  * "min" then minimise the aic
#' @param how_cents How the centroids will be initialised. If "rand (default)"
#'   then the coordinates are created using a uniform distribution on the range
#'   of each axis. If "points" then the centroid coordinates are assign to
#'   nclust random points from the dataspace.
#' @param clust_prob_on Bool switch. If TRUE (default) then cluster
#'   probability if used to weight the cluster scores.
#' @param norm_typ The type of norm to be used for distance calculations.
#'   The default is the Froebenius norm "F".
#' @param threshold The threshold for distance between cluster centres
#'   for clusters to be considered converged.
#' @param narm Bool to indicate how to handle NaNs. If TRUE (default)
#'   NaNs are ignored in calculations.
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
cluster_kmeans <- function(data_list, iter_traits,
                           nclust = 10,
                           max_dist = 10.0,
                           space_typ = "regular",
                           clust_typ = "basic",
                           how_cents = "rand",
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
  #crop_se_snp_list <- find_low_percent_col_mat(se_mat,
  #                                           # col_percent = 0.1,
  #                                           # na_rm = narm)
  crop_se <- se_mat #se_mat[crop_se_snp_list, ]
  crop_se <- crop_se[stats::complete.cases(crop_se), ]
  crop_snp_list <- rownames(crop_se)

  # Crop the data to focus on the snps with the lowest standard error
  b_mat_crop <- b_mat[crop_snp_list, ]
  if (grepl("min", clust_typ)) {
    # Initial cluster dataframe
    clust_re_list <- lapply(1:(nclust + 1), km_nan,
                            b_mat = b_mat_crop,
                            clust_threshold = threshold,
                            norm_typ = norm_typ,
                            na_rm = narm,
                            prob_on = clust_prob_on,
                            how_cents = how_cents)
    ic_list <- lapply(clust_re_list, find_all_ic,
                      num_axis = ncol(b_mat))
    ic_df <- Reduce(rbind, ic_list)
    # Find the number of centres that minimizes the AIC
    min_cents <- ic_df$ncents[which.min(ic_df$aic)]
    centroids_df <- clust_re_list[[min_cents]]$centres
    clust_dist_df <- clust_re_list[[min_cents]]$clust_dist
    # Use the number of centres to locate corresponding clusters since there
    # maybe variations due to machine precision in the aic values.
    cluster_df <- ic_df[which(ic_df$ncents == min_cents), ]
    rownames(cluster_df) <- NULL
    cluster_df <- tibble::column_to_rownames(cluster_df, var = "snp_id")
    clust_out <- list("clusters" = cluster_df,
                      "centres" = centroids_df,
                      "clust_dist" = clust_dist_df)
  } else if (grepl("basic", clust_typ)) {
    # Initial cluster dataframe
    clust_out <- km_nan(b_mat_crop,
                        nclust = nclust,
                        clust_threshold = threshold,
                        norm_typ = norm_typ,
                        prob_on = clust_prob_on,
                        na_rm = narm,
                        how_cents = how_cents)
  } else if (grepl("mrclust", clust_typ)) {
    col1 <- colnames(b_mat)[1]
    col2 <- colnames(b_mat)[2]
    bx <- b_mat_crop[, col1]
    by <- b_mat_crop[, col2]
    bxse <- crop_se[, col1]
    byse <- crop_se[, col2]
    ratio_est <- by / bx
    ratio_est_se <- byse / abs(bx)
    res_em <- mrclust::mr_clust_em(theta = ratio_est,
                                   theta_se = ratio_est_se,
                                   bx = bx,
                                   by = by,
                                   bxse = bxse,
                                   byse = byse,
                                   obs_names = crop_snp_list,
                                   junk_sd = 0.5)
    print(res_em)
  }
  # cluster number identification for each observation
  nan_snp_list <- lapply(setdiff(rownames(b_mat), crop_snp_list),
                         find_closest_clust_snp,
                         b_mat = b_mat,
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
  clust_out$clusters <- rbind(clust_out$clusters, nan_cluster_df)
  clust_out$clust_dist <- rbind(clust_out$clust_dist, nan_clust_dist_df)
  #DEBUG START
  clust_out$clusters <- tibble::rownames_to_column(clust_out$clusters, "snp_id")
  c1 <- colnames(b_mat_crop)[1]
  c2 <- colnames(b_mat_crop)[2]
  num_axis <- ncol(b_mat_crop)
  clust_out$clusters["num_axis"] <- num_axis
  if (space_typ == "angle") {
    plot_clust_angle_scatter_test(clust_out$clusters, b_mat,
                                  se_mat,
                                  iter_traits,
                                  c1,
                                  c2,
                                  num_axis = num_axis)
  }
  clust_out$clusters <- tibble::column_to_rownames(clust_out$clusters, "snp_id")
  #DEBUG END
  # ADDFEATURE - Assign junk clusters.
  return(clust_out)
}