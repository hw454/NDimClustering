#' Apply kmeans clustering ignoring NaN terms
#'
#' @description
#' Use common axes for distance calculations to igore NaNs and
#' retain data.
#'
#' @param b_mat Matrix of the data, each row is a snp and the columns
#'   are the traits.
#' @param nclust The number of clusters to allocate. default\:5
#' @param iter_max The maximum number of iterations to reach cluster
#'   convergence. default\:300
#' @param clust_threshold The threshold for how close cluster
#'   centres need to be to be considered converged. default\:1e-5
#' @param norm_typ The type of norm to be used in distance calculation.
#'   The default is the Froebenius norm "F".
#' @param na_rm Bool to indicate how to deal with NaNs. If TRUE (default)
#'   then NaNs are ignored.
#' @param prob_on Bool to indicate whether cluster probability is to be
#'   assign. default\:TRUE
#'
#' @details
#' The cluster centres are randomly assigned using \link{make_rand_cent}.
#' The cluster centres are check and reassigned using \link{check_clust_cent},
#' When the clusters are converged or the maximum number of iterations are
#' used then the clusters are return in the dataframe "cluster_df".
#' "clust_out" is the list containing the "clusters_df" dataframe
#' labelled "cluster" of the cluster membership with columns\:
#'   * "clust_num" the number of clusters
#'   * "clust_dist" the distance from the snp to the cluster centre
#'   (or distance between angles)
#'   * "clust_prob" probability the snp is in the cluster. Calculated
#'   using \link{calc_clust_prob}.
#' "clust_out" also contains the "centroids_df" dataframe labelled "centres"
#' whose columns are the traits and rows are the cluster numbers. Additional
#' column "thresh_check" contains Bool indicating whether that cluster
#' converged.
#'
#' @return clust_out
#'
#' @export
km_nan <- function(b_mat,
                  nclust = 5,
                  iter_max = 300,
                  clust_threshold = 1e-5,
                  norm_typ = "F",
                  na_rm = TRUE,
                  prob_on = TRUE){
  set.seed(123)
  snp_list <- rownames(b_mat)
  # Generate data frame with max and min data.
  # Min and Max are columns. Rows for each column of b_mat
  min_max_df <- data.frame(
    row.names = colnames(b_mat),
    min = apply(b_mat, 2, min, na.rm = na_rm),
    max = apply(b_mat, 2, max, na.rm = na_rm)
  )
  min_max_df <- na.omit(min_max_df)
  # Randomly assign central coords per cluster.
  centroid_list <- lapply(rownames(min_max_df), make_rand_cent,
                          n_cents = nclust, min_max_df = min_max_df)
  centroids_df <- Reduce(cbind, centroid_list)
  # Randomly assign cluster to snps
  nsnps <- length(snp_list)
  clust_samp <- replicate(nsnps, sample(1:nclust, 1))
  cluster_df <- data.frame(
    row.names = snp_list,
    clust_num = clust_samp
  )
  cluster_df["clust_prob"] <- numeric()
  clust_dist_list <- lapply(snp_list, calc_dist_cent,
                          b_mat = b_mat,
                          cluster_df = cluster_df,
                          centroids_df = centroids_df,
                          norm_typ = norm_typ)
  clust_dist_df <- Reduce(rbind, clust_dist_list)
  cluster_df <- cbind(cluster_df, clust_dist_df)
  for (iter in 1:iter_max){
    # For each SNP find the cluster with the closest centre.
    snp_clust_list <- lapply(snp_list, find_closest_clust_snp,
                            b_mat = b_mat,
                            cluster_df = cluster_df,
                            centroids_df = centroids_df,
                            norm_typ = norm_typ)
    # Combine the list of dataframes into one dataframe.
    # Override Cluster_df with the new assignment
    cluster_df <- Reduce(rbind, snp_clust_list)
    # Recompute the centroids based on the average of the clusters.
    # Check if the previous centres differ from the cluster means.
    thresh_list <- lapply(rownames(centroids), check_clust_cent,
                        clustnum_df = cluster_df$clust_num,
                        b_mat = b_mat, centroids_df = centroids,
                        na_rm = na_rm, norm_typ = norm_typ,
                        clust_threshold = clust_threshold)#,
                        #axes_nonnan=axes_nonnan)
    thresh_check_df <- Reduce(rbind, thresh_list)
    # Are all the new centroids within the threshold of the previous
    if (all(thresh_check_df$thresh_check)) {
      print(paste("Clusters converged", iter))
      break
    }
    centroids_df <- thresh_check_df[,
                  !names(thresh_check_df) %in% c("thresh_check")]
    # Are all the new centroids within the threshold of the previous
  }
  if (prob_on) {
    cluster_df$clust_prob <- clust_prob_calc(cluster_df$clust_dist)
  }
 clust_out <- list("clusters" = cluster_df, "centres" = centroids_df)
 return(clust_out)
}