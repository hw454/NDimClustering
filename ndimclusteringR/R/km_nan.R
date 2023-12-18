#' Apply kmeans clustering ignoring NaN terms
#'
#' @description
#' Use common axes for distance calculations to igore NaNs and
#' retain data.
#'
#' @param b_mat Matrix of the data, each row is a snp and the columns
#'   are the traits.
#' @param p_mat matrix of p-values corresponding to the scores in b_mat
#' @param nclust The number of clusters to allocate. default\:5
#' @param iter_max The maximum number of iterations to reach cluster
#'   convergence. default\:300
#' @param clust_threshold The threshold for how close cluster
#'   centres need to be to be considered converged. default\:1e-5
#' @param bin_p_clust Bool to indicate whether cluster probability is to be
#'   assign. default\:TRUE
#' @param how_cents How the centroids will be initialised. If "rand (default)"
#'   then the coordinates are created using a uniform distribution on the range
#'   of each axis. If "points" then the centroid coordinates are assign to
#'   nclust random points from the dataspace.
#'
#' @details
#' The cluster centres are randomly assigned using [make_rand_cent].
#' The cluster centres are check and reassigned using [check_clust_cent],
#' When the clusters are converged or the maximum number of iterations are
#' used then the clusters are return in the dataframe "cluster_df".
#' "clust_out" is the list containing the "clusters_df" dataframe
#' labelled "cluster" of the cluster membership with columns\:
#'   * "clust_num" the number of clusters
#'   * "clust_dist" the distance from the snp to the cluster centre
#'   (or distance between angles)
#'   * "clust_prob" probability the snp is in the cluster. Calculated
#'   using [calc_clust_prob].
#' "clust_out" also contains the "centroids_df" dataframe labelled "centres"
#' whose columns are the traits and rows are the cluster numbers. Additional
#' column "thresh_check" contains Bool indicating whether that cluster
#' converged.
#'
#' @return clust_out
#'
#' @export
km_nan <- function(b_mat, p_mat,
  nclust = 5, iter_max = 300, clust_threshold = 1e-5,
  bin_p_clust = TRUE, how_cents = "point"
) {
  snp_list <- rownames(b_mat)
  # Generate data frame with max and min data.
  # Min and Max are columns. Rows for each column of b_mat
  if (grepl("rand", how_cents)) {
    min_max_df <- data.frame(row.names = colnames(b_mat),
      min = apply(b_mat, 2, min),
      max = apply(b_mat, 2, max)
    )
    min_max_df <- stats::na.omit(min_max_df)
    # Randomly assign centroids coords per cluster.
    centroid_list <- lapply(colnames(b_mat), make_rand_cent,
                            n_cents = nclust, min_max_df = min_max_df)
    centroids_df <- Reduce(cbind, centroid_list)
  } else if (grepl("point", how_cents)) {
    if (nclust == 1) {
      b_mat_tmp <- as.data.frame(b_mat)
      init_cent_id <- sample(snp_list, nclust, replace = FALSE)
      centroids_df <- b_mat_tmp[init_cent_id, ]
      rownames(centroids_df) <- 1
    } else {
      init_cent_id <- sample(snp_list, nclust, replace = FALSE)
      centroids_df <- as.data.frame(b_mat[init_cent_id, ])
      rownames(centroids_df) <- seq_len(nclust)
    }
    colnames(centroids_df) <- colnames(b_mat)
  }
  # Randomly assign cluster to snps
  nsnps <- length(snp_list)
  clust_samp <- sample(rownames(centroids_df), nsnps, replace = TRUE)
  cluster_df <- data.frame(
    row.names = snp_list,
    clust_num = clust_samp
  )
  cluster_df["clust_prob"] <- numeric()
  clust_dist_mem_list <- lapply(snp_list, calc_member_dist_cent,
                                b_mat = b_mat,
                                cluster_df = cluster_df,
                                centroids_df = centroids_df)
  clust_dist_mem_df <- Reduce(rbind, clust_dist_mem_list)
  cluster_df <- cbind(rn = rownames(cluster_df),
                      cluster_df,
                      clust_dist_mem_df,
                      row.names = NULL)
  cluster_df <- tibble::column_to_rownames(cluster_df, var = "rn")
  for (iter in 1:iter_max){
    # For each SNP find the cluster with the closest centre.
    snp_clust_list <- lapply(snp_list, find_closest_clust_snp,
                             b_mat = b_mat,
                             cluster_df = cluster_df,
                             centroids_df = centroids_df)
    # Combine the list of dataframes into one dataframe.
    # Override Cluster_df with the new assignment
    df_cols <- function(df, col) {
      df[[col]]
    }
    cluster_df_list <- lapply(snp_clust_list,
                              df_cols,
                              col = "clusters")
    cluster_df <- Reduce(rbind, cluster_df_list)
    clust_dist_df_list <- lapply(snp_clust_list,
                                 df_cols,
                                 col = "clust_dist")
    clust_dist_df <- Reduce(rbind, clust_dist_df_list)
    # Recompute the centroids based on the average of the clusters.
    # Check if the previous centres differ from the cluster means.
    thresh_list <- lapply(rownames(centroids_df), check_clust_cent,
                          clustnum_df = cluster_df$clust_num,
                          b_mat = b_mat,
                          centroids_df = centroids_df,
                          p_mat = p_mat,
                          clust_threshold = clust_threshold,
                          bin_p_clust = bin_p_clust)
    thresh_check_df <- Reduce(rbind, thresh_list)
    # Are all the new centroids within the threshold of the previous
    if (all(thresh_check_df$thresh_check)) {
      print(paste("Clusters converged", iter))
      break
    }
    centroids_df <- thresh_check_df[, !names(thresh_check_df)
                                    %in% c("thresh_check")]
    # Are all the new centroids within the threshold of the previous
  }
  # Calculate the probability each snp is in the cluster
  cluster_df$clust_prob <- calc_clust_prob(cluster_df)
  clust_out <- list("clusters" = cluster_df, "centres" = centroids_df,
                    "clust_dist" = clust_dist_df)
  return(clust_out)
}