#' Find the cluster whose centroid is closest to snp_id
#'
#' @description
#'   Find the distance between the snp and each cluster centroid.
#'   If the distance is less then the currently assigned cluster then
#'   reassign the cluster number and cluster distance.
#'
#' @param snp_id the label for the snp
#' @param b_mat the matrix of the data. The rows correpond to the snps.
#' @param cluster_df dataframe of cluster membership. Each row is a snp,
#'   it has columns for\:
#'    * "clust_num" the number of the allocated cluster.
#'    * "clust_dist" the distane to the centre of the cluster.
#'    * "clust_prob" probability the snp is in the assigned cluster.
#' @param centroids_df dataframe of the centroid co-ordinates. Each row
#'   corresponds to a cluster number and the columns are the trait axes.
#' @param max_dist The maximum distance between any two points. Used
#'   for initialisation if snp is not already in the cluster assignment.
#'
#' @details
#'   The dataframe "snp_clust_df"  has a row labelled by the "snp_id".
#'   The columns are\:
#'   * "clust_num" The number of the closest cluster
#'   * "clust_dist" The distance from the snp to the centroid of that cluster.
#'   * "clust_prob" probability the snp is in the cluster it's been assigned.
#'   The dataframe "snp_clust_dist_df"  has a row labelled by the "snp_id".
#'   The columns are the clust numbers, each value is the distance from the
#'   snp to that cluster.
#'
#'   out_list = ("clusters" = snp_clust_df, "clust_dist" = snp_clust_dist_df)
#'
#' @return out_list
#'
#' @family distance_functions
#'
#' @export
find_closest_clust_snp <- function(snp_id, b_mat, cluster_df, centroids_df,
  max_dist = 10.0
) {
  # Initialise cluster distance dataframe
  snp_clust_dist_df <- data.frame(
    row.names = snp_id
  )
  p_cols <- lapply(rownames(centroids_df),
                   make_clust_col,
                   rows = row.names(snp_clust_dist_df))
  np_df <- Reduce(cbind, p_cols)
  snp_clust_dist_df <- cbind(
    rn = rownames(snp_clust_dist_df),
    snp_clust_dist_df,
    np_df,
    row.names = NULL
  )
  snp_clust_dist_df <- tibble::column_to_rownames(
    snp_clust_dist_df,
    var = "rn"
  )
  snp_score <- b_mat[snp_id, ]
  # Initialise distance and cluster number
  if (!(snp_id %in% row.names(cluster_df))) {
    snp_dist <- max_dist
    snp_clust_num <- 0
  } else {
    snp_dist <- cluster_df[snp_id, "clust_dist"]
    snp_clust_num <- cluster_df[snp_id, "clust_num"]
  }
  # Find the distance between the snp and all clusters. Store in
  # cluster distance dataframe and return the distance and cluster
  # number for the closest cluster.
  for (c_num in rownames(centroids_df)) {
    dist <- norm(data.matrix(
      stats::na.omit(centroids_df[c_num, ] - snp_score)
    ))
    c_col <- make_clust_col_name(c_num)
    snp_clust_dist_df[snp_id, c_col] <- dist
    if (dist < snp_dist) {
      snp_dist <- dist
      snp_clust_num <- c_num
    }
  }
  snp_cluster_df <- data.frame(
    row.names = snp_id,
    clust_num = snp_clust_num,
    clust_dist = snp_dist,
    clust_prob = 1.0
  )
  out_list <- list("clusters" = snp_cluster_df,
                   "clust_dist" = snp_clust_dist_df)
  return(out_list)
}