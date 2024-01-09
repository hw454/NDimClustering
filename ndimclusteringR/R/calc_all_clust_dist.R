#' Find the table of distances between point and cluster centroid.
#'
#' @description
#'   Find the distance between the snp and each cluster centroid.
#'
#' @param b_mat the matrix of the data. The rows correpond to the snps.
#' @param centroids_df dataframe of the centroid co-ordinates. Each row
#'   corresponds to a cluster number and the columns are the trait axes.
#'
#' @details
#'   The dataframe "snp_clust_dist_df"  has a row labelled by the "snp_id".
#'   The columns are the clust numbers, each value is the distance from the
#'   snp to that cluster.
#'
#' @return snp_clust_dist_df
#'
#' @family distance_functions
#'
#' @export
calc_all_clust_dist <- function(b_mat, centroids_df
) {
  calc_snp_all_clust_dist <- function(snp_id) {
    # Initialise cluster distance dataframe
    snp_clust_dist_df <- data.frame(
      row.names = snp_id
    )
    p_cols <- lapply(rownames(centroids_df),
      make_clust_col,
      rows = rownames(snp_clust_dist_df)
    )
    np_df <- Reduce(cbind, p_cols)
    snp_clust_dist_df <- cbind(
      snp_clust_dist_df,
      np_df
    )
    snp_score <- b_mat[snp_id, ]
    # Find the distance between the snp and all clusters. Store in
    # cluster distance dataframe.
    for (c_num in rownames(centroids_df)) {
      dist <- norm(data.matrix(
        stats::na.omit(centroids_df[c_num, ] - snp_score)
      ))
      c_col <- make_clust_col_name(c_num)
      snp_clust_dist_df[snp_id, c_col] <- dist
    }
    return(snp_clust_dist_df)
  }
  clust_dist_df_list <- lapply(rownames(b_mat),
    calc_snp_all_clust_dist
  )
  clust_dist_df <- Reduce(rbind, clust_dist_df_list)
  return(clust_dist_df)
}