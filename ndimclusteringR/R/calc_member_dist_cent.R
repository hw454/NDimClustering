#' Calculate the distance between the snps and their allocated cluster centroid
#'
#' @param snp_id the label for the snp
#' @param b_mat the data matrix used for the coordinate of the snp.
#'   Rows are snps, columns are the trait axes.
#' @param cluster_df The dataframe with the snp cluster memberships.
#'   Columns are\:
#'    * "clust_num" the cluster number.
#'    * "clust_dist" distance from the snp to the cluster centroid.
#'    * "clust_prob" probability the snp is correctly allocated to the cluster.
#'   Rows are each snp.
#' @param centroids_df The dataframe of the centroid co-ordinates for each
#'   cluster. The rows are each cluster number and the columns are the traits.
#' @param norm_typ The type of norm to use in the distance calculation. The
#'   default is the Froebenius norm "F".
#'
#' @description Calculate the distance between the point in "b_mat"
#'   corresponding to "snp_id" and the cluster centroid allocated to "snp_id".
#'   This is stored in a dataframe "snp_clust_df" which has rowname "snp_id",
#'   and column "clust_dist", with the calculated distance.
#'
#' @return snp_clust_df
#'
#' @family distance_functions
#' @family centroid_functions
#' @family clustering_components
#'
#' @export
calc_member_dist_cent <- function(snp_id, b_mat, cluster_df, centroids_df,
  norm_typ = "F"
) {
  c_num <- cluster_df[snp_id, "clust_num"]
  cent <- centroids_df[c_num, ]
  snp_score <- b_mat[snp_id, ]
  c_dist <- norm(data.matrix(stats::na.omit(cent - snp_score)), norm_typ)
  snp_clust_df <- data.frame(row.names = snp_id,
    "clust_dist" = c_dist
  )
  return(snp_clust_df)
}