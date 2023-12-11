#' Calculate the distance between the snps and each cluster centroid
#'
#' @param b_mat the data matrix used for the coordinate of the snp.
#'   Rows are snps, columns are the trait axes.
#' @param centroids_df The dataframe of the centroid co-ordinates for each
#'   cluster. The rows are each cluster number and the columns are the traits.
#' @param norm_typ The type of norm to use in the distance calculation. The
#'   default is the Froebenius norm "F".
#'
#' @description Calculate the distance between the point in "b_mat"
#'   corresponding to "snp_id" and each cluster centroid.
#'   This is stored in a dataframe "snp_clust_df" which has rowname "snp_id",
#'   and column "clust_dist", with the calculated distance.
#'
#' @return snp_clust_df
#'
#' @family distance_functions
#' @family centroid functions
#' @family clustering_components
#'
#' @export
calc_clust_dist <- function(b_mat, centroids_df,
                           norm_typ = "F") {
   clust_nums <- row.names(centroids_df)
   clust_dist_lists <- lapply(clust_nums, calc_clust_dist_col,
        centroids_df = centroids_df,
        data_mat = b_mat)
   clust_dist_df <- Reduce(cbind, clust_dist_lists)
  return(clust_dist_df)
}