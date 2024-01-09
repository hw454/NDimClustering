#' Calculate the distance between the snps and their allocated cluster centroid
#'
#' @param snp_id the label for the snp
#' @param data_mat the data matrix used for the coordinate of the snp.
#'   Rows are snps, columns are the trait axes.
#' @param cent The centre for the cluster
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
#' @family centroid functions
#' @family clustering_components
#'
#' @export
calc_snp_cent_dist <- function(snp_id, data_mat, cent,
                           norm_typ = "F") {
   snp_score <- data_mat[snp_id, ]
   c_dist <- norm(data.matrix(stats::na.omit(cent - snp_score)), norm_typ)
   snp_clust_df <- data.frame(row.names = snp_id,
                       "clust_dist" = c_dist)
  return(snp_clust_df)
}