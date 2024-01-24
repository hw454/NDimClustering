#' Check convergence for the cluster centroid
#'
#' @param clustnum_df The cluster membership dataframe. The rows are the snps,
#'  The columns are "clust_num".
#' @param b_mat The data for each snp. Each row is the snps coordinate.
#' @param p_mat Matrix of p-values corresponding to b_mat
#' @param bin_p_clust Bool switch 1 if pvalues should be used to weight averages
#' of snps
#'
#' @details
#'   If the cluster is not empty then the new centres are found using
#'     the means of the snps in the cluster.
#'   If the bin_p_clust then weight the average by the p-values.
#'
#' @return out_centroid_df
#'
#' @family clustering_components
#' @family centroid_functions
#' @family dbscan_functions
#'
#' @export
calc_clust_cent <- function(clustnum_df, b_mat, p_mat, bin_p_clust = TRUE
) {
  clust_average <- function(c_num, clustnum_df,
    b_mat, p_mat, bin_p_clust
  ) {
    # Function for calculating the centroid for a given cluster.
    clust_snps <- which(clustnum_df$clust_num == c_num)
    n_terms <- length(clust_snps)
    if (n_terms == 0) {
      cent_df <- data.frame(col.names = colnames(b_mat))
    } else if (n_terms == 1) {
      cent_df <- as.data.frame(t(b_mat[clust_snps, ]))
      rownames(cent_df) <- c_num
    } else {
      b_clust <- b_mat[clust_snps, ]
      if (bin_p_clust) {
        p_clust <- - log10(p_mat[clust_snps, ])
        tot <- colSums(p_clust)
        scores <- b_clust * -log10(p_clust)
        cent <- colSums(scores) / tot
      } else {
        cent <- colMeans(b_clust)
      }
      cent_df <- as.data.frame(t(cent))
      rownames(cent_df) <- c_num
    }
    return(cent_df)
  }
  clust_nums <- unique(clustnum_df$clust_num)
  cent_list <- lapply(clust_nums,
    clust_average,
    clustnum_df = clustnum_df,
    b_mat = b_mat,
    p_mat = p_mat,
    bin_p_clust = bin_p_clust
  )
  centroids_df <- Reduce(rbind, cent_list)
  return(centroids_df)
}