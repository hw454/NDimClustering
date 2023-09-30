#' Calculate the total sum of the cluster distances divided by the variance
#'
#' @param data matrix of data
#' @param group_col The column indicate terms to go together
#'   default\:"clust_num"
#' @param dist_col The name of the column containing the term distances.
#'   default\:"clust_dist"
#'
#' @description
#'   Calculate the sum of the squares for each cluster group using
#'   [calc_sum_sq_clusts]
#'
#' @return The sum of the clusters sum of squares/ divided by the variance
#'   of the distance data.
#'
#' @family cluster_properties
#' @family ic_functions
#'
#' @export
calc_tot_within <- function(data,
                          group_col = "clust_num",
                          dist_col = "clust_dist") {
  ll <- lapply(unique(data[, group_col]),
                    calc_sum_sq_clusts,
                    data = data,
                    group_col = group_col,
                    dist_col = dist_col)
  sig <- as.numeric((stats::var(data[dist_col])))
  l <- Reduce(sum, ll) / sig
  return(l)
}