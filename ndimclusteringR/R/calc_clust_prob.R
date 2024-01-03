#' Calculate the probability the snp is in the cluster
#'
#' @description
#'   \deqn{p=\frac{1.0}{1.0+d}}
#'
#' @param cluster_df The cluster dataframe.
#'
#' @return p
#'
#' @family probability_functions
#' @family clustering_components
#'
#' @export
calc_clust_prob <- function(cluster_df) {
  d <- cluster_df$clust_dist
  prob <- 1.0 / (1.0 + d)
  return(prob)
}