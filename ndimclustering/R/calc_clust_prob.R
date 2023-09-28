#' Calculate the probability the snp is in the cluster
#'
#' @description
#'   \deqn{p=\frac{1.0}{1.0+d}}
#'
#' @param d The snp distance to the cluster centre
#'
#' @return p
#'
#' @export
calc_clust_prob <- function(d) {
  dist <- 1.0 / (1.0 + d)
  return(dist)
}