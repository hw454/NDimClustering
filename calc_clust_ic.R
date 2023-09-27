#' Find the BIC and AIC for the clusters.
#' @param group_col The name of the column to group terms on
#' @param dist_col The naame of the column to use distances from.
#' @param num_axis The number of axis i.e the dimensionality
#' @description
#' Calculate the sum of squares using \link{calc_tot_withins}
#' and \link{calc_sum_sq_clusts}, this is "l".
#' "k" is the number of cluster groups.
#' "d" is the dimensionality = num_axis
#' "n" is the number of points
#' \deqn{
#' aic = d + 2*k*l}
#' \deqn{
#' bic = d+log(n)*k*l}
#' out = list("aic" = aic, "bic" = bic)
#' @return out
#' @export
calc_clust_ic <- function(clust, group_col, dist_c, num_axis) {
  # Number of estimated parameters
  k <- length(unique(clust$clust_num))
  # Number of data points
  n <- nrow(clust)
  # Dimension of each data point
  d <- num_axis
  # Log-Likelihood term
  l <- calc_tot_withins(clust, group_col, dist_col)
  aic_n <- d + 2 * k * l
  bic_n <- d + log(n) * k * l
  ic_out <- list("aic" = aic_n, "bic" = bic_n)
  return(ic_out)
}