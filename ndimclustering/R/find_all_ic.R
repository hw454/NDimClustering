#' For each cluster find the aic and bic
#'
#' @description Using [calc_clust_ic]
#'   for each cluster find the aic and bic and add this to the
#'   cluster results dataframe.
#'
#' @param clust_re cluster results dataframe found using [km_nan]
#' @param num_axis the number of axis
#'
#' @details
#'   The aic and bic are found for each cluster using [calc_clust_ic]
#'   this is added to clust_re\$clusters and clust_re\$clusters is output.
#'   Ncents is also added as an axis to clust_re\$clusters.
#'
#' @return clust_re\$clusters
#'
#' @export
find_all_ic <- function(clust_re, num_axis) {
  set.seed(240) # setting seed
  ncents <- length(unique(clust_re$clusters$clust_num))
  ic <- calc_clust_ic(clust_re,
                group_col = "clust_num",
                dist_col = "clust_dist",
                num_axis = num_axis)
  clust_re$clusters <- dplyr::mutate(clust_re, "aic" = ic$aic)
  clust_re$clusters <- dplyr::mutate(clust_re, "bic" = ic$bic)
  clust_re$clusters <- dplyr::mutate(clust_re, "ncents" = ncents)
  return(clust_re$clusters)
}