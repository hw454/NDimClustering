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
#' @family ic_functions
#'
#' @export
find_all_ic <- function(clust_re, num_axis) {
  set.seed(240) # setting seed
  # The row names will get overridden in rbind so assign rownames to column
  # then reassign when ncents has been chosen.
  clust_re$clusters <- tibble::rownames_to_column(clust_re$clusters, "snp_id")
  ic <- calc_clust_ic(clust_re$clusters,
                group_col = "clust_num",
                dist_col = "clust_dist",
                num_axis = num_axis)
  clust_re$clusters <- dplyr::mutate(clust_re$clusters, "aic" = ic$aic)
  clust_re$clusters <- dplyr::mutate(clust_re$clusters, "bic" = ic$bic)
  return(clust_re$clusters)
}