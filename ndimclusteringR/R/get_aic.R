#' Get the set of clusters which minimizes aic
#'
#' @param clust_out list of cluster dataframes for number of clusters
#'
#' @export
#' @family cluster_properties
#' @family cluster_functions
get_aic <- function(clust_out) {
  cluster_df <- clust_out$clusters
  num_axis <- cluster_df[1, "num_axis"]
  calc_sum_sq <- function(num, data){
    sq <- data[data[, "clust_num"] == num, "clust_dist"]**2
    s <- sum(sq)
    return(s)
  }
  k <- length(unique(cluster_df$clust_num))
  # Dimension of each data point
  d <- num_axis
  # Log-Likelihood term
  ll <- lapply(unique(cluster_df[, "clust_num"]),
    calc_sum_sq,
    data = clust_out
  )
  sig <- as.numeric((stats::var(cluster_df["clust_dist"])))
  l <- Reduce(sum, ll) / sig
  aic_n <- d + 2 * k * l
  cluster_df <- dplyr::mutate(cluster_df, "aic" = aic_n)
  return(cluster_df)
}