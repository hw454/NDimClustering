#' Function for making the label for a cluster column
#'
#' @param i the cluster number
#'
#' @export
#' @family cluster_labels
#' @family clustering_components
make_clust_col_name <- function(i) {
  return(paste0("clust_", i))
}