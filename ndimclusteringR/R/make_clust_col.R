#' Function for making the columns for clusters
#'
#' @param i the cluster number
#'
#' @export
#'
#' @family cluster_labels
#' @family cluster_properties
make_clust_col <- function(i, rows = character(), col_val = NA) {
  # Create a dataframe with column Pi
  cname <- make_clust_col_name(i)
  out <- data.frame(row.names = rows)
  if (length(rows) > 0) {
    out[cname] <- col_val
  } else {
    out[cname] <- numeric()
  }
  return(out)
}