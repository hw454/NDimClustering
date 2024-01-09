#' Create a dataframe with column clust_i.
#'
#' @details i is the input number
#'
#' @param i - The number for the column label
#' @param rows - The rows for the dataframe,
#'  defalut\: character() for an empty dataframe
#' @param col_val Value for column, default\: NA
#'
#' @return dataframe(Pi = numeric())
#'
#' @export
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