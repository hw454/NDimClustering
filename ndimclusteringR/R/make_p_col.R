#' Create a dataframe with column Pi.
#'
#' @details i is the input number
#'
#' @param i - The number for the column label
#'
#' @return dataframe(Pi = numeric())
#'
#' @export
make_p_col <- function(i, rows = character()) {
  # Create a dataframe with column Pi
  cname <- paste0("P", i)
  out <- data.frame(row.names = rows)
  if (length(rows) > 0) {
    out[cname] <- NA
  } else {
    out[cname] <- numeric()
  }
  return(out)
}