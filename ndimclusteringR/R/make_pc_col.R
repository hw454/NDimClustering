#' Create a dataframe with column Pi.
#'
#' @details i is the input number
#'
#' @param i - The number for the column label
#' @param rows - The rows for the dataframe,
#'  defalut\: character() for an empty dataframe
#'
#' @return dataframe(Pi = numeric())
#'
#' @export
make_pc_col <- function(i, rows = character()) {
  # Create a dataframe with column Pi
  cname <- paste0("PC", i)
  out <- data.frame(row.names = rows)
  if (length(rows) > 0) {
    out[cname] <- NA
  } else {
    out[cname] <- numeric()
  }
  return(out)
}