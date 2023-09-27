#' Create a dataframe with column Pi.
#'
#' @details i is the input number
#'
#' @param i - The number for the column label
#'
#' @return dataframe(Pi = numeric())
#'
#' @export
make_p_col <- function(i) {
  #' Create a dataframe with column Pi
  cname <- paste0("P", i)
  out <- data.frame(col = numeric())
  colnames(out) <- c(cname)
  return(out)
}