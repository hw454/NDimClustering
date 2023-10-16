#' Create a dataframe with column clust_i.
#'
#' @details i is the input number
#'
#' @param i - The number for the column label
#'
#' @return dataframe(Pi = numeric())
#'
#' @export
make_p_col <- function(i) {
  # Create a dataframe with column Pi
  cname <- paste0("clust_", i)
  out <- data.frame(col = numeric())
  colnames(out) <- c(cname)
  return(out)
}