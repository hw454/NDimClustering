#' Sum the squares of the distances for a group of terms
#'
#' @description
#'   for all the terms in group_col with value num sum the
#'   square of the dist_col values
#'
#' @param num The number to group the terms on
#' @param data The data to take the values from
#' @param dist_col The name of the column to get the values from
#' @param group_col The name of the column to group terms on
#'
#' @return sum of the squares
#'
#' @family cluster_properties
#' @family ic_functions
#'
#' @export
calc_sum_sq_clusts <- function(num, data, dist_col, group_col) {
  sq <- data[data[, group_col] == num, dist_col]**2
  s <- sum(sq)
  return(s)
}
