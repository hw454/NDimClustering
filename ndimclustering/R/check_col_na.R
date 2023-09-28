#' Check the percentage of NaNs in column
#'
#' @param percent - The percentage of the column which must be not NaN
#' @param b_col - The column vector from the matrix
#'
#' @return 0 or 1 - 1 if `b_col` has an acceptable number of non_NaN, else 0.
#'
#' @export
check_col_na <- function(b_col, percent = 0.95) {
  n_accept <- length(b_col) * (1 - percent)
  narows <- which(is.na(b_col))
  if (length(narows) > n_accept) {
    # All rows are NaN so trait will be removed from trait
    return(1)
  } else {
    return(0)
    }
}