#' Test whether the function for checking the NaNs in a column works.
#' @param b_mat data matrix
#' @param nan_col The column name to check for NaNs. default\:"30600_irnt"
#' Should return 0
#' @return 1 if succesful, 0 if not
#' @export
test_all_na <- function(b_mat, nan_col = "30600_irnt") {
  test <- check_col_na(b_mat[, nan_col])
  print("test")
  if (test) {
    return(1)
  } else {
    return(0)
  }
}