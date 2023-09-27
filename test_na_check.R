#' Check the `na_col_check` function
#' @param b_mat - Data matrix with trait labelled columns
#' @param nan_col - Column label for an invalid column
#' @return 0,1 - 1 if work, 0 if not
#' @export
test_na_check <- function(b_df, nan_col = "30600_irnt") {
  #' Test whether the function for checking the NaNs in a column works.
  test <- na_col_check(b_df[, nan_col])
  print("test")
  if (test) {
    return(1)
  } else {
    return(0)
  }
}