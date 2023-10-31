#' Test the `check_col_na` function
#'
#' @return
#'
#' @export
#'
#' @family tests
test_na_check <- function() {
  # Test whether the function for checking the NaNs in a column works.
  df <- data.frame(row.names = 1:20)
  df["A"] <- 1
  df["B"] <- -1
  df["C"] <- NA
  nan_col <- "C"
  not_nan_col <- "B"
  print("... Testing the NaN columns")
  test1 <- check_col_na(df[, nan_col])
  testit::assert("Check_col_na didn't fail NaN col", {
    test1 == 1
  })
  print("... Testing a non-NaN column")
  test2 <- check_col_na(df[, not_nan_col])
  testit::assert("Check_col_na fail none NaN col", {
    test2 == 0
  })
}