#' Add "np" columns to a dataframe.
#'
#' @description Create an empty dataframe.
#' Check that the column names after add_np_cols(5)
#' are (P1,P2,P3,P4,P5)
#'
#' @family tests
#'
#' @export
test_add_np_cols_empty <- function() {
  df <- data.frame(row.names = character())
  df <- add_np_cols(df, 5)
  expect_cols <- c("P1", "P2", "P3", "P4", "P5")
  print("... Testing empty case")
  testit::assert("Add np has found the wrong column names when empty",
    all(colnames(df) == expect_cols))
}