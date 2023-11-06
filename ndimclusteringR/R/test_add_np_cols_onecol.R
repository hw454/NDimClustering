#' Add "np" columns to a dataframe.
#'
#' @description Create an empty dataframe with column A.
#' Check that the column names after add_np_cols(5)
#' are (A, P1,P2,P3,P4,P5)
#'
#' @family test
#'
#' @export
test_add_np_cols_onecol <- function() {
  df <- data.frame(row.names = character(), A = integer())
  df <- add_np_cols(df, 5)
  expect_cols <- c("A", "P1", "P2", "P3", "P4", "P5")
  print("... Testing one col case")
  testit::assert("Add np has found the wrong column names with one col",
    all(colnames(df) == expect_cols))
}