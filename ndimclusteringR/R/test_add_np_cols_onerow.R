#' Add "np" columns to a dataframe.
#'
#' @description Create an empty dataframe.
#' Check that the column names after add_npcols(5)
#' are (P1,P2,P3,P4,P5)
#'
#' @family tests
#'
#' @export
test_add_np_cols_onerow <- function() {
  df <- data.frame(row.names = "A")
  df <- add_np_cols(df, 5)
  expect_cols <- c("P1", "P2", "P3", "P4", "P5")
  print("... Testing one row, no col")
  testit::assert("Add np has found the wrong column names with one row",
    colnames(df) == expect_cols)
  testit::assert("Add np has found the wrong row names with one row",
    rownames(df) == c("A"))
}