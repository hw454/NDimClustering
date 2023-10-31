#' Add "np" columns to a dataframe.
#'
#' Create an empty dataframe with column A.
#' Check that the column names after add_nclust_cols(5)
#' are (A, C1,C2,C3,C4,C5)
#'
#' @return
#'
#' @family test
#'
#' @export
test_add_nclust_cols_onecol <- function() {
  df <- data.frame(row.names = character(), A = integer())
  df <- add_nclust_cols(df, 5)
  expect_cols <- c("A", "C1", "C2", "C3", "C4", "C5")
  print("... Testing one col case")
  testit::assert("Add clust has found the wrong column names with one col",
    all(colnames(df) == expect_cols))
  return()
}