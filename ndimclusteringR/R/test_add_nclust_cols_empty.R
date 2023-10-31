#' Add "np" columns to a dataframe.
#'
#' Create an empty dataframe.
#' Check that the column names after add_nclust_cols(5)
#' are (C1,C2,C3,C4,C5)
#'
#' @return
#'
#' @family tests
#'
#' @export
test_add_nclust_cols_empty <- function() {
  df <- data.frame(row.names = character())
  df <- add_nclust_cols(df, 5)
  expect_cols <- c("C1", "C2", "C3", "C4", "C5")
  print("... Testing empty case")
  testit::assert("Add clust has found the wrong column names when empty",
    all(colnames(df) == expect_cols))
  return()
}