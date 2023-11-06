#' Add "np" columns to a dataframe.
#'
#' @description Create an empty dataframe.
#' Check that the column names after add_nclust_cols(5)
#' are (C1,C2,C3,C4,C5)
#'
#' @family tests
#'
#' @export
test_add_nclust_cols_onerow <- function() {
  df <- data.frame(row.names = "A")
  df <- add_nclust_cols(df, 5)
  expect_cols <- c("C1", "C2", "C3", "C4", "C5")
  print("... Testing one row, no col")
  testit::assert("Add clust has found the wrong column names with one row",
    colnames(df) == expect_cols)
  testit::assert("Add clust has found the wrong row names with one row",
    rownames(df) == c("A"))
}