#' Add "np" columns to a dataframe.
#'
#' @description Create an empty dataframe.
#' Check that the column names after add_nclust_cols(5)
#' are (C1,C2,C3,C4,C5)
#'
#' @family tests
#'
#' @export
test_add_nclust_cols_empty <- function() {
  df <- data.frame(row.names = character())
  df <- add_nclust_cols(df, 5)
  c_1 <- make_clust_col_name(1)
  c_2 <- make_clust_col_name(2)
  c_3 <- make_clust_col_name(3)
  c_4 <- make_clust_col_name(4)
  c_5 <- make_clust_col_name(5)
  expect_cols <- c(c_1, c_2, c_3, c_4, c_5)
  print("... Testing empty case")
  testit::assert("Add clust has found the wrong column names when empty",
    all(colnames(df) == expect_cols))
}