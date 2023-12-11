#' Add "np" columns to a dataframe.
#'
#' @description Create an empty dataframe with column A.
#' Check that the column names after add_nclust_cols(5)
#' are (A, C1,C2,C3,C4,C5)
#'
#' @family dataframe_editors
#'
#' @export
test_add_nclust_cols_onecolrow <- function() {
  df <- data.frame(row.names = "A", A = NA)
  df <- add_nclust_cols(df, 5)
  c_1 <- make_clust_col_name(1)
  c_2 <- make_clust_col_name(2)
  c_3 <- make_clust_col_name(3)
  c_4 <- make_clust_col_name(4)
  c_5 <- make_clust_col_name(5)
  expect_cols <- c("A", c_1, c_2, c_3, c_4, c_5)
  print("... Testing one col one row")
  testit::assert("Add clust has found the wrong column names with one col",
                 all(colnames(df) == expect_cols))
  testit::assert("Add clust has found the wrong row names when empty",
                 rownames(df) == c("A"))
}