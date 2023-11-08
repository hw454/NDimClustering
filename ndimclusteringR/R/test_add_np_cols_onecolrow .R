#' Add "np" columns to a dataframe.
#'
#' @description Create an empty dataframe with column A.
#' Check that the column names after add_nclust_cols(5)
#' are (A, P1,P2,P3,P4,P5)
#'
#' @family dataframe_editors
#'
#' @export
test_add_np_cols_onecolrow <- function() {
  df <- data.frame(row.names = "A", A = NA)
  print(colnames(df))
  df <- add_np_cols(df, 5)
  pc_1 <- make_pc_col_name(1)
  pc_2 <- make_pc_col_name(2)
  pc_3 <- make_pc_col_name(3)
  pc_4 <- make_pc_col_name(4)
  pc_5 <- make_pc_col_name(5)
  expect_cols <- c("A", pc_1, pc_2, pc_3, pc_4, pc_5)
  print(colnames(df))
  print(expect_cols)
  print("... Testing one col one row")
  testit::assert("Add np has found the wrong column names with one col",
    all(colnames(df) == expect_cols))
  testit::assert("Add np has found the wrong row names when empty",
    rownames(df) == c("A"))
}