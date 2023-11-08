#' Add "np" columns to a dataframe.
#'
#' @description Create an empty dataframe.
#' Check that the column names after add_nclust_cols(5)
#' are (C1,C2,C3,C4,C5)
#'
#' @family tests
#'
#' @export
test_make_clust_col <- function() {
  df <- data.frame(row.names = character())
  ex_rows <- c("A", "B", "C")
  c_num <- as.integer(5)
  df <- make_clust_col(c_num, rows = ex_rows)
  df2 <- make_clust_col(c_num)
  expect_cols <- c(make_clust_col_name(c_num))
  print("... test case with rows")
  testit::assert("make_clust_col has found the wrong column names",
    all(colnames(df) == expect_cols))
  testit::assert("make_clust_col has found the wrong row names",
    all(rownames(df) == ex_rows))
  print("... test empty case")
  testit::assert("make_clust_col has found the wrong row names when empty",
    all(length(rownames(df2)) == 0))
}