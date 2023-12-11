#' Test add a column with value to a list of dataframes
#'
#' @export
test_add_col_df_list <- function() {
  dummy_df0 <- data.frame(row.names = character(), A = integer())
  dummy_df1 <- data.frame(row.names = "row_1", A = 0.0)
  dummy_df_list <- list("df_0" = dummy_df0, "df_1" = dummy_df1)
  col <- "B"
  col_val <- 0.0
  expect_cols <- c("A", col)
  dummy_df_list <- add_col_df_list(dummy_df_list, col, col_val)
  df_0 <- dummy_df_list[["df_0"]]
  df_1 <- dummy_df_list[["df_1"]]
  print("... Testing type list")
  testit::assert("... ... Add col df_list has wrong type",
                 is.list(dummy_df_list))
  print("... Testing items")
  testit::assert("... ...Wrong type for dataframe",
                 is.data.frame(df_0))
  testit::assert("... ...Wrong number of rows for first dataframe",
                 all(nrow(df_0) == nrow(dummy_df0)))
  testit::assert("... ...Wrong number of rows for second dataframe",
                 all(nrow(df_1) == nrow(dummy_df1)))
  testit::assert("... ...Wrong column labels first dataframe",
                 all(colnames(df_0) == expect_cols))
  testit::assert("... ...Wrong column labels second dataframe",
                 all(colnames(df_1) == expect_cols))
}