#' Test add a column with value to a dataframe
#' @export
test_add_col_df <- function(){
    dummy_df0 <- data.frame(row.names = character(), A = integer())
    col <- "B"
    col_val <- 0.0
    expect_cols <- c("A", col)
    dummy_df0 <- add_col_df(dummy_df0, col, col_val)
    print("... Testing empty case")
    testit::assert("Add col empty case has wrong cols",
      all(colnames(dummy_df0) == expect_cols))
    dummy_df1 <- data.frame(row.names = "row_1", A = 0.0)
    col <- "B"
    col_val <- 0.0
    expect_cols <- c("A", col)
    dummy_df1 <- add_col_df(dummy_df1, col, col_val)
    print("... Testing empty case")
    testit::assert("Add col empty case has wrong cols",
      all(colnames(dummy_df1) == expect_cols))
    testit::assert("Add col one row has wrong number of rows.",
      all(nrow(dummy_df1) == 1))
    return()
}