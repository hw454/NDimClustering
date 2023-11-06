#' Add a column with value to dataframes
#'
#' @param df Dataframes
#' @param col Label for column name
#' @param col_val value for column
#'
#' @return list of dataframes
#'
#' @export
add_col_df <- function(df, col, col_val = NULL) {
    if (nrow(df) <= 0){
        df[col] <- numeric()
    } else {
        df[col] <- col_val
    }
    return(df)
}