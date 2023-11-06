#' Add a column with value to all dataframes in a list
#'
#' @param df_list list of dataframes
#' @param col Label for column name
#' @param col_val value for column
#'
#' @return list of dataframes
#'
#' @export
add_col_df_list <- function(df_list, col, col_val) {
    df_list_out <- lapply(df_list, add_col_df,
                      col = col,
                      col_val = col_val)
    df_list_out <- setNames(df_list_out, names(df_list))
    return(df_list)
}