#' Filter the dataframe and get the columns names that aren't in "ignore_cols"
#'
#' @description
#'  * Filter the dataframe "c_scores" by the column "filter_col"
#'   with value "filter_val".
#' * Once filtered removed any all NaN columns.
#' * Extract the column names which aren't in "ignore_cols" to col_list
#'
#' @param c_scores dataframe of interest
#' @param filter_col string label for column to filter by
#' @param filter_val value to filter "filter_col" with
#' @param ignore_cols list of column names to not include in output.
#'
#' @return col_list
#'
#' @export
crop_col_names <- function(c_scores, filter_col, filter_val, ignore_cols) {
    c_term <- c_scores[c_scores[filter_col] == filter_val, ]
    c_term <- c_term[, colSums(!is.na(c_term)) != 0]
    col_list <- colnames(c_term[!(colnames(c_term) %in% ignore_cols)])
    return(col_list)

}