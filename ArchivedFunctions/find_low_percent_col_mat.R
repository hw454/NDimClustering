#' Find the rownames corresponding to the lowest `col_percent` of each column
#'
#' @param data_mat the data matrix
#' @param col_percent the threshold percentage
#' @param na_rm whether to remove NaNs or not
#'
#' @export
find_low_percent_col_mat <- function(data_mat, col_percent, na_rm = TRUE) {
   sub_rows_list <- lapply(colnames(data_mat), find_low_percent_col,
                           data_mat = data_mat,
                           col_percent = col_percent,
                           na_rm = na_rm)
   sub_rows <- Reduce(rbind, sub_rows_list)
   sub_rows <- as.vector(sub_rows)
   sub_rows <- unique(sub_rows)
   return(sub_rows)
}
