#' Find rownames corresponding to the lowest `col_percent` of the column `col`
#'
#' @param data_mat the data matrix
#' @param col the column label
#' @param col_percent the threshold percentage
#' @param na_rm whether to remove NaNs or not
#'
#' @export
find_low_percent_col <- function(data_mat, col,
                                col_percent = 0.1,
                                na_rm = TRUE) {
    ten_p <- quantile(col, prob = col_percent, na.rm = na_rm)
    data_sub <- data_mat[, data_mat[col] < ten_p]
    sub_rows <- rownames(data_sub)
    return(sub_rows)
}