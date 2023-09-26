#' Crop matrix using column names
#' @param mat matrix to be cropped
#' @param num_rows the number of rows to keep
#' @param col_names the column names to include.
#' @description Return the matrix with row 1:num_rows
#' and column names in col_names
#' @export
crop_mat_colnames <- function(mat, num_rows, col_names) {
  mat_out <- mat[1:num_rows,
                 which(colnames(mat) %in% col_names)]
  return(mat_out)
}