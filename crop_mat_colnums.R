#' Return the matrix with rows 1: num_rows and columns col0:col1
#' @param mat the matrix to be cropped
#' @param num_rows the number of rows to keep.
#' @param col0 the first col to keep
#' @param col1 the last col to keep
#' @return cropped matrix
#' @export
crop_mat_colnums <- function(mat, num_rows, col0, col1) {
  mat_out <- mat[1:num_rows, col0:col1]
  return(mat_out)
}