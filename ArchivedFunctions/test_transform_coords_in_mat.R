#' Test the function `transform_coords_in_mat`
#'
#' @description
#' * Check the types
#' * Check the row labels
#' * Check the column labels
#' * Check the dimensions
#'
#' @return vec_mat
#'
#' @export
test_transform_coords_in_mat <- function() {
  nt <- 20
  nr <- 100
  np <- 3
  p_mat <- matrix(runif(nt * nr, 0, 1), nrow = nr)
  col_nms <- 1:nt
  row_nms <- 1:nr
  colnames(p_mat) <- col_nms
  rownames(p_mat) <- row_nms
  expec_cols <- c("PC1", "PC2", "PC3")
  expec_rows <- row_nms
  t_mat <- matrix(runif(nt * np, 0, 1), nrow = nt)
  pt_mat <- transform_coords_in_mat(p_mat, t_mat)
  testit::assert("... Transform is not a matrix", {
    is.matrix(pt_mat)
  })
  testit::assert("... Wrong number of rows", {
    nrow(pt_mat) == nr
  })
  testit::assert("... Wrong number of columns", {
    ncol(pt_mat) == np
  })
  testit::assert("... Wrong column labels", {
    colnames(pt_mat) == expec_cols
  })
  testit::assert("... Wrong row labels", {
    rownames(pt_mat) == expec_rows
  })

}
