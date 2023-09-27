#' Rescale the columns of the matrix "mat"
#'
#' @description
#'   Rescale the columns of the matrix `mat` so that the values lie
#'   between -1 and 1 with the with the same distribution.
#'
#' @param mat The matrix whose columns are to be rescaled
#' @param narm Bool to indicate whether NaNs should be ignored in calculations.
#'   default \: TRUE
#'
#' @details Rescale of each column is found using \link{calc_col_scale} then
#'   column binding the results using \link[SparkR]{cbind} into "out_mat".
#'
#' @return out_mat
#'
#' @export
calc_scale_mat <- function(mat, narm = TRUE) {
  # Compute the Sample Mean with na_rm
  xbar <- apply(mat, 2, mean, na.rm = narm)
  # Compute the Sample SE with na_rm
  se <- apply(mat, 2, sd, na.rm = narm)
  out_list <- lapply(colnames(mat), calc_col_scale,
  data = mat, mu = xbar, se = se)
  out_mat <- Reduce(cbind, out_list)
  return(out_mat)
}