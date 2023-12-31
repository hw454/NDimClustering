#' Scale the values in the column of a matrix around the mean
#'
#' @description Check if SE is 0, if 0 then all points are set to 0.
#'   If not 0 then rescale using (x-mu)/se
#'
#' @param col Column to be scaled
#' @param data Matrix of data
#' @param mu Mean of the column
#' @param se Standard error of the column
#'
#' @details
#'   if (se==0){data_out = 0}
#'   else{ data_out = (data[,col]-mu)/se}
#'
#' @return data_out
#'
#' @family pca_functions
#'
#' @export
calc_col_scale <- function(col, data, mu, se) {
  if (is.na(se[col])) {
    # If SE is NULL
   return(data[, col] - mu[col])

  } else if (se[col] == 0) {
    # If SE is 0 or NULL
   return(data[, col] - mu[col])
  } else {
    # SE is defined
    return((data[, col] - mu[col]) / se[col])
  }
}