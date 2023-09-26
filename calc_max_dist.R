#' Calculate the maximum distance between any two points.
#' @param score_mat data matrix
#' @param norm_typ type of data to use for distance
#' @param na_rm bool switch on whether to remove NaNs in min and max.
#' @description
#' Points are given by rows of the dataframe.
#' Find the maximum distance between any two points in "score_mat".
#' Distance is given by the norm "norm_typ"
#' @details
#' To get an upper and lower bound for the points find the max in each column
#' and the min in each column. The distance between these is an upper bound
#' for the maximum distance between the points in the matrix
#' @return max_dist
#' @export
max_dist_calc <- function(score_mat, norm_typ = "F", na_rm = TRUE) {
  # Find the max distance based on range on each axis.
  max_p <- apply(score_mat, 2, max, na.rm = na_rm)
  min_p <- apply(score_mat, 2, min, na.rm = na_rm)
  max_dist <- norm(as.matrix(max_p - min_p), norm_typ)
return(max_dist)
}