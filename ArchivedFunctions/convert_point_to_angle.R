#' Convert points to angle values
#'
#' @description Convert the point co-ordinates to values corresponding
#'   to angles to the unit vectors for each axis.
#'
#' @param col The column on the unit matrix to find the angle to
#' @param unit_mat the matrix of unit vectors
#' @param p The point to find the angle for.
#'
#' @details
#'   \deqn{
#'   \theta = \acos{\frac{p \cdot u}{||u|| ||p||}}.
#'   }
#'   If \eqn{\theta>\pi/2} then reassign to the accute angle of the extended
#'   vector.
#'   \deqn{\theta = \pi - \theta}
#'
#' @return \eqn{\theta}
#'
#' @export
convert_point_to_angle <- function(col, unit_mat, p) {
  # Convert the scores on the axis for the point
  # p to the angle to the unit vec at unit_mat[col,]
  unit_vec <- as.matrix(unit_mat[col, ])
  dot_prod <- t(p) %*% unit_vec
  norm_x <- norm(p, type = "F")
  norm_y <- norm(unit_vec, type = "F")
  theta <- acos(dot_prod / (norm_x * norm_y))
  if (theta > pi * 0.5) {
    theta <- pi - theta
  }
  return(as.numeric(theta))
}