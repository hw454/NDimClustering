#' Convert spatial co-ordinate to vector of angles
#' @description
#' Convert a point vector to a vector of values corresponding
#' to the angles to the unit axis.
#' @param p The point label to be converted
#' @param mat The matrix containing the co-ordinates for the points.
#' @param unit_mat The matrix of unit vectors in columns.
#' @details
#' Calculate the angles using \link{convert_point_to_angle} for
#' each column in unit_mat.
#' Combine the angles into a vector "ang_vec"
#' @return ang_vec
#' @export
convert_point_to_angles_all <- function(p, mat, unit_mat) {
  vec <- as.matrix(mat[p, ])
  nc <- ncol(mat)
  ang_list <- lapply(1:nc, convert_point_to_angle,
                    unit_mat = unit_mat,
                    p = vec)
  ang_vec <- Reduce(cbind, ang_list)
  return(ang_vec)
}