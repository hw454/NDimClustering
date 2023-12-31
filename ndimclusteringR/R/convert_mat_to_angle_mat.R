#' Convert a matrix of spatial data into a matrix of angles
#'
#' @param mat the matrix of spatial data, each row corresponds to a point.
#'
#' @description The angles correspond to the angles to the unit vector
#'   on each axis.
#'   The angles are calculated for each point (row) in mat
#'   using [convert_point_to_angles_all] which uses
#'   [convert_point_to_angle] for each column.
#'   The angles for each point are then row bound using "rbind".
#'   into "ang_mat".
#'
#' @return ang_mat
#'
#' @family clustering_components
#'
#' @export
convert_mat_to_angle_mat <- function(mat) {
  nc <- ncol(mat)
  nr <- nrow(mat)
  unit_mat <- diag(nc)
  p_list <- lapply(1:nr, convert_point_to_angles_all,
                   mat = mat,
                   unit_mat = unit_mat)
  ang_mat <- Reduce(rbind, p_list)
  ang_mat <- as.matrix(ang_mat)
  rownames(ang_mat) <- rownames(mat)
  colnames(ang_mat) <- colnames(mat)
  return(ang_mat)
}