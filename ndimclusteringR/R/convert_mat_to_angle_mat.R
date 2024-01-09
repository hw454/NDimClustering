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
#' @family preconditioning_functions
#'
#' @export
convert_mat_to_angle_mat <- function(mat) {
  nc <- ncol(mat)
  nang <- nc - 1
  nr <- nrow(mat)
  unit_mat <- diag(nc)
  # Find the angles to all axis except 1 for a point p
  convert_point_to_angles_all <- function(p, mat, unit_mat) {
    vec <- as.matrix(mat[p, ])
    ang_list <- lapply(1:nang, convert_point_to_angle,
      unit_mat = unit_mat,
      p = vec
    )
    ang_vec <- Reduce(cbind, ang_list)
    return(ang_vec)
  }
  # Iterate through each point in the matrix calculating the angles
  p_list <- lapply(1:nr, convert_point_to_angles_all,
                   mat = mat,
                   unit_mat = unit_mat)
  # Convert list of rows to matrix
  ang_mat <- as.matrix(Reduce(rbind, p_list))
  rownames(ang_mat) <- rownames(mat)
  colnames(ang_mat) <- paste0("Ang_to_", colnames(mat)[1:nang])
  return(ang_mat)
}