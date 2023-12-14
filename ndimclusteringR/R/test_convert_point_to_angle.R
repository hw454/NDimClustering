#' Test the `crop_matrices` function
#'
#' @export
#'
#' @family tests
test_convert_point_to_angle <- function() {
  # Test whether the function for checking the NaNs in a column works.
  dummy_angles <- c(0, pi / 2, pi / 3, pi / 4, pi)
  dummy_traits <- c("T1", "T2")
  nc <- length(dummy_traits)
  dummy_points_t1 <- Reduce(rbind, lapply(dummy_angles, cos))
  dummy_points_t2 <- Reduce(rbind, lapply(dummy_angles, sin))
  dummy_points <- as.matrix(cbind(dummy_points_t1, dummy_points_t2))
  unit_mat <- diag(nc)
  expec_nang <- ncol(dummy_points) - 1
  ang <- convert_point_to_angle(col = 1,
    unit_mat = unit_mat,
    p = as.matrix(dummy_points[1, ])
  )
  expec_ang <- 0
  testit::assert("Angle 0 not found", {
    ang == expec_ang
  })
  ang <- convert_point_to_angle(col = 1,
    unit_mat = unit_mat,
    p = as.matrix(dummy_points[2, ])
  )
  expec_ang <- pi / 2
  err <- abs(ang - expec_ang)
  testit::assert("Angle pi/2 not found", {
    err < 1e-8
  })
  ang <- convert_point_to_angle(col = 1,
    unit_mat = unit_mat,
    p = as.matrix(dummy_points[3, ])
  )
  expec_ang <- pi / 3
  err <- abs(ang - expec_ang)
  testit::assert("Angle pi/3 not found", {
    err < 1e-8
  })
  ang <- convert_point_to_angle(col = 1,
    unit_mat = unit_mat,
    p = as.matrix(dummy_points[4, ])
  )
  expec_ang <- pi / 4
  err <- abs(ang - expec_ang)
  testit::assert("Angle pi / 4 not found", {
    err < 1e-8
  })
  ang <- convert_point_to_angle(col = 1,
    unit_mat = unit_mat,
    p = as.matrix(dummy_points[5, ])
  )
  expec_ang <- pi
  err <- abs(ang - expec_ang)
  testit::assert("Angle pi not found", {
    err < 1e-8
  })
  print("... Test the angle conversions on the full matrix")
  ang_mat <- convert_mat_to_angle_mat(
    dummy_points
  )
  err <- ang_mat - dummy_angles
  testit::assert("Angle 0 not found", {
    err[1, ] < 1e-8
  })
  testit::assert("Angle pi / 2 not found", {
    err[2, ] < 1e-8
  })
  testit::assert("Angle pi/ 3 not found", {
    err[3, ] < 1e-8
  })
  testit::assert("Angle pi/4 not found", {
    err[4, ] < 1e-8
  })
  testit::assert("Angle pi not found", {
    err[5, ] < 1e-8
  })
  testit::assert("Ang_mat has wrong number of columns", {
    ncol(ang_mat) == expec_nang
  })
  testit::assert("Ang_mat has wrong number of rows", {
    nrow(ang_mat) == nrow(dummy_points)
  })
}