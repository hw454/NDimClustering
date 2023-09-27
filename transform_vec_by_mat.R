  #' Transform a vector onto a new coordinate system using transform matrix
  #' @description
  #' transform the vector p_mat[r,] to out_vec on a co-ordinate
  #'   system whose basis is the columns of t_mat.
  #' @param p_mat Matrix whose rows are points
  #' @param r The row of matrix corresponding to the point to be
  #'   transformed
  #' @param t_mat The transform matrix whose rows correspond to the
  #' columns of p_mat
  #' @description
  #' \deqn{v=p_mat[r,] * t_mat}
  #' @details
  #' Any NaNs in "p_mat[r,]" are set to 0 before transform. Due to
  #' dimensionality reduction it's not possible to maintain NaN tracking.
  #' @return v
  #' @export
transform_vec_by_mat <- function(p_mat, r, t_mat) {
  q <- p_mat[r, ]
  q[is.na(q)] <- 0.0
  out_vec <- q %*% t_mat
  rownames(out_vec) <- c(r)
  colnames(out_vec) <- lapply(1:dim(t_mat)[2], pc_name)
  return(out_vec)
}