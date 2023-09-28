  #' Compute the np principal components of b_mat.
  #'
  #' @param b_mat - Data matrix, columns are axis, rows are points.
  #' @param pval_mat - P-values associated with each point.
  #' @param se_mat - Standard errors of the data points.
  #' @param np - Integer, number of prinipal components, default \:3
  #' @param narm - Bool to indicate whether remove NaN values in
  #'   sample variables. default\: TRUE
  #'
  #' @description
  #'   * Find the np largest Eigen values and eigen vectors.
  #'   * Combine the vectors with column binding to form the transform
  #'   matrix.
  #'   * Transform the matrices into the princpal component space.
  #'   * Combine the matrices into "out_list"
  #'   \code{
  #'   out_list   <- list("transform" = e_mat,
  #'                  "beta" = b_pc_mat,
  #'                   "pval" = pval_pc_mat,
  #'                   "se" = se_pc_mat)
  #'                   }
  #'
  #' @return out_list
  #'
  #' @export
find_principal_components <- function(b_mat, pval_mat, se_mat,
                np = 3, narm = TRUE) {
  # Normalise each Column
  std_b_mat <- scale_mat(b_mat, narm)
  # Compute Correlation Matrix
  c_mat <- stats::cor(std_b_mat, method = "pearson", use = "pairwise.complete.obs")
  # Find Eigen Vectors and Eigen Values
  e_mat <- find_np_eigen_mat(c_mat, np = np)
  # Ensure the transform rownames correspond to the axes in the original data.
  row.names(e_mat) <- colnames(b_mat)
  # Map Scores onto the space of the components
  # represented by np largest Eigen values.
  b_pc_mat   <- transform_coords_in_mat(b_mat, e_mat)
  pval_pc_mat <- transform_coords_in_mat(pval_mat, e_mat)
  se_pc_mat <- transform_coords_in_mat(se_mat, e_mat)
  out_list   <- list("transform" = e_mat,
                     "beta" = b_pc_mat,
                     "pval" = pval_pc_mat,
                     "se" = se_pc_mat)
  return(out_list)
}
