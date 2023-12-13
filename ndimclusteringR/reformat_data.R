#' Reformat the data into the structure desired for clustering.
#'
#' @param data_matrices The list of matrices to be transformed.
#'   The data governing any transformation is in the "beta" matrix
#' @param bin_angles binary variable indicating whether to convert to
#'  angle data or not. 1 is angles, 0 is not.
#' @param pca_type Method for principal components to use. Options are
#'  "manual", "prcomp" or "none"
#'
#' @description
#'   1. If the angles switch is set convert the data to the angles.
#'     The is only done for the main beta data matrix.
#'   2. Find the principal components for the beta matrix. Use the method set
#'     in pca_type, if none then skip and return the data_matrices.
#'
#' @return data_matrices List of matrices
#'
#' @export
#' @family transform_data

reformat_data <- function(data_matrices,
  bin_angles = 1, pca_type = "none", np = 1
) {
  # If the angles switch is set convert the data to the angles.
  # The is only done for the main beta data matrix.
  if (bin_angles) {
    data_matrices$beta <- convert_mat_to_angle_mat(data_matrices$beta)
  }
  # If pca_type is "none" then return data
  # If pca_type is "prcomp" then find the pca vectors using prcomp and transform
  # all the data.
  # If pca_type is "manual" then find the pca vectors using the manual pca
  # algorithm. This method is not currently supported.
  if (pca_type == "none") {
    nc <- ncol(data_matrices$beta)
    pca_list <- list("beta_pc" = data_matrices$beta,
                     "pval_pc" = data_matrices$pval,
                     "se_pc" = data_matrices$se,
                     "transform" = diag(nc))
    data_matrices <- append(data_matrices, pca_list)
  } else if (pca_type == "prcomp") {
    pca_beta <- stats::prcomp(data_matrices$beta,
      center = TRUE,
      scale = TRUE,
      rank = np
    )
    t_mat <- pca_beta$rotation
    b_pc_mat <- pca_beta$x
    p_pc_mat <- transform_coords_in_mat(
      data_matrices$pval, t_mat
    )
    se_pc_mat <- transform_coords_in_mat(
      data_matrices$se, t_mat
    )
  }
  # Store the matrices for result output
  pca_list <- list("beta_pc" = b_pc_mat,
                   "pval_pc" = p_pc_mat,
                   "se_pc" = se_pc_mat,
                   "transform" = t_mat)
  data_matrices <- append(data_matrices, pca_list)
  return(data_matrices)
}