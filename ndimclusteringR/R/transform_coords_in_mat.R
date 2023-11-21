#' Transform the co-ordinates given by the rows in p_mat
#' to a coordinate system with basis given by the columns of t_mat
#'
#' @param p_mat the Matrix whos rows are points to be transformed
#' @param t_mat the transform matrix, rows correspond to the columns
#'   of p_mat. Columns are the number of dimensioned beinng reduced to.
#'
#' @description
#' For each row in p_mat multiply by t_mat. Bind the output vectors into
#' "vec_mat" using "rbind".
#'
#' @return vec_mat
#'
#' @export
transform_coords_in_mat <- function(p_mat, t_mat) {
  vec_list <- lapply(row.names(p_mat), transform_vec_by_mat,
                     p_mat = p_mat, t_mat = t_mat)
  vec_mat <- Reduce(rbind, vec_list)
  colnames(vec_mat) <- lapply(seq_len(ncol(vec_mat)),
                              make_pc_col_name)
  rownames(vec_mat) <- rownames(p_mat)
  return(vec_mat)
}
