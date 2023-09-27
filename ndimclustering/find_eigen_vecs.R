  #' Find the eigen values and vectors for the matrix mat.
  #'
  #' @description
  #'   Take the np largest eigen values and their corresponding eigen vectors.
  #'   Pairs with NaN correlation have no correlating scores and we can
  #'   therefore set their correlation to 0.
  #'
  #' @param mat the matrix to find the eigen values and vectors for.
  #' @param np the number of eigen vectors to keep.
  #'
  #' @details
  #'   "em_mat" is the matrix formed by column binding the eigen 
  #'   vectors with "np" largest eigen values.
  #'
  #' @return e_mat
find_np_eigen_mat <- function(mat, np) {
  mat[is.na(mat)] <- 0.0
  ev <- eigen(mat)
  # The vectors matrix contains the unit eigen vectors in the columns
  nend <- min(np, dim(mat)[1])
  vecs <- ev$vectors[, 1:nend]
  return(vecs)
}