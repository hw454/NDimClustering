pca <- function(b_df, pval_df, se_df, 
                np = 3, narm = TRUE) {
  #' Compute the np principal components of b_df.
  #' Inputs: b_df - Data matrix, columns are axis, rows are points.
  #' pval_df - P-values associated with each point, these need to be
  #'    mapped to transformed values also.
  #' np - Integer, number of prinipal components to find,
  #'    given by the np largest Eigen values.
  #' narm - Bool to indicate whether remove NaN values in sample variables.

  # Normalise each Column
  std_b_df <- scale_mat(b_df, narm)
  # Compute Correlation Matrix
  c_mat <- cor(std_b_df, method = "pearson", use = "pairwise.complete.obs")
  # Find Eigen Vectors and Eigen Values
  e_mat <- find_np_eigen_mat(c_mat, np = np)
  # Map Scores onto the space of the components
  # represented by np largest Eigen values.
  b_pc_mat   <- transform_coords(b_df, e_mat)
  pval_pc_mat <- transform_coords(pval_df, e_mat)
  se_pc_mat <- transform_coords(se_df,e_mat)
  out_list   <- list("transform" = e_mat, 
                     "beta" = b_pc_mat, 
                     "pval" = pval_pc_mat,
                     "se" = se_pc_mat)
  return(out_list)
}
transform_coords <- function(p_mat, t_mat) {
  #' Transform the co-ordinates given by the rows
  #' in p_mat to a coordinate system
  #' with basis given by the columns of t_mat
  vec_list <- lapply(row.names(p_mat), trans_vec,
                   p_mat = p_mat, t_mat = t_mat)
  vec_df <- Reduce(rbind, vec_list)
  return(vec_df)
}

trans_vec <- function(p_mat, r, t_mat) {
  #' Transform the vector p_mat[r] to out_vec
  #' on a co-ordinate system whose basis is
  #' the columns of t_mat
  out_vec <- p_mat[r,] %*% t_mat
  rownames(out_vec) <- c(r)
  colnames(out_vec) <- lapply(1:dim(t_mat)[2],pc_name)
  return(out_vec)
}

pc_name <- function(i) {
  return(paste0("P", i))
}

find_np_eigen_mat <- function(mat, np) {
  #' Find the eigen values and vectors for the matrix mat.
  #' Take the np largest eigen values and their corresponding eigen vectors.
  #' Pairs with NaN correlation have no correlating scores and we can therefore 
  #' set their correlation to 0.
  mat[is.na(mat)] <- 0.0
  ev <- eigen(mat)
  # The vectors matrix contains the unit eigen vectors in the columns
  nend <- np %>% min(dim(mat)[1])
  vecs <- ev$vectors[, 1:nend]
  return(vecs)
}

scale_mat <- function(mat, narm = TRUE) {
  #' Rescale the columns of the matrix `mat` so that the values lie 
  #' between -1 and 1 with the with the same distribution.
  #' narm indicates whether NaNs should be removed in the sample calculations.
  # Compute the Sample Mean with na_rm
  xbar <- apply(mat, 2, mean, na.rm = narm)
  # Compute the Sample SE with na_rm
  se <- apply(mat, 2, sd, na.rm = narm)
  out_list <- lapply(colnames(mat), scale_nan_col,
  data = mat, mu = xbar, se = se)
  out_mat <- Reduce(cbind, out_list)
  return(out_mat)
}
scale_nan_col <- function(col, data, mu, se) {
  # Check if SE is 0, if 0 then all points are set to 0.
  # If not 0 then rescale using (x-mu)/se
  return((data[, col] - mu) / se)
}
test_scale <- function(b_df, se_df, p_df) {
  b_scal <- scale_mat(b_df, TRUE)
  print(b_scal)
}
test_pca <- function(b_df, pval_df,
np = 3, narm = TRUE) {
  pca(b_df, pval_df, np, narm)
}
# test_pca(unstd_beta_df[, 1:10], pval_df[, 1:10])
