#' Crop all of the matrices in mat_list
#'
#' @description include only the first n_rows
#'   and columns between n_col0 and n_col1,
#'   ensure that out_pheno and exp_pheno are included.
#'
#' @details
#'   1. Get the columns of each matrix between n_col0 and n_col1.
#'     Crop to n_rows at this step at the same time
#'   2. Check if the exposure and outcome columns are in the n_col range.
#'   3. If not then:
#'      a. Get the matrix with just the exposure and outcome columns
#'      using `crop_mat_colnames`. Crop to n_rows at this step at the same time
#'   4. Merge the two sets to create a matrix with columns in the n_col
#'   range and the exposure and outcome columns.
#'
#' @param mat_list The list of matrices that need to be cropped.
#' @param trait_df The dataframe of the traits.
#' @param out_pheno The outcome phenotype label
#' @param exp_pheno The ecposure phenotype label.
#' @param n_rows The number of rows to crop by.
#' @param n_col0 The first column to include
#' @param n_col1 The last column to include.
#'
#' @return List of cropped matrices.
#'
#' @family preconditioning_functions
#'
#' @export
crop_matrices <- function(mat_list, trait_df,
  out_pheno, exp_pheno, n_rows, n_col0, n_col1
) {
  mat_names <- names(mat_list)
  mat_nrows <- nrow(mat_list[[mat_names[1]]])
  if (n_rows > mat_nrows) {
    n_rows <- mat_nrows
  }
  mat_crop <- lapply(mat_list, crop_mat_colnums,
    num_rows = n_rows,
    col0 = n_col0,
    col1 = n_col1
  )
  names(mat_crop) <- names(mat_list)
  out_check <-  out_pheno %in% trait_df$phenotype[n_col0:n_col1]
  exp_check <- exp_pheno %in% trait_df$phenotype[n_col0:n_col1]
  if (!out_check || !exp_check) {
    mat_out_exp <- lapply(mat_list, crop_mat_colnames,
      num_rows = n_rows,
      col_names = c(out_pheno, exp_pheno)
    )
    names(mat_out_exp) <- names(mat_list)
    # Combine matrices and remove duplicates
    filt_mat_list <- lapply(names(mat_list),
      join_mat_on_name,
      mat_list_a = mat_out_exp,
      mat_list_b = mat_crop
    )
    names(filt_mat_list) <- names(mat_list)

    trait_list <- colnames(mat_crop[[mat_names[1]]])
    trait_df$phenotype <- trait_list
    # Collect the matrices into one object
    data_matrices <- append(filt_mat_list, list("trait_info" = trait_df))
  } else {
    trait_list <- colnames(mat_crop[[mat_names[1]]])
    data_matrices <- append(mat_crop, list("trait_info" = trait_list))
  }
  return(data_matrices)
}
