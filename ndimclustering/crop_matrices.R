#' Crop all of the matrices in mat_list
#'
#' @description include only the first n_rows
#'   and columns between n_col0 and n_col1,
#'   ensure that out_pheno and exp_pheno are included.
#'
#' @param mat_list The list of matrices that need to be cropped.
#' @param trait_df The dataframe of the traits.
#' @param out_pheon The outcome phenotype label
#' @param exp_pheno The ecposure phenotype label.
#' @param n_rows The number of rows to crop by.
#' @param n_col0 The first column to include
#' @param n_col1 The last column to include.
#'
#' @return List of cropped matrices.
#'
#' @export
crop_mat_list <- function(mat_list, trait_df,
                        out_pheno, exp_pheno,
                        n_rows, n_col0, n_col1) {
  mat_out <- lapply(mat_list, crop_mat_colnames,
                    num_rows = n_rows,
                    col_names = c(out_pheno, exp_pheno))
  names(mat_out) <- names(mat_list)
  trait_out <- trait_info[trait_info$phenotype %in% c(out_pheno, exp_pheno)]

  mat_crop <- lapply(mat_list, crop_mat_colnums,
                          num_rows = n_rows,
                          col0 = n_col0,
                          col1 = n_col1)
  names(mat_crop) <- names(mat_list)

  filt_mat_list <- lapply(names(mat_list),
                          join_out_mat_on_name,
                          mat_out_list = mat_out,
                          mat_crop_list = mat_crop
  )
  names(filt_mat_list) <- names(mat_list)

 trait_info   <- rbind(trait_info[n_col0:n_col1], trait_out)

  # Collect the matrices into one object
  data_matrices <- append(filt_mat_list, list("trait_info" = trait_info))
  return(data_matrices)
}
