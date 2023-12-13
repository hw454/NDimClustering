#' read in data from csv files as matrices.
#'
#' @description The directory must contain files with the names\:
#'   * "pval_df.csv"
#'   * "unstdBeta_df.csv"
#'   * "unstdSE_df.csv"
#'   * "trait_info_nfil.csv"
#'  If test parameters are set the data will be cropped to run
#' the functions on a smaller dataset.
#'
#' @param data_dir the directory location of the data
#' @param test Switch 0,1. If 1 then crop the data with the
#'   test variables. default\:0
#' @param num_rows If test==1 then this is the number of rows
#'   of the data that will be used. default\:0
#' @param num_trait0 if test==1 then this is the first column
#'   that will be used. default\:0
#' @param num_trait1 if test==1 then this is the last column
#'   that will be used. default\:0
#'
#' @details
#' rows labeled by SNP_id and columns labelled by trait info.
#' The entries at each position correspond to the values for beta, se, t, and p
#'
#' @export
setup_matrices <- function(data_dir,
  test = 0, num_rows = 0, num_trait0 = 0, num_trait1 = 0
) {
  be_name <- paste0(data_dir, "unstdBeta_df.csv")
  se_name <- paste0(data_dir, "unstdSE_df.csv")
  pv_name <- paste0(data_dir, "pval_df.csv")
  trait_name <- paste0(data_dir, "trait_info_nfil.csv")
  beta_mat <- as.matrix(data.table::fread(be_name), rownames = 1)
  se_mat   <- as.matrix(data.table::fread(se_name), rownames = 1)
  pv_mat   <- as.matrix(data.table::fread(pv_name), rownames = 1)
  trait_info   <- as.data.frame(data.table::fread(trait_name, quote = ""))
  print(trait_info)

  # Find the label for the Outcome trait and the first Exposure trait
  out_row <- which(trait_info$pheno_category == "Outcome")
  out_pheno <- trait_info$phenotype[out_row]
  exp_row <- which(trait_info$pheno_category == "Exposure")[1]
  exp_pheno <- trait_info$phenotype[exp_row]

  # Crop data for testing
  if (test) {
    mat_list <- list("beta" = beta_mat,
                     "pval" = pv_mat,
                     "se" = se_mat)
    data_matrices <- crop_matrices(mat_list = mat_list,
      trait_df = trait_info,
      out_pheno = out_pheno,
      exp_pheno = exp_pheno,
      n_rows = num_rows,
      n_col0 = num_trait0,
      n_col1 = num_trait1
    )
  } else {
    # Collect the matrices into one object
    data_matrices <- list("beta" = beta_mat,
                          "pval" = pv_mat,
                          "se" = se_mat,
                          "trait_info" = trait_info)
  }
  return(data_matrices)
}