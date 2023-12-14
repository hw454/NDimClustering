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
  test = 0, num_rows = 5, num_trait0 = 0, num_trait1 = 5
) {
  be_name <- paste0(data_dir, "unstdBeta_df.csv")
  se_name <- paste0(data_dir, "unstdSE_df.csv")
  pv_name <- paste0(data_dir, "pval_df.csv")
  trait_name <- paste0(data_dir, "trait_info_nfil.csv")
  trait_info   <- data.table::fread(trait_name, quote = "")
  beta_df <- as.data.frame(read.csv(be_name)
  )
  se_df   <- as.data.frame(read.csv(se_name)
  )
  pv_df   <- as.data.frame(read.csv(pv_name)
  )
  # write_csv in R will assign X to the row.names. IF the data was
  # created elsewhere then the row.names will have no column label
  # instead. The unlabelled column is the default for loading the row names.
  if ("X" %in% colnames(beta_df)){
    beta_df <- tibble::column_to_rownames(beta_df, var = "X")
    beta_mat <- as.matrix(beta_df)
  } else {
    beta_mat <- as.matrix(beta_df)
  }
  if ("X" %in% colnames(se_df)){
    se_df <- tibble::column_to_rownames(se_df, var = "X")
    se_mat <- as.data.frame(se_df)
  } else {
    se_mat <- as.data.frame(se_df)
  }
  if ("X" %in% colnames(pv_df)){
    pv_df <- tibble::column_to_rownames(pv_df, var = "X")
    pv_mat <- as.matrix(pv_df)
  } else {
    se_mat <- as.data.frame(se_df)
  }
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