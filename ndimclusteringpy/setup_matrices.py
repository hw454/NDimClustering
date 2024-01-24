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
#'
#' @family setup_functions
import pandas as pd
def setup_matrices (data_dir,
  test = 0, num_rows = 5, num_trait0 = 0, num_trait1 = 5
):
  be_name = data_dir + "unstdBeta_df.csv"
  se_name = data_dir + "unstdSE_df.csv"
  pv_name = data_dir + "pval_df.csv"
  trait_name = data_dir + "trait_info_nfil.csv"
  # Read the data from the files
  trait_info  = pd.read_csv(trait_name,
    sep = ",",
    index_col = 0
  )
  beta_df = pd.read.csv(be_name)
  se_df   = pd.read.csv(se_name)
  pv_df   = pd.read.csv(pv_name)

  # Find the label for the Outcome trait and the first Exposure trait
  out_row = trait_info.loc[trait_info["pheno_category"] == "Outcome"]
  out_pheno = out_row.phenotype
  exp_row = trait_info.loc[trait_info["pheno_category"] ==  "Exposure"][0]
  exp_pheno = out_row.phenotype

  # Crop data for testing
  if (test):
    mat_list = {"beta" : beta_df,
                "pval" : pv_df,
                "se" : se_df}
    data_matrices = crop_matrices(
      mat_list = mat_list,
      trait_df = trait_info,
      out_pheno = out_pheno,
      exp_pheno = exp_pheno,
      n_rows = num_rows,
      n_col0 = num_trait0,
      n_col1 = num_trait1
    )
  else:
    # Collect the matrices into one object
    data_matrices = {"beta" : beta_df,
                     "pval" : pv_df,
                     "se" : se_df,
                     "trait_info" : trait_info}
  return(data_matrices)