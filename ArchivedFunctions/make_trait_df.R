  #' Create a dataframe for traits
  #'
  #' @description
  #' The dataframe "trait_df" has the columns "label" and "a_ind" the index
  #' position for the trait. Only include traits whose column
  #' in data_mat meets the na_percent threshold.
  #'
  #' @param pheno_list - list of all traits
  #' @param data_mat - matrix of the data.
  #' @param na_percent - the percent of a data column that is
  #' required to be not NaN.
  #'
  #' @details
  #' Check traits using the function [check_trait].
  #' trait_df is the "rbind" of all the valid traits.
  #' trait_df has columns:
  #' * "label" - The labels for the traits.
  #' Also the column names of the data matrix.
  #' * "a_ind" - The position in the trait list.
  #'
  #' @return trait_df
  #'
  #' @export
make_trait_df <- function(pheno_list, data_mat, na_percent) {
  trait_df_list <- lapply(pheno_list, check_trait,
                          pheno_list = pheno_list,
                          data_mat = data_mat,
                          na_percent = na_percent)
  trait_df <- Reduce(rbind, trait_df_list) # Combine all traits into one df.
  trait_df <- dplyr::distinct(trait_df)  # Remove duplicate traits.
  return(trait_df)
}