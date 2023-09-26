na_col_check <- function(b_col, percent = 0.95) {
  #' Check if the number of NaNs in a clolumn is more than 1-percent.
  n_accept <- length(b_col) * (1 - percent)
  narows <- which(is.na(b_col))
  if (length(narows) > n_accept) {
    # All rows are NaN so trait will be removed from trait
    return(1)
  } else {
    return(0)
    }
}

check_trait <- function(a, pheno_list, data_mat, na_percent) {
    # Add the trait to the trait dataframe
    allna <- na_col_check(data_mat[, a], na_percent)
    if (!allna) {
      a_ind <- which(pheno_list == a)[1]
        trait_single_df <- data.frame(
                                   label = a,
                                   axes_ind = a_ind)
    } else {
      trait_single_df <- data.frame(
                                    label = character(),
                                    axes_ind = integer()
                                    )
    }
    return(trait_single_df)
}

test_all_na <- function(b_df, nan_col = "30600_irnt") {
  #' Test whether the function for checking the NaNs in a column works.
  test <- na_col_check(b_df[, nan_col])
  print("test")
  if (test) {
    return(1)
  } else {
    return(0)
  }
}