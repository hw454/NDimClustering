#' Check is a trait is valid
#'
#' @param a - The label for the trait.
#' @param pheno_list - The list of all traits.
#' @param data_mat - The matrix containing the data of interest.
#'   With columns labelled by trait labels.
#' @param na_percent - percentage of a column which is required to be NOT NaN.
#'
#' @description
#'   if trait passes [na_col_check]:
#'     trait_single_df = data.frame(label = a,
#'                              axes_ind = position in pheno_list)
#'   else: trait_single_df = data.frame()
#'
#' @return trait_single_df
#'
#' @family checking_functions
#'
#' @export
check_trait <- function(a, pheno_list, data_mat, na_percent) {
    # Add the trait to the trait dataframe
    allna <- check_col_na(data_mat[, a], na_percent)
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