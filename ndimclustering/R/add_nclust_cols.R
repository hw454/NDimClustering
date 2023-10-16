#' Add "np" columns to a dataframe.
#'
#' @description Representing the clusters a column for the clust distance
#'   is added for each cluster.
#'   The label for each column is given by the function [make_clust_cols].
#'
#' @param df the dataframe to add columns to
#' @param np the number of columns to add.
#'
#' @return df
#'
#' @family dataframe_editors
#'
#' @export
add_np_cols <- function(df, nr) {
  # Create dataframe with columns with labels Pi
  # for i in 1 to np to the dataframe.
  p_cols <- lapply(1:nr, make_clust_col)
  np_df <- Reduce(cbind, p_cols)
  full_df <- cbind(df, np_df)
  return(full_df)
}