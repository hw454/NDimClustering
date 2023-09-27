#' Create a set of randoms values
#' @description
#' The points are selected using a uniform distribution
#' between the maximum and minimum values on the axis a.
#' @param a the trait to create random values for.
#' @param min_max_df the dataframe of minimum and maximum
#'   values for each trait. The rows are trait labelled and
#'   the columns are "min" and "max"
#' @param n_cents the number of random values to create. default\:5
#' @details
#' Create a dataframe "cent_df" with column label "a" and the values are the
#' random values created using \link[stats]{runif}.
#' @return cent_df
#' @export
rand_cent <- function(a, min_max_df, n_cents = 5) {
  cent <- runif(n_cents, min_max_df[a, "min"], min_max_df[a, "max"])
  cent_df <- data.frame(a = cent)
  colnames(cent_df) <- c(a)
  return(cent_df)
}