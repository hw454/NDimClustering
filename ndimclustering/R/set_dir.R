#' Set the directory using the base directory and iteration vairables
#'
#' @param res_dir0 - string for the base directory.
#' @param iter_traits - The dataframe containing the iteration variables.
#'
#' @description
#' Create the string for the directory from the iteration variables using
#' the function [method_str].
#' Combine the base directory and iteration directory.
#' Create the directory if it doesn't exist.
#'
#' @return the directory path as a string
#'
#' @export
set_directory <- function(res_dir0, iter_traits) {
  # Set the directory for the results using the base directory
  # and the iteration parameters
  res_dir <- paste0(res_dir0, make_path_label_str(iter_traits), "/")
  # Create results directory if it doesn't exist
  if (!dir.exists(res_dir)) {
    dir.create(res_dir)
  }
  return(res_dir)
}