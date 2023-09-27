#' Take the matrices with matching label name from each list
#' and column bind them
#'
#' @param name The name of the matrix in each list.
#' @param mat_list_a first matrix list, includes a matrix labelled name
#' @param mat_list_b second matrix list, includes a matrix labelled name
#'
#' @return matrix combine the matrix named name in each list. Bound on columns.
#'
#' @export
join_mat_on_name <- function(name, mat_list_a, mat_list_b) {
  out <- cbind(mat_list_a[name], mat_list_b[name])
  return(out)
}