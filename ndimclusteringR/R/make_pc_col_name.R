#' make_pc_col_name
#'
#' @description Function for creating the cluster label.
#'
#' @param i integer for PC
#'
#' @return "PC"c_num
#'
#' @family label_functions
#'
#' @export
make_pc_col_name <- function(i) {
return(paste0("PC", i))
}