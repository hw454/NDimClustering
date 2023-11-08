#' make_clust_col_name
#'
#' @description Function for creating the cluster label.
#'
#' @param c_num integer for cluster number
#'
#' @return "clust_"c_num
#'
#' @family label_functions
#'
#' @export
make_clust_col_name <- function(c_num) {
return(paste0("clust_", c_num))
}