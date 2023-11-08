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
testit::assert("Error in `make clust col name`.
               Cluster number should be an integer",
               typeof(c_num) == integer())
return(paste0("clust_", c_num))
}