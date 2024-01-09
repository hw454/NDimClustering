#' Setup your data for running clustering functions
#'
#' @description The function is written to use create variables that will be
#' needed globally for running the package functions.
#'
#' @param threshmul The multiplier for the thresholding of the cluster
#'   score difference. default\:5
#' @param clust_threshold The threshold for cluster centre difference to
#'   determine whether clusters have converged.
#' @param na_percent The percentage of a column that must be not NaN.
#'   default\:0.8
#' @param nr The maximum number of clusters. default\:10
#' @param np The number of principal components. default\:np
#' @param clust_norm The type of norm used for finding the
#'   distance between points. default is the Froebenius norm "F".
#' @param thresh_norm The type of norm used for finding the
#'   distance between scores. default is the Froebenius norm "F".
#'
#' @details
#' Variables are created and declared globally.
#' thresholds <- ("threshmul" = threshmul,
#'                  "diff" = 1e-5,
#'                 "clust" = clust_threshold)
#' na_handling <- ("narm" = TRUE, "percent" = na_percent)
#' nums <- ("nr" = nr, "np" = np, "max_dist" =  1.0)
#' norm_typs <- ("clust" =  clust_norm, "thresh_norm" =  thresh_norm)
#'
#' @return empty
#'
#' @export
setup_algorithm_data <- function(threshmul = 5.0,
                               clust_threshold = 1e-5,
                               na_percent = 0.8,
                               nr = 10,
                               np = 3,
                               clust_norm = "F",
                               thresh_norm = "F") {
# Contain the variables into lists, threshold$diff and
# nums$max_dist will both be updated later.
thresholds <<- list("threshmul" = threshmul,
                  "diff" = 1e-5,
                  "clust" = clust_threshold)
na_handling <<- list("narm" = TRUE, "percent" = na_percent)
nums <<- list("nr" = nr, "np" = np, "max_dist" =  1.0)
norm_typs <<- list("clust" =  clust_norm, "thresh_norm" =  thresh_norm)
return()
}
