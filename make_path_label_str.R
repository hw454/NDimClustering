#' Create the string representing the variations in the iterations.
#' @description
#'  This string will have no space so it is suitable for file and
#'  directory names.
#' @param iter_traits - Data frame containing the iteration variables
#'   "bp_on", "clust_prob_on",  "clust_typ"
#' @details
#' The string will start with the clustering type, then the bp switch,
#' then the cluster probability switch.
#' @return label_str
#' @export
make_path_label_str <- function(iter_traits) {
  #' Create the string that describes the type of method for filenames
  if (iter_traits$bp_on) {
    bp_str <- "_bpON"
  } else {
    bp_str <- "_bpOFF"
  }
  if (iter_traits$clust_prob_on) {
    clust_prob_str <- "_clustprobON"
  } else {
    clust_prob_str <- "_clustprobOFF"
  }
  return(paste0(iter_traits$clust_typ, bp_str, clust_prob_str))
}