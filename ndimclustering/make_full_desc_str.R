#' Create a sentence string describing the iteration.
#'
#' @description For titles and plot captions create a string with
#'  correct sentence structure which decribes the options in the iteration.
#'
#' @param iter_traits - Data frame containing the iteration variables
#'   "bp_on", "clust_prob_on",  "clust_typ"
#'
#' @details
#'  The string will start with the clustering type, then the bp switch,
#'  then the cluster probability switch.
#'
#' @return description_str
#'
#' @export
create_full_desc_str <- function(iter_traits) {
  #' Create the string that describes the type of method for titles and captions
  if (iter_traits$bp_on) {
    bp_str <- "bp on"
  } else {
    bp_str <- "bp off"
  }
  if (iter_traits$clust_prob_on) {
    clust_prob_str <- "ClustProb on"
  } else {
    clust_prob_str <- "ClustProb off"
  }
  return(paste(iter_traits$clust_typ, "and", bp_str, "and", clust_prob_str))
}