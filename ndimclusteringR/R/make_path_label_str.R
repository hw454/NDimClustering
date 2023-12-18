#' Create the string representing the variations in the iterations.
#'
#' @description
#'  This string will have no space so it is suitable for file and
#'  directory names.
#'
#' @param iter_traits - Data frame containing the iteration variables
#'   "bp_on", "clust_prob_on",  "clust_typ"
#'
#' @details
#' The string will start with the clustering type, then the bp switch,
#' then the cluster probability switch.
#'
#' @return label_str
#'
#' @export
make_path_label_str <- function(iter_traits) {
  # Create the string that describes the type of method for filenames
  if (iter_traits$bin_p_clust) {
    bp_str <- "_bpON"
  } else {
    bp_str <- "_bpOFF"
  }
  if (iter_traits$bin_angles) {
    angle_str <- "_anglesON"
  } else {
    angle_str <- "_anglesOFF"
  }
  dir_name <- paste0(iter_traits$clust_typ,
    bp_str,
    angle_str,
    iter_traits$how_cents,
    "_", iter_traits$n_pcs, "pcs",
    "_", iter_traits$nclust, "clusts/"
  )
  return(dir_name)
}