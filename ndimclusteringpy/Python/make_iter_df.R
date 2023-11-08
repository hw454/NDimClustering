#' Iteration setup function.
#'
#' @description
#' Create a dataframe whose rows correspond to iterations of the
#' [full_cluster_and_plot] function. Each row indicates a different
#' set of input terms.
#'
#' @param clust_typ_list - List of the types of clustering to run on
#' @param bp_on_list - List of the variations of whether to include the
#'   probability accoisated with the scores.
#' @param clust_prob_on_list - List of the variations of whether to
#'   include the probability of whether a point is in a cluster based on
#'   it's distance to the centre in cluster scoring.
#' @param ndim_typ - List of options for whether to iterate through
#'   terms trait_list or run the program on them all.
#'
#' @details
#' iter_df = data.frame() with columns: "bp_on", "clust_prob_on"
#' "ndim_typ", "clust_typ".
#' Each row is numbered and the dataframe contains all combinations
#' from the input lists.
#'
#' @return iter_df
#'
#' @export
make_iter_df <- function(clust_typ_list,
                        bp_on_list,
                        clust_prob_on_list,
                        ndim_typ) {
  iter_df_full <- data.frame(row.names = integer(),
                             "bp_on" = logical(),
                             "clust_prob_on" = logical(),
                             "clust_typ" = character(),
                             "ndim_typ" = character())
  for (clust_typ_str in clust_typ_list) {
    for (bp_on in bp_on_list) {
      for (clust_prob_on in clust_prob_on_list) {
        for (ndim in ndim_typ) {
          iter_traits <- data.frame(
            "bp_on" = bp_on,
            "clust_prob_on" = clust_prob_on,
            "clust_typ" = clust_typ_str,
            "ndim_typ" = ndim)
          iter_df_full <- rbind(iter_df_full, iter_traits)
        }
      }
    }
  }
  return(iter_df_full)
}