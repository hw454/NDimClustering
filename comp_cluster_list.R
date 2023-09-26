#' Compare the cluster cn1 to all the clusters in the cluster list.
#' Outpout data frame of the difference in cluster scores between each pair.
#' -------------------------------------------------------------------------
#' @ param cn1 - The number of the cluster being compared.
#' @ param c_nums - The list of all the cluster numbers.
#' @ param c_score0 - The dataframe of cluster scores.
#' @ param axis - The traits/ columns the clusters are being compared on.
#' @ param clust_norm - The type of norm to be used when compareing the
#' difference accross multiple axis. Default is "F" which is the
#' Froebenius norm.
#' ----------------------------------------------------------------------
#' As each difference is found the pair is added to the dataframe `diff_df`.
#' This dataframe has columns "c_num1", "c_num2", "diff"
#' @ return diff_df
#' @ export
compare_oneclust_tolist <- function(cn1, c_nums, c_score0, axis,
                                    clust_norm = "F") {
  #FIXME remove cn1 from c_nums each time.
  diff_score_list <- lapply(c_nums, clust_score_diff,
                            cn1 = cn1,
                            c_score0 = c_score0,
                            axis = axis,
                            norm_typ = clust_norm)
  diff_score_list <- diff_score_list[!sapply(diff_score_list, is.null)]
  diff_scores <- Reduce(rbind, diff_score_list)
  return(diff_scores)
}