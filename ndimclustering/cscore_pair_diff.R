  #' Find the score difference between two clusters.
  #'
  #' @param cn1 - First cluster number
  #' @param cn2 - Second cluster number
  #' @param axis - The traits/ columns to be included.
  #' @param norm_typ - The type of norm to be used when
  #'   calculating the difference accross multiple distance.
  #'   The default is the Froebenius norm "F".
  #'
  #' @description
  #'   score_diff is the dataframe with terms:
  #'   * "c_num1" - the first cluster number
  #'   * "c_num2" - the second cluster number
  #'   * "diff" - The norm of the difference between the
  #'     two cluster scores accross the axis.
  #'
  #' @details
  #' If the difference is NaN then score_diff is returned
  #' with the same columns but empty rows.
  #'
  #' @return score_diff
  #'
  #' @export
cscore_pair_diff <- function(cn1, cn2, c_score0, axis, norm_typ = "F") {
  cs1 <- as.matrix(c_score0[c_score0$clust_num == cn1, axis])
  cs2 <- as.matrix(c_score0[c_score0$clust_num == cn2, axis])
  # If no members to cluster with score on that axis it will be NaN
  if (!all(is.na(cs1)) && !all(is.na(cs2)) && nrow(cs1) && nrow(cs2)) {
    metric_score <- dist_funcs::clust_metric(cs1, cs2, norm_typ)
    if (is.na(metric_score)) {
      return()
    } else {
      score_diff <- data.frame(
        c_num1 = cn1,
        c_num2 = cn2,
        diff = metric_score
      )
    }
  } else {
    score_diff <- data.frame(
        c_num1 = integer(),
        c_num2 = integer(),
        diff = numeric()
      )
  }
  return(score_diff)
}