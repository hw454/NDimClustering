#' Calculate the distance between each snp and every other snp.
#' @description Using \link{calc_col_dist} and \link{calc_pair_dist}
#' find the distance between each snp and every other snp.
#' @param score_mat the data matrix. Each row is a snp and the values
#'   are the co-ordinates of the snp accross the trait axis.
#' @param norm_typ the type of norm to use in the distance calculations
#'   when comparig the snps. Default is the Froebenius norm "F".
#' @details "dist_df" is a dataframe with columns\:
#'   * "snp1" snp label
#'   * "snp2" snp label for the compared snp
#'   * "dist" the distance between snp1 and snp2
#' @return dist_df
#' @export
calc_all_dist_pairs <- function(score_mat, norm_typ = "F") {
  # Find the distance between all pairs of points.
  dist_df <- data.frame(
    snp1 = character(),
    snp2 = character(),
    dist = numeric()
  )
  dist_list <- lapply(rownames(score_mat), calc_col_dist,
                      score_mat = score_mat,
                      norm_typ = norm_typ)
  dist_df <- Reduce(rbind, dist_list)
return(dist_df)
}