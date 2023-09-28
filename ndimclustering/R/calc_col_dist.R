#' Find the distance between a snp and all other snps.
#'
#' @param score_mat Matrix of scores for each snp
#' @param snp1 the snp name to compare all others to
#' @param norm_typ the type of norm to use in distance calculations.
#'   default is the Froebenius norm "F".
#'
#' @description
#'   Use [calc_pair_dist_df] to find the distance between pairs of snps.
#'   The results from this for each pair is bound on the rows using
#'   [SparkR::rbind] into the datframe "dist_df".
#'
#' @return dist_df
#'
#' @export
dist_col_calc <- function(score_mat, snp1, norm_typ = "F") {
  dist_list <- lapply(rownames(score_mat), calc_pair_dist_df,
                      score_mat = score_mat,
                      snp1 = snp1,
                      norm_typ = norm_typ)
  dist_df <- Reduce(rbind, dist_list)
  return(dist_df)
}