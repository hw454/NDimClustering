#' Calculate the distance between a pair of snps and return dataframe
#'
#' @param score_mat the data matrix. Rows correspond to snps.
#' @param snp1 the first snp name to use.
#' @param snp2 the second snp name to use
#' @param norm_typ the type of norm to use in the distance calculation.
#'   The default is the Froebenius norm "F".
#'
#' @description
#'   The distance between the snps is given by the norm of the difference
#'   between their co-ordinates. With the co-ordinates being the corresponding
#'   rows in the "score_mat". The type of norm used given by the "norm_typ".
#'   NaNs are ignored so snps are only compared on common axis.
#'   The output dataframe "dist_df" is 1 row with
#'     * "snp1" the name of the first snp
#'     * "snp2" the name of the second snp
#'     * "dist" the distance between the snps.
#'
#' @return dist_df
#'
#' @family distance_functions
#'
#' @export
calc_pair_dist_df <- function(score_mat, snp1, snp2, norm_typ = "F") {
  # Find the metric distance between the points given by snp1
  # and snp2 in score_df using the norm_typ metric.
  x1 <- score_mat[rownames(score_mat) == snp1]
  x2 <- score_mat[rownames(score_mat) == snp2]
  xp <- data.matrix(stats::na.omit(x1 - x2))
  d <- norm(xp, norm_typ)
  dist_df <- data.frame(snp1 = snp1,
                      snp2 = snp2,
                      dist = d)
  return(dist_df)
}