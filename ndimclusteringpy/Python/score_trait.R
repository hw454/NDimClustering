#' Get the score for a cluster on the trait axis.
#'
#' @param beta_mat - Matrix of the data. Rows-snps, columns-traits.
#' @param pval_mat - Probabilities associated with beta_mat
#' @param snp_list - The list of snps in this cluster
#' @param a - The trait axis
#' @param c_id - The id for the cluster score.
#' @param cprob_mat - The probability the snp is in the cluster.
#' @param bp_on - Boolean switch, if TRUE (default)
#'   use the pval_mat to weight the snp contribution to the cluster score.
#'
#' @description
#' b_sub=Get the data at the snp_rows and trait column of beta_mat
#' c_sub=Get the data at the snp_rows and trait column of cprob_mat
#' Omit any NaNs.
#' if bp_on:
#'  p_sub = Get the data at the snp_rows and trait column of pval_mat
#'  out = b_sub*p_sub*c_sub/ (p_sub*c_sub)
#' else:
#'  out = b_sub * c_sub / (sum(c_sub))
#'
#' @return out
#'
#' @family scoring_functions
#'
#' @export
score_trait <- function(beta_mat, pval_mat, snp_list, a, c_id,
                        cprob_mat, bp_on = TRUE) {
  b_sub <- stats::na.omit(beta_mat[snp_list, a])
  c_sub <- stats::na.omit(cprob_mat[snp_list, ])
  if (bp_on) {
    p_sub <- stats::na.omit(pval_mat[snp_list, a])
    snp_assoc <- abs(b_sub * p_sub * c_sub)
    total_probs <- sum(p_sub * c_sub)
  } else {
    snp_assoc <- abs(b_sub * c_sub)
    total_probs <- sum(c_sub)
  }
  axis_snp_assoc <- sum(snp_assoc) / total_probs
  trait_score_df <- data.frame(row.names =  c_id,
                               a = axis_snp_assoc)
  colnames(trait_score_df) <- c(a)
  return(trait_score_df)
}