#' Find the score for each trait for a given cluster.
#' @ param c_num - cluster number
#' @ param beta_mat - data matrix
#' @ param cluster_df - The cluster membership dataframe.
#' row names are the SNPs
#' Columns are:
#' "clust_num" the number of the cluster the snp has been assigned to.
#' "clust_prob" the probability the snp is in the cluster.
#' "clust_dist" the distance to the cluster centre.
#' @ param num_axis - The number of traits in use.
#' @ param pval_mat - The matrix of probabilities for the data in beta_mat.
#' @ param bp_on - Boolean, if TRUE the probabilities in pval_mat are used
#' in the cluster scoring. Default is TRUE
#' ------------------------------------------------------------------------
#' Create a dataframe of the score for each trait using \link{score_trait}.
#' Column bind together into cscore
#' Add cluster number to column in cscore.
#' @ return cscore
#' @ export
score_cluster <- function(c_num, beta_mat, clusters_df, num_axis,
                          pval_mat, bp_on = TRUE) {
  traits <- colnames(beta_mat)
  c_id <- paste0("na", num_axis, "_cn", c_num)
  snp_list <- rownames(clusters_df[clusters_df$clust_num == c_num, ])
  cprob_mat <- data.matrix(clusters_df[snp_list, "clust_prob"], rownames = TRUE)
  trait_score_df_list <- lapply(traits, score_trait,
                                c_id = c_id,
                                snp_list = snp_list,
                                beta_mat = beta_mat,
                                pval_mat = pval_mat,
                                cprob_mat = cprob_mat,
                                bp_on = bp_on)
  c_score <- Reduce(cbind, trait_score_df_list)
  c_score["clust_num"] <- c_num
  return(c_score)
}
