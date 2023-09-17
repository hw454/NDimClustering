clust_score <- function(clusters_df, beta_mat, pval_mat,
                        bp_on = TRUE,
                        clust_prob_on = TRUE,
                        num_axis = ncol(beta_df)) {
  #' The dataframe of clusters is used alongside the beta and pvalue
  #' for each SNP to score the clusters on their association with
  #' each trait. Information for the traits is stored in aim_df
  clust_nums <- unique(clusters_df$clust_num)
  clust_scores_list <- lapply(clust_nums, score_cluster,
    beta_mat = beta_mat, clusters_df = clusters_df,
    pval_mat = pval_mat, bp_on = bp_on, num_axis = num_axis)
  clust_scores <- Reduce(rbind, clust_scores_list)
  return(clust_scores)
}
score_cluster <- function(c_num, beta_mat, clusters_df, num_axis,
                          pval_mat, bp_on) {
  traits <- colnames(beta_mat)
  c_id <- paste0("na", num_axis, "_cn", c_num)
  snp_list <- rownames(clusters_df[clusters_df$clust_num == c_num, ])
  clust_probs <- clusters_df[snp_list, "clust_prob"]
  trait_score_df_list <- lapply(traits, trait_score,
                                c_id = c_id,
                                snp_list = snp_list,
                                beta_mat = beta_mat,
                                pval_mat = pval_mat,
                                clust_probs = clust_probs,
                                bp_on = bp_on)
  c_score0 <- Reduce(cbind, trait_score_df_list)
  c_score0["clust_num"] <- c_num
  return(c_score0)
}

trait_score <- function(beta_mat, pval_mat, snp_list, a, c_id,
                        clust_probs = TRUE, bp_on = TRUE) {
  #' Get the score for all points in the cluster weighting by the
  #' pvalues of the association and their weighting within the cluster
  # Get the row indices for the traits in the cluster
  b_sub <- na.omit(subset(beta_mat, row.names = snp_list, select = a))
  c_sub <- na.omit(clust_probs)
  if (bp_on) {
    p_sub <- na.omit(subset(pval_mat, row.names = snp_list, select = a))
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