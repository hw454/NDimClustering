clust_score <- function(clusters_df, beta_df, pval_df,
bp_on = TRUE, clust_prob_on = TRUE,   num_axis = ncol(beta_df)) {
  #' The dataframe of clusters is used alongside the beta and pvalue
  #' for each SNP to score the clusters on their association with
  #' each trait. Information for the traits is stored in aim_df
  clust_nums <- unique(clusters_df$clust_num)
  clust_scores_list <- lapply(clust_nums, score_cluster,
    beta_df = beta_df, clusters_df = clusters_df,
    pval_df = pval_df, bp_on = bp_on, num_axis = num_axis)
  clust_scores <- Reduce(rbind, clust_scores_list)
  return(clust_scores)
}
score_cluster <- function(c_num, beta_df, clusters_df, num_axis,
                          pval_df, bp_on) {
  traits <- colnames(beta_df)
  c_id <- paste0("na", num_axis, "_cn", c_num)
  snp_list <- rownames(clusters_df[clusters_df$clust_num == c_num, ])
  clust_probs <- clusters_df[snp_list, "clust_prob"]
  trait_score_df_list <- lapply(traits, axis_score,
                                c_id = c_id, 
                                snp_list = snp_list, 
                                beta_df = beta_df,
                                pval_df = pval_df, 
                                clust_probs = clust_probs, 
                                bp_on = bp_on)
  c_score0 <- Reduce(cbind, trait_score_df_list)
  c_score0["clust_num"] <- c_num
  return(c_score0)
}

axis_score <- function(beta_df, pval_df, snp_list, a, clust_probs, c_id,
bp_on = TRUE) {
    #' Get the score for all points in the cluster weighting by the
  #' pvalues of the association and their weighting within the cluster
  # Get the row indices for the traits in the cluster
  b_rows <- which(rownames(beta_df) %in% snp_list)
  b_sub <- na.omit(beta_df[snp_list, a])
  c_sub <- na.omit(clust_probs)
  if (bp_on) {
    p_sub <- na.omit(pval_df[snp_list, a])
    snp_assoc <- abs(b_sub * p_sub * clust_probs)
    total_probs <- sum(p_sub * clust_probs)
  } else {
    snp_assoc <- abs(b_sub * clust_probs)
    total_probs <- sum(clust_probs)
  }
  axis_snp_assoc <- sum(snp_assoc) / total_probs
  trait_score_df <- data.frame(row.names =  c_id, 
                               a = axis_snp_assoc)
  colnames(trait_score_df) <- c(a)
  return(trait_score_df)
}