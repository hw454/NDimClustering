clust_score <- function(clusters_df, beta_df, pval_df, aim_df, 
bp_on = TRUE, clust_prob_on = TRUE) {
  #' The dataframe of clusters is used alongside the beta and pvalue
  #' for each SNP to score the clusters on their association with
  #' each trait. Information for the traits is stored in aim_df
  traits     <- colnames(beta_df)
  num_axis   <- length(traits)
  clust_nums <- unique(clusters_df$cluster)
  for (i in clust_nums) {
    # Store the cluster scores for each cluster number in c_score0
    c_score0 <- data.frame(
      id = character(),
      num_axis = integer(),
      clust_num  = integer(),
      clust_size = integer()
    )
    # Cluster id
    c_id <- paste0("na", num_axis, "_cn", i)
    # Get the associated SNPs with the cluster
    SNP_list <- names(clusters_df$cluster[clusters_df$cluster==i]) # nolint
    # ? How to deal with multiple SNP entries in data input.
    # Check the number of associated SNPS
    l1 <- length(SNP_list)
    c_score0 <- c_score0 %>% add_row(num_axis = num_axis, clust_num = i,
                                     clust_size = l1)
    for (a in traits){
      # Add the cluster score for each axis to c_score0
      # Track the number of terms in each cluster.
      if (clust_prob_on) {
        clust_probs <- clusters_df$clust_prob[SNP_list]
        names(clust_probs) <- SNP_list
      } else {
        clust_probs <- clusters_df$clust_prob[SNP_list]
        names(clust_probs) <- SNP_list
        clust_probs[SNP_list] <- 1.0
      }
      ax_score <- axis_score(beta_df, pval_df, SNP_list,
      aim_df, a, clust_probs, bp_on)
      # Account for the probability of being in the cluster
      c_score0["id"] <- c_id
      if (!(a %in% colnames(c_score0))) {
        c_score0[a] <- ax_score
      } else {
        cscore_col <- c_score0[c_score0$id == c_id, a]
        c_score0[c_score0$id == c_id, a] <- cscore_col + ax_score
      }
    }
    #' Add the cluster scores for that cluster number to the dataframe
    #' for them all. The full dataframe is created in the first look.
    #' This is done instead of initialisation due to columns being
    #' set in the loop.
    if (i == clust_nums[1]) {
      clust_scores <- c_score0
    } else {
    clust_scores <- rbind(c_score0, clust_scores)
    }
  }
  return(clust_scores)
}

clust_pca_score <- function(clusters_df, b_pc_mat, p_pc_mat,
bp_on = TRUE, clust_prob_on = TRUE, num_axis = 0) {
  #' The dataframe of clusters is used alongside the beta and
  #' pvalue for each SNP to score the clusters on their
  #' association with each trait. Information for the traits
  #' is stored in aim_df
  clust_nums <- unique(clusters_df$cluster)
  # Number of PCs
  np <- dim(b_pc_mat)[2]
  for (i in clust_nums) {
    # Store the cluster scores for each cluster number in c_score0
    c_score0 <- data.frame(
      id = character(),
      num_axis = integer(),
      clust_num  = integer(),
      clust_size = integer()
    )
    # Add columns for the PCs
    c_score0 <- add_np_pca_col(c_score0, np)
    # Cluster id
    c_id <- paste0("na", num_axis, "_cn", i)
    # Get the associated SNPs with the cluster
    snp_list <- names(clusters_df$cluster[clusters_df$cluster == i])
    # ? How to deal with multiple SNP entries in data input.
    # Check the number of associated SNPS
    l1 <- length(snp_list)
    c_score0 <- c_score0 %>% add_row(num_axis = num_axis, clust_num = i,
                                     clust_size = l1)
    for (c in colnames(b_pc_mat)){
      # Add the cluster score for each axis to c_score0
      # Track the number of terms in each cluster.
      if (clust_prob_on) {
        clust_probs <- clusters_df$clust_prob[snp_list]
        names(clust_probs) <- snp_list
      } else {
        clust_probs <- clusters_df$clust_prob[snp_list]
        names(clust_probs) <- snp_list
        clust_probs[snp_list] <- 1.0
      }
      ax_score <- axis_score(b_pc_mat, p_pc_mat, snp_list,
      c, clust_probs, bp_on)
      # Account for the probability of being in the cluster
      c_score0["id"] <- c_id
      if (!(c %in% colnames(c_score0))) {
        c_score0[c] <- ax_score
      } else {
        c_score_col <- c_score0[c_score0$id == c_id, c]
        c_score0[c_score0$id == c_id, c] <- c_score_col + ax_score
      }
    }
    # Add the cluster scores for that cluster number
    # to the dataframe for them all
    # The full dataframe is created in the first look. This is done instead of
    # initialisation due to columns being set in the loop.
    if (i == clust_nums[1]) {
      clust_scores <- c_score0
    } else {
      clust_scores <- rbind(c_score0, clust_scores)
      }
  }
  return(clust_scores)
}


axis_score <- function(beta_df, pval_df, snp_list, a, clust_probs,
bp_on = TRUE) {
  # Get the row indices for the traits in the cluster
  b_rows <- which(rownames(beta_df) %in% snp_list)
  # Get the pathway colocization score from the input data for the SNP
  if (bp_on) {
    bp <- beta_df[b_rows, a] * pval_df[b_rows, a]
    snp_assoc <- abs(bp * clust_probs[snp_list])
    # beta_df is a matrix with rows labelled by SNP_id
    # and columns labelled by trait label.
  } else {
    snp_assoc <- abs(beta_df[b_rows, a] * clust_probs[snp_list])
    # beta_df is a matrix with rows labelled by SNP_id
    # and columns labelled by trait label.
  }
    return(mean(snp_assoc))
}