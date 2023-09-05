clust_pca_compare <- function(unstd_beta_df, unstd_se_df, pval_df, # nolint
                         diff_threshold, thresh_norm, clust_threshold,
                         clus_norm, np = 3, nr = 10, which_clust = "basic",
                         bp_on = TRUE, clust_prob_on = TRUE, narm = TRUE) {
  #' Iterate through the columns in axes and clusters the data.
  #' If there is a distinct difference between two clusters exit.
  aim_df <- data.frame(
    label = character(),
    axes_ind = integer(),
    b_df_ind = integer(),
    nSNPs = integer()
  )
  # Data frame for recording the cluster scores.
  c_scores <- data.frame(
    id = character(),
    num_axis = integer(),
    clust_num = integer(),
    clust_size = integer(),
    clust_prob = numeric(),
    NA_Count = integer()
  )
  # Add np columns for each PC
  c_scores <- add_np_pca_cols(c_scores, np)
  # Initialise with outcome
  out_col <- which(colnames(unstd_beta_df) == out_pheno)[1]
  aim_df <- data.frame(label = OUT_pheno,
                       axes_ind = which(trait_info$phenotype == out_pheno)[1],
                       b_df_ind = out_col,
                       nSNPs = sum(!is.na(unstd_beta_df[, out_pheno]))[1]
  )
  max_df <- data.frame(num_axis = integer(),
                       cn1 = integer(),
                       cn2 = integer(),
                       max_diff = integer())
  for (ai in 1:length(trait_info$phenotype)){
    # Update the traits forming the axis
    a <- trait_info$phenotype[ai]
    # Add the trait to the trait dataframe
    covered <- (a %in% aim_df$label)
    if (!covered) {
      # If the axis is already included then add row is not run.
      aim_df <- aim_df_add_a(aim_df, a, trait_info$phenotype,
                             unstd_beta_df)
      allna <- all_na_check(unstd_beta_df)
      if (allna) {
        #' If the trait column was removed from the beta_df
        #' during na removal then remove the trait from the
        #' trait data frame and move to the next trait.
        aim_df <- aim_df[aim_df$label != a, ]
      } else {
        # If the trait is not all NaN then run clustering.
        print(paste0("PCA on ", length(aim_df$label), " axes"))
        print(paste0("New axis ", a))
              # Num axis
        num_axis <- length(aim_df$label)
        b_iter_df <- unstd_beta_df[, aim_df$label]
        pval_iter_df <- pval_df[, aim_df$label]
        # Cluster the data on these axes
        # Get the data upto this axis
        b_df <- unstd_beta_df[,aim_df$label]
        se_df <- unstd_se_df[,aim_df$label]
        p_df <- pval_df[,aim_df$label]
        pca_list   <- pca(b_df, p_df, se_df, np, narm)
        b_pc_mat    <- pca_list$beta
        p_pc_mat    <- pca_list$pval
        se_pc_mat   <- pca_list$se
        t_mat       <- pca_list$transform
        # Get column names for PCs
        pc_cols <- colnames(b_pc_mat)
        # Cluster the data on these axes
        allna <- all_na_check(b_df)
        # Cluster the data on these axes
        if (which_clust == "min") {
          cluster_df <- cluster_kmeans_min(b_pc_mat, nr)
        } else {
          cluster_df <- cluster_kmeans_basic(b_pc_mat, nr,
          clust_threshold, clust_norm)
          }
        # Find the set of cluster numbers
        c_nums <- unique(cluster_df$cluster)
        # Score the clustered data based on affiliations with axes.
        # Find the score for each PC
        c_score0 <- clust_score(cluster_df, beta_df = b_pc_mat, pval_df = p_pc_mat,
                                bp_on = bp_on, clust_prob_on = clust_prob_on, 
                                num_axis = num_axis)
        # Iterate through each cluster and compare across the others to find if
        # any pair have a distinct difference.
        # Iterate through each cluster and compare across the others to find if
        # any pair have a distinct difference.
        n2 <- length(c_nums)
        n2 <- n2 - 1
        diff_score_list <- lapply(c_nums, compare_oneclust_tolist,
                                  c_nums = c_nums, c_score0 = c_score0,
                                  axis = colnames(b_pc_mat))
        diff_score_list <- diff_score_list[! sapply(diff_score_list, is.null)]
        diff_scores <- Reduce(rbind, diff_score_list)
        row <- which.max(diff_scores$diff)
        max_df0 <- diff_scores[row,]
        max_df <-rbind(max_df0,max_df)
        if (max(diff_scores$diff) > diff_threshold) {
          print(paste("Threshold met on outcome", a))
          return(c_scores)
        }
        c_score0["num_axis"] <- num_axis
        c_scores <- rbind(c_scores, c_score0)
      }
    }
  }
  print("None of the outcomes clustering met the thresholding test. ")
  out_list <- list("clust_scores" = c_scores,
                   "max_diff" = max_df)
  return(out_list)
}

clust_compare <- function(unstd_beta_df, unstd_se_df, pval_df, tstat_df,
                         threshold, thresh_norm, clust_threshold, nr,
                         clust_norm, which_clust = "basic", bp_on = TRUE,
                         clust_prob_on = TRUE) {

  #' Iterate through the columns in axes and clusters the data.
  #' If there is a distinct difference between two clusters exit.

  aim_df <- data.frame(
    label = character(),
    axes_ind = integer(),
    b_df_ind = integer(),
    nSNPs = integer()
  )
  # Data frame for recording the cluster scores.
  c_scores <- data.frame(
    id = character(),
    num_axis = integer(),
    clust_num = integer(),
    clust_size = integer(),
    clust_prob = numeric(),
    NA_Count = integer()
  )
  # Initialise with outcome
  b_col <- which(colnames(unstd_beta_df) == OUT_pheno)[1]
  aim_df <- data.frame(label = OUT_pheno,
                       axes_ind = which(trait_info$phenotype == OUT_pheno)[1],
                       b_df_ind = b_col,
                       nSNPs = sum(!is.na(unstdBeta_df[, OUT_pheno]))[1]
                       )
  for (a in trait_info$phenotype) {
    # Add the trait to the trait dataframe
    covered <- (a %in% aim_df$label)
    if (!covered) {
      # If the axis is already included then add row is not run.
      aim_df <- aim_df_add_a(aim_df, a, trait_info$phenotype,
                 unstd_beta_df)
      # Cluster the data on these axes
      allna <- all_na_check(unstd_beta_df)
      if (allna) {
        #' If the trait column was removed from the beta_df
        #' during na removal then remove the trait from the
        #' trait data frame and move to the next trait.
        aim_df <- aim_df[aim_df$label != a, ]
      } else {
        # If the trait is not all NaN then run clustering.
        print(paste0("Cluster on ", length(aim_df$label), " axes"))
        print(paste0("New axis ", a))
        num_axis <- length(aim_df$label)
        b_iter_df <- unstd_beta_df[, aim_df$label]
        pval_iter_df <- pval_df[, aim_df$label]
        # Cluster the data on these axes
        if (which_clust == "min") {
          cluster_df <- cluster_kmeans_min(b_iter_df, aim_df, nr)
          } else {
            cluster_df <- cluster_kmeans_basic(pval_iter_df, aim_df,
            nr, clust_threshold, clust_norm)
            }
        # Find the set of cluster numbers
        c_nums <- unique(cluster_df$cluster)
        # Score the clustered data based on affilliations with axes.
        # Find the score for this axis
        c_score0 <- clust_score(cluster_df,
                                b_iter_df, pval_iter_df,
                                bp_on = bp_on, clust_prob_on = clust_prob_on,
                                num_axis = num_axis) # Score across all axis
        # Initialise new column for new axis
        c_scores[a] <- numeric()
        # Iterate through each cluster and compare across the others to find if
        # any pair have a distinct difference.
        n2 <- length(c_nums)
        n2 <- n2 - 1
        diff_score_list <- lapply(c_nums, compare_oneclust_tolist,
                                  c_nums = c_nums, c_score0 = c_score0,
                                  axis = aim_df$label)
        diff_score_list <- diff_score_list[! sapply(diff_score_list, is.null)]
        diff_scores <- Reduce(rbind, diff_score_list)
        if (max(diff_scores$diff) > Threshold) {
          print(paste("Threshold met on outcome", a))
          return(c_scores)
        }
        c_score0["num_axis"] <- num_axis
        c_scores <- rbind(c_scores, c_score0)
      }
    }
  }
  print("None of the outcomes clustering met the thresholding test. ")
  return(c_scores)
}

compare_oneclust_tolist <- function(cn1, c_nums, c_score0, axis) {
  diff_score_list <- lapply(c_nums, clust_score_diff,
                            cn1 = cn1,
                            cscore0 =  c_score0,
                            axis = axis)
  diff_score_list <- diff_score_list[!sapply(diff_score_list, is.null)]
  diff_scores <- Reduce(rbind, diff_score_list)
  return(diff_scores)
}

clust_score_diff <- function(cn1, cn2, c_score0, axis) {
  #' Find the score difference between two clusters.
  #' If NaN return Null
  cs1 <- as.matrix(c_score0[c_score0$clust_num == cn1, axis])
  cs2 <- as.matrix(c_score0[c_score0$clust_num == cn2, axis])
  # If no members to cluster with score on that axis it will be NaN
  if (!all(is.na(cs1)) && !all(is.na(cs2)) && nrow(cs1) && nrow(cs2)) {
    metric_score <- clust_metric(cs1, cs2, thresh_norm)
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
    return()
  }
  return(score_diff)
}
add_np_pca_cols <- function(df, np) {
  p_cols <- lapply(1:np, p_col_df)
  np_df <- Reduce(cbind, p_cols)
  full_df <- cbind(df, np_df)
  return(full_df)
}

p_col_df <- function(i) {
  cname <- paste0("P", i)
  out <- data.frame(col = numeric())
  colnames(out) <- c(cname)
  return(out)
}