clust_pca_compare <- function(unstd_beta_df, unstd_se_df, pval_df,
                         axes, threshold, thresh_norm, clust_threshold,
                         clus_norm, np = 3, nr = 10, which_clust = "basic",
                         bp_on = TRUE, clust_prob_on = TRUE) {
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
  aim_df <- data.frame(label = OUT_pheno,
                       axes_ind = which(trait_info$phenotype ==OUT_pheno)[1],
                       b_df_ind = which(colnames(unstd_beta_df) == OUT_pheno)[1],
                       nSNPs = sum(!is.na(unstd_beta_df[, OUT_pheno]))[1]
  )
  for (ai in 1:length(trait_info$phenotype)){
    # Update the traits forming the axis
    a <- trait_info$phenotype[ai]
    # Add the trait to the trait dataframe
    covered <- (a %in% aim_df$label)
    if (!covered){
      # If the axis is already included then add row is not run.
      aim_df <- aim_df_add_a(aim_df, a, trait_info$phenotype,
                             unstd_beta_df)
      if (allna) {
        #' If the trait column was removed from the beta_df during na removal then 
        #' remove the trait from the trait data frame and move to the next trait.
        aim_df <- aim_df[aim_df$label != a, ]
      }
      else{
        # If the trait is not all NaN then run clustering.
        print(paste0("PCA on ", length(aim_df$label), " axes"))
        print(paste0("New axis ", a))
              # Num axis
        num_axis <- length(aim_df$label)
        # Get the data upto this axis
        b_df <- unstd_beta_df[aim_df$label]
        se_df <- unstd_se_df[aim_df$label]
        p_df <- pval_df[aim_df$label]
        pca_list   <- pca(b_df, p_df, np, narm)
        b_pc_mat    <- pca_list$beta
        p_pc_mat    <- pca_list$pval
        t_mat       <- pca_list$Transform
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
        c_score0 <- clust_pca_score(cluster_df, b_pc_mat, p_pc_mat,
        bp_on, clust_prob_on)
        # Iterate through each cluster and compare across the others to find if 
        # any pair have a distinct difference.
        N2 <- length(c_nums)
        N1 <- N2-1
        for (i in 1:N1){
          for (j in (i+1):N2){
            ci <- c_nums[i]
            cj <- c_nums[j]
            cs1 <- as.matrix(c_score0[c_score0$clust_num == ci, pc_cols])
            cs2 <- as.matrix(c_score0[c_score0$clust_num == cj, pc_cols])
            metric_score <- clust_metric(cs1, cs2, thresh_norm)
            if (is.na(metric_score)) {
            } else {
              if (metric_score>threshold) {
                print(paste("Threshold met on outcome", a))
                c_scores <- rbind(c_scores, c_score0)
                return(c_scores)
              }
            }
          }
        }
      }
      c_scores <- rbind(c_scores, c_score0)
    }
  }
  print("None of the outcomes clustering met the thresholding test. ")
  return(c_scores)
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
#test2 <- testna(unstdBeta_df,trait_info$phenotype)
# test<-clust_compare(stdBeta_df ,SNP_ind,axes,threshold,thresh_norm, clust_norm)
#test<-clust_compare(unstdBeta_df,unstdSE_df,pval_df,tstat_df,
# trait_axes,threshold,thresh_norm,clust_norm)