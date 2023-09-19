clust_pca_compare_iterative <- function(data_matrices,
                         out_pheno,
                         thresholds = list("threshmul" = 5,
                                            "diff" = 1e-5,
                                            "clust" = 1e-5),
                         na_handling = list("narm" = TRUE, "percent" = 0.95),
                         iter_traits = list(iter = 0,
                                            "bp_on" = TRUE,
                                            "clust_prob_on" = TRUE,
                                            "clust_typ" = "basic"),
                         norm_typs = list("clust" = "F", "thresh" = "F"),
                         nums = list("max_dist" = 1, "np" = 3, "nr" = 5)
                        ) {
  #' Iterate through the traits in trait_info contained in the
  #' data_matrices list. Find the principal components.
  #' Transform the data onto these.
  #' Cluster the transformed data at each iteration.
  #' If there is a distinct difference between two clusters exit.
  trait_df <- data.frame(
    label = character(),
    axes_ind = integer(),
    b_df_ind = integer(),
    nSNPs = integer()
  )
  # Data frame for recording the cluster scores.
  c_scores <- data.frame(
    num_axis = integer(),
    clust_num = integer()
  )
  # Add np columns for each PC
  c_scores <- add_np_pca_cols(c_scores, nums$np)
  # Initialise with outcome
  trait_df <- data.frame(label = out_pheno,
                       axes_ind = which(data_matrices$trait_info$phenotype == out_pheno)[1] #nolint
                       )
  max_df <- data.frame(num_axis = integer(),
                       cn1 = integer(),
                       cn2 = integer(),
                       max_diff = integer())
  cluster_df <- data.frame("snp_id" = character(),
                          "clust_num" = integer(),
                          "clust_prob" = numeric(),
                          "clust_dist" = numeric(),
                          "num_axis" = integer())
  cluster_df <- add_np_pca_cols(cluster_df, nums$np)
  df_list <- list("clust_scores" = c_scores,
                  "max_diff" = max_df,
                  "e_list" = list(),
                  "trait" = list(),
                  "cluster_items" = cluster_df)
  for (a in data_matrices$trait_info$phenotype){
    # Update the traits forming the axis
    # Add the trait to the trait dataframe
    covered <- (a %in% trait_df$label)
    allna <- na_col_check(data_matrices$b_df[, a], na_handling$percent)
    if (!covered && !allna) {
      # Add trait to trait dataframe
      trait_df <- trait_df_add_a(trait_df,
                                   a = a,
                                   axes = data_matrices$trait_info$phenotype)
      # Num axis
      num_axis <- length(trait_df$label)
      trait_lab <- paste0("traits_", num_axis)
      trait_list <- list(trait_df)
      names(trait_list) <- trait_lab
      df_list$trait <- append(df_list$trait, trait_list)
      # If the trait is not all NaN then run clustering.
      print(paste("PCA on", num_axis, "axes"))
      print(paste("New axis", a))
      # Get the data upto this axis
      b_iter_df <- data_matrices$beta[, trait_df$label]
      se_iter_df <- data_matrices$se[, trait_df$label]
      pval_iter_df <- data_matrices$pval[, trait_df$label]
      # Cluster the data on these axes
      pca_list   <- pca(b_iter_df, pval_iter_df, se_iter_df,
                        nums$np, na_handling$narm)
      b_pc_mat    <- pca_list$beta
      p_pc_mat <- pca_list$pval
      t_mat       <- pca_list$transform
      # Store the matrices for result output
      df_list$e_list <- append(df_list$e_list, list(a = t_mat))
      # Get column names for PCs
      pc_cols <- colnames(b_pc_mat)
      # Before clustering find the maximum distance and
      # distance variance on the transformed axes.
      nums$max_dist <- max_dist_calc(b_pc_mat,
                                norm_typ = norm_typs$clust,
                                na_rm = na_handling$narm)
      dist_mat <- setup_dist(b_pc_mat, norm_typ = norm_typs$clust)
      thresholds$diff <- thresholds$threshmul * var(dist_mat$dist,
                                                    na.rm = na_handling$narm)
      # Cluster the data on these axes
      if (grepl("angle", iter_traits$clust_typ)){
          st <- "angle"
      } else {
          st <- "regular"
       }
      if (grepl("min", iter_traits$clust_typ)) {
        cluster_df <- cluster_kmeans_min(b_pc_mat,
                                         nums$nr,
                                         nums$max_dist,
                                         space_typ = st,
                                         clust_prob_on = iter_traits$clust_prob_on, # nolint
                                         norm_typ = norm_typs$clust,
                                         threshold = thresholds$clust,
                                         narm = na_handling$narm)
      } else if (grepl("basic", iter_traits$clust_typ)) {
        cluster_df <- cluster_kmeans_basic(b_pc_mat,
                                           nums$nr,
                                           nums$max_dist,
                                           space_typ = st,
                                           clust_prob_on = iter_traits$clust_prob_on, # nolint
                                           threshold = thresholds$clust,
                                           norm_typ = norm_typs$clust,
                                           narm = na_handling$narm)
      }
      cluster_df["num_axis"] <- num_axis
      for (col in pc_cols) {
        cluster_df[col] <- b_pc_mat[col]
      }
      cluster_df["snp_id"] <- row.names(cluster_df)
      df_list$cluster_items <- rbind(df_list$cluster_items,
                                    cluster_df)
      # Find the set of cluster numbers
      c_nums <- unique(cluster_df$clust_num)
      # Score the clustered data based on affiliations with axes.
      # Find the score for each PC
      c_score0 <- clust_score(cluster_df,
                              beta_mat = b_pc_mat,
                              pval_mat = p_pc_mat,
                              bp_on = iter_traits$bp_on,
                              clust_prob_on = iter_traits$clust_prob_on,
                              num_axis = num_axis)
      # Iterate through each cluster and compare across the others to find if
      # any pair have a distinct difference.
      diff_score_list <- lapply(c_nums, compare_oneclust_tolist,
                                c_nums = c_nums,
                                c_score0 = c_score0,
                                axis = pc_cols,
                                clust_norm = norm_typs$clust)
      diff_score_list <- diff_score_list[!sapply(diff_score_list, is.null)]
      diff_scores <- Reduce(rbind, diff_score_list)
      # Find the pair with maximum cluster score diff and store in max_df0
      row <- which.max(diff_scores$diff)
      max_df0 <- diff_scores[row, ]
      max_df0["num_axis"] <- num_axis
      df_list$max_diff <- rbind(df_list$max_diff, max_df0)
      c_score0["num_axis"] <- num_axis
      df_list$clust_scores <- rbind(df_list$clust_scores, c_score0)
      # Plot the scatter for this iteraion
      clust_scatter(cluster_df, b_pc_mat,
                    pca_list$se,
                    iter_traits,
                    num_axis)
      print("scatter plot done")
      # Plot the transform heatmap.
      pca_list$transform %>% plot_transform_heatmap(iter_traits,
                                      num_axis)
      if (max_df0$diff > thresholds$diff && num_axis > 3) {
        print(paste("Threshold met on outcome", a))
        return(df_list)
      }
    }
  }
  print("None of the outcomes clustering met the thresholding test. ")
  return(df_list)
}

clust_pca_compare_all <- function(data_matrices,
                         out_pheno,
                         thresholds = list("threshmul" = 5,
                                            "diff" = 1e-5,
                                            "clust" = 1e-5),
                         na_handling = list("narm" = TRUE, "percent" = 0.95),
                         iter_traits = list(iter = 0,
                                            "bp_on" = TRUE,
                                            "clust_prob_on" = TRUE,
                                            "clust_typ" = "basic"),
                         norm_typs = list("clust" = "F", "thresh" = "F"),
                         nums = list("max_dist" = 1, "np" = 3, "nr" = 5)
                        ) {
  #' Iterate through the traits in trait_info contained in the
  #' data_matrices list. Find the principal components.
  #' Transform the data onto these.
  #' Cluster the transformed data at each iteration.
  #' If there is a distinct difference between two clusters exit.
  # Data frame for recording the cluster scores.
  c_scores <- data.frame(
    num_axis = integer(),
    clust_num = integer()
  )
  # Add np columns for each PC
  c_scores <- add_np_pca_cols(c_scores, nums$np)
  # Initialise with outcome
  trait_df <- data.frame(label = out_pheno,
                       axes_ind = which(data_matrices$trait_info$phenotype == out_pheno)[1] #nolint
                       )
  max_df <- data.frame(num_axis = integer(),
                       cn1 = integer(),
                       cn2 = integer(),
                       diff = integer())
  cluster_df <- data.frame("snp_id" = character(),
                          "clust_num" = integer(),
                          "clust_prob" = numeric(),
                          "clust_dist" = numeric(),
                          "num_axis" = integer())
  cluster_df <- add_np_pca_cols(cluster_df, nums$np)
  df_list <- list("clust_scores" = c_scores,
                  "max_diff" = max_df,
                  "e_mat" = data_matrices$beta,
                  "b_pc" = data_matrices$beta,
                  "se_pc" = data_matrices$se,
                  "trait" = data.frame,
                  "cluster_items" = cluster_df)
  # Fill trait_df with valid traits before running main program.
  trait_df <- make_trait_df(pheno_list = data_matrices$trait_info$phenotype,
                            data_mat = data_matrices$b_df,
                            na_percent = na_handling$percent
                            )
  # Num axis
  num_axis <- length(trait_df$label)
  df_list$trait <- trait_df
  # If the trait is not all NaN then run clustering.
  print(paste("PCA on", num_axis, "axes"))
  # Get the data upto this axis
  b_iter_df <- data_matrices$beta[, trait_df$label]
  se_iter_df <- data_matrices$se[, trait_df$label]
  pval_iter_df <- data_matrices$pval[, trait_df$label]
  # Cluster the data on these axes
  pca_list   <- pca(b_iter_df, pval_iter_df, se_iter_df,
                          nums$np, na_handling$narm)
  b_pc_mat    <- pca_list$beta
  p_pc_mat <- pca_list$pval
  t_mat       <- pca_list$transform
  # Store the matrices for result output
  df_list$e_mat <- t_mat
  df_list$b_pc <- b_pc_mat
  df_list$se_pc <- pca_list$se
  # Get column names for PCs
  pc_cols <- colnames(b_pc_mat)
  # Before clustering find the maximum distance and
  # distance variance on the transformed axes.
  nums$max_dist <- max_dist_calc(b_pc_mat,
                                  norm_typ = norm_typs$clust,
                                  na_rm = na_handling$narm)
  dist_mat <- setup_dist(b_pc_mat, norm_typ = norm_typs$clust)
  thresholds$diff <- thresholds$threshmul * var(dist_mat$dist,
                                                na.rm = na_handling$narm)
  # Cluster the data on these axes
  if (grepl("angle",iter_traits$clust_typ, fixed = TRUE)) {
    st <- "angle"
  } else {
    st <- "regular"
  }
  if (grepl("min", iter_traits$clust_typ)) {
    cluster_df <- cluster_kmeans_min(b_pc_mat,
                                      nums$nr,
                                      nums$max_dist,
                                      space_typ = st,
                                      clust_prob_on = iter_traits$clust_prob_on, # nolint
                                      norm_typ = norm_typs$clust,
                                      threshold = thresholds$clust,
                                      narm = na_handling$narm)
  } else if (grepl("basic", iter_traits$clust_typ)) {
    cluster_df <- cluster_kmeans_basic(b_pc_mat,
                                        nums$nr,
                                        nums$max_dist,
                                        space_typ = st,
                                        clust_prob_on = iter_traits$clust_prob_on, # nolint
                                        threshold = thresholds$clust,
                                        norm_typ = norm_typs$clust,
                                        narm = na_handling$narm)
  }
  cluster_df$num_axis <- num_axis

  df_list$clust_items <- rbind(df_list$clust_items, cluster_df)
  # Find the set of cluster numbers
  c_nums <- unique(cluster_df$clust_num)
  # Score the clustered data based on affiliations with axes.
  # Find the score for each PC
  c_score0 <- clust_score(cluster_df,
                                beta_mat = b_pc_mat,
                                pval_mat = p_pc_mat,
                                bp_on = iter_traits$bp_on,
                                clust_prob_on = iter_traits$clust_prob_on,
                                num_axis = num_axis)
  # Iterate through each cluster and compare across the others to find if
  # any pair have a distinct difference.
  diff_score_list <- lapply(c_nums, compare_oneclust_tolist,
          c_nums = c_nums,
          c_score0 = c_score0,
          axis = pc_cols,
          clust_norm = norm_typs$clust
        )
  diff_score_list <- diff_score_list[!sapply(diff_score_list, is.null)]
  diff_scores <- Reduce(rbind, diff_score_list)
  # Find the pair with maximum cluster score diff and store in max_df0
  row <- which.max(diff_scores$diff)
  max_df0 <- diff_scores[row, ]
  max_df0["num_axis"] <- num_axis
  df_list$max_diff <- rbind(df_list$max_diff, max_df0)
  c_score0["num_axis"] <- num_axis
  df_list$clust_scores <- rbind(df_list$clust_scores, c_score0)
  print("None of the outcomes clustering met the thresholding test. ")
  return(df_list)
}

make_trait_df <- function(pheno_list, data_mat, na_percent) {
  #' Create a dataframe with the trait labels and the index
  #' position for the trait. Only include traits whose column
  #'  in data_mat meets the na_percent threshold.
  trait_df_list <- lapply(pheno_list, check_trait,
                          pheno_list = pheno_list,
                          data_mat = data_mat,
                          na_percent = na_percent)
  trait_df <- Reduce(rbind, trait_df_list) # Combine all traits into one df.
  trait_df <- trait_df %>% distinct()  # Remove duplicate traits.
  return(trait_df)
}

check_trait <- function(a, pheno_list, data_mat, na_percent) {
    # Add the trait to the trait dataframe
    allna <- na_col_check(data_mat[, a], na_handling$percent)
    if (!allna) {
      a_ind <- which(pheno_list == a)[1]
        trait_single_df <- data.frame(
                                   label = a,
                                   axes_ind = a_ind)
    } else {
      trait_single_df <- data.frame(
                                    label = character(),
                                    axes_ind = integer()
                                    )
    }
    return(trait_single_df)
}


na_col_check <- function(b_col, percent = 0.95) {
  #' Check if the number of NaNs in a clolumn is more than 1-percent.
  n_accept <- length(b_col) * (1 - percent)
  narows <- which(b_col %>% is.na())
  if (length(narows) > n_accept) {
    # All rows are NaN so trait will be removed from trait
    return(1)
  } else {
    return(0)
    }
}

compare_oneclust_tolist <- function(cn1, c_nums, c_score0, axis,
                                    clust_norm = "F") {
  #' Compare all clusters in the list c_nums to cn1 and store these differences
  diff_score_list <- lapply(c_nums, clust_score_diff,
                            cn1 = cn1,
                            c_score0 = c_score0,
                            axis = axis,
                            norm_typ = clust_norm)
  diff_score_list <- diff_score_list[!sapply(diff_score_list, is.null)]
  diff_scores <- Reduce(rbind, diff_score_list)
  return(diff_scores)
}

clust_score_diff <- function(cn1, cn2, c_score0, axis, norm_typ = "F") {
  #' Find the score difference between two clusters.
  #' If NaN return Null
  cs1 <- as.matrix(c_score0[c_score0$clust_num == cn1, axis])
  cs2 <- as.matrix(c_score0[c_score0$clust_num == cn2, axis])
  # If no members to cluster with score on that axis it will be NaN
  if (!all(is.na(cs1)) && !all(is.na(cs2)) && nrow(cs1) && nrow(cs2)) {
    metric_score <- clust_metric(cs1, cs2, norm_typ)
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
  #' Create dataframe with columns with labels Pi
  #' for i in 1 to np to the dataframe.
  p_cols <- lapply(1:np, p_col_df)
  np_df <- Reduce(cbind, p_cols)
  full_df <- cbind(df, np_df)
  return(full_df)
}

p_col_df <- function(i) {
  #' Create a dataframe with column Pi
  cname <- paste0("P", i)
  out <- data.frame(col = numeric())
  colnames(out) <- c(cname)
  return(out)
}

pc_name_list <- function(np) {
  p_cols <- lapply(1:np, p_name)
  return(p_cols)
}

p_name <- function(i) {
  #' Create a dataframe with column Pi
  cname <- paste0("P", i)
  return(c_name)
}