cluster_and_plot <- function(data_matrices,
                            out_pheno,
                            iter = 1,
                            iter_df = data.frame(
                                        "iter" = 0,
                                        "bp_on" = FALSE,
                                        "clust_prob_on" = FALSE,
                                        "clust_typ" = "basic"),
                            thresholds = list("diff" = 1e-5, "clust" = 1e-5),
                            na_handling = list("narm" = TRUE, "percent" = 0.95),
                            norm_typs = list("clust" = "F", "thresh" = "F"),
                            nums = list("max_dist" = 1, "np" = 3, "nr" = 5)
) {
  iter_traits <- iter_df[iter, ]
  res_dir <- set_directory(res_dir0, iter_traits)
  # Find the distances between all points to initialise the threshold
  # for cluster difference.
  dist_df <- setup_dist(data_matrices$beta, norm_typs$clust)
  thresholds$diff <- thresholds$diff_mul * var(dist_df$dist, na.rm = TRUE)
  max_dist <- max(dist_df$dist, na.rm = na_handling$nam)
  nums$max_dist <- max_dist
  print("Begining algorithm for inputs")
  print(iter_traits)
  out <- clust_pca_compare(data_matrices = data_matrices,
                          out_pheno = out_pheno,
                          na_handling = na_handling,
                          iter_traits = iter_traits,
                          norm_typs = norm_typs,
                          nums = nums
  )
  print("Clust done")
  max_diff_df <- out$max_diff
  c_scores <- out$clust_scores
  print("Diff done")
  c_scores %>% plot_trait_heatmap(iter_traits)
  print("Heatmap plot done")
  max_diff_df %>% plot_max_diff(iter_traits)
  print("Diff plot done")
  out <- list("iter_df" = iter_df,
              "c_scores" = c_scores,
              "max_df" = max_diff_df)
  return(out)
}
setup_dist <- function(score_df, norm_typ = "F") {
  # Setup distances
  dist_df <- data.frame(
    snp1 = character(),
    snp2 = character(),
    dist = numeric()
  )
  dist_list <- lapply(rownames(score_df), dist_col_calc,
                      score_df = score_df,
                      norm_typ = norm_typ)
  dist_df <- Reduce(rbind, dist_list)
return(dist_df)
}

dist_col_calc <- function(score_df, snp1, norm_typ = "F") {
  dist_list <- lapply(rownames(score_df), pair_dist_calc,
                      score_df = score_df,
                      snp1 = snp1,
                      norm_typ = norm_typ)
  dist_df <- Reduce(rbind, dist_list)
  return(dist_df)
}

pair_dist_calc <- function(score_df, snp1, snp2, norm_typ = "F") {
  #' Find the metric distance between the points given by snp1
  #' and snp2 in score_df using the norm_typ metric.
  x1 <- score_df[rownames(score_df) == snp1]
  x2 <- score_df[rownames(score_df) == snp2]
  xp <- data.matrix(na.omit(x1 - x2))
  d <- norm(xp, norm_typ)
  dist_df <- data.frame(snp1 = snp1,
                      snp2 = snp2,
                      dist = d)
  return(dist_df)
}

set_directory <- function(res_dir0, iter_traits) {
  #' Set the directory for the results using the base directory
  #' and the iteration parameters
  res_dir <- paste0(res_dir0, method_str(iter_traits), "/")
  # Create results directory if it doesn't exist
  if (!dir.exists(res_dir)) {
    dir.create(res_dir)
  }
  return(res_dir)
}

method_str <- function(iter_traits) {
  #' Create the string that describes the type of method for filenames
  if (iter_traits$bp_on) {
    bp_str <- "_bpON"
  } else {
    bp_str <- "_bpOFF"
  }
  if (iter_traits$clust_prob_on) {
    clust_prob_str <- "_clustprobON"
  } else {
    clust_prob_str <- "_clustprobOFF"
  }
  return(paste0(iter_traits$clust_typ_str, bp_str, clust_prob_str))
}

desc_str <- function(iter_traits) {
  #' Create the string that describes the type of method for titles and captions
  if (iter_traits$bp_on) {
    bp_str <- "bp on"
  } else {
    bp_str <- "bp off"
  }
  if (iter_traits$clust_prob_on) {
    clust_prob_str <- "ClustProb on"
  } else {
    clust_prob_str <- "ClustProb off"
  }
  return(paste(iter_traits$clust_typ_str, "and", bp_str, "and", clust_prob_str))
}

make_iter_df <- function(clust_typ_list, bp_on_list, clust_prob_on_list) {
  #' Create a dataframe whose rows correspond to iterations of the
  #' ClusterAndPlot programe. Each row indicates a different
  #' set of input terms. 
  iter_df_full <- data.frame(row.names = integer(),
                             "bp_on" = logical(),
                             "clust_prob_on" = logical(),
                             "clust_typ" = character())
  for (clust_typ_str in clust_typ_list) {
    for (bp_on in bp_on_list) {
      for (clust_prob_on in clust_prob_on_list) {
        iter_traits <- data.frame(
          "bp_on" = bp_on,
          "clust_prob_on" = clust_prob_on,
          "clust_typ" = clust_typ_str)
        iter_df_full <- rbind(iter_df_full,iter_traits)
      }
    }
  }
  return(iter_df_full)
}

crop_mat_list <- function(mat_list, trait_df,
                        out_pheno, exp_pheno,
                        n_rows, n_col0, n_col1) {
  #' Crop all of the matrices in mat_list to include only the
  #' first n_rows and columns between n_col0 and n_col1,
  #' ensure that out_pheno and exp_pheno are included.
  mat_out <- lapply(mat_list, crop_mat_colnames,
                    num_rows = n_rows,
                    col_names = c(out_pheno, exp_pheno))
  names(mat_out) <- names(mat_list)
  trait_out <- trait_info[trait_info$phenotype %in% c(out_pheno, exp_pheno)]

  mat_crop <- lapply(mat_list, crop_mat_colnums,
                          num_rows = n_rows,
                          col0 = n_col0,
                          col1 = n_col1)
  names(mat_crop) <- names(mat_list)

  filt_mat_list <- lapply(names(mat_list),
                          join_out_mat_on_name,
                          mat_out_list = mat_out,
                          mat_crop_list = mat_crop
  )
  names(filt_mat_list) <- names(mat_list)

 trait_info   <- rbind(trait_info[n_col0:n_col1], trait_out)

  # Collect the matrices into one object
  data_matrices <- append(filt_mat_list, list("trait_info" = trait_info))
  return(data_matrices)
}

join_out_mat_on_name <- function(name, mat_out_list, mat_crop_list) {
  #' Take the matrices with matching label name from each list
  #' and column bind them
  out <- cbind(mat_out_list[name], mat_crop_list[name])
  return(out)
}

crop_mat_colnames <- function(mat, num_rows, col_names) {
  #' Return the matrix with row 1:num_rows and column names in col_names
  mat_out <- mat[1:num_rows,
                 which(colnames(mat) %in% col_names)]
  return(mat_out)
}
crop_mat_colnums <- function(mat, num_rows, col0, col1) {
  #' Return the matrix with rows 1: num_rows and columns col0:col1
  mat_out <- mat[1:num_rows, col0:col1]
  return(mat_out)
}