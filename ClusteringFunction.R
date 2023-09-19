### 3) Clustering
# Fitting K-Means clustering Model
find_ic <- function(clust, group_col, dist_c, num_axis) {
  # Find the BIC and AIC for the clusters.
  # Number of estimated parameters
  k <- length(unique(clust$clust_num))
  # Number of data points
  n <- nrow(clust)
  # Dimension of each data point
  d <- num_axis
  # Log-Likelihood term
  l <- tot_withins(clust, group_col, dist_c)
  aic_n <- d + 2 * k * l
  bic_n <- d + log(n) * k * l
  ic_out <- list("aic" = aic_n, "bic" = bic_n)
  return(ic_out)
}

tot_withins <- function(data,
                          group_col = "clust_num",
                          dist_c = "clust_dist") {
  ll <- lapply(unique(data[group_col]),
                    sum_sq_clusts,
                    data = data,
                    dist_c = dist_c)
  sig <- as.numeric((var(data[group_col])))
  l <- Reduce(sum, ll) / sig
  return(l)
}

sum_sq_clusts <- function(num, data, dist_c) {
  sq <- data[dist_c]**2
  s <- sum(sq)
  return(s)
}

kmeans_ic <- function(cents, data_mat,
                      clust_thres = 1e-5,
                      clust_norm = "F",
                      na_rm = TRUE,
                      clust_prob_on = TRUE) {
  set.seed(240) # setting seed
  clust_re <- kmeans_skip_nan(data_mat,
                              centers = cents,
                              clust_threshold = clust_thres,
                              norm_typ = "F",
                              na_rm = FALSE,
                              prob_on = clust_prob_on)
  ic <- find_ic(clust_re,
                group_col = "clust_num",
                dist_c = "clust_dist",
                num_axis = ncol(data_mat))
  clust_re <- clust_re %>% mutate("aic" = ic$aic)
  clust_re <- clust_re %>% mutate("bic" = ic$bic)
  clust_re <- clust_re %>% mutate("ncents" = cents)
  return(clust_re)
}
closest_clust <- function(snp_id, b_df, clust_re,
                          max_dist = 10.0, norm_typ = "F") {
  #' Find the cluster whose centre is closest to the point with id snp_id.
  #' Used for points with NaNs that are not included in the
  #' initial clustering algorithm
  snp_dist <- max_dist
  clust_num <- clust_re$clusters[snp_id]
  c_nums <- clust_re$cluster %>% unique()
  for (c_num in c_nums){
    cent <- clust_re$centers[c_num]
    snp_clust_dist <- norm(cent - b_df[snp_id, ], norm_typ)
    clust_re$clust_dist[snp_id] <- snp_clust_dist
    if (snp_clust_dist < snp_dist) {
      # Assign cluster number when distance is less than previous
      snp_dist <- snp_clust_dist
      clust_num <- c_num
    }
  }
  out_df <- data.frame(
    row.names = snp_id,
    clust_num = clust_num,
    clust_dist = snp_dist,
    clust_prob = 1.0
  )
  return(out_df)
}

clust_metric <- function(cs1, cs2, norm_typ) {
  if (length(cs1) <= 1) {
    return(abs(cs1 - cs2))
  } else {
    clustid_1 <- which(!is.na(cs1))
    clustid_2 <- which(!is.na(cs2))
    clust_ids <- intersect(clustid_1, clustid_2)
    if (length(clust_ids) == 0) {
      return(NaN)
      } else if (length(clust_ids) == 1) {
        return(abs(cs1[clust_ids] - cs2[clust_ids]))
        } else {
      return(norm(as.matrix(cs1[clust_ids] - cs2[clust_ids]), norm_typ))
    }
  }
}


trait_df_add_a <- function(trait_df, a, axes) {
  #' Add the trait a to the dataframe of traits.
  #' Add the traits label, it's index in the trait list and it's index in the
  #' beta dataframe.
  if (a %in% trait_df$label) {
    return(trait_df)
  } else {
    trait_df0 <- data.frame(label = a,
                            axes_ind = which(axes == a)[1]
                            )
    trait_df <- rbind(trait_df, trait_df0)
    return(trait_df)
  }
}

#- Main clustering functions

cluster_kmeans_basic <- function(b_df,
                                nr = 10,
                                max_dist = 10.0,
                                space_typ = "regular",
                                clust_prob_on = TRUE,
                                norm_typ = "F",
                                threshold = 1e-5,
                                narm = TRUE) {
  #' Using the association scores for each SNP accross traits cluster the traits
  #' using kmeans. Return the cluster setup which minimises AIC.

  # Columns
  ax <- colnames(b_df)
  # Filter NaNs before clustering
  if (space_typ == "angle") {
    # For each point in b_df_comp convert the score to the angle
    # between the vectors to the origin and the unit vectors on the axis.
    b_df_clust <- mat_to_angle_mat(b_df)
  } else {
    b_df_clust <- b_df
  }

  b_df_comp <- b_df_clust[complete.cases(b_df_clust), ]

  # Initial cluster dataframe
  clust_re <- kmeans_skip_nan(b_df_comp,
                              centers = nr,
                              clust_threshold = threshold,
                              norm_typ = norm_typ,
                              prob_on = clust_prob_on,
                              na_rm = narm)
  # cluster number identification for each observation
  snp_cluster_list <- lapply(setdiff(names(b_df_clust), names(b_df_comp)),
                            closest_clust,
                            b_df = b_df_clust,
                            clust_re = clust_re,
                            max_dist = max_dist,
                            norm_typ = norm_typ)
  nan_cluster_df <- Reduce(rbind, snp_cluster_list)
  if (clust_prob_on) {
    nan_cluster_df$clust_prob <- nan_cluster_df$clust_dist %>% clust_prob_calc()
  }
  clust_re <- rbind(clust_re, nan_cluster_df)
  return(clust_re)
}

cluster_kmeans_min <- function(b_df,
                              nr = 10,
                              max_dist = 10.0,
                              space_typ = "regular",
                              clust_prob_on = TRUE,
                              norm_typ = "F", threshold = 1e-5, narm = TRUE) {
  #' Using the association scores for each SNP accross traits cluster the traits
  #' using kmeans. Return the cluster setup which minimises AIC.

  if (space_typ == "angle") {
    # For each point in b_df_comp convert the score to the angle
    # between the vectors to the origin and the unit vectors on the axis.
    b_df_clust <- mat_to_angle_mat(b_df)
  } else {
    b_df_clust <- b_df
  }

  # Filter complete cases
  b_df_comp <- b_df_clust[complete.cases(b_df_clust), ]

  # Initial cluster dataframe
  ic_list <- lapply(2:(nr + 1), kmeans_ic,
                    data_mat = b_df_comp,
                    clust_thres = threshold,
                    clust_norm = norm_typ,
                    na_rm = narm,
                    clust_prob_on = clust_prob_on)
  ic_df <- Reduce(rbind, ic_list)
  # Find the number of centres that minimizes the AIC
  min_cents <- ic_df$ncents[which.min(ic_df$aic)]
  # Use the number of centres to locate corresponding clusters since there
  # maybe variations due to machine precision in the aic values.
  clust_re_min_aic <- ic_df[which(ic_df$ncents == min_cents), ]
  # Assign terms with NaNs to nearest cluster
  nan_clusts <- lapply(setdiff(rownames(b_df_clust), rownames(b_df_comp)),
                       closest_clust,
                       b_df = b_df_clust,
                       clust_re = clust_re_min_aic,
                       max_dist = max_dist,
                       norm_typ = norm_typ)
  nan_clusts_df <- Reduce(rbind, nan_clusts)
  clust_out <- rbind(clust_re_min_aic, nan_clusts_df)
  if (clust_prob_on) {
    clust_out$clust_prob <- clust_out$clust_dist %>% clust_prob_calc()
  }
  return(clust_out)
}

mat_to_angle_mat <- function(mat) {
  nc <- ncol(mat)
  nr <- nrow(mat)
  unit_mat <- diag(nc)
  p_list <- lapply(1:nr, point_to_angles,
                  mat = mat,
                  unit_mat = unit_mat)
  out_mat <- Reduce(rbind, p_list)
  out_mat <- as.matrix(out_mat)
  row.names(out_mat) <- row.names(mat)
  colnames(out_mat) <- colnames(mat)
  return(out_mat)
}

make_unit_vec <- function(i, tot_n) {
  u <- integer(tot_n)
  u[i] <- 1
  return(u)
}
point_to_angles <- function(p, mat, unit_mat) {
  vec <- as.matrix(mat[p, ])
  nc <- ncol(mat)
  ang_list <- lapply(1:nc, point_to_angle,
                    unit_mat = unit_mat,
                    p = vec)
  ang_row <- Reduce(cbind, ang_list)
  return(ang_row)
}
point_to_angle <- function(col, unit_mat, p) {
  unit_vec <- as.matrix(unit_mat[col, ])
  dot_prod <- t(p) %*% unit_vec
  norm_x <- norm(p, type = "F")
  norm_y <- norm(unit_vec, type = "F")
  theta <- acos(dot_prod / (norm_x * norm_y))
  return(as.numeric(theta))
}