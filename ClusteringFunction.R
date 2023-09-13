### 3) Clustering
# Fitting K-Means clustering Model
find_ic <- function(fit) {
  # Find the BIC and AIC for the clusters. 
  m <- ncol(fit$centers)
  n <- length(fit$cluster)
  k <- nrow(fit$centers)
  d <- fit$tot.withinss
  return(data.frame(aic = d + 2 * m * k,
                    bic = d + log(n) * m * k))
}

kmeans_ic <- function(cents, data_mat, 
                      clust_thres = 1e-5, 
                      clust_norm = "F", 
                      na_rm = TRUE, 
                      clust_prob_on = TRUE){
  set.seed(240) # setting seed
  clust_re <- kmeans_skip_nan(data_mat, 
                              centers = cents,
                              clust_threshold = clust_thres,
                              norm_typ = "F", 
                              na_rm = FALSE, 
                              prob_on = clust_prob_on) 
  ic <- find_ic(clust_re)
  clust_re["aic"] <- ic$aic
  clust_re["bic"] <- ic$bic
  clust_re["ncents"] <- cents
  return(ic)
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
cluster_kmeans_basic <- function(b_df, nr, max_dist, clust_prob_on = TRUE,
norm_typ = "F", threshold = 1e-5, narm = TRUE) {
  #' Using the association scores for each SNP accross traits cluster the traits
  #' using kmeans. Return the cluster setup which minimises AIC.

  # Columns
  ax <- colnames(b_df)
  # Filter NaNs before clustering
  b_df_comp <- b_df[complete.cases(b_df),]
  # Initial cluster dataframe
  clust_re <- kmeans_skip_nan(b_df_comp,
                              centers = nr, iter.max = 300,
                              clust_threshold = threshold, norm_typ = norm_typ,
                              prob_on = clust_prob_on, na_rm = narm)
  # cluster number identification for each observation
  snp_cluster_list <- lapply(setdiff(names(b_df), names(b_df_comp)),
                            closest_clust, b_df = b_df, clust_re = clust_re,
                            max_dist = max_dist, norm_typ = norm_typ)
  nan_cluster_df <- Reduce(rbind, snp_cluster_list)
  if (clust_prob_on) {
    nan_cluster_df$clust_prob <- nan_cluster_df$clust_dist %>% clust_prob_calc()
  }
  clust_re <- rbind(clust_re, nan_cluster_df)
  return(clust_re)
}

cluster_kmeans_min <- function(b_df, nr = 10, max_dist = 10.0, clust_prob_on = TRUE,
                               norm_typ = "F", threshold = 1e-5, narm = TRUE) {
  #' Using the association scores for each SNP accross traits cluster the traits
  #' using kmeans. Return the cluster setup which minimises AIC.

  # Columns
  ax <- colnames(b_df)
  # Filter complete cases
  b_df_comp <- b_df[complete.cases(b_df),]

  # Initial cluster dataframe
  ic_list <- lapply(2:(nr+1), kmeans_ic, 
                    data_mat = b_df_comp,
                    clust_thres = threshold, 
                    clust_norm = norm_typ, 
                    na_rm = narm, 
                    clust_prob_on = clust_prob_on)
  ic_df <- Reduce(rbind,ic_list)
  # Find the number of centres that minimizes the AIC
  min_cents <- ic_df$ncents[which.min(ic_df$aic),]
  # Use the number of centres to locate corresponding clusters since there 
  # maybe variations due to machine precision in the aic values.
  clust_re_min_aic <- ic_df[which(ic_df$ncents == min_cents)] 
  # cluster number identification for each observation
  #FIXME dist now calculated and stored within clustering algorithm
  #snp_dist_list <- lapply(rownames(b_df_comp),
   ##                       clust_dist_calc,
    #                      clust_re = clust_re_min_aic,
     #                     b_df = b_df,
      #                    norm_typ = norm_typ)
  #clust_re_min_aic <- Reduce(rbind, snp_dist_list)
  # Assign terms with NaNs to nearest cluster
  nan_clusts <- lapply(setdiff(rownames(b_df), rownames(b_df_comp)),
                       closest_clust,
                       b_df = b_df, 
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

na_col_check <- function (b_col, percent = 0.95) {
  #' Check if the entire row is NaN
  n_accept <- length(b_col)*(1-percent)
  narows <- which(b_col %>% is.na())
  if (length(narows) > n_accept) {
    # All rows are NaN so trait will be removed from trait
    return(1)
  } else {
    return(0)
    }
}

test_all_na <- function(b_df, nan_col = "30600_irnt"){
  test <- na_col_check(b_df[, nan_col])
  print('test')
  if (test) {
    return(1)
  } else {
    return(0)
  }
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
clust_dist_calc <- function(snp, b_df, clust_re, norm_typ = "F"){
  #' Using the clustering result find the corresponding cluster centre and
  #' the distance to the centre from the point.
  clust_num <- clust_re$cluster[snp]
  clust_cent <- clust_re$centers
  snp_clust_dist <- clust_metric(clust_cent, b_df[snp, ], norm_typ)
  out <- data.frame(
    row.names = snp,
    clust_num = clust_num,
    clust_dist = snp_clust_dist,
    clust_prob = 1.0
  )
  return(out)
}

clust_prob_calc <- function(clust_df, b_df, dist_typ) {
  #' For each SNP calculate the distance to the centre of the cluster
  #' distance calculated based of distance type dist_typ

  # iterate through snp
  # find clust centre by clust number
  for (SNP_clust in clust_df$id) {
    point0 <- b_df[SNP_clust]
  }
  # find centre-snp distance
  # Invert distance for prob
}

trait_df_add_a <- function(trait_df, a, axes){
  #' Add the trait a to the dataframe of traits.
  #' Add the traits label, it's index in the trait list and it's index in the
  #' beta dataframe.
  if (a %in% trait_df$label) {
    return(trait_df)
  } else {
    trait_df <- trait_df %>% add_row(label = a, "axes_ind" = which(axes == a)[1])
    return(trait_df)
  }
}
testna <- function(unstd_beta_df, trait_axes) {
  aim_df <- data.frame(
    label = character(),
    axes_ind = integer(),
    b_df_ind = integer(),
    nSNPs = integer()
  )
  # Initialise with outcome
  aim_df <- aim_df_add_a(aim_df,trait_info$phenotype[125], trait_info$phenotype,
                         unstd_beta_df)
  b2_df <- remove_na_from_row(unstd_beta_df, aim_df)
  return(1)
}
aim_df <- data.frame(
  label = character(),
  axes_ind = integer(),
  b_df_ind = integer(),
  nSNPs = integer()
)
# Initialise with outcome
# aim_df <- aim_df_add_a(aim_df,OUT_pheno,trait_info$phenotype,
#                        unstdBeta_df)
# aim_df <- aim_df_add_a(aim_df,trait_info$phenotype[2],trait_info$phenotype,
#                        unstdBeta_df)
# nr=3
# test <- cluster_kmeans(unstdBeta_df,aim_df,nr)
