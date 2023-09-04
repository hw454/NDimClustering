### 3) Clustering
# Fitting K-Means clustering Model
kmeans_ic <- function(fit) {
  #https://stackoverflow.com/questions/15839774/how-to-calculate-bic-for-k-means-clustering-in-r # nolint
  m <- ncol(fit$centers)
  n <- length(fit$cluster)
  k <- nrow(fit$centers)
  D <- fit$tot.withinss
  return(data.frame(aic = D + 2 * m * k,
                    bic = D + log(n) * m * k))
}
closest_clust <- function(snp_id, b_df, clust_re, max_dist, norm_typ = "F") {
  #' Find the cluster whose centre is closest to the point with id snp_id.
  #' Used for points with NaNs that are not included in the
  #' initial clustering algorithm
  snp_dist <- max_dist
  for (cent in clust_re$centers){
    snp_clust_dist <- norm(cent - b_df[snp_id, ], norm_typ)
    clust_re$clust_dist[snp_id] <- snp_clust_dist
    if (snp_clust_dist < snp_dist) {
      # Assign cluster number when distance is less than previous
      snp_dist <- snp_clust_dist
      clust_num <- which(clust_re$centers == cent)
    }
  }
  out_df <- data.frame(
    row.names = snp_id,
    clust_num = clust_num,
    clust_dist = clust_dist,
    clust_prob = 1.0
  )
  return(out_df)
}
cluster_kmeans_basic <- function(b_df, nr, max_dist, clust_prob_on = TRUE,
norm_typ = "F", threshold = 1e-5) {
  #' Using the association scores for each SNP accross traits cluster the traits
  #' using kmeans. Return the cluster setup which minimises AIC.

  # Columns
  ax <- colnames(b_df)
  # Filter NaNs before clustering
  b_df_comp <- b_df[complete.cases(b_df)]
  # Initial cluster dataframe
  set.seed(240) # setting seed
  clust_re <- kmeans_skip_nan(b_df_comp,
                              centers = nr, nstart = (nr+1), iter.max = 300,
                              threshold = threshold, norm_typ = norm_typ)
  # cluster number identification for each observation
  snp_cluster_list <- lapply(setdiff(names(beta_df), names(b_df_comp)),
                            closest_clust, b_df = b_df, clust_re = clust_re,
                            max_dist = max_dist, norm_typ = norm_typ)
  nan_cluster_df <- Reduce(rbind, snp_cluster_list)
  if (clust_prob_on){
    nan_cluster_df$clust_prob <- nan_cluster_df$clust_dist %>% clust_prob_calc()
  }
  clust_re <- clust_re(rbind, nan_cluster_df)
  return(clust_re)
}
cluster_kmeans_min <- function(beta_df, nr, max_dist, clust_prob_on = TRUE
norm_typ = "F", threshold = 1e-5) {
  #' Using the association scores for each SNP accross traits cluster the traits
  #' using kmeans. Return the cluster setup which minimises AIC.

  # Columns
  ax <- colnames(b_df)
  # Filter complete cases
  b_df_comp <- b_df[complete.cases(b_df)]

  # Initial cluster dataframe
  cluster_list <- list()
  ic_df <- data.frame(matrix(data = NA, nrow = nr, ncol = 1))
  colnames(ic_df) <- c("aic")
  ic_df$nCluster <- 2:(nr + 1)
  for (i in 2:(nr + 1)) {
    set.seed(240) # setting seed
    clust_re <- kmeans(eff_dfs,
                       centers = i, nstart = (nr + 1), iter.max = 300)
    cluster_list[[length(cluster_list) + 1]] <- clust_re
    ic <- kmeans_ic(clust_re)
    ic_df[(i - 1), 1] <- ic$aic
  }
  # cluster number identification for each observation
  clust_re_min_aic <- cluster_list[[which(ic_df$aic == min(ic_df$aic))]]
  # cluster number identification for each observation
  for (snp in names(clust_re_min_aic$cluster)){
    clust_num <- clust_re_min_aic$cluster[snp]
    clust_cent <- clust_re_min_aic$centers
    snp_clust_dist <- norm(clust_cent - b_df[snp, ], "F")
    clust_re_min_aic$clust_dist[snp] <- snp_clust_dist
  }
  # Assign NaNs to nearest cluster
  nan_clusts <- lapply(setdiff(names(b_df),names(b_df_comp)),
                      b_df = b_df, clust_re = clust_re_min_aic,
                      max_dist = max_dist, norm_typ = norm_typ)
  nan_clusts_df <- Reduce(rbind,nan_clusts)
  names(clust_re_min_aic$clust_dist) <- names(clust_re_min_aic$cluster)
  clust_out <- data.frame(
    row.names = names(clust_re_min_aic$cluster),
    clust_num = clust_re_min_aic$cluster,
    clust_dist = clust_re_min_aic$clust_dist,
    clust_prob = clust_re_min_aic$clust_prob
  )
  clust_out <- rbind(clust_out, nan_clusts_df)
  if (clust_prob_on){
    clust_out$clust_prob <- clust_out$clust_dist %>% clust_prob_calc()
  }
  return(clust_out)
}

all_na_check <- function (b_df) {
  #' Check if the entire row is NaN
  narows <- which(b_df[, dim(b_df)[2]] %>% is.na())
  if (length(narows) == dim(b_df)[1]) {
    # All rows are NaN so trait will be removed from trait
    return(1)
  } else {
    return(0)
    }
}

clust_metric <- function (cs1, cs2, norm_typ) {
  if (length(cs1) <= 1) {
    return(abs(cs1 - cs2))
  } else {
    clustid_1 <- which(!is.na(cs1))
    clustid_2 <- which(!is.na(cs2))
    clust_ids <- intersect(clustid_1, clustid_2)
    if (length(clust_ids) == 0) {
      return( NaN)
      } else if (length(clust_ids) == 1){
        return(abs(cs1[clust_ids] - cs2[clust_ids]))
        } else {
      return(norm(as.matrix(cs1[clust_ids] - cs2[clust_ids]), norm_typ))
    }
  }
}
get_dist <- function(a, center, dist_typ) {
  if (dist_typ == "p_euclid") {
    # calculate the euclidean distance to the centre point
    dist <- norm(a - b, "F")
  } else if (dist_typ=="line_euclid"){
    # Calculate the distance to the line given by center
    dist <- 1
  } else {
    dist <- 1
    }
  return(dist)
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

aim_df_add_a <- function(aim_df, a, axes, b_df){
  #' Add the trait a to the dataframe of traits.
  #' Add the traits label, it's index in the trait list and it's index in the 
  #' beta dataframe.
  if (a %in% aim_df$label) {
    return(aim_df)
    } else {
    aim_df <- aim_df %>% add_row(label= a,"axes_ind" = which(axes == a)[1],
    "b_df_ind" = 0, "nSNPs" = 0)
    aim_df[aim_df$label == a, "b_df_ind"] <- which(colnames(b_df) == a)[1]
    aim_df[aim_df$label == a, "nSNPs"] <- sum(!is.na(b_df[, a]))[1]
    return(aim_df)}
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

#nr=10
#trait_axes_ind=which(colnames(stdBeta_df_noEXP) %in% trait_axes)
#load(paste0(res_dir,"QCdata_",EXP_pheno,"ClusterIter",1,".Rdata"))
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
