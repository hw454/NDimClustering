kmeans_skip_nan <- function(b_mat, centers = nr,
                            iter_max = 300, clust_threshold = 1e-5,
                            norm_typ = "F", na_rm = FALSE, prob_on = TRUE) {
  set.seed(123)
  snp_list <- rownames(b_mat)
  # Generate data frame with max and min data.
  # Min and Max are columns. Rows for each column of b_mat
  min_max_df <- data.frame(
    row.names = colnames(b_mat),
    min = apply(b_mat, 2, min, na.rm = na_rm),
    max = apply(b_mat, 2, max, na.rm = na_rm)
  )
  min_max_df <- min_max_df %>% na.omit()
  # Randomly assign central coords per cluster.
  centroid_list <- lapply(rownames(min_max_df), rand_cent,
                          n_cents = centers, min_max_df = min_max_df)
  centroids <- Reduce(cbind, centroid_list)
  # Randomly assign cluster to snps

  nsnps <- length(snp_list)
  clust_samp <- replicate(nsnps, sample(1:centers, 1))
  cluster_df <- data.frame(
    row.names = snp_list,
    clust_num = clust_samp
  )
  cluster_df["clust_prob"] <- numeric()
  clust_dist_list <- lapply(snp_list, cent_dist_calc,
                          b_mat = b_mat, cluster_df = cluster_df,
                          centroids = centroids, norm_typ = norm_typ)
  clust_dist_df <- Reduce(rbind, clust_dist_list)
  cluster_df <- cbind(cluster_df, clust_dist_df)
  for (iter in 1:iter_max){
    # For each SNP find the cluster with the closest centre.
    snp_clust_list <- lapply(snp_list, snp_closest_clust,
                              b_mat = b_mat, cluster_df = cluster_df,
                             centroids = centroids, norm_typ = norm_typ)
    # Combine the list of dataframes into one dataframe.
    # Override Cluster_df with the new assignment
    cluster_df <- Reduce(rbind, snp_clust_list)
    # Recompute the centroids based on the average of the clusters.
    # Check if the previous centres differ from the cluster means.
    #print(clust_threshold)
    thresh_list <- lapply(rownames(centroids), clust_cent_check,
                        iter = iter, cluster_df = cluster_df,
                        b_dfs = b_mat, centroids = centroids,
                        na_rm = na_rm, norm_typ = norm_typ,
                        clust_threshold = clust_threshold)#,
                        #axes_nonnan=axes_nonnan)
    thresh_check_df <- Reduce(rbind, thresh_list)
    # Are all the new centroids within the threshold of the previous
    if (all(thresh_check_df$thresh_check)) {
      print(paste("Clusters converged", iter))
      break
    }
    centroids <- thresh_check_df[,
                  !names(thresh_check_df) %in% c("thresh_check")]
    # Are all the new centroids within the threshold of the previous
  }
  if (prob_on) {
    cluster_df$clust_prob <- cluster_df$clust_dist %>% clust_prob_calc()
  }
 return(cluster_df)
}

snp_closest_clust <- function(snp_id, b_mat, cluster_df,
                              centroids, norm_typ) {
  #' For each cluster check whether the centre is closer than the currently
  #' assigned cluster
  snp_score <- b_mat[snp_id, ]
  snp_dist <- cluster_df[snp_id, "clust_dist"]
  snp_clust_num <- cluster_df[snp_id, "clust_num"]
  for (c_num in rownames(centroids)) {
    dist <- norm(data.matrix(
              na.omit(centroids[c_num, ] - snp_score)), norm_typ)
     if (dist < snp_dist) {
       snp_dist <- dist
        snp_clust_num <- c_num
     }
   }
   snp_cluster_df <- data.frame(
     row.names = snp_id,
     clust_dist = snp_dist,
     clust_num = snp_clust_num,
     clust_prob = 1.0
   )
  return(snp_cluster_df)
}
cent_dist_calc <- function(snp_id, b_mat, cluster_df, centroids, norm_typ) {
   snp_score <- b_mat[snp_id, ]
   c_num <- cluster_df[snp_id, "clust_num"]
   cent <- centroids[c_num, ]
   c_dist <- norm(data.matrix(na.omit(cent - snp_score)), norm_typ)
   clust_df <- data.frame(row.names = snp_id,
                       "clust_dist" = c_dist)
  return(clust_df)
}

clust_prob_calc <- function(d) {
  dist <- 1.0 / (1.0 + d)
  return(dist)
}

clust_cent_check <- function(c_num, iter, cluster_df, b_dfs, centroids,
                            na_rm = TRUE, norm_typ = "F",
                            clust_threshold = 1e-5) {
  sub_snp_list <- which(cluster_df$clust_num == c_num)
  snp_scores <- b_dfs[sub_snp_list, ]
  nterms <- length(sub_snp_list)
  centroidscheck <- centroids # Initialise the centroids checking dataframe
  # Empty clusters not changed
  if (nterms) {
    #' Compute the new centers based on the mean of the snps in the clusters
    #' If there's only one term then don't use mean
    #' else: Mean for each column gives the value for the centre on each axis
    if (nterms == 1) {
      centroidscheck[c_num, ] <- snp_scores
    } else {
      centroidscheck[c_num, ] <- colMeans(snp_scores, na.rm = na_rm)
    }
    # Calculate how much the centroid has moved.
    centroiddiff <- data.matrix(
                                na.omit(centroidscheck[c_num, ]
                                - centroids[c_num, ]))
    centroidchange <- norm(centroiddiff, norm_typ)
    # Check if the diff between new and old centres is below the threshold
  } else {
    centroidchange <- 0
  }
  centroid_df <- data.frame(
    row.names = c_num,
    thresh_check = (centroidchange < clust_threshold && iter > 1)
  )
  if (centroid_df[c_num, "thresh_check"]) {
    centroid_df <- cbind(centroid_df, centroids[c_num, ])
  } else {
     centroid_df <- cbind(centroid_df, centroidscheck[c_num, ])
  }
  return(centroid_df)
 }

rand_cent <- function(a, min_max_df, n_cents) {
  cent <- runif(n_cents, min_max_df[a, "min"], min_max_df[a, "max"])
  out_cent <- data.frame(a = cent)
  colnames(out_cent) <- c(a)
  return(out_cent)
}