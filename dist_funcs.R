#- Matrix distance functions
max_dist_calc <- function(score_mat, norm_typ = "F", na_rm = TRUE) {
  # Find the max distance based on range on each axis.
  max_p <- apply(score_mat, 2, max, na.rm = na_rm)
  min_p <- apply(score_mat, 2, min, na.rm = na_rm)
  max_dist <- norm(as.matrix(max_p - min_p), norm_typ)
return(max_dist)
}

setup_dist <- function(score_df, norm_typ = "F") {
  # Find the distance between all pairs of points.
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