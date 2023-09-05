setup_dist <- function(score_df,norm_typ="F"){
  # Setup distances
  dist_df = data.frame(
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
                      snp1=snp1,
                      norm_typ = norm_typ)
  dist_df <- Reduce(rbind,dist_list)
  return(dist_df)
}

pair_dist_calc <- function(score_df, snp1, snp2, norm_typ = "F") {
  x1 <- score_df[rownames(score_df) == snp1]
  x2 <- score_df[rownames(score_df) == snp2]
  xp <- data.matrix(na.omit(x1-x2))
  d <- norm(xp,norm_typ)
  dist_df <- data.frame(snp1 = snp1,
                      snp2 = snp2,
                      dist=d)
  return(dist_df)
}