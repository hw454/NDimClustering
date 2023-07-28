setup_dist <- function(score_df,norm_typ="F"){
  # Setup distances
  dist_df=data.frame(
    snp1=character(),
    snp2=character(),
    dist=numeric()
  )
  for (snp1 in rownames(score_df)){
    for (snp2 in rownames(score_df)){
      x1=unstdBeta_df[rownames(unstdBeta_df)==snp1]
      x2=unstdBeta_df[rownames(unstdBeta_df)==snp2]
      xp=data.matrix(na.omit(x1-x2))
      d=norm(xp,norm_typ)
      dist_df<-dist_df %>% add_row(snp1=snp1,snp2=snp2,
                           dist=d)
    }
  }
return(dist_df)
}