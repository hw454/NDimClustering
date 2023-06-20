clust_compare <-function(clusters_df,beta_df,axes,threshold){
  N<-length(cluster_df$cluster)
  c_scores<- data.frame(
    clust_num=1:Km,
    clust_size=0,
    NA_Count=0,
  )
  i=0
  for (a in axes){
    clust_scores=clust_score(cluster_df,beta_df)
    if (i==0){
      c_scores <- clust_scores
    }
    else{
      c_scores <- merge(c_scores,clust_scores, by="clust_num")
    }
    for (i in 1:N){
      cs1=c_scores[c_scores$clust_num==i,a]
      for (j in i:N){
        cs2=c_scores[c_scores$clust_num==j,a]
        if (abs(cs1-cs2)>threshold){
          print("Test based on outcome",a)
          return(1)
          }
      }
    }
  }
  print("None of the outcomes clustering met the thresholding test. ")
  return(0)
}



