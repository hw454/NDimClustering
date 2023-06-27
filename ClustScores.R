clust_score <- function(clusters_df,beta_df,aim,aim_ind){
  Km <- length(unique(clusters_df$cluster))
  SM <- length(clusters_df$clusters)
  clust_scores <- data.frame(
    clust_num  =unique(clusters_df$cluster),
    clust_size =0,
    NA_Count   =0
  )
  for (ai in 1:length(aim)){
    a=aim[ai]
    a_ind=aim_ind[ai]
    clust_scores[a]=0.0
    for (i in unique(clusters_df$cluster)) {
      # Extract the data from the clustering for each snp
      # No probability output for prob in cluster 
      # psnp <- row$probability # Probability of snp in cluster nc
      # Get the pathway colocization score from the input data for the SNP
      snp_assoc<-beta_df[ i,a_ind] # beta_df is a matrix with rows labelled by SNP_id and columns labelled by trait label.
      # ? How to deal with multiple SNP entries in data input.
      # Retrieve non-NAN values. 
      na_l <- length(snp_assoc[is.na(snp_assoc)])
      snp_path1 <- snp_assoc[!is.na(snp_assoc)]
      # Check length of remaining values-> n=1, store, -> n>1, store average, -> n=0 pass.
      l1 <- length(snp_assoc)
      # Track the number of terms in each cluster.
      clust_scores[clust_scores$clust_num==i,'clust_size'] <- clust_scores[clust_scores$clust_num==i,'clust_size']+l1
      clust_scores[clust_scores$clust_num==i,'NA_Count']   <- clust_scores[clust_scores$clust_num==i,'NA_Count']+na_l
    
      if (l1==0 ){# Do nothing and pass to the next snp
      }
      else if(l1==1){
        # Add the association score to the cluster pathway score. Weight by the probability.
        clust_scores[clust_scores$clust_num==i,a] <- clust_scores[ clust_scores$clust_num==i,a]+snp_assoc
      }
      else{
        # Average remaining values
        # Add the pathway1 coloc score to the cluster pathway score. Weight by the probability.
        clust_scores[clust_scores$clust_num==i,a]<-clust_scores[clust_scores$clust_num==i,a]+sum(snp_assoc)
      }
    }
  }
  return(clust_scores)
}
