clust_score <- function(clusters_df,beta_df,b_lab){
  Km <- length(unique(clusters_df$clusters))-1
  SM <- length(clusters_df$clusters)
  clust_scores <- data.frame(
    clust_num=1:Km,
    clust_size=0,
    NA_Count=0,
    b_lab=0.0
  )
  for (i in 1:SM) {
    # Extract the data from the clustering for each snp
    nc  <- clusters_df$cluster[clusters_df$cluster==i] # Cluster Number
    # No probability output for prob in cluster 
    # psnp <- row$probability # Probability of snp in cluster nc
    # Get the pathway colocization score from the input data for the SNP
    snp_assoc<-beta_df[ beta_df$SNP==i,b_lab]
    # ? How to deal with multiple SNP entries in data input.
    # Retrieve non-NAN values. 
    na_l <- length(snp_assoc[is.na(snp_assoc)])
    snp_path1 <- snp_assoc[!is.na(snp_assoc)]
    # Check length of remaining values-> n=1, store, -> n>1, store average, -> n=0 pass.
    l1 <- length(snp_assoc)
    # Track the number of terms in each cluster.
    clust_scores[clust_scores$clust_num==nc,'clust_size'] <- clust_scores[clust_scores$clust_num==nc,'clust_size']+l1
    clust_scores[clust_scores$clust_num==nc,'NA_Count']   <- clust_scores[clust_scores$clust_num==nc,'NA_Count']+na_l
    
    if (l1==0 ){# Do nothing and pass to the next snp
    }
    else if(l1==1){
      # Add the association score to the cluster pathway score. Weight by the probability.
      clust_scores[clust_scores$clust_num==nc,b_lab] <- clust_scores[ clust_scores$clust_num==nc,b_lab]+snp_assoc
    }
    else{
      # Average remaining values
      # Add the pathway1 coloc score to the cluster pathway score. Weight by the probability.
      clust_scores[clust_scores$clust_num==nc,b_lab]<-clust_scores[clust_scores$clust_num==nc,p1_lab]+sum(snp_assoc)}
  }
  return(clust_scores)
}