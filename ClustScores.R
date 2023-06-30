clust_score <- function(clusters_df,beta_df,aim_df){
  Km <- length(unique(clusters_df$cluster))
  SM <- length(clusters_df$clusters)
  num_axis=length(aim_df$label)
  for (i in unique(clusters_df$cluster)) {
    c_score0<- data.frame(
      id =character(),
      num_axis=integer(),
      clust_num  =integer(),
      clust_size =integer()
    )
    # Cluster id
    c_id <- paste0('na',num_axis,'_cn',i)
    # Get the associated SNPs with the cluster
    SNP_list=names(cluster_df$cluster[cluster_df$cluster==i])
    # ? How to deal with multiple SNP entries in data input.
    # Check the number of associated SNPS
    l1 <- length(SNP_list)
    c_score0 <- c_score0 %>% add_row(num_axis=num_axis,clust_num=i,
                                     clust_size=l1)
    for (ai in 1:num_axis){
      a=aim_df$label[ai]
      # Track the number of terms in each cluster.
      c_score<-axis_score(beta_df,SNP_list,aim_df)
      c_score0['id']<-c_id
      if (!(a %in% colnames(c_score0))){
        c_score0[a]=c_score
      }
      else{
        c_score0[c_score0$id==c_id,a]<- c_score0[c_score0$id==c_id,a]+c_score
      }
    }
    if (i==unique(clusters_df$cluster)[1]){
      clust_scores<-c_score0
    } else{ clust_scores <-rbind(c_score0,clust_scores) }
  }
  return(clust_scores)
}

axis_score <- function(beta_df,SNP_list,aim_df){
  # Get the column index for the trait.
  b_cols <-aim_df[aim_df$label==a,'b_df_ind']
  # Get the row indices for the traits in the cluster
  b_rows<-which(rownames(beta_df) %in% SNP_list)
  # Get the pathway colocization score from the input data for the SNP
  snp_assoc<-beta_df[b_rows,b_cols] # beta_df is a matrix with rows labelled by SNP_id and columns labelled by trait label.
  return(mean(snp_assoc))
}