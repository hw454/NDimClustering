clust_score <- function(clusters_df,beta_df,pval_df,aim_df,
                        bp_on=TRUE,clust_prob_on=TRUE){
  #' The dataframe of clusters is used alongside the beta and pvalue for each SNP
  #' to score the clusters on their association with each trait.
  #' Information for the traits is stored in aim_df
  #Km <- length(unique(clusters_df$cluster))
  #SM <- length(clusters_df$clusters)
  traits <- aim_df$label
  num_axis<-length(traits)
  clust_nums=unique(clusters_df$cluster)
  clust_scores<- data.frame(
    id =character(),
    num_axis=integer(),
    clust_num  =integer(),
    clust_size =integer()
  )
  for (i in clust_nums) {
    # Store the cluster scores for each cluster number in c_score0
    c_score0<- data.frame(
      id =character(),
      num_axis=integer(),
      clust_num  =integer(),
      clust_size =integer()
    )
    # Cluster id
    c_id <- paste0('na',num_axis,'_cn',i)
    # Get the associated SNPs with the cluster
    #print('inside scores')
    #print(clusters_df)
    SNP_list=rownames(clusters_df[clusters_df$cluster==i,])
    # ? How to deal with multiple SNP entries in data input.
    # Check the number of associated SNPS
    l1 <- length(SNP_list)
    c_score0 <- c_score0 %>% add_row(num_axis=num_axis,clust_num=i,
                                     clust_size=l1)
    c_score0['id']<-c_id
    for (a in traits){
      # Add the cluster score for each axis to c_score0
      #a=aim_df$label[ai]
      # Track the number of terms in each cluster.
      if (clust_prob_on){
        clust_probs=subset(subset(clusters_df,rownames(clusters_df) %in% SNP_list),select='clust_prob')
      }
      else{
        clust_probs=subset(subset(clusters_df,rownames(clusters_df) %in% SNP_list),select='clust_prob')
        clust_probs[SNP_list,'clust_prob']=1.0
      }
      ax_score<-axis_score(beta_df,pval_df,SNP_list,aim_df,a,clust_probs,bp_on)
      # Ensure that a is a column 
      if (!(a %in% colnames(clust_scores))){
        clust_scores[a]=numeric()
      }
      # Ensure that a is a column in c_score0
      if (!(a %in% colnames(c_score0))){
        c_score0[a]=numeric()
        c_score0[a]=ax_score
      }
      else{
        c_score0[c_score0$id==c_id,a]<- c_score0[c_score0$id==c_id,a]+ax_score
      }
    }
    # Add the cluster scores for that cluster number to the dataframe for them all
    # The trait columns are added within the loop if they are not already there 
    # (i.e the first time they are checked)
    clust_scores <-rbind(c_score0,clust_scores) 
  }
  clust_scores <- clust_scores %>% column_to_rownames('id')
  return(clust_scores)
}

axis_score <- function(beta_df,pval_df,SNP_list,aim_df,a,clust_probs,bp_on=TRUE){
  # Get the column index for the trait.
  b_col_names <-aim_df$label==a
  # Get the pathway colocization score from the input data for the SNP
  b_sub=beta_df[SNP_list,b_col_names]
  clust_prob<-clust_probs[SNP_list,'clust_prob']
  if (bp_on){
    p_sub=pval_df[SNP_list,b_col_names]
    snp_assoc<-abs(b_sub*p_sub*clust_prob)# beta_df is a matrix with rows labelled by SNP_id and columns labelled by trait label.
    total_probs<-sum(p_sub*clust_prob) # If beta probas are used then averaging needs to be over total p-values
  }
  else{
    snp_assoc<-abs(b_sub*clust_prob) # beta_df is a matrix with rows labelled by SNP_id and columns labelled by trait label.
    total_probs<-sum(clust_prob) # If beta probs aren't used than use the total umber of entries on axis
  }
  axis_snp_assoc<-snp_assoc/total_probs
  return(mean(snp_assoc))
}