clust_score <- function(clusters_df,beta_df,pval_df,aim_df,
                        bp_on=TRUE,clust_prob_on=TRUE){
  #' The dataframe of clusters is used alongside the beta and pvalue for each SNP
  #' to score the clusters on their association with each trait.
  #' Information for the traits is stored in aim_df
  #Km <- length(unique(clusters_df$cluster))
  #SM <- length(clusters_df$clusters)
  clust_nums=unique(clusters_df$cluster)
  clust_scores_list<-lapply(clust_nums,score_cluster,
                            beta_df=beta_df,clusters_df=clusters_df,
                            pval_df=pval_df,aim_df=aim_df,
                            bp_on=bp_on)
  clust_scores<-Reduce(rbind,clust_scores_list)
  return(clust_scores)
}

score_cluster <- function(c_num,beta_df,clusters_df,num_axis,
                          pval_df,aim_df,bp_on){
  # Get the list of traits and the number of them.
  traits <- aim_df$label
  num_axis<-length(traits)
  # Cluster id
  c_id <- paste0('na',num_axis,'_cn',c_num)
  # Get the associated SNPs with the cluster
  # Check the number of associated SNPS
  SNP_list=rownames(clusters_df[clusters_df$cluster==c_num,])
  l1 <- length(SNP_list)
  clust_probs<-clusters_df[SNP_list,'clust_prob']
  # Create a dataframe for each trait score for the cluster then combine.
  trait_score_df_list<-lapply(traits,axis_score,
                              c_id=c_id,SNP_list,b_df=beta_df,
                              pval_df=pval_df,aim_df=aim_df,
                              clust_probs=clust_probs,bp_on=bp_on)
  c_score0=Reduce(cbind,trait_score_df_list)
  c_score0['clust_num']<-c_num
  return(c_score0)
}

axis_score <- function(a,c_id,SNP_list,b_df,pval_df,aim_df,clust_probs,bp_on=TRUE){
  # Get the column index for the trait.
  b_col_name <-aim_df$label==a
  # Get the pathway colocization score from the input data for the SNP
  b_sub=na.omit(b_df[SNP_list,b_col_name])
  clust_probs<-na.omit(clust_probs)
  if (bp_on){
    p_sub=na.omit(pval_df[SNP_list,b_col_name])
    snp_assoc<-abs(b_sub*p_sub*clust_probs) # beta_df is a matrix with rows labelled by SNP_id and columns labelled by trait label.
    total_probs<-sum(p_sub*clust_probs) # If beta probas are used then averaging needs to be over total p-values
  }
  else{
    snp_assoc<-abs(b_sub*clust_probs) # beta_df is a matrix with rows labelled by SNP_id and columns labelled by trait label.
    total_probs<-sum(clust_probs) # If beta probs aren't used than use the total umber of entries on axis
  }
  axis_snp_assoc<-sum(snp_assoc)/total_probs
  trait_score_df=data.frame(
     row.names=c_id
  )
  trait_score_df[a]=axis_snp_assoc
  return(trait_score_df)
}