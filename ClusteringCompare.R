clust_compare <-function(unstdBeta_df,unstdSE_df,pval_df,tstat_df,
                         axes,threshold,thresh_norm,clus_norm){
  #' Iterate through the columns in axes and clusters the data. 
  #' If there is a distinct difference between two clusters exit.
  
  aim_df <- data.frame(
    label =character(),
    axes_ind = integer(),
    b_df_ind = integer(),
    nSNPs = integer()
  )
  # Data frame for recording the cluster scores.
  c_scores<- data.frame(
    id=character(),
    num_axis=integer(),
    clust_num =integer(),
    clust_size=integer(),
    NA_Count=integer()
  )
  # Initialise with outcome
  aim_df <- aim_df_add_a(aim_df,OUT_pheno,trait_info$phenotype,
               unstdBeta_df)
  #FIXME
  for (ai in 1:length(trait_info$phenotype)){
    # Update the traits forming the axis
    a=trait_info$phenotype[ai]
    # Add the trait to the trait dataframe
    aim_df <- aim_df_add_a(aim_df,a,trait_info$phenotype,
                 unstdBeta_df)
    nr=10.0#/length(aim_df$label) # FIXME - Max number of clusters.
    # Cluster the data on these axes
    unstdBeta_df <- remove_na_from_row(unstdBeta_df,aim_df)
    allna <- all_na_check(unstdBeta_df,aim_df)
    if (allna){
      #' If the trait column was removed from the beta_df during na removal then 
      #' remove the trait from the trait data frame and move to the next trait.
      aim_df <- aim_df[aim_df$label!=a,]
    }
    else{
      print(paste0('Cluster on ',length(aim_df$label),' axes'))
      print(paste0('New axis ',a))
      # Cluster the data on these axes
      cluster_df=cluster_kmeans(unstdBeta_df,aim_df,nr) 
      # Find the set of cluster numbers
      c_nums<-unique(cluster_df$cluster) 
      # Score the clustered data based on affilliations with axes.
      # Find the score for this axis
      c_score0 <- clust_score(cluster_df,
                            unstdBeta_df,aim_df) # Score across all axis
      # Initialise new column for new axis
      c_scores[a]=numeric()
      # Iterate through each cluster and compare across the others to find if 
      # any pair have a distinct difference.
      N2=length(c_nums)
      N1=N2-1
      for (i in 1:N1){
        for (j in (i+1):N2){
          ci=c_nums[i]
          cj=c_nums[j]
          cs1<-as.matrix(c_score0[c_score0$clust_num==i,aim_df$label])
          cs2<-as.matrix(c_score0[c_score0$clust_num==j,aim_df$label])
          if (clust_metric(cs1,cs2,thresh_norm)>(threshold)){
            print("Threshold met on outcome")
            print(a)
            return(c_scores)
          }
        }
      }
      c_scores <- rbind(c_scores,c_score0)
    }
  }
  print("None of the outcomes clustering met the thresholding test. ")
  return(c_scores)
}
#test2 <- testna(unstdBeta_df,trait_info$phenotype)
# test<-clust_compare(stdBeta_df ,SNP_ind,axes,threshold,thresh_norm, clust_norm)
#test<-clust_compare(unstdBeta_df,unstdSE_df,pval_df,tstat_df,
# trait_axes,threshold,thresh_norm,clust_norm)