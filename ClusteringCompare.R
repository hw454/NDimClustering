clust_compare <-function(unstdBeta_df,unstdSE_df,pval_df,tstat_df,
                         axes,threshold,thresh_norm,clust_threshold,clus_norm,
                         which_clust='basic',bp_on=TRUE,clust_prob_on=TRUE,
                         max_dist=0.0,nr=5){
  #' Iterate through the columns in axes and clusters the data. 
  #' If there is a distinct difference between two clusters exit.
  
  # Data frame for recording the cluster scores.
  c_scores<- data.frame(
    id=character(),
    num_axis=integer(),
    clust_num =integer(),
    clust_size=integer(),
    clust_prob=numeric(),
    NA_Count=integer()
  )
  # Initialise with outcome
  aim_df <- data.frame(label=OUT_pheno,
                       axes_ind=which(trait_info$phenotype==OUT_pheno)[1],
                       b_df_ind= which(colnames(unstdBeta_df)== OUT_pheno)[1],
                       nSNPs=sum(!is.na(unstdBeta_df[,OUT_pheno]))[1]
                       )
  c_scores[OUT_pheno]=numeric()
  for (ai in 1:length(trait_info$phenotype)){
    # Update the traits forming the axis
    a=trait_info$phenotype[ai]
    # Add the trait to the trait dataframe
    covered <- (a %in% aim_df$label)
    if (!covered){
      #print(covered)
      # If the axis is already included then add row is not run.
      aim_df <- aim_df_add_a(aim_df,a,trait_info$phenotype,
                 unstdBeta_df)
      # Cluster the data on these axes
      #unstdBeta_df <- remove_na_from_row(unstdBeta_df,aim_df)
      allna <- all_na_check(unstdBeta_df,aim_df)
      if (allna){
        #' If the trait column was removed from the beta_df during na removal then 
        #' remove the trait from the trait data frame and move to the next trait.
        aim_df <- aim_df[aim_df$label!=a,]
      }
      else{
        # If the trait is not all NaN then run clustering.
        print(paste0('Cluster on ',length(aim_df$label),' axes'))
        print(paste0('New axis ',a))
        # Cluster the data on these axes
        if (which_clust=='min'){cluster_df=cluster_kmeans_min(unstdBeta_df,aim_df,nr,max_dist)}
        else{cluster_df=cluster_kmeans_basic(unstdBeta_df,aim_df,nr,max_dist,clust_norm,clust_threshold)}
        # Find the set of cluster numbers
        c_nums<-unique(cluster_df$cluster) 
        # Score the clustered data based on affilliations with axes.
        # Find the score for this axis
        c_score0 <- clust_score(cluster_df,unstdBeta_df,pval_df,aim_df,
                                bp_on,clust_prob_on) # Score across all axis
        # Initialise new column for new axis
        # Create column for merge and fill.
        c_scores[a]=numeric()
        # Iterate through each cluster and compare across the others to find if 
        # any pair have a distinct difference.
        N2=length(c_nums)
        N1=N2-1
        #FIXME find all diffs then take max
        diff_scores_list <- lapply(c_nums,compare_oneclust_tolist,
                                   c_nums=c_nums,
                                   c_score0=c_score0,
                                   aim_df=aim_df)
        diff_score_list<-diff_score_list[!sapply(diff_score_list,is.null)]
        diff_scores<-Reduce(rbind,diff_score_list)
        row=which(diff_scores,)
      }
      c_score0['num_axis']=length(aim_df$label)
      c_scores <- rbind(c_scores,c_score0)   
    }
    
  }
  print("None of the outcomes clustering met the thresholding test. ")
  return(c_scores)
}

compare_oneclust_tolist<-function(cn1,c_nums,c_score0,aim_df){
  diff_score_list<-lapply(c_nums,clust_score_diff,
                          cn1=cn1,
                          c_score0=c_score0,
                          aim_df)
  diff_score_list<-diff_score_list[!sapply(diff_score_list,is.null)]
  diff_scores<-Reduce(rbind,diff_score_list)
  return(diff_scores)
}

clust_score_diff <- function(cn1,cn2,c_score0,aim_df){
  #' Find the score difference between two clusters
  #' If NaN return Null
  cs1<-as.matrix(c_score0[c_score0$clust_num==cn1,aim_df$label])
  cs2<-as.matrix(c_score0[c_score0$clust_num==cn2,aim_df$label])
  # IF there are no members to a cluster it's score will be NaN
  if (!all(is.na(cs1)) && !all(is.na(cs2)) && nrow(cs1) && nrow(cs2)==0){
    metric_score=clust_metric(cs1,cs2,thresh_norm)
    if (is.na(metric_score)){return()}
    else{
      score_diff=data.frame(
        c_num1=cn1,
        c_num2=cn2,
        diff=metric_score
      )
    }
  }
  else{return()}
}
#test2 <- testna(unstdBeta_df,trait_info$phenotype)
# test<-clust_compare(stdBeta_df ,SNP_ind,axes,threshold,thresh_norm, clust_norm)
#test<-clust_compare(unstdBeta_df,unstdSE_df,pval_df,tstat_df,
# trait_axes,threshold,thresh_norm,clust_norm)