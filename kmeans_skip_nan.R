kmeans_skip_nan <- function(eff_dfs, centers = nr, nstart = (nr+1), iter.max = 300,
                            clust_threshold=1e-5,norm_typ="F",na_rm=TRUE){
  set.seed(123)
  centroids=data.frame(cluster=1:centers
  )
  axes_nonnan=c()
  for(a in colnames(eff_dfs)){
    mina=min(eff_dfs[,a],na.rm=na_rm)
    maxa=max(eff_dfs[,a],na.rm=na_rm)
    if (!is.na(mina) & !is.na(maxa)){
      axes_nonnan<-axes_nonnan %>% append(a)
      centps=runif(centers,mina,maxa)
      centroids[a]=centps
    }
  }
  centroids<- centroids%>% column_to_rownames('cluster')
  max_dist=max(eff_dfs,na.rm=na_rm)-min(eff_dfs,na.rm=na_rm)
  # Randomly assign all SNPS a cluster number before iteration
  cluster_df=data.frame(
    snp_id=character(),
    clusters=integer(),
    clust_dist=numeric(),
    clust_prob=numeric()
  )
  #FIXME
  # Do the random assignment in dataframe setup not in loop
  for(snp in rownames(eff_dfs)){
    row=which(cluster_df$snp_id==snp)
    snp_score=eff_dfs[snp,]
    c_num=sample(1:centers,1)
    cent=centroids[c_num,]
    c_dist<-norm(data.matrix(na.omit(cent-snp_score)),norm_typ)
    c_prob<-clust_prob_calc(c_dist)
    cluster_df<-cluster_df %>% add_row(snp_id=snp,clusters=c_num,
                                       clust_dist=c_dist,clust_prob=c_prob)
  }
  # Repeat for less cluster numbers and choose the min result
  for (iter in 1:iter.max){
    #print('in iter')
    #print(cluster_df)
    # Initialise clusters for next iteration.
    for(clust_snp_id in rownames(eff_dfs)){
      snp_score=eff_dfs[clust_snp_id,]
      snp_dist=cluster_df[cluster_df$snp_id==clust_snp_id,'clust_dist']
      snp_clust_num=cluster_df[cluster_df$snp_id==clust_snp_id,'clusters']
      for(clust_num in 1:centers){
        cent=centroids[clust_num,]
        dist=norm(data.matrix(na.omit(cent-snp_score)),norm_typ)
        if(dist<snp_dist){
          snp_dist=dist
          snp_clust_num=clust_num
        }
      }
      row=which(cluster_df$snp_id==clust_snp_id)
      cluster_df[row,'clust_dist']<-snp_dist
      cluster_df[row,'clusters']<-snp_clust_num
      cluster_df[row,'clust_prob']<-clust_prob_calc(snp_dist)
      #print('After iter')
      #print(cluster_df)
    }
    snp_list<-rownames(eff_dfs)
    # Recompute the centroids based on the average of the clusters.
    centroidscheck<-centroids
    thresh_check=c()
    for (i in 1:length(centroids)){
      snp_list=cluster_df[cluster_df$clusters==i,'snp_id']
      snp_scores=eff_dfs[snp_list,]
      if (length(snp_list)){
        # Empty clusters not changed
        for(a in axes_nonnan){
          centroidscheck[i,a]=mean(eff_dfs[snp_list,a],na.rm=na_rm)
        }
        centroiddiff=data.matrix(na.omit(centroidscheck[i,axes_nonnan]-centroids[i,axes_nonnan]))
        centroidchange=norm(centroiddiff,norm_typ)
        if(centroidchange<clust_threshold & iter>1){
          thresh_check<-append(thresh_check,TRUE)
        }
        else{
          thresh_check<-append(thresh_check,FALSE)
        }
      }
    }
    if(all(thresh_check)){
      print('Clusters converged')
      print(iter)
      break
    }
    else{
      centroids=centroidscheck
    }
    # Are all the new centroids within the threshold of the previous
    
  }
 cluster_df<-column_to_rownames(cluster_df,'snp_id')
 #print('end clust')
 #print(cluster_df)
 return(cluster_df)
}

clust_prob_calc <-function(d){
  dist=1.0/(1.0+d)
  return(dist)
}