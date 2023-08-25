snp_closest_clust<-function(snp_id,eff_dfs,cluster_df,centroids,norm_typ){
  # For each cluster check whether the centre is closer than the currently assigned cluster.
  snp_score=eff_dfs[snp_id,]
  snp_dist=cluster_df[snp_id,'clust_dist']
  snp_clust_num=cluster_df[snp_id,'clusters']
  for(clust_num in rownames(centroids)){
    dist=norm(data.matrix(na.omit(centroids[clust_num,]-snp_score)),norm_typ)
    if(dist<snp_dist){
      snp_dist=dist
      snp_clust_num=clust_num
    }
  }
  snp_cluster_df=data.frame(
    row.names = snp_id,
    clust_dist=snp_dist,
    clusters=snp_clust_num,
    clust_prob=1.0
  )
  return(snp_cluster_df)
}

cent_dist_calc <-function(snp_id,eff_dfs,cluster_df,centroids,norm_typ){
  snp_score=eff_dfs[snp_id,]
  c_num=cluster_df[snp_id,'clusters']
  cent=centroids[c_num,]
  c_dist<-norm(data.matrix(na.omit(cent-snp_score)),norm_typ)
  clust_df=data.frame(row.names = snp_id,
                      'clust_dist'=c_dist)
  return(clust_df)
}
clust_prob_calc <-function(d){
  dist=1.0/(1.0+d)
  return(dist)
}

clust_cent_check<-function(c_num,iter,cluster_df,eff_dfs,centroids,na_rm=TRUE,norm_typ="F"){
  sub_snp_list=which(cluster_df$clusters==c_num)
  snp_scores=eff_dfs[sub_snp_list,]
  nterms=length(sub_snp_list)
  centroidscheck<-centroids # Initialise the centroids checking dataframe
  # Empty clusters not changed
  if (nterms){
    # Compute the new centers based on the mean of the snps in the clusters
    # If there's only one term then don't use mean
    # else: Mean for each column gives the value for the centre on each axis
    if (nterms==1){centroidscheck[c_num,]<-snp_scores} 
    else{centroidscheck[c_num,]<-colMeans(snp_scores,na.rm=na_rm)}
    # Calculate how much the centroid has moved. 
    centroiddiff=data.matrix(na.omit(centroidscheck[c_num,]-centroids[c_num,]))
    centroidchange=norm(centroiddiff,norm_typ)
    # Check if the difference between the new and old centres has gone below the threshold
  }
  else{
    centroidchange=0
  }
  centroid_df=data.frame(
    row.names=c_num,
    thresh_check=(centroidchange<clust_threshold & iter>1)
  )
  if(centroid_df[c_num,'thresh_check']){
    centroid_df<-cbind(centroid_df,centroids[c_num,])
  }
  else{
    centroid_df<-cbind(centroid_df,centroidscheck[c_num,])
  }
  return(centroid_df)
}

rand_cent<-function(a,min_max_df,n_cents){
  cent<-runif(n_cents,min_max_df[a,'min'],min_max_df[a,'max'])
  out_cent=data.frame(a=cent)
  colnames(out_cent)=c(a)
  return(out_cent)
}

kmeans_skip_nan <- function(eff_dfs, centers = nr, nstart = (nr+1), iter.max = 300,
                            clust_threshold=1e-5,norm_typ="F",na_rm=TRUE,prob_on=TRUE){
  set.seed(123)
  # Find the range for each trait in the association matrix
  #min_mat=apply(eff_dfs,2,min,na.rm=na_rm)
  max_mat=apply(eff_dfs,2,max,na.rm=na_rm)
  min_max_df=data.frame(
    row.names=names(eff_dfs),
    min=apply(eff_dfs,2,min,na.rm=na_rm),
    max=apply(eff_dfs,2,max,na.rm=na_rm)
  )
  min_max_df<-min_max_df %>% na.omit()
  # Store the non-NaN traits
  #axes_nonnan=names(!is.na(min_mat))
  # Create a centre point for each cluster. Randomly select within the range for each trait.
  #centroids=data.frame(row.names=1:centers)
  centroid_list=lapply(rownames(min_max_df),rand_cent,
                       n_cents=centers,min_max_df=min_max_df)
  centroids<-Reduce(cbind,centroid_list)
  SNP_list<-rownames(eff_dfs)
  nSNPs=length(SNP_list)
  clust_samp<-replicate(nSNPs,sample(1:centers,1))
  # Create cluster dataframe for snps
  cluster_df=data.frame(
    row.names = rownames(eff_dfs),
    clusters=clust_samp
  )
  cluster_df['clust_prob']=numeric()
  # Do the random assignment in dataframe setup not in loop
  clust_dist_list<-lapply(SNP_list,cent_dist_calc,
                          eff_dfs=eff_dfs,cluster_df=cluster_df,
                          centroids=centroids,norm_typ=norm_typ)
  clust_dist_df<-Reduce(rbind,clust_dist_list)
  cluster_df<-cbind(cluster_df,clust_dist_df)
  # Repeat for less cluster numbers and choose the min result
  # Recluster until the cluster centres converge or iter.max reached
  for (iter in 1:iter.max){
    # For each SNP find the cluster with the closest centre.
    snp_clust_list <- lapply(SNP_list, snp_closest_clust,
                             eff_dfs=eff_dfs,cluster_df=cluster_df,
                             centroids=centroids,norm_typ=norm_typ)#,axes_nonnan=axes_nonnan)
    # Combine the list of dataframes into one dataframe.
    # Override Cluster_df with the new assignment
    cluster_df<-Reduce(rbind,snp_clust_list)
    # Recompute the centroids based on the average of the clusters.
    # Check if the previous centres differ from the cluster means.
    thresh_list<-lapply(rownames(centroids),clust_cent_check,
                        iter=iter,cluster_df=cluster_df,
                        eff_dfs=eff_dfs,centroids=centroids,
                        na_rm=na_rm,norm_typ=norm_typ)#,
                        #axes_nonnan=axes_nonnan)
    thresh_check_df<-Reduce(rbind,thresh_list)
    # Are all the new centroids within the threshold of the previous
    thresh_t=sum(thresh_check_df$thresh_check==FALSE)
    if(all(thresh_check_df$thresh_check)){
      print('Clusters converged')
      print(iter)
      break
    }
    centroids<-thresh_check_df[,!names(thresh_check_df) %in% c('thresh_check')]
  }
 if (prob_on){cluster_df$clust_prob<-cluster_df$clust_dist %>% clust_prob_calc()
 }
 return(cluster_df)
}
