kmeans_skip_nan <- function(eff_dfs, centers = i, nstart = (nr+1), iter.max = 300,
                            threshold=1e-5,norm_typ="F",na_rm=FALSE){
  set.seed(123)
  axes=names(eff_dfs)
  centroids=data.frame(
  )
  for(a in axes){
    mina=min(eff_dfs[a],na.rm=na_rm)
    maxa=max(eff_dfs[a],na.rm=na_rm)
    centps=runif(i,mina,maxa)
    centroids[a]=centps
  }
  max_dist=max(eff_dfs,na.rm=na_rm)-min(eff_dfs,na.rm=na_rm)
  clusters=c()
  clust_dist=c()
  clust_prob=c()
  for (iter in 1:iter.max){
    for(snp_score in eff_dfs){
      snp_dist=max_dist
      snp_id=which(eff_dfs==snp_dist)
      snp_clust_num=0
      for(cent in centroids){
        clust_num=which(centroids==cent)
        dist=norm(na.omit(cent-snp_score),norm_typ)
        if(dist<snp_dist){
          snp_dist=dist
          snp_clust_num=clust_num
        }
      }
      clusters[snp_id]=snp_clust_num
      clust_dist[snp_id]=snp_dist
      clust_prob[snp_id]=1/(1+snp_dist)
    }
    # Recompute the centroids based on the average of the clusters.
    centroidscheck=data.frame()
    thresh_check=c()
    for (i in 1:length(centroids)){
      snp_list=which(clusters==i)
      snp_scores=eff_dfs[snp_list]
      for(a in axes){
        centroidscheck[i,a]=mean(snp_scores[a],na.rm=na_rm)
      }
      print(centroidscheck)
      centroiddiff=na.omit(centroidscheck[i]-centroids[i])
      centroidchange=norm(centroiddiff,norm_type)
      if(centroidchange<threshold){
        thresh_check<-append(thresh_check,TRUE)
      }
      else{
        thresh_check<-append(thresh_check,FALSE)
      }
    }
    if(all(thresh_check)){
      break
    }
    # Are all the new centroids within the threshold of the previous
    
  }
 print('Maximum number of iterations reached without converging. ')
 cluster_df=data.frame(
   clusters=clusters,
   clust_dist=clust_dist,
   clust_prob=clust_prob
 )
 return(cluster_df)
}