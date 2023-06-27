clust_compare <-function(unstdBeta_df,unstdSE_df,pval_df,tstat_df,
                         axes,threshold,thresh_norm,clus_norm){
  aN=length(axes)
  aim=c(OUT_pheno)
  aim_df <- data.frame(
    label =aim[1],
    axes_ind = which(axes==OUT_pheno),
    b_df_ind = which(axes==OUT_pheno)
  )
  for (ai in 1:aN){
    # Update the traits forming the axis
    a=axes[ai]
    if (a %in% aim){# a is already stored in aim.
      }
    else {
      nr=30.0/length(aim)
      # Filter the data based on the current axes
      trait_info_iter=trait_info[trait_info$phenotype %in% aim]
      aim <- append(aim,a)
      QCFilter_Iter(unstdBeta_df,unstdSE_df,pval_df,tstat_df,trait_info,fpaths_fil_dir,ai)
      # Load the filtered data
      load(paste0(res_dir,"QCdata_",EXP_pheno,"ClusterIter",ai,".Rdata"))
      # Cluster the data on these axes
      if (a %in% colnames(stdBeta_df_noEXP)){
        if (a %in% aim_df$label){
          aim_df[aim_df$label==a,b_df_ind]=which(colnames(stdBeta_df_noEXP)== a)
        }
        else{
          aim_df_0 <- data.frame(
          label=a,
          axes_ind=which(axes==a),
          b_df_ind=which(colnames(stdBeta_df_noEXP)== a)
          )
          aim_df <- aim_df %>% rows_insert(aim_df_0)
        }
        cluster_df=cluster_kmeans(stdBeta_df_noEXP,sus_SNP_ind,aim_df$b_df_ind,nr) 
        c_nums<-unique(cluster_df$cluster)
        Nsnps=length(cluster_df)
        # Score the clustered data based on affilliations with axes.
        c_scores<- data.frame(
          clust_num =c_nums,
          clust_size=0,
          NA_Count=0
        )
        c_score0 <- clust_score(cluster_df,stdBeta_df_noEXP,aim_df$label,aim_df$b_df_ind) # Score accross all axis
        print(c_score0)
        if (ai==0){
          c_scores <- c_score0
        }
        else{
          c_scores <- merge(c_scores,c_score0, by="clust_num")
        }
        i=0
        N2=length(c_nums)
        N1=N2-1
        for (i in 1:N1){
          for (j in (i+1):N2){
            ci=c_nums[i]
            cj=c_nums[j]
            cs1<-as.matrix(c_scores[c_scores$clust_num==i,aim_df$label])
            cs2<-as.matrix(c_scores[c_scores$clust_num==j,aim_df$label])
            print(cs1)
            print(cs2)
            print(threshold*Nsnps)
            print(norm(cs1-cs2,thresh_norm))
            print(clust_metric(cs1,cs2,thresh_norm))
            if (clust_metric(cs1,cs2,thresh_norm)>(threshold*Nsnps)){
              print("Threshold met on outcome")
              print(a)
              return(c_scores)
            }
          }
        }
      }
      else if(a %in% aim_df$label){
        aim_df <- aim_df[!aim_df$label==a]
      }
    }
  }
  print("None of the outcomes clustering met the thresholding test. ")
  return(c_scores)
}
clust_metric <- function(cs1,cs2,norm_typ){
  if (length(cs1)<=1){
    return(abs(cs1-cs2))
  }
  else{return(norm(cs1-cs2,norm_typ))}
}

#test<-clust_compare(stdBeta_df ,SNP_ind,axes,threshold,thresh_norm, clust_norm)
test<-clust_compare(unstdBeta_df,unstdSE_df,pval_df,tstat_df,
                    trait_axes,threshold,thresh_norm,clust_norm)
