### 3) Clustering

# cochran's Qtest >>>>
q.meta.test = function(b_x, se_x){  #https://cjvanlissa.github.io/Doing-Meta-Analysis-in-R/heterogeneity-statistics.html
  se_x[which(se_x==0)]=1  ### check!
  w0 = se_x^-2
  w = w0/sum(w0)
  
  meta.bet_x = sum(b_x*w)
  meta.se_x = (sum(w^2*se_x^2))^0.5
  Q = sum( (b_x-meta.bet_x)^2 * w0 )
  df = length(b_x)-1
  Q.pval = 1 - pchisq(Q,df)
  return(list("Q"=Q, "Q.pval"=Q.pval))
}
# Fitting K-Means clustering Model 
kmeansIC = function(fit){
  #https://stackoverflow.com/questions/15839774/how-to-calculate-bic-for-k-means-clustering-in-r
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(data.frame(AIC = D + 2*m*k,
                    BIC = D + log(n)*m*k))
}

cluster_kmeans <- function(beta_df,aim_df,nr){
  # Set up for k-means clustering
  b_df <- beta_df[,aim_df$b_df_ind] 
  # Delete Any rows with NaNs
  # Test function: more axes should yeild less or equal data points. Point has full co-ordinate accros axes
  
  eff_df <- abs(b_df)
  if (length(aim_df$label)<=1){eff_dfs<-t(eff_df)}
  else {eff_dfs <- t(scale(t(eff_df)))}
  
  # Initial cluster dataframe
  cluster_list = list()
  IC_df = data.frame(matrix(data=NA, nrow = nr, ncol=1))
  colnames(IC_df) = c("AIC")
  IC_df$nCluster = 2:(nr+1)
  for(i in 2:(nr+1)){
    set.seed(240) # setting seed
    kmeans.re <- kmeans(eff_dfs, centers = i, nstart = (nr+1), iter.max = 300)
    cluster_list[[length(cluster_list)+1]] = kmeans.re
    IC = kmeansIC(kmeans.re)
    IC_df[(i-1),1] = IC$AIC
  }
  # cluster number identification for each observation
  kmeans.minAIC = cluster_list[[which(IC_df$AIC==min(IC_df$AIC))]]
  nClust.AIC = max(kmeans.minAIC$cluster)
  
  # plot AIC
  pdf(file=paste0(res_dir,"AIC-ncluster_",EXP_pheno,".pdf"), width=5, height = 7)
  plot(IC_df$nCluster,IC_df$AIC, xlab="nCluster", ylab="AIC")
  abline(v=nClust.AIC,col="red")
  dev.off()
  
  # assigning SNPs to clusters
  AICclusters_rsid = list()
  for(i in 1:nClust.AIC){
    AICclusters_rsid[[i]] = names(kmeans.minAIC$cluster)[which(kmeans.minAIC$cluster==i)]
  }
  print(paste0("Number of SNPs in each AIC grouped clusters: ", paste0(lengths(AICclusters_rsid), collapse = ", ")))
  AICclusters_rsid_df = t(plyr::ldply(AICclusters_rsid, rbind))
  write.csv(AICclusters_rsid_df, paste0(res_dir,"AICclusters_rsid_",EXP_pheno,".csv"), row.names = FALSE)
  
  return(kmeans.minAIC)
}
all_na_check <- function(b_df,aim_df){
  #' Check if the entire row is NaN
  nend=length(aim_df$b_df_ind)
  narows=which(is.na(b_df[,aim_df$b_df_ind[nend]]))
  if (length(narows)==dim(b_df)[1]){
    # All rows are NaN so trait will be removed from trait
    return(1)
  }
  else {return(0)}
}
remove_na_from_row <- function(b_df,aim_df){
  #' Remove rows which contain NA values but only from the columns aim_df$b_df_ind.
  #' If any of the columns contain only NaNs then return the dataframe with the column removed. 
  nend=length(aim_df$b_df_ind)
  narows=which(is.na(b_df[,aim_df$b_df_ind[nend]]))
  if (length(narows)==0){# No NaN rows
  }
  else if (all_na_check(b_df,aim_df)){
    # All rows are NaN so trait will be removed from trait
  }
  else{
    # Remove the NaN rows from the dataframe
    b_df <- b_df[-c(narows)]
  }
  return(b_df)
}

clust_metric <- function(cs1,cs2,norm_typ){
  if (length(cs1)<=1){
    return(abs(cs1-cs2))
  }
  else{
    clustid_1=which(!is.na(cs1))
    clustid_2=which(!is.na(cs2))
    clust_ids=intersect(clustid_1,clustid_2)
    if (length(clust_ids)==0){return(NaN)}
    else if(length(clust_ids)==1){return(abs(cs1[clust_ids]-cs2[clust_ids]))}
    else{
      return(norm(as.matrix(cs1[clust_ids]-cs2[clust_ids]),norm_typ))
    }
  }
}

aim_df_add_a <- function(aim_df,a,axes,b_df){
  #' Add the trait a to the dataframe of traits.
  #' Add the traits label, it's index in the trait list and it's index in the 
  #' beta dataframe.
  if (a %in% aim_df$label){
    return(aim_df)}
  else{
  aim_df <- aim_df %>% add_row(label= a,'axes_ind'=which(axes==a)[1],'b_df_ind'=0)
  aim_df[aim_df$label==a,'b_df_ind'] <- which(colnames(b_df)== a)[1]
  aim_df[aim_df$label==a,'nSNPs'] <- sum(!is.na(b_df[,a]))[1]
  return(aim_df)}
}
testna <- function(unstdBeta_df,trait_axes){
  aim_df <- data.frame(
    label =character(),
    axes_ind = integer(),
    b_df_ind = integer(),
    nSNPs = integer()
  )
  # Initialise with outcome
  aim_df <- aim_df_add_a(aim_df,trait_info$phenotype[125],trait_info$phenotype,
                         unstdBeta_df)
  b2_df <- remove_na_from_row(unstdBeta_df,aim_df)
  return(1)
}

#nr=10
#trait_axes_ind=which(colnames(stdBeta_df_noEXP) %in% trait_axes)
#load(paste0(res_dir,"QCdata_",EXP_pheno,"ClusterIter",1,".Rdata"))
#test <- cluster_kmeans(unstdBeta_df,aim_df,nr)
