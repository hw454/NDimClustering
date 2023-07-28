### 3) Clustering
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
cluster_kmeans_basic <- function(beta_df,aim_df,nr,norm_typ="F",threshold=1e-5){
  #' Using the association scores for each SNP accross traits cluster the traits
  #' using kmeans. Return the cluster setup which minimises AIC.
  # Set up for k-means clustering
  b_df <- beta_df[,aim_df$b_df_ind] 
  # Filter NaNs before clustering
  b_df_comp <- b_df[complete.cases(b_df)]
  print(b_df_comp)
  # Test function: more axes should yeild less or equal data points. Point has full co-ordinate accros axes
  
  eff_df <- abs(b_df_comp)
  if (length(aim_df$label)<=1){eff_dfs<-t(eff_df)}
  else {eff_dfs <- t(scale(t(eff_df)))}
  # Initial cluster dataframe
  print(eff_dfs)
  set.seed(240) # setting seed
  clust_re <- kmeans_skip_nan(eff_dfs, centers = nr, nstart = (nr+1), iter.max = 300,threshold=threshold,norm_typ=norm_typ)
  # cluster number identification for each observation
  for (snp in names(clust_re$cluster)){
    clust_num<-clust_re$cluster[snp]
    clust_cent<-clust_re$centers
    snp_clust_dist<-norm(clust_cent-b_df[snp,],"F")
    clust_re$clust_dist[snp]=snp_clust_dist
  }
  # Assign NaNs to nearest cluster
  for (snp in (names(beta_df) & !names(b_df_comp))){
    snp_dist=max(dist_df)
    for (cent in clust_re$centers){
      snp_clust_dist<-norm(cent-b_df[snp,],"F")
      clust_re$clust_dist[snp]=snp_clust_dist
      if (snp_clust_dist<snp_dist){
        snp_dist<-snp_clust_dist
        clust_num<-which(clust_re$centers==cent)
      }
    }
    clust_re$cluster[snp]=clust_num
    clust_re$clust_dist[snp]=snp_dist
  }  
  clust_re$clust_prob<-1/(1+clust_re$clust_dist)
  names(clust_re$clust_prob)<-names(clust_re$cluster)
  return(clust_re)
}
cluster_kmeans_min <- function(beta_df,aim_df,nr){
  #' Using the association scores for each SNP accross traits cluster the traits
  #' using kmeans. Return the cluster setup which minimises AIC.
  # Set up for k-means clustering
  b_df <- beta_df[,aim_df$b_df_ind] 
  # Delete Any rows with NaNs
  # Test function: more axes should yeild less or equal data points. Point has full co-ordinate accros axes
  # Filter complete cases
  b_df_comp <- b_df[complete.cases(b_df)]
  eff_df <- abs(b_df)
  if (length(aim_df$label)<=1){eff_dfs<-t(eff_df_comp)}
  else {eff_dfs <- t(scale(t(eff_df)))}
  
  # Initial cluster dataframe
  cluster_list = list()
  IC_df = data.frame(matrix(data=NA, nrow = nr, ncol=1))
  colnames(IC_df) = c("AIC")
  IC_df$nCluster = 2:(nr+1)
  for(i in 2:(nr+1)){
    set.seed(240) # setting seed
    clust_re <- kmeans(eff_dfs, centers = i, nstart = (nr+1), iter.max = 300)
    cluster_list[[length(cluster_list)+1]] = clust_re
    IC = kmeansIC(clust_re)
    IC_df[(i-1),1] = IC$AIC
  }
  # cluster number identification for each observation
  clust_re_minAIC = cluster_list[[which(IC_df$AIC==min(IC_df$AIC))]]
  # cluster number identification for each observation
  for (snp in names(clust_re_minAIC$cluster)){
    clust_num<-clust_re_minAIC$cluster[snp]
    clust_cent<-clust_re_minAIC$centers
    snp_clust_dist<-norm(clust_cent-b_df[snp,],"F")
    clust_re_minAIC$clust_dist[snp]=snp_clust_dist
  }
  # Assign NaNs to nearest cluster
  for (snp in (names(beta_df) & !names(b_df_comp))){
    snp_dist=max(dist_df)
    for (cent in clust_re_minAIC$centers){
      snp_clust_dist<-norm(cent-b_df[snp,],"F")
      clust_re_minAIC$clust_dist[snp]=snp_clust_dist
      if (snp_clust_dist<snp_dist){
        snp_dist<-snp_clust_dist
        clust_num<-which(clust_re_minAIC$centers==cent)
      }
    }
    clust_re_minAIC$cluster[snp]=clust_num
    clust_re_minAIC$clust_dist[snp]=snp_dist
  }  
  names(clust_re_minAIC$clust_dist)<-names(clust_re_minAIC$cluster)
  clust_re_minAIC$clust_prob<-1/(1+clust_re_minAIC$clust_dist)
  return(clust_re_minAIC)
}

all_na_check <- function(b_df,aim_df){
  #' Check if the entire row is NaN
  nend=length(aim_df$b_df_ind)[1]
  narows=which(b_df[,aim_df$b_df_ind[nend]] %>% is.na())
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
get_dist <- function(a,center,dist_typ){
  if (dist_typ=='p_euclid'){
    # calculate the euclidean distance to the centre point
    dist<-norm(a-b,"F")
  }
  else if (dist_typ=='line_euclid'){
    # Calculate the distance to the line given by center
    dist<-1
  }
  else{dist<-1}
  return(dist)
}
clust_prob <- function(clust_df,b_df,dist_typ){
  #' For each SNP calculate the distance to the centre of the cluster
  #' distance calculated based of distance type dist_typ
   
  # iterate through snp
  # find clust centre by clust number
  for (SNP_clust in clust_df$id){
    point0=b_df[SNP_clust]
    
  }
  # find centre-snp distance
  # Invert distance for prob
}

aim_df_add_a <- function(aim_df,a,axes,b_df){
  #' Add the trait a to the dataframe of traits.
  #' Add the traits label, it's index in the trait list and it's index in the 
  #' beta dataframe.
  if (a %in% aim_df$label){
    return(aim_df)}
  else{
    aim_df <- aim_df %>% add_row(label= a,'axes_ind'=which(axes==a)[1],'b_df_ind'=0,'nSNPs'=0)
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
aim_df <- data.frame(
  label =character(),
  axes_ind = integer(),
  b_df_ind = integer(),
  nSNPs = integer()
)
# Initialise with outcome
# aim_df <- aim_df_add_a(aim_df,OUT_pheno,trait_info$phenotype,
#                        unstdBeta_df)
# aim_df <- aim_df_add_a(aim_df,trait_info$phenotype[2],trait_info$phenotype,
#                        unstdBeta_df)
# nr=3
# test <- cluster_kmeans(unstdBeta_df,aim_df,nr)
