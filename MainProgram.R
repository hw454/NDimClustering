############ 
# Template file for n-dimensional clustering
# Take m SNPs and n axis.
###########
# Inputs:
# theta=nxnxm array of ratios of beta values
# theta_se=nxnxm array of ratios of SE values
# b= nxm array of beta values 
# bse = nxm array of beta_se values
# pval =nxm array of p-val for b
# tstat = nxm array of t-stat for b

# 0 Setup the packages and programs

# 1 Use OpenGWAS to load csv for first pair. Exposure and Outcome. 
#source("QC_filtering.R")

# 2 Run MRanalysis and save results to theta, theta_se, b, bse, pval, tstat
# For the worked-example this is already done
# You may also load MR results from openGWAS data

# 3 Cluster based on results.
#source("Clustering.R") # This will move into loop

# 3 Test clustering on threshold
# -> Fail - Return to 1 with extra dimension
# -> Pass - Exit with results. 

#source("QC_filtering_iterative.R")
source("ClusteringFunction.R")
source("ClustScores.R")
source("NClust_Plots.R")
source("ClusteringCompare.R") # This call must be last

EXP_pheno = "21001"
data_dir = "./working-example/data/"                   # Location of the data directory
OUT_pheno = "845"
res_dir0 = "../NDimClustResults/working-example/"   

clust_typ_str1='basic'
bp_on1=TRUE
clust_typ_str2='min'
bp_on2=FALSE
clust_typ_list=c(clust_typ_str1,clust_typ_str2)
bp_on_list=c(bp_on1,bp_on2)

iter_df<-data.frame(
  index=integer(),
  clust_typ=character(),
  bp_on=logical()
)
iter=1
for (clust_typ_str in clust_typ_list){
  for (bp_on in bp_on_list){
    print(iter)
    if (bp_on){bp_str<-'_bpON'}
    else{bp_str<-'_bpOFF'}
    res_dir<- paste0(res_dir0,clust_typ_str,bp_str,'/')
    source("SetupNDimClust.R")
    iter_df<- iter_df %>% add_row('index'=iter,'clust_typ'=clust_typ_str,'bp_on'=bp_on)
    test1<-clust_compare(unstdBeta_df,unstdSE_df,pval_df,tstat_df,
                    trait_info$phenotype,threshold,thresh_norm,clust_norm,
                    clust_typ_str,bp_on)
    print('Clust done')
    max_diff_df1<- test1 %>% create_max_diff(thresh_norm)
    print('Diff done')
    test1 %>% plot_trait_heatmap(clust_typ_str,bp_on)
    print('Heatmap plot done')
    max_diff_df1 %>% plot_max_diff(clust_typ_str,bp_on)
    print('Diff plot done')
    if (iter==1){
      clust_out<-test1
      clust_out['input_iter']<-iter
      max_diff_df1['input_iter']<-iter
      max_diff_df0<-max_diff_df1
    }
    else{
      test1['input_iter']<-iter
      clust_out<-rbind(clust_out,test1)
      max_diff_df1['input_iter']<-iter
      max_diff_df0<- rbind(max_diff_df0,max_diff_df1)
    }
    iter<-iter+1
  }
}

# MAKE A LIST VERSION
res_dir<- paste0(res_dir0,'/')
# Plot both max_diffs on the same plot
plot_max_diff_list(max_diff_df0,iter_df)

# EXP_pheno = "21001"
# data_dir = "../NDimClustInputs/"                   # Location of the data directory
# OUT_pheno = "845"
# res_dir = "../NDimClustResults/OpenGWASdata/"   
# test1<-clust_compare(unstdBeta_df,unstdSE_df,pval_df,tstat_df,
#                      trait_info$phenotype,threshold,thresh_norm,clust_norm,'basic')
# max_diff_df1<- test1 %>% create_max_diff(thresh_norm)
# test1 %>% plot_trait_heatmap()
# max_diff_df1 %>% plot_max_diff()
# 
# res_dir = "../NDimClustResults/working-exampleMinClust/"   
# test2<-clust_compare(unstdBeta_df,unstdSE_df,pval_df,tstat_df,
#                      trait_info$phenotype,threshold,thresh_norm,clust_norm,'min')
# test2 %>% plot_trait_heatmap()
# max_diff_df2<- test2 %>% create_max_diff(thresh_norm)
# max_diff_df2 %>% plot_max_diff()
# 
# # Plot both max_diffs on the same plot
# max_diff_df1 %>% plot_max_diff_both(max_diff_df2)
