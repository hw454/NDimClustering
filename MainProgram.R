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
source("SetupNDimClust.R")

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
test<-clust_compare(unstdBeta_df,unstdSE_df,pval_df,tstat_df,
                    trait_info$phenotype,threshold,thresh_norm,clust_norm)
test %>% plot_max_diff(thresh_norm)
test %>% plot_trait_heatmap()
