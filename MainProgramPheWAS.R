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

stdir=getwd()
setwd(paste0(stdir,"/NDimClustering/"))
# 1 Use OpenGWAS to load csv for first pair. Exposure and Outcome. 
source("QC_filtering.R")

# 2 Run MRanalysis and save results to theta, theta_se, b, bse, pval, tstat

# 3 Cluster based on results.
# CLUSTERS from PheWAS k-means clustering.
setwd(paste0(stdir,"/NDimClustering/"))
source("Clustering.R")

# 3 Test clustering on threshold
# -> Fail - Return to 1 with extra dimension
# -> Pass - Exit with results. 
threshold=0.2
axes=unique(trait_info$phenotype)
source("ClusteringCompare.R")
test<-clust_compare(kmeanns.minAIC,stdBeta_df,axes,threshold)