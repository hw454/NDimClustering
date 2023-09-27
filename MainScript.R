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

library('ndimclustering')
# 1 Load csv

# 2 Compute principal components
# 3 Cluster based on results.
# This will move into loop

# 3 Test clustering on threshold
# -> Fail - Return to 1 with extra dimension
# -> Pass - Exit with results.

# Location of the data directory
data_dir <- "../NDimClustInputs/BP_Cancer/"
res_dir0 <- "../NDimClustResults/BP_BreastCancer/"
if (!file.exists(res_dir0)) {
  dir.create(file.path(getwd(), res_dir0))
}

# Versions to be iterated through.
# - type of clustering
clust_typ_str1 <- "basic_angle"
clust_typ_str2 <- "min_angle"
clust_typ_list <- c(clust_typ_str2, clust_typ_str1)
# - probability used from beta value
bp_on1 <- TRUE
bp_on2 <- FALSE
bp_on_list <- c(bp_on1, bp_on2)
# - probability used for cluster scoring
clust_prob_on1 <- TRUE
clust_prob_on2 <- FALSE
clust_prob_on_list <- c(clust_prob_on1, clust_prob_on2)
# - Iterate through all change or increment
ndim_typ <- "iterative"

source("SetupNDimClust.R")


# Initialise the dataframe for storing the run details of each iteration
iter_df <- make_iter_df(clust_typ_list,
                            bp_on_list,
                            clust_prob_on_list,
                            ndim_typ)

niter <- dim(iter_df)[1]

# For each each set of input parameters run the full clustering program
res_out <- lapply(1:niter, cluster_and_plot,
                  data_matrices = data_matrices,
                  iter_df = iter_df,
                  out_pheno = out_pheno,
                  thresholds = thresholds,
                  na_handling = na_handling,
                  norm_typs = norm_typs,
                  nums = nums)