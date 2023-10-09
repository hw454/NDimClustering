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
# 1 Setup the directory and algorithm variables
# 2 Setup the data matrices loaded from the defined directory
# 3 For each set of traits run the following\:
# 3a Find the principal components.
# 3b If requrired convert data to angles.
# 3c Cluster the data.
# 3d Score the clusters on their association with PCs.
# 3e Plot the clusters
# 3f Plot the cluster scores
# 3g Plot the PC trait scores

devtools::install("ndimclustering")
library("ndimclustering")

# USER INPUTS
# Location of the data directory
data_dir <- "../NDimClustInputs/BMI_SBP/"
res_dir0 <- "../NDimClustResults/BMI_SBP/"

# Versions to be iterated through.
# - type of clustering
clust_typ_str1 <- "basic_angle"
clust_typ_str2 <- "min_angle"
# - probability used from beta value
bp_on1 <- TRUE
bp_on2 <- FALSE
# - probability used for cluster scoring
clust_prob_on1 <- TRUE
clust_prob_on2 <- FALSE
# - Iterate through all change or increment
ndim_typ <- "all"

# Create diretory for results if it doesn't exist
if (!file.exists(res_dir0)) {
  dir.create(file.path(res_dir0))
}
# Group terms
clust_typ_list <- c(clust_typ_str2, clust_typ_str1)
bp_on_list <- c(bp_on1, bp_on2)
clust_prob_on_list <- c(clust_prob_on1, clust_prob_on2)

source("SetupNDimClust.R")


# Initialise the dataframe for storing the run details of each iteration
iter_df <- make_iter_df(clust_typ_list,
                            bp_on_list,
                            clust_prob_on_list,
                            ndim_typ)

niter <- dim(iter_df)[1]

# For each each set of input parameters run the full clustering program
res_out <- lapply(1:niter, full_cluster_and_plot,
                  res_dir0 = res_dir0,
                  data_matrices = data_matrices,
                  iter_df = iter_df,
                  out_pheno = out_pheno,
                  thresholds = thresholds,
                  na_handling = na_handling,
                  norm_typs = norm_typs,
                  nums = nums)