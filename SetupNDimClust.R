# 0 Setup
# - Install packages required for NDim clustering.
# - Load the data to be clustered.

#library(data.table)
library(dplyr)
#library(ggplot2)
#library(ggrepel)
#library(tidyr)
#library(TwoSampleMR)
#library(tidyverse)

#' Data inputs needed: unstanderdised association matrix between IV and traits,
#' SE matrix of association matrix, t-stat matrix of association matrix,
#' P-value matrix of association matrix,
#' trait info (trait, description, sample size)

# variable set-up
threshmul <- 5.0
clust_threshold <- 1e-5
thresh_norm <- "F"
clust_norm <- "F"
nr <- 10 # Number of clusters
np <- 2 # Number of PCs
na_percent <- 0.25 # The percentage of a column that can acceptably be not NaN

# Testing dimensions
test <- 0 # testing switch
num_trait0 <- 1
num_trait1 <- 410
num_rows <- 50

# Files containing data
hail_gcorr_dir <- paste0(data_dir, "Hail_AllxAll.csv")
fpaths_fil_dir <- paste0(data_dir, "fpaths_fil_nfil.txt")

# Contain the variables into lists, threshold$diff and
# nums$max_dist will both be updated later.
thresholds <- list("diff_mul" = threshmul,
                       "diff" = 1e-5,
                       "clust" = clust_threshold)
na_handling <- list("narm" = TRUE, "percent" = na_percent)
nums <- list("nr" = nr, "np" = np, "max_dist" =  1.0)
norm_typs <- list("clust" =  clust_norm, "thresh_norm" =  thresh_norm)

#' read in data as matrices.
#' rows labeled by SNP_id and columns labelled by trait info.
#' The entries at each position correspond to the values for beta, se, t, and p
unstd_beta_df <- as.matrix(data.table::fread(
                        paste0(data_dir, "unstdBeta_df.csv")),
                        rownames = 1)
unstd_se_df   <- as.matrix(data.table::fread(
                        paste0(data_dir, "unstdSE_df.csv")),
                        rownames = 1)
tstat_df     <- as.matrix(data.table::fread(paste0(data_dir, "tstat_df.csv")),
                        rownames = 1)
pval_df      <- as.matrix(data.table::fread(paste0(data_dir, "pval_df.csv")),
                        rownames = 1)
trait_info   <- data.table::fread(paste0(data_dir, "trait_info_nfil.csv"))

# Find the label for the Outcome trait and the first Exposure trait
row <- which(trait_info$pheno_category == "Outcome")
out_pheno <- trait_info$phenotype[row]
row <- which(trait_info$pheno_category == "Exposure")[1]
exp_pheno <- trait_info$phenotype[row]

# Crop data for testing
if (test) {
data_matrics <- crop_data(mat_list = list("beta" = unstd_beta_df,
                               "pval" = pval_df,
                               "se" = unstd_se_df),
                          trait_df = trait_info,
                          out_pheno = out_pheno,
                          exp_pheno = exp_pheno,
                          n_rows = num_rows,
                          n_col0 = num_trait0,
                          n_col1 = num_trait1)
} else {
  # Collect the matrices into one object
  data_matrices <- list("beta" = unstd_beta_df,
                        "pval" = pval_df,
                        "se" = unstd_se_df,
                        "trait_info" = trait_info)
}
