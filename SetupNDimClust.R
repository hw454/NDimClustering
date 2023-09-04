# 0 Setup
# - Install packages required for NDim clustering.
# - Load the data to be clustered.

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(TwoSampleMR)
library(tidyverse)

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
np <- 3 # Number of PCs

# Testing dimensions
test <- 0 # testing switch
num_trait0 <- 400
num_trait1 <- 410
num_rows <- 50


# Files containing data
hail_gcorr_dir <- paste0(data_dir, "Hail_AllxAll.csv")
fpaths_fil_dir <- paste0(data_dir, "fpaths_fil_nfil.txt")

# Create results directory if it doesn't exist
if (!dir.exists(res_dir)) {
  dir.create(res_dir)
}

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

row <- which(trait_info$pheno_category == "Outcome")
out_pheno <- trait_info$phenotype[row]
row <- which(trait_info$pheno_category == "Exposure")[1]
exp_pheno <- trait_info$phenotype[row]

# Crop data for testing
if (test) {
  b_out <- unstd_beta_df[which(
    colnames(unstd_beta_df) %in% c(out_pheno, exp_pheno))]
  se_out <- unstd_se_df[which(
    colnames(unstd_se_df) %in% c(out_pheno, exp_pheno))]
  p_out <- pval_df[which(colnames(pval_df) %in% c(out_pheno, exp_pheno))]
  t_out <- tstat_df[which(colnames(tstat_df) %in% c(out_pheno, exp_pheno))]
  trait_out <- trait_info[trait_info$phenotype %in% c(out_pheno, exp_pheno)]

  unstd_beta_df <- cbind(unstd_beta_df[1:num_rows, num_trait0:num_trait1],
   b_out)
  unstd_se_df <- cbind(unstd_se_df[1:num_rows, num_trait0:num_trait1], se_out)
  tstat_df     <- cbind(tstat_df[1:num_rows, num_trait0:num_trait1], t_out)
  pval_df      <- cbind(pval_df[1:num_rows, num_trait0:num_trait1], p_out)
  trait_info   <- rbind(trait_info[num_trait0:num_trait1], trait_out)

  colnames(unstd_beta_df)[dim(unstdBeta_df)[2]] <- out_pheno
  colnames(unstd_se_df)[dim(unstdSE_df)[2]] <- out_pheno
  colnames(tstat_df)[dim(tstat_df)[2]] <- out_pheno
  colnames(pval_df)[dim(pval_df)[2]] <- out_pheno
  }
