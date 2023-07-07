# 0 Setup 
# - Install packages required for NDim clustering. 
# - Load the data to be clustered. 

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)
#library(xlsx)
library(TwoSampleMR)
library(tidyverse)

## having obtained the top SNPs associated with BMI (genome-wide signficance), and ran phewas on them to get their effects on all UKBB traits we have.
## Data inputs needed: unstanderdised association matrix between IV and traits, SE matrix of association matrix, t-stat matrix of association matrix, P-value matrix of association matrix, trait info (trait, description, sample size)

# variable set-up
EXP_pheno = "21001"
data_dir = "./working-example/data/"                   # Location of the data directory
res_dir = "../NDimClustResults/working-exampleMinClust/"                       # Location of the results directory
OUT_pheno = "845"
threshold=5.0 #FIXME threshold should be based on dataspread
thresh_norm="M"
clust_norm="M"

# Testing dimensions
test=0 # testing swicth
num_trait0=400
num_trait1=410
num_rows=100

# Files containing data
hail_gcorr_dir = paste0(data_dir,"Hail_AllxAll.csv")   # Location of the data files
fpaths_fil_dir = paste0(data_dir,"fpaths_fil_nfil.txt")

# Create results directory if it doesn't exist
if (!dir.exists(res_dir)){dir.create(res_dir)}

# read in data as matrices. 
# rows labeled by SNP_id and columns labelled by trait info. The entries at each position correspond to the values for beta, se, t, and p
unstdBeta_df = as.matrix(data.table::fread(paste0(data_dir,"unstdBeta_df.csv")), rownames=1)
unstdSE_df   = as.matrix(data.table::fread(paste0(data_dir,"unstdSE_df.csv"))  , rownames=1)
tstat_df     = as.matrix(data.table::fread(paste0(data_dir,"tstat_df.csv"))    , rownames=1)
pval_df      = as.matrix(data.table::fread(paste0(data_dir,"pval_df.csv"))     , rownames=1)
trait_info   = data.table::fread(paste0(data_dir,"trait_info_nfil.csv"))

# Crop data for testing
if (test){
  b_out     <- unstdBeta_df[which(colnames(unstdBeta_df) %in% c(OUT_pheno,EXP_pheno))]
  se_out    <-unstdSE_df[   which(colnames(unstdSE_df)   %in% c(OUT_pheno,EXP_pheno))]
  p_out     <- pval_df[     which(colnames(pval_df)      %in% c(OUT_pheno,EXP_pheno))]
  t_out     <- tstat_df[    which(colnames(tstat_df)     %in% c(OUT_pheno,EXP_pheno))]
  trait_out <-trait_info[   trait_info$phenotype         %in% c(OUT_pheno,EXP_pheno)]
  
  unstdBeta_df <- cbind(unstdBeta_df[1:num_rows,num_trait0:num_trait1],b_out)
  unstdSE_df   <- cbind(unstdSE_df[  1:num_rows,num_trait0:num_trait1],se_out)
  tstat_df     <- cbind(tstat_df[    1:num_rows,num_trait0:num_trait1],t_out)
  pval_df      <- cbind(pval_df[     1:num_rows,num_trait0:num_trait1],p_out)
  trait_info   <- rbind(trait_info[             num_trait0:num_trait1],trait_out)
  
  colnames(unstdBeta_df)[dim(unstdBeta_df)[2]] <- OUT_pheno
  colnames(unstdSE_df)[  dim(unstdSE_df  )[2]] <- OUT_pheno
  colnames(tstat_df)[    dim(tstat_df    )[2]] <- OUT_pheno
  colnames(pval_df)[     dim(pval_df     )[2]] <- OUT_pheno
  }
