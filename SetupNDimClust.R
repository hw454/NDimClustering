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

## having obtained the top SNPs associated with BMI (genome-wide signficance), and ran phewas on them to get their effects on all UKBB traits we have.
## Data inputs needed: unstanderdised association matrix between IV and traits, SE matrix of association matrix, t-stat matrix of association matrix, P-value matrix of association matrix, trait info (trait, description, sample size)

# variable set-up
EXP_pheno = "21001"
data_dir = "./working-example/data/"                   # Location of the data directory
res_dir = "../NDimClustResults/"                       # Location of the results directory
exp_gcorr_thresh = 0.75
OUT_pheno = "845"
pheno_irnt = TRUE  #maybe if true change all grepl to paste0(EXP_pheno,"_irnt")
threshold=0.5
thresh_norm="M"
clust_norm="M"

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
