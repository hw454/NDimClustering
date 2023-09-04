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

# 1 Load csv

# 2 Compute principal components
# 3 Cluster based on results.
# This will move into loop

# 3 Test clustering on threshold
# -> Fail - Return to 1 with extra dimension
# -> Pass - Exit with results.

source("./ClusteringFunction.R")
source("./ClustScores.R")
source("./NClust_Plots.R")
source("./CalcDist.R")
source("./kmeans_skip_nan.R")
source("./ClusteringCompare.R")
source("./PrincipalComponentAnalysis.R")

# Location of the data directory
data_dir <- "./working-example/data/"
res_dir0 <- "../NDimClustResults/working-example/"
if (!file.exists(res_dir0)) {
  dir.create(file.path(getwd(), res_dir0))
}

# Versions to be iterated through.
# - type of clustering
clust_typ_str1 <- "basic"
clust_typ_str2 <- "min"
clust_typ_list <- c(clust_typ_str1, clust_typ_str2)
# - probability used from beta value
bp_on1 <- TRUE
bp_on2 <- FALSE
bp_on_list <- c(bp_on1, bp_on2)
# - probability used for cluster scoring
clust_prob_on1 <- TRUE
clust_prob_on2 <- FALSE
clust_prob_on_list <- c(clust_prob_on1, clust_prob_on2)

# Initialise the dataframe for storing the run details of each iteration
iter_df <- data.frame(
  index = integer(),
  clust_typ = character(),
  bp_on = logical(),
  clust_prob_on = logical()
)
iter <- 1
for (clust_typ_str in clust_typ_list) {
  for (bp_on in bp_on_list) {
    for (clust_prob_on in clust_prob_on_list) {
      if (bp_on) {
        bp_str <- "_bpON"
      } else {
        bp_str <- "_bpOFF"
      }
      if (clust_prob_on) {
        clust_prob_str <- "_clustprobON"
      } else {
        clust_prob_str <- "_clustprobOFF"
      }
      res_dir <- paste0(res_dir0, clust_typ_str, bp_str, clust_prob_str, "/")
      source("SetupNDimClust.R")
      # Find the distances between all points to initialise the threshold
      # for cluster difference.
      dist_df <- setup_dist(unstdBeta_df, clust_norm)
      diff_threshold <- threshmul * var(dist_df$dist, na.rm = TRUE)
      iter_df <- iter_df %>% add_row("index" = iter,
      "clust_typ" = clust_typ_str, "bp_on" = bp_on,
      "clust_prob_on" = clust_prob_on)
      print("Begining algorithm for inputs")
      print(iter_df[length(iter_df)])
      test1 <- clust_pca_compare(unstd_beta_df,unstd_se_df, pval_df,
                                 diff_threshold = diff_threshold, 
                                 thresh_norm = thresh_norm, 
                                 clust_threshold = clust_threshold,
                                 clust_norm, np = np, nr = np, which_clust = clust_typ_str,
                                 bp_on = bp_on, clust_prob_on = clust_prob_on, narm = TRUE)
      print("Clust done")
      max_diff_df1 <- test1 %>% create_max_diff(thresh_norm)
      print("Diff done")
      test1 %>% plot_trait_heatmap(clust_typ_str, bp_on, clust_prob_on)
      print("Heatmap plot done")
      max_diff_df1 %>% plot_max_diff(clust_typ_str, bp_on, clust_prob_on)
      print("Diff plot done")
      if (iter == 1) {
        clust_out <- test1
        clust_out["input_iter"] <- iter
        max_diff_df1["input_iter"] <- iter
        max_diff_df0 <- max_diff_df1
      } else {
        test1["input_iter"] <- iter
        clust_out <- rbind(clust_out, test1)
        max_diff_df1["input_iter"] <- iter
        max_diff_df0 <- rbind(max_diff_df0, max_diff_df1)
      }
      iter <- iter + 1
    }
  }
}

# MAKE A LIST VERSION
res_dir <- paste0(res_dir0, "/")
# Plot both max_diffs on the same plot
plot_max_diff_list(max_diff_df0, iter_df)
