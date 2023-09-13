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
clust_typ_list <- c(clust_typ_str2, clust_typ_str1)
# - probability used from beta value
bp_on1 <- TRUE
bp_on2 <- FALSE
bp_on_list <- c(bp_on1, bp_on2)
# - probability used for cluster scoring
clust_prob_on1 <- TRUE
clust_prob_on2 <- FALSE
clust_prob_on_list <- c(clust_prob_on1, clust_prob_on2)

set_directory <- function(res_dir0,iter_traits) {
  res_dir <- paste0(res_dir0, method_str(iter_traits ), "/")
  return(res_dir)
}

method_str <- function(iter_traits) {
  if (iter_traits$bp_on) {
    bp_str <- "_bpON"
  } else {
    bp_str <- "_bpOFF"
  }
  if (iter_traits$clust_prob_on) {
    clust_prob_str <- "_clustprobON"
  } else {
    clust_prob_str <- "_clustprobOFF"
  }
  return(paste0(clust_typ_str, bp_str, clust_prob_str))
}

desc_str <- function(iter_traits){
  if (iter_traits$bp_on) {
    bp_str <- "bp on"
  } else {
    bp_str <- "bp off"
  }
  if (iter_traits$clust_prob_on) {
    clust_prob_str <- "ClustProb on"
  } else {
    clust_prob_str <- "ClustProb off"
  }
  return(paste(clust_typ_str, "and", bp_str, "and", clust_prob_str))
}

make_iter_df <- function(clust_typ_list,bp_on_list,clust_prob_on_list) {
  iter_df_full <- data.frame(row.names= integer(),
                             "bp_on" = logical(),
                             "clust_prob_on" = logical(),
                             "clust_typ" = character())
  for (clust_typ_str in clust_typ_list) {
    for (bp_on in bp_on_list) {
      for (clust_prob_on in clust_prob_on_list) {
        iter_traits <- data.frame(
          "bp_on" = bp_on,
          "clust_prob_on" = clust_prob_on,
          "clust_typ" = clust_typ_str)
        iter_df_full <-rbind(iter_df_full,iter_traits)
      }
    }
  }
  return(iter_df_full)
}


full_prog <- function(iter, iter_df){
  iter_traits <- iter_df[iter,]
  res_dir <- set_directory(res_dir0, iter_traits)
  source("SetupNDimClust.R")
  # Find the distances between all points to initialise the threshold
  # for cluster difference.
  dist_df <- setup_dist(unstd_beta_df, norm_typs$clust)
  threshold$diff <- threshold$diff_mul * var(dist_df$dist, na.rm = TRUE)
  max_dist <- max(dist_df$dist,na.rm = na_rm)
  nums$max_dist <- max_dist
  print("Begining algorithm for inputs")
  print(iter_traits)
  out <- clust_pca_compare(data_matrices = data_matrices,
                           thresholds = thresholds,
                           na_handling = na_handling,
                           iter_traits = iter_traits,
                           norm_typs = norm_typs,
                           nums = nums
                          )
  print("Clust done")
  max_diff_df <- out$max_diff
  c_scores <- out$clust_scores
  #max_diff_df1 <- test1 %>% create_max_diff(thresh_norm)
  print("Diff done")
  c_scores %>% plot_trait_heatmap(iter_traits)
  print("Heatmap plot done")
  max_diff_df %>% plot_max_diff(iter_traits)
  print("Diff plot done")
  out <- list("iter_df" = iter_df,"c_scores" = c_scores,"max_df" = max_diff_df)
  return(out)
}

# Initialise the dataframe for storing the run details of each iteration
iter_df <- make_iter_df(clust_typ_list,bp_on_list,clust_prob_on_list)

niter <- dim(iter_df)[1]

res_out <- lapply(1:niter, full_prog,
                  iter_df = iter_df)

# MAKE A LIST VERSION
#res_dir <- paste0(res_dir0, "/")
# Plot both max_diffs on the same plot
#plot_max_diff_list(max_diff_df0, iter_df)
