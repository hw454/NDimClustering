#' Function for scoring clusters
#'
#' @param clusters_df dataframe of snps and their cluster numbers,
#'   cluster distance and cluster probability.
#' @param beta_mat matrix with coloumns corrsponding to components and rows
#'   corresponding to snps/ items.
#' @param pval_mat matrix with coloumns corrsponding to components and rows
#'   corresponding to snps/ items.
#' @param bp_on Bool switch, if TRUE (default) use pvals in cluster scores.
#' @param clust_prob_on Bool switch, if TRUE (default) weight the cluster score
#'   by the probability the snp is in the cluster.
#' @param num_axis The number of trait axis in use. default\: ncol("beta_mat")
#'
#' @description Score each cluster with \link{score_cluster}
#' Score each trait within each cluster with \link{score_trait}
#' Combine all the scores into a data frame with rows given by the cluster id
#' and columns are the trait, the values are the score.
#' There is also an additional column for the cluster number and number of axis.
#'
#' @return score dataframe
#'
#' @export
score_all_clusters <- function(clusters_df, beta_mat, pval_mat,
                        bp_on = TRUE,
                        clust_prob_on = TRUE,
                        num_axis = ncol(beta_mat)) {
  # The dataframe of clusters is used alongside the beta and pvalue
  # for each SNP to score the clusters on their association with
  # each trait. Information for the traits is stored in aim_df
  clust_nums <- unique(clusters_df$clust_num)
  clust_scores_list <- lapply(clust_nums, score_cluster,
    beta_mat = beta_mat, clusters_df = clusters_df,
    pval_mat = pval_mat, bp_on = bp_on, num_axis = num_axis)
  clust_scores <- Reduce(rbind, clust_scores_list)
  return(clust_scores)
}