% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score_all_clusters.R
\name{score_all_clusters}
\alias{score_all_clusters}
\title{Function for scoring clusters}
\usage{
score_all_clusters(
  clusters_df,
  beta_mat,
  pval_mat,
  bp_on = TRUE,
  clust_prob_on = TRUE,
  num_axis = ncol(beta_mat)
)
}
\arguments{
\item{clusters_df}{dataframe of snps and their cluster numbers,
cluster distance and cluster probability.}

\item{beta_mat}{matrix with coloumns corrsponding to components and rows
corresponding to snps/ items.}

\item{pval_mat}{matrix with coloumns corrsponding to components and rows
corresponding to snps/ items.}

\item{bp_on}{Bool switch, if TRUE (default) use pvals in cluster scores.}

\item{clust_prob_on}{Bool switch, if TRUE (default) weight the cluster score
by the probability the snp is in the cluster.}

\item{num_axis}{The number of trait axis in use. default\: ncol("beta_mat")}
}
\value{
score dataframe
}
\description{
Score each cluster with \link{score_cluster}
Score each trait within each cluster with \link{score_trait}
Combine all the scores into a data frame with rows given by the cluster id
and columns are the trait, the values are the score.
There is also an additional column for the cluster number and number of axis.
}
\seealso{
Other scoring_functions: 
\code{\link{score_cluster}()},
\code{\link{score_trait}()}
}
\concept{scoring_functions}
