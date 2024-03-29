% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_all_clust_dist.R
\name{calc_all_clust_dist}
\alias{calc_all_clust_dist}
\title{Find the table of distances between point and cluster centroid.}
\usage{
calc_all_clust_dist(b_mat, centroids_df)
}
\arguments{
\item{b_mat}{the matrix of the data. The rows correpond to the snps.}

\item{centroids_df}{dataframe of the centroid co-ordinates. Each row
corresponds to a cluster number and the columns are the trait axes.}
}
\value{
snp_clust_dist_df
}
\description{
Find the distance between the snp and each cluster centroid.
}
\details{
The dataframe "snp_clust_dist_df"  has a row labelled by the "snp_id".
  The columns are the clust numbers, each value is the distance from the
  snp to that cluster.
}
\seealso{
Other distance_functions: 
\code{\link{calc_member_dist_cent}()},
\code{\link{find_closest_clust_snp}()}
}
\concept{distance_functions}
