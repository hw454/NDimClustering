% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_tot_within.R
\name{calc_tot_within}
\alias{calc_tot_within}
\title{Calculate the total sum of the cluster distances divided by the variance}
\usage{
calc_tot_within(data, group_col = "clust_num", dist_col = "clust_dist")
}
\arguments{
\item{data}{matrix of data}

\item{group_col}{The column indicate terms to go together
default\:"clust_num"}

\item{dist_col}{The name of the column containing the term distances.
default\:"clust_dist"}
}
\value{
The sum of the clusters sum of squares/ divided by the variance
  of the distance data.
}
\description{
Calculate the sum of the squares for each cluster group using
  [calc_sum_sq_clusts]
}
\seealso{
Other cluster_properties: 
\code{\link{calc_clust_ic}()},
\code{\link{calc_clust_prob}()},
\code{\link{calc_sum_sq_clusts}()}

Other ic_functions: 
\code{\link{calc_sum_sq_clusts}()},
\code{\link{find_all_ic}()}
}
\concept{cluster_properties}
\concept{ic_functions}
