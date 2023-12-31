% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_clust_prob.R
\name{calc_clust_prob}
\alias{calc_clust_prob}
\title{Calculate the probability the snp is in the cluster}
\usage{
calc_clust_prob(d)
}
\arguments{
\item{d}{The snp distance to the cluster centre}
}
\value{
p
}
\description{
\deqn{p=\frac{1.0}{1.0+d}}
}
\seealso{
Other cluster_properties: 
\code{\link{calc_clust_ic}()},
\code{\link{calc_sum_sq_clusts}()},
\code{\link{calc_tot_within}()}
}
\concept{cluster_properties}
\concept{probability_functions}
