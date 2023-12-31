% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_max_dist.R
\name{max_dist_calc}
\alias{max_dist_calc}
\title{Calculate the maximum distance between any two points.}
\usage{
max_dist_calc(score_mat, norm_typ = "F", na_rm = TRUE)
}
\arguments{
\item{score_mat}{data matrix}

\item{norm_typ}{type of data to use for distance}

\item{na_rm}{bool switch on whether to remove NaNs in min and max.}
}
\value{
max_dist
}
\description{
Points are given by rows of the dataframe.
  Find the maximum distance between any two points in "score_mat".
  Distance is given by the norm "norm_typ"
}
\details{
To get an upper and lower bound for the points find the max in each column
  and the min in each column. The distance between these is an upper bound
  for the maximum distance between the points in the matrix
}
\seealso{
Other distance_functions: 
\code{\link{calc_clust_dist_col}()},
\code{\link{calc_clust_dist}()},
\code{\link{calc_member_dist_cent}()},
\code{\link{calc_pair_dist_df}()},
\code{\link{calc_snp_cent_dist}()}
}
\concept{distance_functions}
