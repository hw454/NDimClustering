% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cscore_pair_diff.R
\name{cscore_pair_diff}
\alias{cscore_pair_diff}
\title{Find the score difference between two clusters.}
\usage{
cscore_pair_diff(cn1, cn2, c_score0, axis, norm_typ = "F")
}
\arguments{
\item{cn1}{- First cluster number}

\item{cn2}{- Second cluster number}

\item{c_score0}{- Dataframe for storing cluster scores within}

\item{axis}{- The traits/ columns to be included.}

\item{norm_typ}{- The type of norm to be used when
calculating the difference accross multiple distance.
The default is the Froebenius norm "F".}
}
\value{
score_diff
}
\description{
score_diff is the dataframe with terms:
  * "c_num1" - the first cluster number
  * "c_num2" - the second cluster number
  * "diff" - The norm of the difference between the
    two cluster scores accross the axis.
}
\details{
If the difference is NaN then score_diff is returned
with the same columns but empty rows.
}
