% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_clust_scatter_rgb_test.R
\name{plot_clust_scatter_rgb_test}
\alias{plot_clust_scatter_rgb_test}
\title{Plot the scatter plot or snps and principal components.}
\usage{
plot_clust_scatter_rgb_test(
  clust_dist_df,
  b_mat,
  se_mat,
  iter_traits,
  num_axis = 1,
  pw = 8,
  ph = 4
)
}
\arguments{
\item{clust_dist_df}{distances between each snp and cluster centres}

\item{b_mat}{The matrix of the score data}

\item{se_mat}{The matrix of the standarad errors associated with the scores.}

\item{iter_traits}{The iteration variables for the type of iteration.}

\item{num_axis}{The number of trait axis, default\:0}

\item{pw}{The plot width, default\:8}

\item{ph}{The plot heigh, default\:4}
}
\description{
Colour by cluster number.
The x and y axis are the assocation with the principal component
The width and height of the errorbars are given by the standard error.
}
