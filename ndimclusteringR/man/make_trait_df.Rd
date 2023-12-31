% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_trait_df.R
\name{make_trait_df}
\alias{make_trait_df}
\title{Create a dataframe for traits}
\usage{
make_trait_df(pheno_list, data_mat, na_percent)
}
\arguments{
\item{pheno_list}{- list of all traits}

\item{data_mat}{- matrix of the data.}

\item{na_percent}{- the percent of a data column that is
required to be not NaN.}
}
\value{
trait_df
}
\description{
The dataframe "trait_df" has the columns "label" and "a_ind" the index
position for the trait. Only include traits whose column
in data_mat meets the na_percent threshold.
}
\details{
Check traits using the function [check_trait].
trait_df is the "rbind" of all the valid traits.
trait_df has columns:
* "label" - The labels for the traits.
Also the column names of the data matrix.
* "a_ind" - The position in the trait list.
}
