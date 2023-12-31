% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_km_nan.R
\name{test_km_nan}
\alias{test_km_nan}
\title{Test kmeans clustering ignoring NaN terms}
\usage{
test_km_nan()
}
\description{
Use common axes for distance calculations to igore NaNs and
retain data.
}
\details{
The cluster centres are randomly assigned using [make_rand_cent].
The cluster centres are check and reassigned using [check_clust_cent],
When the clusters are converged or the maximum number of iterations are
used then the clusters are return in the dataframe "cluster_df".
"clust_out" is the list containing the "clusters_df" dataframe
labelled "cluster" of the cluster membership with columns\:
  * "clust_num" the number of clusters
  * "clust_dist" the distance from the snp to the cluster centre
  (or distance between angles)
  * "clust_prob" probability the snp is in the cluster. Calculated
  using [calc_clust_prob].
"clust_out" also contains the "centroids_df" dataframe labelled "centres"
whose columns are the traits and rows are the cluster numbers. Additional
column "thresh_check" contains Bool indicating whether that cluster
converged.
}
\seealso{
Other tests: 
\code{\link{test_add_nclust_cols_empty}()},
\code{\link{test_add_nclust_cols_onerow}()},
\code{\link{test_add_np_cols_empty}()},
\code{\link{test_add_np_cols_onerow}()},
\code{\link{test_calc_clust_dist_col}()},
\code{\link{test_calc_clust_dist}()},
\code{\link{test_calc_member_dist_cent}()},
\code{\link{test_find_aic}()},
\code{\link{test_find_low_percent_col_mat}()},
\code{\link{test_find_low_percent_col}()},
\code{\link{test_make_clust_col}()},
\code{\link{test_na_check}()}
}
\concept{tests}
