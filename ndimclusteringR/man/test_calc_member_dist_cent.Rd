% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_calc_member_dist_cent.R
\name{test_calc_member_dist_cent}
\alias{test_calc_member_dist_cent}
\title{Test `calc_member_dist_cent` the function which finds the distance
from a point to the centre of a cluster numbered c_num}
\usage{
test_calc_member_dist_cent()
}
\description{
* Check output is dataframe
 * Check the number of rows is same data_mat
 * Check there is a column for each cluster
 * Check the column names matches the cluster number labels
}
\seealso{
Other tests: 
\code{\link{test_calc_clust_cent}()},
\code{\link{test_cluster_dbscan}()},
\code{\link{test_cluster_kmeans_min}()},
\code{\link{test_cluster_kmeans}()},
\code{\link{test_convert_point_to_angle}()},
\code{\link{test_crop_matrices}()},
\code{\link{test_get_aic}()},
\code{\link{test_join_mat_on_name}()},
\code{\link{test_km_nan}()},
\code{\link{test_reformat_data}()},
\code{\link{test_rescale_by_end_col}()}
}
\concept{tests}
