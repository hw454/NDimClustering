% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_kmeans.R
\name{cluster_kmeans}
\alias{cluster_kmeans}
\title{Cluster the data using kmeans then minimising aic.}
\usage{
cluster_kmeans(
  data_matrices,
  nclust = 10,
  max_dist = 10,
  how_cents = "rand",
  bin_p_clust = TRUE,
  threshold = 1e-05
)
}
\arguments{
\item{data_matrices}{List of the matrices of data.
Conatins the\:
  * "beta" initial matrix of assocaition data
  * "se" initial matrix of standard error data
  * "pval" initial matrix of p-value data
  * "trait_info" list of informations about the traits
  * "beta_pc" reformatted matrix of assocaition data
  * "se" reformatted matrix of standard error data
  * "pval" reformatted matrix of p-value data
  * "transform" matrix used for the pca transformation
  * Rows correspond to snps, columns to traits}

\item{nclust}{The maximum number of clusters to consider. default\:10}

\item{max_dist}{The maximum distance between any two points.}

\item{how_cents}{How the centroids will be initialised. If "rand (default)"
then the coordinates are created using a uniform distribution on the range
of each axis. If "points" then the centroid coordinates are assign to
nclust random points from the dataspace.}

\item{bin_p_clust}{Bool switch. If TRUE (default) then cluster
probability if used to weight the cluster scores.}

\item{threshold}{The threshold for distance between cluster centres
for clusters to be considered converged.}
}
\value{
clusters_df
}
\description{
1. Data is cropped to the complete cases.
2. The points are clustered using a k-means clustering [km_nan].
3. The data removed in the crop is assigned to the closest cluster.
4. The "clusters_df" dataframe of the cluster membership with columns\:
    * "clust_num" the number of clusters
    * "clust_dist" the distance from the snp to the cluster centre
      (or distance between angles)
    * "clust_prob" probability the snp is in the cluster. Calculated
      using [calc_clust_prob].
  clust_out = [clusters, clust_dist, centres].
  clust_out$clusters is cluster_df.
  The clust_dist dataframe has columns corresponding to the clusters,
    the rows are the datapoints and each value is the distance from that
    point to the cluster centre.
  The centres dataframe has columns corresponding to the axis and the rows
    corresponding to the clusters. Each row is the co-ordinate of the cluster
    centroid.
}
\seealso{
Other cluster_functions: 
\code{\link{cluster_dbscan}()},
\code{\link{cluster_kmeans_min}()},
\code{\link{get_aic}()},
\code{\link{km_nan}()}

Other k_means: 
\code{\link{cluster_kmeans_min}()},
\code{\link{km_nan}()},
\code{\link{make_rand_cent}()}
}
\concept{cluster_functions}
\concept{k_means}
