import pandas as pd
""" Check convergence for the cluster centroid

 :param c_num: The number of the cluster to be checked
 :param clustnum_df: The cluster membership dataframe. The rows are the snps,
  The columns are "clust_num".
 :param b_mat: The data for each snp. Each row is the snps coordinate.
 :param centroids_df The dataframe of the centroids for each cluster on
   the previous iteration. The rows are the cluster numbers and the
   columns are the centroid value on that trait axis.
 :param na_rm: Boolean. If TRUE (default) then remove NaNs in means.
 :param norm_typ: String indicating the type of norm to be used when
   calculating the distance between the centroids. Default is the Froebenius
   norm "fro".
 :param clust_threshold: The threshold for how close the clusters should
   be to be considered converged. default \:1e-5

   The check can only be made if there are snps which are members of the
   cluster.
   If the cluster is not empty then the new centres are found using
     the means of the snps in the cluster.
   If the new cluster centroid is within "clust_threshold" of the previous
     centroid then "thresh_check"=TRUE.
   If thresh_check then create a dataframe "out_centroid_df" with the
     previous centroid and the column "thesh_check" with TRUE in it.
   else create a dataframe "out_centroid_df" of the new centroids and
     the column "thresh_check" with FALSE in it.

 :return: out_centroid_df
 :rtype: dataframe
"""
def check_clust_cent (c_num, clustnum_df, b_mat, centroids_df,
                            na_rm = TRUE, norm_typ = "fro",
                            clust_threshold = 1e-5) :
  sub_snp_list = clustnum_df[clustnum_df == c_num].index
  snp_scores = b_mat[sub_snp_list]
  nterms = len(sub_snp_list)
  # Initialise the centroids checking dataframe
  new_centroids_df = centroids_df
  # Empty clusters not changed
  if (nterms):
    # Compute the new centers based on the mean of the snps in the clusters
    # If there's only one term then don't use mean
    # else: Mean for each column gives the value for the centre on each axis
    if (nterms == 1) :
      new_centroids_df[c_num] = snp_scores
    else :
      new_centroids_df[c_num] = snp_scores.mean(axis=1, skipna = na_rm)
    # Calculate how much the centroid has moved.
    centroiddiff = new_centroids_df[c_num]- centroids_df[c_num]
    centroidchange <- norm(centroiddiff, norm_typ)
    # Check if the diff between new and old centres is below the threshold
  else :
    centroidchange = 0

  check_df = pd.Dataframe(
    index = c_num,
    data = {
      "thresh_check" : (centroidchange < clust_threshold)
    }
  )
  if (check_df[c_num, "thresh_check"]) :
    out_centroid_df = pd.concat([check_df, centroids_df[c_num]], axis =1 )
  else :
    out_centroid_df = pd.concat([check_df, new_centroids_df[c_num]], axis = 1)
  return(out_centroid_df)