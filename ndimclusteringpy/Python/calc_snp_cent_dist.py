import numpy as np
import pandas as pd
""" Calculate the distance between the snps and their allocated cluster centroid

 :param snp_id: the label for the snp
 :param b_mat: the data matrix used for the coordinate of the snp.
   Rows are snps, columns are the trait axes.
 :param cluster_df: The dataframe with the snp cluster memberships.
   Columns are\:
    * "clust_num" the cluster number.
    * "clust_dist" distance from the snp to the cluster centroid.
    * "clust_prob" probability the snp is correctly allocated to the cluster.
   Rows are each snp.
 :param centroids_df: The dataframe of the centroid co-ordinates for each
   cluster. The rows are each cluster number and the columns are the traits.
 :param norm_typ: The type of norm to use in the distance calculation. The
   default is the Froebenius norm "fro".
 :param narm: Bool, whether to skip NaN in norm, default is "TRUE".

  Calculate the distance between the point in "b_mat"
   corresponding to "snp_id" and the cluster centroid allocated to "snp_id".
   This is stored in a dataframe "snp_clust_df" which has rowname "snp_id",
   and column "clust_dist", with the calculated distance.

 :return: snp_clust_df
 :rtype: dataframe"""
def calc_snp_cent_dist (snp_id, data_mat, cent, norm_typ = "fro", narm = True):
   snp_score = data_mat[snp_id, :]
   c_dist = np.linalg.norm(cent - snp_score, ord = norm_typ, skipna = narm)
   snp_clust_df <- pd.Dataframe(index = snp_id,
                       data={"clust_dist" : c_dist})
   return(snp_clust_df)