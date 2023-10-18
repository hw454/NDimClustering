import pandas as pd

""" Calculate the distance between the snps and each cluster centroid

 :param b_mat: the data matrix used for the coordinate of the snp.
   Rows are snps, columns are the trait axes.
 :param centroids_df: The dataframe of the centroid co-ordinates for each
   cluster. The rows are each cluster number and the columns are the traits.
 :param norm_typ: The type of norm to use in the distance calculation. The
   default is the Froebenius norm "fro".

 Calculate the distance between the point in "b_mat"
   corresponding to "snp_id" and the cluster centroid allocated to "snp_id".
   This is stored in a dataframe "snp_clust_df" which has rowname "snp_id",
   and column "clust_dist", with the calculated distance.

 :return: snp_clust_df"""

def calc_clust_dist (b_mat, centroids_df,
                           norm_typ = "fro") :
   clust_nums = centroids_df.index()
   clust_dist_lists =[calc_clust_dist_col(c,
                                          centroids_df = centroids_df, 
                                          data_mat = b_mat)
                                          for c in clust_nums]
   clust_dist_df = pd.concat(clust_dist_lists, axis = 1)
   return(clust_dist_df)
