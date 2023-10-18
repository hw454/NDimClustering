""" Calculate the distance between the snps and the cluster centroid

 :param c_num: The clust number to calculate the distance to
 :param b_mat: the data matrix used for the coordinate of the snp.
   Rows are snps, columns are the trait axes.
 :param centroids_df: The dataframe of the centroid co-ordinates for each
   cluster. The rows are each cluster number and the columns are the traits.
 :param norm_typ: The type of norm to use in the distance calculation. The
   default is the Froebenius norm "fro".

 Calculate the distance between each point in "b_mat"
   and the cluster centroid allocated to c_num.
   This is stored in a dataframe "clust_dist_df" which has rownames for
   the snps "snp_id", and columns "clust_i", for each column.

 :return: clust_dist_df
 :rtype: dataframe
"""
def calc_clust_dist_col (c_num, data_mat, centroids_df,
                           norm_typ = "F"):
   snp_list = data_mat.index
   cent = centroids_df[c_num, ]
   clust_dist_lists =[calc_snp_cent_dist(snp,
                                        cent = cent,
                                        data_mat = data_mat) 
                                        for snp in snp_list]
   clust_dist_col_df = pd.concat[clust_dist_lists]  
   clust_dist_col_df.columns = ["clust_", c_num]
   return(clust_dist_col_df)
