""" Add "np" columns to a dataframe.

  Representing the clusters a column for the clust distance
   is added for each cluster.
   The label for each column is given by the function [make_clust_cols].

 :param df: the dataframe to add columns to
 :param np: the number of columns to add.

 :return: df
 :rtype: dataframe
"""
def add_nclust_cols(df,nr):
  # Create dataframe with columns with labels Pi
  # for i in 1 to np to the dataframe.
  p_cols = [make_clust_col(i) for i in range(1,nr)]
  np_df = pd.concat(p_cols, axis =1)
  full_df= pd.concat([df, np_df], axis = 1)
  return(full_df)