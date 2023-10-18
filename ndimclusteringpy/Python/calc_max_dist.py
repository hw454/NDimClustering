import numpy as np
""" Calculate the maximum distance between any two points.

 :param score_mat: data matrix or scores
 :param norm_typ: type of data to use for distance. Default is the Frobenius norm "fro".
 :param na_rm: bool switch on whether to remove NaNs in min and max. Default is "TRUE".

   Points are given by rows of the dataframe.
   Find the maximum distance between any two points in "score_mat".
   Distance is given by the norm "norm_typ"

   To get an upper and lower bound for the points find the max in each column
   and the min in each column. The distance between these is an upper bound
   for the maximum distance between the points in the matrix

 :return: max_dist"""
def max_dist_calc (score_mat, norm_typ = "fro", na_rm = TRUE) :
  # Find the max distance based on range on each axis.
  max_p = score_mat.max(axis=1, skipna = na_rm)
  min_p = score_mat.min(axis=1, skipna = na_rm)
  max_dist == np.linalg.norm(max_p - min_p, ord = norm_typ)
return(max_dist)
