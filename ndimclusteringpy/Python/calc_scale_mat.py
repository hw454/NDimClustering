import pandas as pd
""" Rescale the columns of the matrix "mat"

   Rescale the columns of the matrix `mat` so that the values lie
   between -1 and 1 with the with the same distribution.

 :param mat: The matrix whose columns are to be rescaled
 :param narm: Bool to indicate whether NaNs should be ignored in calculations.
   default "TRUE"

  Rescale of each column is found using [calc_col_scale] then
   column binding the results using cbind into "out_mat".

 :return: out_mat"""
def calc_scale_mat (mat, narm = TRUE) :
  # Compute the Sample Mean with na_rm
  xbar = mat.mean(axis=1, skipna = narm)
  # Compute the Sample SE with na_rm
  se = mat.std(axis=1, skipna = narm)
  out_list = [calc_col_scale(c, data = mat, mu = xbar, se = se) for c in mat.columns]
  out_mat = pd.concat(out_list, axis = 1)
  return(out_mat)