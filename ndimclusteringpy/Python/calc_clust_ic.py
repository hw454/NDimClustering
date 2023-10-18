""" Find the BIC and AIC for the clusters.

 :param clust: The cluster membership dataframe.
 :param group_col: The name of the column to group terms on
 :param dist_col: The naame of the column to use distances from.
 :param num_axis: The number of axis i.e the dimensionality

   Calculate the sum of squares using [calc_tot_withins]
   and [calc_sum_sq_clusts], this is "l".
   "k" is the number of cluster groups.
   "d" is the dimensionality = num_axis
   "n" is the number of points
   :math:`aic = d + 2*k*l`
   :math:`bic = d+log(n)*k*l`
   out = {"aic" : aic, "bic" : bic}.

 :return: out
 :rtype: dictionary of 2 floats.
"""
def calc_clust_ic (clust, group_col, dist_col, num_axis) :
  # Number of estimated parameters
  k = len(clust.clust_num.uniqe())
  # Number of data points
  n = len(clust.axes[0])
  # Dimension of each data point
  d = num_axis
  # Log-Likelihood term
  l = calc_tot_within(clust, group_col = group_col, dist_col = dist_col)
  aic_n = d + 2 * k * l
  bic_n = d + log(n) * k * l
  ic_out = {"aic" : aic_n, "bic" : bic_n}
  return(ic_out)