""" Calculate the total sum of the cluster distances divided by the variance

 :param data: matrix of data
 :param group_col: The column indicate terms to go together
   default "clust_num"
 :param dist_col: The name of the column containing the term distances.
   default "clust_dist"

   Calculate the sum of the squares for each cluster group using
   :func:`calc_sum_sq_clusts`

 :return: The sum of the clusters sum of squares/ divided by the variance
   of the distance data."""
def calc_tot_within (data,group_col = "clust_num", dist_col = "clust_dist"):
  ll =[ calc_sum_sq_clusts(clust,
                    data = data,
                    group_col = group_col,
                    dist_col = dist_col) 
                    for clust in unique(data[:, group_col])]
  sig = data[dist_col].var()
  l = sum(ll) / sig
  return(l)