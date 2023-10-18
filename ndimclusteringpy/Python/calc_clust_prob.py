""" Calculate the probability the snp is in the cluster

 :math:`p=\frac{1.0}{1.0+d}`

 :param d: The snp distance to the cluster centre

 :return: p
 :retype: number or array same dimension as d"""
def calc_clust_prob (d) :
  dist = 1.0 / (1.0 + d)
  return(dist)