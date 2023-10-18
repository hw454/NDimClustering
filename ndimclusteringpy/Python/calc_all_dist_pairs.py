""" Calculate the distance between each snp and every other snp.

 Using [calc_col_dist] and [calc_pair_dist]
 to find the distance between each snp and every other snp.

 :param score_mat: the data matrix. Each row is a snp and the values
    are the co-ordinates of the snp accross the trait axis.
 :param norm_typ: the type of norm to use in the distance calculations
    when comparig the snps. Default is the Froebenius norm "fro".

 "dist_df" is a dataframe with columns\:
   * "snp1" snp label
   * "snp2" snp label for the compared snp
   * "dist" the distance between snp1 and snp2

 :return dist_df:
 :rtype: dataframe
"""
def calc_all_dist_pairs (score_mat, norm_typ = "fro"):
  # Find the distance between all pairs of points.
  dist_df = pd.Dataframe({
    "snp1" : pd.Series(dtype = "str"),
    "snp2" :  pd.Series(dtype = "str"),
    "dist" :  pd.Series(dtype = "float")}
  )
  dist_list = [calc_col_dist(row,
                      score_mat = score_mat,
                      norm_typ = norm_typ) for row in score_mat.index]

  dist_df <- pd.concat(dist_list)
return(dist_df)