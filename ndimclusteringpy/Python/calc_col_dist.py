import pandas as pd
""" Find the distance between a snp and all other snps.

 :param score_mat: Matrix of scores for each snp
 :param snp1: the snp name to compare all others to
 :param norm_typ: the type of norm to use in distance calculations.
   default is the Froebenius norm "fro".

   Use [calc_pair_dist_df] to find the distance between pairs of snps.
   The results from this for each pair is bound on the rows using
   "rbind" into the datframe "dist_df".

 :return: dist_df
"""
def calc_col_dist (score_mat, snp1, norm_typ = "fro") :
  dist_list = [calc_pair_dist_df(row,
                                  score_mat = score_mat,
                                  snp1 = snp1,
                                  norm_typ = norm_typ) for row in score_mat.index]
  dist_df = pd.concat(dist_list)
  return(dist_df)