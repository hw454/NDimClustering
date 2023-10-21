""" Find the principal components (PCs) and cluster the data.
 Scores the clusters on the PCs

 :param data_matrices: A list containing
   * "beta" matrix corresponding to the data
   * "pval" probabilities associated with "beta" values.
   * "se" standard errors associated with "beta" values.
   * "trait_info" dataframe of all the traits.
 :param out_pheno: Label for the outcome phenotype
 :param thresholds: List for the  threshold related variables.
   * "threshmul" is multiplied by the data variance for cluster
   difference threshold.
   * "diff" is the the cluster difference threshold once found.
   * "clust" is the required distance between cluster centres for a
   cluster to be considered converged.
 :param na_handling: List of the terms related to handling NaNs
   * "narm" - TRUE (default) if NaNs are to be removed in calculations.
   * "percent" - percentage of a column which must be not NaN. default\:0.95
 :param iter_traits: list of terms indicating the type of program.
   * "iter" - default\:0
   * "bp_on" - default\:TRUE switch if pval scores is to be used.
   * "clust_prob_on" - default\:TRUE switch to use prob of being in cluster
   * "clust_typ" - default="basic", the clustering method to use.
 :param norm_typs: List of norm types
   * "clust" - default="F"
   * "thresh" - default="F"
 :param nums: List of important numbers.
   * "max_dist" - Maximum distance between points, default=1
   * "np" - Number of principal components
   * "nc" - Number of clusters.

   out_dict is the dictionary of the clusters, the principal component matrices,
   cluster scores, and the maximum difference between clusters.
   Generated using :func:`clust_pca_compare_single`.
   
 :return: out_dict"""
def clust_pca_all (data_matrices,out_pheno,
                   thresholds = {"threshmul" : 5, "diff" : 1e-5,"clust" : 1e-5},
                   iter_traits = {"iter" : 0, "bp_on" : True, "clust_prob_on" : True, "clust_typ" : "basic"},
                   na_handling = {"narm" : True, "percent" : 0.95},
                   norm_typs = {"clust" : "fro", "thresh" : "fro"},
                   nums = {"max_dist" : 1, "np" : 3, "nr" : 5} ):
  # Data frame for recording the cluster scores.
  c_scores = pd.DataFrame( columns = ["num_axis", "clust_num" ])
  # Add np columns for each PC
  c_scores = add_np_cols(c_scores, nums["np"])
  # Initialise with outcome
  trait_ind = np.where(data_matrices["trait_info"]["phenotype"] == out_pheno)[1] 
  trait_df = pd.DataFrame(data = {
                       "num_axis" : 1,
                       "label" : out_pheno,
                       "axes_ind" : trait_ind
                       })
  max_df = pd.DataFrame(columns = ["num_axis", "cn1", "cn2", "diff"])
  cluster_df = pd.DataFrame(columns = ["snp_id", "clust_num", "clust_prob", "clust_dist", "num_axis"])
  cluster_df = add_np_cols(cluster_df, nums["np"])
  clust_dist_df = pd.DataFrame(columns = ["snp_id" ])
  clust_dist_df = add_nclust_cols(clust_dist_df, nums["nr"])
  df_dict = {"clust_scores" : c_scores,
              "max_diff" : max_df,
              "e_list" : [],
              "b_pc_list" : [],
              "se_pc_list" : [],
              "trait" : trait_df,
              "clust_items" : cluster_df,
              "clust_membership" : clust_dist_df}
  # Fill trait_df with valid traits before running main program.
  trait_df = make_trait_df(pheno_list = data_matrices["trait_info"]["phenotype"],
                           data_mat = data_matrices["b_df"],
                           na_percent = na_handling["percent"])
  num_axis = len(trait_df)
  trait_df["num_axis"] = num_axis
  df_dict["trait"] = trait_df
  out_dict = clust_pca_compare_single(df_list = df_list,
                                      iter_traits = iter_traits,
                                      num_axis = num_axis,
                                      data_matrices = data_matrices,
                                      na_handling = na_handling,
                                      thresholds = thresholds,
                                      norm_typs = norm_typs,
                                      nums = nums)
  return(out_dict)