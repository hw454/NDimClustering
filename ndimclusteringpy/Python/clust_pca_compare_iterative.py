""" Cluster data on the principal components and score the clusters on
 association. Compare the cluster scores and iteratively add a
 trait each time.

 :param data_matrices: A list containing the matrices
   * "beta" dataframe corresponding to the association data.
   * "pval" the dataframe of probabilities associated with each score.
   * "se" dataframe of standard errors associated with each score.
   * "trait_info" which is the dataframe of all the traits.
 :param out_pheno: Label for the outcome phenotype
 :param thresholds: Dictionary for the three threshold related variables.
   * "threshmul" is multiplied by the data variance for cluster
   difference threshold.
   * "diff" is the the cluster difference threshold once found.
   * "clust" is the required distance between cluster centres for a
   cluster to be considered converged.
 :param na_handling: Dictionary of the terms related to handling NaNs
   * "narm" - TRUE (default) if NaNs are to be removed in calculations.
   * "percent" - percentage of a column which must be not NaN. default\:0.95
 :param iter_traits: dictionary of terms indicating the type of program.
   * "iter" - integer, (default 0)
   * "bp_on" - TRUE (default) if probability of scores is to be used.
   * "clust_prob_on" - TRUE (default) if prob of being in cluster to be used.
   * "clust_typ" - default="basic", the clustering method to use.
 :param norm_typs: Dictionary of norm types
   * "clust" - default = "fro"
   * "thresh" - default = "fro"
 :param nums: Dictionary of important numbers.
   * "max_dist" - Maximum distance between points, default=1
   * "np" - Number of principal components
   * "nc" - Number of clusters.

   df_list is the list of of the clusters, the principal component matrices,
   cluster scores, and the maximum difference between clusters.
   Each iteration generated using [clust_pca_compare_single].

 :return: df_dict"""
def clust_pca_compare_iterative (data_matrices, out_pheno,
                         thresholds = {"threshmul" : 5, "diff" : 1e-5, "clust" : 1e-5},
                         na_handling = {"narm" : True, "percent" : 0.95},
                         iter_traits = {"iter" : 0, "bp_on" : True, "clust_prob_on" : True, "clust_typ" : "basic"},
                         norm_typs = {"clust" : "fro", "thresh" : "fro"},
                         nums = {"max_dist" : 1, "np" : 3, "nr" : 5}
                        ) :
  # Iterate through the traits in trait_info contained in the
  # data_matrices list. Find the principal components.
  # Transform the data onto these.
  # Cluster the transformed data at each iteration.
  # If there is a distinct difference between two clusters exit.
  trait_df = pd. DataFrame(columns = ["label", "axes_ind"])
  # Data frame for recording the cluster scores.
  c_scores = pd.DataFrame(columns = ["num_axis", "clust_num"])
  # Add np columns for each PC
  c_scores = add_np_cols(c_scores, nums["np"])
  # Initialise with outcome
  ax_ind = which(data_matrices["trait_info"]["phenotype"] == out_pheno)[1]
  trait_df = pd.DataFrame(data = {
                          "label" : out_pheno,
                          "axes_ind" :  ax_ind})
  max_df = pd.DataFrame(columns = ["num_axis", "cn1", "cn2", "max_diff"])
  cluster_df = pd.DataFrame(columns = ["snp_id",
                          "clust_num",
                          "clust_prob",
                          "clust_dist",
                          "num_axis",
                          "ncents"])
  cluster_df = add_np_cols(cluster_df, nums["np"])
  clust_dist_df = pd.DataFrame(columns=["snp_id"])
  clust_dist_df = add_nclust_cols(cluster_df, nums["nr"])
  df_list = {"clust_scores" : c_scores,
             "max_diff" : max_df,
             "e_list" : [],
             "trait" : [],
             "cluster_items" : cluster_df,
             "cluster_membership" : clust_dist_df}
  pheno_list = data_matrices["trait_info"]["phenotype"]
  na = len(pheno_list)
  for ai in range(na):
    # Update the traits forming the axis
    # Add the trait to the trait dataframe
    a = pheno_list[ai]
    print("New trait on axis", a)
    covered = (a in trait_df["label"])
    allna = check_col_na(data_matrices["beta"][a], na_handling["percent"])
    if (not covered and not allna) :
      # Add trait to trait dataframe
      a_ind = np.where(pheno_list == a)[0]
      trait_row = pd.DataFrame(data = {
                              "label" : a,
                              "axes_ind" : a_ind})
      trait_df = pd.concat([trait_df, trait_row])
      df_list["trait"] = pd.concat([df_list["trait"], trait_df])
      df_list = clust_pca_compare_single(df_list,
                                    iter_traits = iter_traits,
                                    num_axis = ai,
                                    data_matrices = data_matrices,
                                    na_handling = na_handling,
                                    thresholds = thresholds,
                                    norm_typs = norm_typs,
                                    nums = nums)
  return(df_list)