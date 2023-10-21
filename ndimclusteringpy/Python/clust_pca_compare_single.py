""" Single iteration for iterating through trait list.

 :param df_dict: the results from previous iterations to append to.
 :param data_matrices: A dictionary containing
   * "beta" main data matrix
   * "pval" probabilities associated with "beta" values
   * "se"  standard errors associated with "beta" values.
   * "trait_info" data frame of all the traits.
 :param num_axis: Number of trait axes in use.
 :param thresholds: Dictionary of threshold related variables.
   * "threshmul" is multiplied by the data variance for cluster
   difference threshold.
   * "diff" is the the cluster difference threshold once found.
   * "clust" is the required distance between cluster centres for a
   cluster to be considered converged.
 :param na_handling: Dictionary of the terms related to handling NaNs
   * "narm" - TRUE (default) if NaNs are to be removed in calculations.
   * "percent" - percentage of a column which must be not NaN. default\:0.95
 :param iter_traits: Dictionary of terms indicating the type of program.
   * "iter" - integer, (default 0)
   * "bp_on" - TRUE (default) if probability of scores is to be used.
   * "clust_prob_on" - TRUE\:default switch for using prob of being in cluster
   * "clust_typ" - default\:"basic", the clustering method to use.
 :param norm_typs: Dictionary of norm types
   * "clust" - default:"fro"
   * "thresh" - default:"fro"
 :param nums: Dictionary of important numbers.
   * "max_dist" - Maximum distance between points, default=1
   * "np" - Number of principal components
   * "nc" - Number of clusters.

   df_list is the list of of the clusters, the principal component matrices,
   cluster scores, and the maximum difference between clusters.
   Each iteration generated using [clust_pca_compare_single].

 :return: df_list"""
def clust_pca_compare_single (df_list, iter_traits, num_axis,
                              data_matrices,
                              na_handling,
                              thresholds,
                              norm_typs,
                              ums) :
  # Extract the trait_df dataframe from df_list
  trait_df = df_list["trait"]
  # If the trait is not all NaN then run clustering.
  print("PCA on", num_axis, "axes")
  # Get the data upto this axis
  b_iter_df = data_matrices["beta"][trait_df["label"]]
  se_iter_df = data_matrices["se"][trait_df["label"]]
  pval_iter_df = data_matrices["pval"][trait_df["label"]]
  # Cluster the data on these axes
  pca_dict = find_principal_components(b_iter_mat, pval_iter_mat, se_iter_mat,
                          nums["np"], na_handling["narm"])
  b_pc_df = pca_list["beta"]
  p_pc_df = pca_list["pval"]
  t_df    = pca_list["transform"]
  # Store the matrices for result output
  df_dict["e_mat"] = t_df
  df_dict["b_pc"] = b_pc_df
  df_dict["se_pc"] = pca_dict["se"]
  # Get column names for PCs
  pc_cols = b_pc_df.columns()
  # Cluster the data on these axes
  if ("angle" in iter_traits["clust_typ"]):
    st = "angle"
  else:
    st = "regular"
  if ("min" in iter_traits["clust_typ"]) :
    cluster_out <- cluster_kmeans_min(pca_dict,
                                      nums["nr"],
                                      space_typ = st,
                                      clust_prob_on = iter_traits["clust_prob_on"],
                                      norm_typ = norm_typs["clust"],
                                      threshold = thresholds["clust"],
                                      narm = na_handling["narm"])
  elif ("basic" in iter_traits["clust_typ"]):
    cluster_out <- cluster_kmeans_basic(pca_list,
                                        nums["nr"],
                                        space_typ = st,
                                        clust_prob_on = iter_traits["clust_prob_on"], # nolint
                                        threshold = thresholds["clust"],
                                        norm_typ = norm_typs["clust"],
                                        narm = na_handling["narm"])
  cluster_df = cluster_out["clusters"]
  centroids_df = cluster_out["centres"]
  # Calculate the distance to all the cluster centres
  clust_dist_df = calc_clust_dist(df_list["b_pc"], centroids_df)
  # Find the set of cluster numbers
  c_nums = cluster_df["clust_num"].unique()
  # Score the clustered data based on affiliations with axes.
  # Find the score for each PC
  c_score0 = score_all_clusters(cluster_df,
                                beta_df = b_pc_df,
                                pval_df = p_pc_df,
                                bp_on = iter_traits["bp_on"],
                                clust_prob_on = iter_traits["clust_prob_on"],
                                num_axis = num_axis)
  # Iterate through each cluster and compare across the others to find if
  # any pair have a distinct difference.
  diff_score_list = [compare_oneclust_tolist(c_n,
          c_nums = c_nums,
          c_score0 = c_score0,
          axis = pc_cols,
          clust_norm = norm_typs["clust"]) for c_n in c_nums]
  diff_score_list = diff_score_list[not diff_score_list.isna()]
  diff_scores = pd.concat(diff_score_list)
  # Find the pair with maximum cluster score diff and store in max_df0
  row = diff_scores["diff"].index(diff_scores["diff"].max())
  max_df0 = diff_scores.loc[row]
  max_df0["num_axis"] = num_axis
  c_score0["num_axis"] = num_axis
  clust_dist_df["num_axis"] = num_axis
  cluster_df["num_axis"] = num_axis
  clust_dist_df = clust_dist_df.rename_axis("snp_id").reset_index()
  cluster_df = cluster_df.rename_axis("snp_id").reset_index()
  df_list["clust_items"] = pd.concat([df_list["clust_items"], cluster_df])
  df_list["max_diff"] = pd.concat([df_list["max_diff"], max_df0])
  df_list["clust_scores"] = pd.concat([df_list["clust_scores"], c_score0])
  df_list["clust_membership"] = pd.concat([df_list["clust_membership"], clust_dist_df],
                                            join = "outer")
  return(df_list)