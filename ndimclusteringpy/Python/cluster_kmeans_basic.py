""" Cluster data with a standard kmeans approach

 :param b_df: Dataframe of data. Rows correspond to snps.
 :param nclust: Number of clusters to allocate
 :param space_typ: Whether to uses angles of spatial co-ordinates.
   If "angle" use angles found using [convert_point_mat_to_angle]
 :param clust_prob_on: Bool switch. If TRUE then calculate the snps
   probability of being in the cluster using [calc_clust_prob]
 :param norm_typ: The type of norm to use in distance calculations.
   The default is the Froebenius norm "F".
 :param threshold: The threshold with which clusters centres must differ
   by for clusters to be considered converged.
 :param narm: Boolean swithc. If TRUE (default) ignore NaNs in calculations.

   If space_typ == "angle" then data is converted to angles.
   The points are clustered using [km_nan] and "nclust" clusters.
   The "cluster_df" dataframe labelled
     "clusters" of the cluster membership with columns\:
     * "clust_num" the number of clusters
     * "clust_dist" the distance from the snp to the cluster centre
       (or distance between angles)
     * "clust_prob" probability the snp is in the cluster. Calculated
       using [calc_clust_prob].

 :return: clusters_df"""
def cluster_kmeans_basic (data_dict,
                          nclust = 10,
                          space_typ = "regular",
                          clust_prob_on = True,
                          norm_typ = "fro",
                          threshold = 1e-5,
                          narm = True):
  # Using the association scores for each SNP accross traits cluster the traits
  # using kmeans. Return the cluster setup which minimises AIC.

  b_df = data_dict["beta"]
  se_df = data_dict["se"]
  # Crop the data to those with the lowest standard error.
  # Add remaining terms to closest cluster once centres have been found.
  crop_lim = [se_df[c].quantile(0.25) for c in se_df.columns()]
  crop_terms = [se_df[c].between(0, crop_lim[i]) for i, c in enumerate(se_df.columns())]
  crop_terms = pd.concat([crop_terms], axis=1)
  crop_se = se_df[crop_terms]
  crop_se = crop_se.dropna(axis=0, how = "any")
  crop_snp_list = crop_se.index

  # Filter NaNs before clustering
  if (space_typ == "angle") :
    # For each point in b_df_comp convert the score to the angle
    # between the vectors to the origin and the unit vectors on the axis.
    b_df_clust = mat_to_angle_mat(b_df)
  else :
    b_df_clust = b_df

  # Crop the data to focus on the snps with the lowest standard error
  b_df_crop = b_df_clust.loc[crop_snp_list]

  # Initial cluster dataframe
  clust_out = km_nan(b_df_crop,
                    nclust = nclust,
                    clust_threshold = threshold,
                    norm_typ = norm_typ,
                    prob_on = clust_prob_on,
                    na_rm = narm)
  # cluster number identification for the snps with higher standard error.
  snp_cluster_list =  [find_closest_clust_snp(row,
                       b_mat = b_mat_clust,
                       cluster_df = clust_out["clusters"],
                       centroids_df = clust_out["centres"],
                       norm_typ = norm_typ) 
                       for row in b_df_clust.index().difference(crop_snp_list.index())]
  nan_cluster_df = pd.concat(snp_cluster_list)
  if (clust_prob_on) :
    nan_cluster_df["clust_prob"] = calc_clust_prob(nan_cluster_df["clust_dist"])
  # ADDFEATURE - Assign junk clusters.
  clust_out["clusters"] = pd.concat([clust_out["clusters"], nan_cluster_df])
  return(clust_out)