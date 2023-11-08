""" Cluster the data using kmeans then minimising aic.

 :param b_df: The association dataframe
 :param nclust: The maximum number of clusters to consider. default=10
 :param max_dist: The maximum distance between any two points.
 :param space_typ: The spatial structure of the data when clustering.
   If "angle" then the data is converted to angles. If "regular" the
   data is left the same.
 :param clust_prob_on: Bool switch. If TRUE (default) then cluster
   probability if used to weight the cluster scores.
 :param norm_typ: The type of norm to be used for distance calculations.
   The default is the Froebenius norm "F".
 :param threshold: The threshold for distance between cluster centres
   for clusters to be considered converged.
 :param narm: Bool to indicate how to handle NaNs. If TRUE (default)
   NaNs are ignored in calculations.

   If space_typ == "angle" then data is converted to angles.
   for i=1:(nr+1)
      The points are clustered using [km_nan] and i clusters
      find the aic for each cluster using [find_ic].
   The set of clusters which minimises the aic is chosen.
   The "clusters_df" dataframe of the cluster membership with columns\:
     * "clust_num" the number of clusters
     * "clust_dist" the distance from the snp to the cluster centre
       (or distance between angles)
     * "clust_prob" probability the snp is in the cluster. Calculated
       using [calc_clust_prob].

 :return: clusters_df"""
def cluster_kmeans_min (data_list,
                              nclust = 10,
                              max_dist = 10.0,
                              space_typ = "regular",
                              clust_prob_on = TRUE,
                              norm_typ = "F",
                              threshold = 1e-5,
                              narm = True) :
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

  if (space_typ == "angle") :
    # For each point in b_df_comp convert the score to the angle
    # between the vectors to the origin and the unit vectors on the axis.
    b_df_clust <- mat_to_angle_mat(b_df)
  else :
    b_df_clust <- b_df

  # Crop the data to focus on the snps with the lowest standard error
  b_df_crop <- b_df_clust.loc[crop_snp_list]

  # Initial cluster dataframe
  clust_re_list = [ km_nan(c_num,
                    b_mat = b_mat_crop,
                    clust_threshold = threshold,
                    norm_typ = norm_typ,
                    na_rm = narm,
                    prob_on = clust_prob_on)
                    for c_num in range(1, nclust+1)]
  ic_list = [find_all_ic(c_l, num_axis = ncol(b_mat)) for c_l in clust_re_list]
  ic_df = pd.concat(ic_list)
  # Find the number of centres that minimizes the AIC
  min_row = ic_df["aic"].index(ic_df["aic"].min())
  min_cents = ic_df["ncents"].loc[min_row]
  centroids_df = clust_re_list[min_cents]["centres"]
  # Use the number of centres to locate corresponding clusters since there
  # maybe variations due to machine precision in the aic values.
  cluster_df = ic_df.loc[min_row]
  cluster_df.set_index("snp_id")
  # cluster number identification for each observation
  snp_cluster_list = [ find_closest_clust_snp(snp,
                            b_mat = b_mat_clust,
                            cluster_df = clust_out["clusters"],
                            centroids_df = clust_out["centres"],
                            norm_typ = norm_typ) 
                            for snp in b_mat_clust.index().difference(crop_snp_list)]
  nan_cluster_df = pd.concat(snp_cluster_list)
  if (clust_prob_on) :
    nan_cluster_df["clust_prob"] = calc_clust_prob(nan_cluster_df["clust_dist"])
  cluster_df = pd.concat([cluster_df, nan_cluster_df])
  # ADDFEATURE - Assign junk clusters.
  clust_out = {"clusters" : cluster_df, "centres": centroids_df}
  return(clust_out)