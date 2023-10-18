import matplotlib.pyplot as mp
""" Plot the scatter plot or snps and principal components.

 Colour by cluster number.
 The x and y axis are the assocation with the principal component
 The width and height of the errorbars are given by the standard error.

 :param clust_dist_df: distances between each snp and cluster centres
 :param b_mat: The matrix of the score data
 :param se_mat: The matrix of the standarad errors associated with the scores.
 :param iter_traits: The iteration variables for the type of iteration.
 :param num_axis: The number of trait axis, default:0
 :param pw: The plot width, default:8
 :param ph: The plot heigh, default:4
 :returns: empty
"""
def rgb_to_hex(r, g, b):
    return "#{:02x}{:02x}{:02x}".format(r, g, b)

def plot_clust_scatter_rgb (clust_dist_df, b_mat,
                          se_mat,
                          iter_traits,
                          num_axis = 1,
                          pw = 8,
                          ph = 4):
  ignore_cols = ["num_axis"]
  crop_clust_dist_df = clust_dist_df[clust_dist_df.num_axis == num_axis, ]
  crop_clust_dist_df.set_index("snp_id")
  full_trait_list = crop_clust_dist_df.columns
  full_trait_list = [f for f in full_trait_list if f not in ignore_cols]
  crop_clust_dist_df = crop_clust_dist_df[: , full_trait_list]
  pnme = str(iter_traits["res_dir"] +
                "clusters_pc_rgb_num_axis" +
                num_axis +
                ".png")
  c1 == b_mat.columns[0]
  c2 = b_mat.columns[1]
  se_max = se_mat.max(axis=1)
  se_min = se_mat.min(axis=1)
  norm_row = (se_mat - se_min) / (se_max - se_min)
  norm_se = norm_row.sum(axis=0)
  alpha_vec = 1.0 / (1.0 + norm_se)
  snp_list = b_mat.index
  max_dist = crop_clust_dist_df.max(axis = 1)
  min_dist = crop_clust_dist_df.min(axis = 1)
  norm_dist_df = (crop_clust_dist_df - min_dist) / (max_dist - min_dist)
  clust_names = norm_dist_df.columns
  colour_vec = reg_to_hex(norm_dist_df[: , clust_names[1]],
                    norm_dist_df[: , clust_names[2]],
                    norm_dist_df[: , clust_names[3]])
  res_df = pd.Dataframe(
    index = snp_list,
    data ={
    "bx" : b_mat[:, c1],
    "by" : b_mat[:, c2],
    "bxse" : se_mat[:, c1],
    "byse" : se_mat[:, c2],
    "cols" : colour_vec,
    "alp" : alpha_vec
    }
  )
  res_df.set_index("rsid")

  ax = res_df.plot.scatter(x = bx, y = by, 
                           xerr = bxse, yerr= byse,
                           c = cols, 
                           title = "Clustered by principal components")
  ax.set_ylabel("Association with PC2")
  ax.set_xlabel("Association with PC1")
  
  res.savefig(pnme)
  return()
