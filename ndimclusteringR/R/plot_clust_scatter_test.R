#' Plot the scatter plot or snps and principal components.
#'
#' @description
#' Colour by cluster number.
#' The x and y axis are the assocation with the principal component
#' The width and height of the errorbars are given by the standard error.
#'
#' @param cluster_df The dataframe of cluster_df membership
#' @param b_mat The matrix of the score data
#' @param iter_traits The iteration variables for the type of iteration.
#' @param c1 The column for the x-axis
#' @param c2 The column for the y-axis
#' @param num_axis The number of trait axis, default\:0
#' @param pw The plot width, default\:8
#' @param ph The plot heigh, default\:4
#' @param save_suffix String to add to plot save name. default \: ""
#'
#' @export
plot_clust_scatter_test <- function(cluster_df, b_mat,
  iter_traits, c1, c2,
  num_axis = 2,
  pw = 8,
  ph = 4,
  save_suffix = ""
) {
  crop_cluster_df <- cluster_df[cluster_df$num_axis == num_axis, ]
  crop_cluster_df <- tibble::column_to_rownames(crop_cluster_df, "snp_id")
  snp_list <- rownames(b_mat)
  clust_nums <- crop_cluster_df[snp_list, "clust_num"]
  res_df <- data.frame(
    row.names = snp_list,
    bx = b_mat[, c1],
    by = b_mat[, c2],
    clust_num = clust_nums
  )
  np <- iter_traits$num_paths + 1
  pnme <- paste0(iter_traits$res_dir,
                 "scatter_withclusts",
                 c1,
                 "_vs_",
                 c2,
                 "_pctype",
                 iter_traits$pc_type,
                 "_numpaths",
                 np,
                 "_clust",
                 iter_traits$clust_typ,
                 "_howcents",
                 iter_traits$how_cents,
                 ".png")
  title_str <- paste("Clusters plotted against the", c1, "and", c2, "traits.")
  caption_str <- paste("Test case with", np,
    "pathways. \n The method for PCA used is",
    iter_traits$pc_type, ".",
    "\n The method for allocating centroids is",
    iter_traits$how_cents,
    "\n The clustering method is", iter_traits$clust_typ
  )
  clust_scatter <- ggplot2::ggplot(data = res_df,
    ggplot2::aes(x = bx, y = by) # nolint: object_usage_linter.
  ) +
    ggplot2::geom_point(ggplot2::aes(
      color = clust_num, # nolint: object_usage_linter
    ),
    shape = 21) + # nolint: object_usage_linter.
    ggplot2::ylab(paste("Association with", c2)) +
    ggplot2::xlab(paste("Association with", c1)) +
    ggplot2::ggtitle(title_str) +
    ggplot2::labs(caption = caption_str)
  ggplot2::ggsave(filename = pnme, width = pw, height = ph)
  print(clust_scatter)
  return(clust_scatter)
}
