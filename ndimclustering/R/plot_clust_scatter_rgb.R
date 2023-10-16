#' Plot the scatter plot or snps and principal components.
#'
#' @description
#' Colour by cluster number.
#' The x and y axis are the assocation with the principal component
#' The width and height of the errorbars are given by the standard error.
#'
#' @param clust_dist_df The dataframe of distances between each snp and the cluster centres
#' @param b_mat The matrix of the score data
#' @param se_mat The matrix of the standarad errors associated with the scores.
#' @param iter_traits The iteration variables for the type of iteration.
#' @param num_axis The number of trait axis, default\:0
#' @param pw The plot width, default\:8
#' @param ph The plot heigh, default\:4
#'
#' @export
plot_clust_scatter_rgb <- function(clust_dist_df, b_mat,
                          se_mat,
                          iter_traits,
                          num_axis = 1,
                          pw = 8,
                          ph = 4) {
  crop_clust_dist_df <- clust_dist_df[clust_dist_df$num_axis==num_axis, ]
  pnme <- paste0(iter_traits$res_dir, "clusters_rgb_num_axis", num_axis, ".png")
  c1 <- colnames(b_mat)[1]
  c2 <- colnames(b_mat)[2]
  se_bar <- colMeans(se_mat)
  norm_se <- rowSums(se_mat / se_bar)
  alpha_vec <- 1.0 / (1.0 + norm_se)
  snp_list <- row.names(b_mat)
  norm_dist <- scale(crop_clust_dist_df)
  clust_names <- colnames(crop_clust_dist_df)
  colour_vec <- rgb(norm_dist[, clust_names[1]],
         norm_dist[, clust_names[2]], norm_dist[, clust_names[3]])
  res_df <- data.frame(
    row.names = snp_list,
    bx = b_mat[, c1],
    by = b_mat[, c2],
    bxse = se_mat[, c1],
    byse = se_mat[, c2],
    cols = colour_vec,
    alp = alpha_vec
  )
  ggplot2::ggplot(data = res_df,
                  ggplot2::aes(x = bx, y = by)) +
  ggplot2::geom_point(ggplot2::aes(
                  colour = cols,
                  size = clust_prob,
                  alpha = alp), shape = 21) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = res_df$bx - 1.96 * res_df$bxse,
                 xmax = res_df$bx + 1.96 * res_df$bxse,
                 colour = cols,
                 alpha = alp), linetype = "solid") +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = res_df$by - 1.96 * res_df$byse,
                 ymax = res_df$by + 1.96 * res_df$byse,
                 colour = cols,
                 alpha = alp), linetype = "solid") +
  ggplot2::ylab("Association with PC2") +
  ggplot2::xlab("Association with PC1") +
  ggplot2::ggtitle("Clustered by principal components") 
  ggplot2::ggsave(filename = pnme, width = pw, height = ph)
  return()
}
