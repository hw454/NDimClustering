#' Plot the scatter plot or snps and principal components.
#'
#' @description
#' Colour by cluster number.
#' The x and y axis are the assocation with the principal component
#' The width and height of the errorbars are given by the standard error.
#'
#' @param cluster_df The dataframe of cluster_df membership
#' @param b_mat The matrix of the score data
#' @param se_mat The matrix of the standarad errors associated with the scores.
#' @param iter_traits The iteration variables for the type of iteration.
#' @param num_axis The number of trait axis, default\:0
#' @param pw The plot width, default\:8
#' @param ph The plot heigh, default\:4
#'
#' @export
plot_clust_scatter <- function(cluster_df, b_mat,
                          se_mat,
                          iter_traits,
                          num_axis = 2,
                          pw = 8,
                          ph = 4) {
  crop_cluster_df <- cluster_df[cluster_df$num_axis == num_axis, ]
  pnme <- paste0(iter_traits$res_dir, "clusters_num_axis", num_axis, ".png")
  c1 <- colnames(b_mat)[1]
  c2 <- colnames(b_mat)[2]
  se_bar <- colMeans(se_mat)
  norm_se <- rowSums(se_mat / se_bar)
  alpha_vec <- 1.0 / (1.0 + norm_se)
  snp_list <- row.names(b_mat)
  res_df <- data.frame(
    row.names = snp_list,
    bx = b_mat[, c1],
    by = b_mat[, c2],
    bxse = se_mat[, c1],
    byse = se_mat[, c2],
    clust_num = cluster_df[snp_list, "clust_num"],
    clust_prob = cluster_df[snp_list, "clust_prob"],
    alp = alpha_vec
  )
  ggplot2::ggplot(data = res_df,
                  ggplot2::aes(x = bx, y = by)) +
  ggplot2::geom_point(ggplot2::aes(colour = clust_num,
                  size = clust_prob,
                  alpha = alp), shape = 21) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = res_df$bx - 1.96 * res_df$bxse,
                 xmax = res_df$bx + 1.96 * res_df$bxse,
                 color = clust_num,
                 alpha = alp), linetype = "solid") +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = res_df$by - 1.96 * res_df$byse,
                 ymax = res_df$by + 1.96 * res_df$byse,
                 color = clust_num,
                 alpha = alp), linetype = "solid") +
  ggplot2::ylab("Association with PC2") +
  ggplot2::xlab("Association with PC1") +
  ggplot2::ggtitle("Clustered by principal components") 
  ggplot2::ggsave(filename = pnme, width = pw, height = ph)
  return()
}
