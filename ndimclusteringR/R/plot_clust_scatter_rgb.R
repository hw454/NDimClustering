#' Plot the scatter plot or snps and principal components.
#'
#' @description
#' Colour by cluster number.
#' The x and y axis are the assocation with the principal component
#' The width and height of the errorbars are given by the standard error.
#'
#' @param clust_dist_df distances between each snp and cluster centres
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
  ignore_cols <- c("num_axis")
  crop_clust_dist_df <- clust_dist_df[clust_dist_df$num_axis == num_axis, ]
  crop_clust_dist_df <- tibble::column_to_rownames(crop_clust_dist_df, "snp_id")
  full_trait_list <- colnames(crop_clust_dist_df)
  full_trait_list <- full_trait_list[!(full_trait_list %in% ignore_cols)]
  crop_clust_dist_df <- crop_clust_dist_df[, full_trait_list]
  pnme <- paste0(iter_traits$res_dir, "clusters_pc_rgb_num_axis", num_axis, ".png")
  c1 <- colnames(b_mat)[1]
  c2 <- colnames(b_mat)[2]
  se_max <- apply(se_mat, 2, max)
  se_min <- apply(se_mat, 2, min)
  norm_se <- rowSums((se_mat - se_min) / (se_max - se_min))
  alpha_vec <- 1.0 / (1.0 + norm_se)
  snp_list <- row.names(b_mat)
  max_dist <- apply(crop_clust_dist_df, 2, max)
  min_dist <- apply(crop_clust_dist_df, 2, min)
  norm_dist_df <- apply(crop_clust_dist_df, 2, scales::rescale)
  clust_names <- colnames(norm_dist_df)
  colour_vec <- paste0("rgb(",norm_dist_df[, clust_names[1]],
                    norm_dist_df[, clust_names[2]],
                    norm_dist_df[, clust_names[3]],")")
  colour_vec <- factor(colour_vec)
  # vector with color values
  my_col_vec <- levels(colour_vec)
  my_col_vec <- sapply(seq_along(my_col_vec),
                  function(i) eval(parse(text = my_col_vec[i])))
  res_df <- data.frame(
    row.names = snp_list,
    bx = b_mat[, c1],
    by = b_mat[, c2],
    bxse = se_mat[, c1],
    byse = se_mat[, c2],
    cols = colour_vec,
    alp = alpha_vec
  )
  print(head(res_df))
  print(head(my_col_vec))
  ggplot2::ggplot(data = res_df,
                  ggplot2::aes(x = bx, y = by)) + # nolint: object_usage_linter.
  ggplot2::geom_point(ggplot2::aes(
                  col = cols # nolint: object_usage_linter.
                 ), shape = 21, show.legend = FALSE) + # nolint: object_usage_linter.
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = res_df$bx - 1.96 * res_df$bxse,
                 xmax = res_df$bx + 1.96 * res_df$bxse,
                 col = cols), linetype = "solid", show.legend = FALSE) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = res_df$by - 1.96 * res_df$byse,
                 ymax = res_df$by + 1.96 * res_df$byse,
                 col = cols), linetype = "solid", show.legend = FALSE) +
  ggplot2::scale_color_manual(values = my_col_vec) +
  ggplot2::ylab("Association with PC2") +
  ggplot2::xlab("Association with PC1") +
  ggplot2::ggtitle("Clustered by principal components")
  ggplot2::ggsave(filename = pnme, width = pw, height = ph)
  return()
}
