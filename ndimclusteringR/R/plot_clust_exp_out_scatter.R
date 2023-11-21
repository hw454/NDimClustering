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
#' @param exp_pheno The label for the exposure phenotype
#' @param out_pheno The label for the output phenotype
#' @param num_axis The number of trait axis, default\:0
#' @param pw The plot width, default\:8
#' @param ph The plot heigh, default\:4
#'
#' @export
plot_clust_exp_out_scatter <- function(cluster_df, b_mat,
                                       se_mat,
                                       iter_traits,
                                       exp_pheno,
                                       out_pheno,
                                       num_axis = 2,
                                       pw = 8,
                                       ph = 4) {
  crop_cluster_df <- cluster_df[cluster_df$num_axis == num_axis, ]
  snp_list <- crop_cluster_df$snp_id
  crop_cluster_df <- tibble::column_to_rownames(crop_cluster_df, "snp_id")

  c1 <- exp_pheno
  c2 <- out_pheno
  se_max <- apply(se_mat[snp_list, ], 2, max)
  se_min <- apply(se_mat[snp_list, ], 2, min)
  norm_se <- rowSums((se_mat[snp_list, ] - se_min) / (se_max - se_min))
  alpha_vec <- 1.0 / (1.0 + norm_se)
  res_df <- data.frame(
    row.names = snp_list,
    bx = b_mat[snp_list, c1],
    by = b_mat[snp_list, c2],
    bxse = se_mat[snp_list, c1],
    byse = se_mat[snp_list, c2],
    clust_num = crop_cluster_df[snp_list, "clust_num"],
    clust_prob = crop_cluster_df[snp_list, "clust_prob"],
    alp = alpha_vec
  )
  # Find the axis limits
  xmin <- min(res_df$bx, na.rm = TRUE)
  xmax <- max(res_df$bx, na.rm = TRUE)
  ymin <- min(res_df$by, na.rm = TRUE)
  ymax <- max(res_df$by, na.rm = TRUE)
  # Create the strings for the filename and labels
  pnme <- paste0(iter_traits$res_dir,
                 "scatter_withclusts",
                 c1,
                 "_vs_",
                 c2,
                 "_naxis",
                 num_axis,
                 ".png")
  title_str <- (paste("Clusters plotted against", c1, "and", c2))
  caption_str <- paste("The clustering type used is", iter_traits$clust_typ)
  ggplot2::ggplot(data = res_df,
                  ggplot2::aes(x = bx, y = by)) + # nolint: object_usage_linter.
    ggplot2::geom_point(ggplot2::aes(
      color = clust_num, # nolint: object_usage_linter.
      size = clust_prob # nolint: object_usage_linter.
    ),
    shape = 1
    ) +
    ggplot2::geom_point(ggplot2::aes(
      color = clust_num, # nolint: object_usage_linter.
      size = clust_prob, # nolint: object_usage_linter.
      alpha = alp # nolint: object_usage_linter.
    ),
    shape = 20
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = res_df$bx - 1.96 * res_df$bxse,
        xmax = res_df$bx + 1.96 * res_df$bxse,
        color = clust_num,
        alpha = alp
      ),
      linetype = "solid"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = res_df$by - 1.96 * res_df$byse,
        ymax = res_df$by + 1.96 * res_df$byse,
        color = clust_num,
        alpha = alp
      ),
      linetype = "solid"
    ) +
    ggplot2::xlim(xmin, xmax) +
    ggplot2::ylim(ymin, ymax) +
    ggplot2::xlab(paste("Association with", c1)) +
    ggplot2::ylab(paste("Association with", c2)) +
    ggplot2::ggtitle(title_str) +
    ggplot2::labs(caption = caption_str)
  ggplot2::ggsave(filename = pnme, width = pw, height = ph)
  return()
}
