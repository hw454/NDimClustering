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
#' @param c1 Label for the data column for x-axis
#' @param c2 Label for the data column for y-axis
#' @param num_axis The number of trait axis, default\:0
#' @param pw The plot width, default\:8
#' @param ph The plot heigh, default\:4
#'
#' @export
plot_clust_scatter_rgb <- function(clust_dist_df, b_mat,
                          se_mat,
                          iter_traits,
                          c1,
                          c2,
                          num_axis = 1,
                          pw = 8,
                          ph = 4) {
  ignore_cols <- c("num_axis", "snp_id")
  # Get the data for the number of axis being plotted
  crop_clust_dist_df <- clust_dist_df[clust_dist_df$num_axis == num_axis, ]
  snp_list <- crop_clust_dist_df$snp_id
  crop_clust_dist_df <- tibble::column_to_rownames(crop_clust_dist_df,
                                                  var = "snp_id")
  # Get the list of trait names
  full_trait_list <- colnames(crop_clust_dist_df)
  full_trait_list <- full_trait_list[!(full_trait_list %in% ignore_cols)]
  crop_clust_dist_df <- crop_clust_dist_df[, full_trait_list]
  # Set the transparency alpha to be higher when the se is lower.
  se_max <- apply(se_mat[snp_list, ], 2, max, na.rm = TRUE)
  se_min <- apply(se_mat[snp_list, ], 2, min, na.rm = TRUE)
  norm_se <- rowSums((se_mat[snp_list, ] - se_min) / (se_max - se_min))
  alpha_vec <- 1.0 / (1.0 + norm_se)
  # Normalise the values in each column then assign to rgb value.
  max_dist <- apply(crop_clust_dist_df[snp_list, ], 2, max, na.rm = TRUE)
  min_dist <- apply(crop_clust_dist_df[snp_list, ], 2, min, na.rm = TRUE)
  norm_dist_df <- (crop_clust_dist_df[snp_list, ] - min_dist) /
                                         (max_dist - min_dist)
  norm_dist_df[norm_dist_df < 0] <- 0.0
  norm_dist_df[norm_dist_df > 1] <- 1.0
  norm_dist_df[is.na(norm_dist_df)] <- 0.0
  clust_names <- colnames(norm_dist_df)
  colour_vec <- paste0("rgb(", norm_dist_df[, clust_names[1]], ",",
                    norm_dist_df[, clust_names[2]], ",",
                    norm_dist_df[, clust_names[3]], ")")
  colour_vec <- factor(colour_vec)
  # vector with color values
  my_col_vec <- levels(colour_vec)
  my_col_vec <- sapply(seq_along(my_col_vec),
                  function(i) eval(parse(text = my_col_vec[i])))
  # Assign the data required for pplotting into a dataframe
  res_df <- data.frame(
    row.names = snp_list,
    bx = b_mat[snp_list, c1],
    by = b_mat[snp_list, c2],
    bxse = se_mat[snp_list, c1],
    byse = se_mat[snp_list, c2],
    cols = colour_vec,
    alp = alpha_vec
  )
  # Create the strings for the filename and labels
  title_str <- paste("Clusters plotted against the", c1, "and", c2, "traits.")
  caption_str <- paste("r score given by weighting to", clust_names[1],
              "\n g score given by weighting to", clust_names[2],
              "\n b score given by weighting to", clust_names[3],
              "\n Clustering method used is", iter_traits$clust_typ)
    # Set the filename
  pnme <- paste0(iter_traits$res_dir,
                "scatter_clustrgb",
                c1,
                "_vs_",
                c2,
                "_naxis",
                num_axis,
                ".png")
  # Find the axis limits
  xmin <- min(res_df$bx, na.rm = TRUE)
  xmax <- max(res_df$bx, na.rm = TRUE)
  ymin <- min(res_df$by, na.rm = TRUE)
  ymax <- max(res_df$by, na.rm = TRUE)
  # Set the main plotting data
  ggplot2::ggplot(data = res_df,
                  ggplot2::aes(x = bx, y = by)) + # nolint: object_usage_linter.
  # Set the colour
  ggplot2::geom_point(ggplot2::aes(
                  col = cols # nolint: object_usage_linter.
                 ),
                 shape = 21,
                 show.legend = FALSE) + # nolint: object_usage_linter.
  # Set the error bars
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = res_df$bx - 1.96 * res_df$bxse,
                 xmax = res_df$bx + 1.96 * res_df$bxse,
                 col = cols), linetype = "solid", show.legend = FALSE) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = res_df$by - 1.96 * res_df$byse,
                 ymax = res_df$by + 1.96 * res_df$byse,
                 col = cols), linetype = "solid", show.legend = FALSE) +
  # Add the colour scale
  ggplot2::scale_color_manual(values = my_col_vec) +
  # Add the labels
  ggplot2::xlim(xmin, xmax) +
  ggplot2::ylim(ymin, ymax) +
  ggplot2::ylab(paste("Association with", c2)) +
  ggplot2::xlab(paste("Association with", c1)) +
  ggplot2::ggtitle(title_str) +
  ggplot2::labs(caption = caption_str)
  # Save
  ggplot2::ggsave(filename = pnme, width = pw, height = ph)
  return()
}
