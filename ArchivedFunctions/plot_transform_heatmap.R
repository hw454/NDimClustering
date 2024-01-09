#' Plot the associations between traits and principal components.
#'
#' @param e_mat the matrix of the PC transform
#' @param iter_traits the iteration details
#' @param num_axis the number of axis.
#' @param pw the width of the plot. default\:16
#' @param ph the height of the plot. default\:4
#'
#' @export
plot_transform_heatmap <- function(e_mat, iter_traits,
  num_axis = 0, pw = 16, ph = 4
) {
  # Plot a heatmap of the scores for each trait for each principal component
  d_str <- make_full_desc_str(iter_traits)
  title_str <- paste("Association score for trait against PC.", d_str)
  pnme <- paste0(iter_traits$res_dir,
    "trait_vs_PC",
    "_num_axis", num_axis, ".png"
  )
  title_iter <- paste0(title_str)
  # Format matrix into long form.
  e_df <- as.data.frame(e_mat)
  e_df <- tibble::rownames_to_column(e_df, "trait")
  e_long_df <- tidyr::pivot_longer(e_df, -c("trait"),
    names_to = "PC", values_to = "association"
  )
  ggplot2::ggplot(data = e_long_df,
    ggplot2::aes(x = PC, # nolint: object_usage_linter.
      y = trait, # nolint: object_usage_linter.
      fill = association # nolint: object_usage_linter.
    )
  ) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::ggtitle(title_iter)
  ggplot2::ggsave(filename = pnme, width = pw, height = ph)
  return()
}