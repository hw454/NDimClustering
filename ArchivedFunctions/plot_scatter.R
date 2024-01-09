#' Plot the scatter plot or snps and principal components.
#'
#' @description
#' Colour by cluster number.
#' The x and y axis are the assocation with the principal component
#' The width and height of the errorbars are given by the standard error.
#'
#' @param b_mat The matrix of the score data
#' @param iter_traits The iteration variables for the type of iteration.
#' @param num_axis The number of trait axis, default\:0
#' @param pw The plot width, default\:8
#' @param ph The plot heigh, default\:4
#'
#' @export
plot_scatter <- function(b_mat, iter_traits,
  num_axis = 2, pw = 8, ph = 4
) {
  c1 <- colnames(b_mat)[1]
  c2 <- colnames(b_mat)[2]
  snp_list <- row.names(b_mat)
  res_df <- data.frame(
    row.names = snp_list,
    bx = b_mat[, c1],
    by = b_mat[, c2]
  )
  xmin <- min(res_df$bx, na.rm = TRUE)
  xmax <- max(res_df$bx, na.rm = TRUE)
  ymin <- min(res_df$by, na.rm = TRUE)
  ymax <- max(res_df$by, na.rm = TRUE)
  pnme <- paste0(iter_traits$res_dir,
    "scatter", c1, "_vs_", c2,
    ".png"
  )
  title_str <- paste("Before clustering, scatter ", c1, "vs", c2)
  ggplot2::ggplot(data = res_df,
    ggplot2::aes(x = bx, y = by) # nolint: object_usage_linter.
  ) +
    ggplot2::geom_point() +
    ggplot2::xlim(xmin, xmax) +
    ggplot2::ylim(ymin, ymax) +
    ggplot2::ylab(paste("Association with", c2)) +
    ggplot2::xlab(paste("Association with", c1)) +
    ggplot2::ggtitle(title_str)
  ggplot2::ggsave(filename = pnme, width = pw, height = ph)
  return()
}
