#' Plot the max difference against the number of traits.
#'
#' @param max_df The dataframe of max difference scores and number of axis.
#' @param iter_traits the iteration options
#' @param pw the width of the plot. default\:16
#' @param ph the height of the plot. default \:4
#'
#' @export
plot_max_diff <- function(max_df, iter_traits, pw = 16, ph = 4) {
  method_str <- make_path_label_str(iter_traits)
  d_str <- make_full_desc_str(iter_traits)
  pnme <- paste0(iter_traits$res_dir,
                 "NumAxis_Vs_MaxScoreDiff",
                 method_str, ".png")
  plot_title <- paste("Number of axis against max cluster difference score. 
  \n Clustering type", d_str)
  lineplot <- ggplot2::ggplot() +
    ggplot2::geom_line(data = max_df,
                       ggplot2::aes(x = "num_axis", y = "diff", group = 1),
                       color = "black") +
    ggplot2::geom_point(data = max_df,
                        ggplot2::aes(x = "num_axis", y = "diff", group = 1),
                        color = "black") +
    ggplot2::ggtitle(plot_title)
  lineplot
  ggplot2::ggsave(filename = pnme, width = pw, height = ph)
  return()
}