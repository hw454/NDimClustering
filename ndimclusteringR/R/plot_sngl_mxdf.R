#' Plot the maximum difference for a single iteration set.
#'
#' @param plot_iter the iteration options row.
#' @param max_df_list the list of max difference dataframes
#'   accross all the iteration options
#' @param iter_df dataframe of iteration options
#' @param plotter Open plot to add to
#'
#' @return plot
#'
#' @export
plot_sngl_mxdf <- function(plot_iter, max_df_list, iter_df, plotter) {
  max_df <- max_df_list[max_df_list$input_iter == plot_iter, ]
  max_df["clust_typ"] <- iter_df$clust_typ[plot_iter]
  max_df["bp_on"] <- iter_df$bp_on[plot_iter]
  max_df["cp_on"] <- iter_df$clust_prob_on[plot_iter]
  plotter +
    ggplot2::geom_line(data = max_df,
              ggplot2::aes(x = "num_axis",
                            y = "diff",
                            col = "clust_typ",
                            linetype = "cp_on")
    ) +
    ggplot2::geom_point(data = max_df,
              ggplot2::aes(x = "num_axis",
                          y = "diff",
                          col = "clust_typ",
                          shape = "bp_on")
    ) +
    ggplot2::labs(x = "Number of iterations",
                           y = "Maximum difference",
                           shape = "bp_on",
                           color = "clust type",
                           linetype = "clust prob on")
  return(plotter)
}