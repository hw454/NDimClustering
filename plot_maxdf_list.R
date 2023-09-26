#' Plot the maximum difference at each iteration for the
#' different input variables.
#' @param max_df dataframe of the maximum differences
#' @param iter_traits the iteration options for each iteration
#' @param pw the width of the plot. default\:4
#' @param ph the height of the plot. default\:4
#' @export
plot_mxdf_list <- function(max_df, iter_traits, pw = 4, ph = 4) {
  pnme <- paste0(iter_traits$res_dir, "NumAxis_Vs_MaxScoreDiff_Compare.png")
  lineplot <- ggplot2::ggplot()
  lineplot <- lineplot +
    ggplot2::geom_line(data = max_df,
              ggplot2::aes(group = "input_iter",
                          x = "num_axis",
                          y = "diff",
                          col = "clust_typ",
                          linetype = "cp_on")
    ) +
    ggplot2::geom_point(data = max_df,
              ggplot2::aes(group = "input_iter",
                        x = "num_axis",
                        y = "diff",
                        col = "clust_typ",
                        shape = "bp_on")
    ) +
    ggplot2::labs(x = "Number of iterations",
                  y = "Maximum difference",
                  shape = "bp_on",
                  color = "clust type",
                  linetype = "clust prob on")

  ggplot2::ggsave(pnme, width = pw, height = ph)
  return(0)
}