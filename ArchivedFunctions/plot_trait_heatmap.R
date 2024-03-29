#' Plot a heatmap of the scores for each trait for each cluster.
#'
#' @param c_scores - Dataframe of clusters and their scores for each axis.
#'   Columns for dataframe are the traits and "clust_num". Rows each cluster
#'   number and number of axis variation. Represented by nai_cnj, where i is
#'   the number of axis and j is the cluster number.
#'   More details in [score_trait], [score_cluster]
#'   and [score_all_clusters].
#' @param iter_traits - The properties of the iteration. List containing\:
#'   * "bp_on",
#'   * "clust_prob_on",
#'   * "resdir",
#'   * "clust_typ",
#'   * "ndim_typ"
#'   Destriptions for these terms in [make_trait_df]
#' @param pw Plot width, default\:16
#' @param ph Plot heigh, default\:4.
#' @param prefix The variable category for y-axis
#'
#' @description
#' For values of num_axis plot the heatmap of traits against cluster score.
#' The limits on the axis are set on the range of the full dataset.
#'
#' @return Final heatmap plot
#'
#' @export
plot_trait_heatmap <- function(c_scores, iter_traits,
  pw = 16, ph = 4, prefix = ""
) {
  ignore_cols <- c("num_axis")
  full_trait_list <- colnames(c_scores)
  full_trait_list <- full_trait_list[!(full_trait_list %in% ignore_cols)]
  full_trait_list <- full_trait_list[full_trait_list != "clust_num"]
  vmin <- min(log(abs(c_scores[full_trait_list])), na.rm = TRUE)
  vmax <- max(log(abs(c_scores[full_trait_list])), na.rm = TRUE)
  break_width <- (vmax - vmin) / 4.0
  colmid <- (vmax + vmin) / 2.0
  print(paste("vmin", vmin))
  print(paste("vmax", vmax))
  print(paste("break_width", break_width))
  d_str <- make_full_desc_str(iter_traits)
  title_str <- paste("Association score for trait against cluster.")
  caption_str <- paste("Cluster type", d_str,
                       "\n The traits are:", full_trait_list)
  for (i in unique(c_scores$num_axis)){
    # Get the traits for this iteration
    caption_str <- paste("Cluster type", d_str,
      " The traits are", prefix,
      "The number of axis are ", i
    )
    trait_list <- crop_col_names(c_scores, "num_axis", i, ignore_cols)
    # FIXME calculate total score for each cluster an add as trait.
    # Extract the association scores for each clust trait pair.
    c_scores_term <- c_scores[c_scores$num_axis == i, ]
    c_scores_term <- c_scores_term[, colnames(c_scores_term) %in% trait_list]
    long_form_df <- tidyr::gather(c_scores_term, "trait", "score", -"clust_num")
    # Use log scale on the association scores
    long_form_df$score <- log(abs(long_form_df$score))
    # Plot the scores against the traits.
    pnme <- paste0(iter_traits$res_dir, prefix, "_vs_clust_naxis", i, ".png")
    title_iter <- paste0(title_str, ". \n Iteration number ", i)
    p <- ggplot2::ggplot(long_form_df,
                      ggplot2::aes(y = trait, #nolint
                                  x = clust_num, #nolint
                                  fill = score)) + #nolint
      ggplot2::theme(axis.text.x =
          ggplot2::element_text(angle = 90, vjust = 0.5)
      ) +
      ggplot2::geom_tile() +
      # The limits are fixed but the colour change points are moving
      # FIXME
      # Mark the two clusters with the highest difference at each step
      ggplot2::scale_fill_gradient2(low = "cyan", high = "blue", mid = "purple",
        midpoint = colmid,
        breaks = seq(vmin, vmax, break_width),
        limits = c(vmin, vmax)
      ) +
      ggplot2::ggtitle(title_iter) +
      ggplot2::labs(caption = caption_str)
    ggplot2::ggsave(filename = pnme, width = pw, height = ph)
  }
  return(p)
}