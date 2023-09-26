plot_trait_heatmap <- function(c_scores, iter_traits, pw = 16, ph = 4) {
  #' Plot a heatmap of the scores for each trait for each cluster.
  ignore_cols <- c("clust_size", "id", "num_axis", "total_score")
  full_trait_list <- colnames(c_scores)
  full_trait_list <- full_trait_list[!(full_trait_list %in% ignore_cols)]
  full_trait_list <- full_trait_list[full_trait_list != "clust_num"]
  vmin <- min(log(abs(c_scores[full_trait_list])), na.rm = TRUE)
  vmax <- max(log(abs(c_scores[full_trait_list])), na.rm = TRUE)
  break_width <- (vmax - vmin) / 4
  colmid <- (vmax + vmin) / 2
  print(paste("vmin", vmin))
  print(paste("vmax", vmax))
  print(paste("break_width", break_width))
  d_str <- str_funcs::desc_str(iter_traits)
  title_str <- paste("Association score for trait against cluster. 
  Cluster type", d_str)
  for (i in unique(c_scores$num_axis)){
    # Get the traits for this iteration
    trait_list <- get_col_list(c_scores, "num_axis", i, ignore_cols)
    # FIXME calculate total score for each cluster an add as trait.
    # Extract the association scores for each clust trait pair.
    c_scores_term <- c_scores[c_scores$num_axis == i, ]
    c_scores_term <- c_scores_term[, colnames(c_scores_term) %in% trait_list]
    long_form_df <- tidyr::gather(c_scores_term, "trait", "score", -"clust_num")
    # Use log scale on the association scores
    long_form_df$score <- log(abs(long_form_df$score))
    # Plot the scores against the traits.
    pnme <- paste0(iter_traits$res_dir, "trait_vs_ClustScores_iter", i, ".png")
    title_iter <- paste0(title_str, ". \n Iteration number ", i)
    ggplot2::ggplot(long_form_df,
                      ggplot2::aes(x = "trait",
                                  y = "clust_num",
                                  fill = "score")) +
      ggplot2::theme(axis.text.x =
                    ggplot2::element_text(angle = 90, vjust = 0.5)) +
      ggplot2::geom_tile() +
      # The limits are fixed but the colour change points are moving
      #FIXME
      # Mark the two clusters with the highest difference at each step
      ggplot2::scale_fill_gradient2(low = "cyan", high = "blue", mid = "purple",
                        na.value = "grey50",
                         midpoint = colmid,
                         breaks = seq(vmin, vmax, break_width),
                         limits = c(vmin, vmax)) +
      ggplot2::ggtitle(title_iter)
    ggplot2::ggsave(filename = pnme, width = pw, height = ph)
  }
  return()
}

plot_transform_heatmap <- function(e_mat, iter_traits,
                                  num_axis = 0, pw = 16, ph = 4) {
  #' Plot a heatmap of the scores for each trait for each principal component
  d_str <- str_funcs::desc_str(iter_traits)
  title_str <- paste("Association score for trait against PC.", d_str)
  pnme <- paste0(iter_traits$res_dir,
                "trait_vs_PC",
                "_num_axis", num_axis, ".png")
  title_iter <- paste0(title_str)
  # Format matrix into long form.
  e_df <- as.data.frame(e_mat)
  e_df <- tibble::rownames_to_column(e_df, "trait")
  e_long_df <- tidyr::pivot_longer(e_df, -c("trait"),
                            names_to = "PC", values_to = "association")
  ggplot2::ggplot(data = e_long_df,
                ggplot2::aes(x = "PC", y = "trait", fill = "association")) +
      ggplot2::geom_raster() +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::ggtitle(title_iter)
  ggplot2::ggsave(filename = pnme, width = pw, height = ph)
  return()
}

get_col_list <- function(df, filter_col, n, ignore_cols = c()) {
  #' Filter the dataframe `df` by the column `filter_col` with value N.
  #' Return the columns names for the columns which are
  #' not all Nan once filtered and are not in `ignore_cols`
  filt_df <- df[df[filter_col] == n, ]
  c_name <- colnames(filt_df)
  keep_cols_list <- lapply(c_name, check_col_nans,
                        ignore_cols = ignore_cols,
                        df = filt_df)
  keep_cols_list <- keep_cols_list[!sapply(keep_cols_list, is.null)]
  return(keep_cols_list)
}

check_col_nans <- function(cn, df, ignore_cols) {
  if (cn %in% ignore_cols) {
    return()
  } else if (!all(is.na(df[cn]))) {
    return(cn)
  } else {
    return()
  }
}

max_diff_single_axis <- function(ni, ignore_cols, c_scores, trait_list,
                                norm_typ = "F") {
  trait_list <- get_col_list(c_scores, "num_axis", ni,
                            ignore_cols = ignore_cols)
  trait_list_no_cnum <- trait_list[trait_list != "clust_num"]
  c_scores_term <- c_scores[c_scores$num_axis == ni, trait_list]
  clust_nums <- unique(c_scores_term$clust_num)
  cs_diff <- 0.0
  c_scores_term["total_score"] <- numeric()
  c1out <- 0
  c2out <- 0
  for (cn1 in clust_nums){
    cs1 <- as.matrix(c_scores_term[c_scores_term == cn1, trait_list_no_cnum])
    tot_score <- norm(as.matrix(cs1), norm_typ)
    c_scores_term[c_scores_term$clust_num == cn1, "total_score"] <- tot_score
    for (cn2 in clust_nums){
      cs2 <- as.matrix(c_scores_term[c_scores_term$clust_num == cn2,
                                    trait_list_no_cnum])
      cs_diff0 <- dist_funcs::clust_metric(cs1, cs2, norm_typ)
      tot_score0 <- norm(as.matrix(cs2), norm_typ)
      c_scores_term[c_scores_term$clust_num == cn2, "total_score"] <- tot_score0
      if (is.na(cs_diff0)) {
      } else if (cs_diff0 > cs_diff) {
        cs_diff <- cs_diff0
        c1out <- cn1
        c2out <- cn2
      }
    }
  }
  out_df <- data.frame(
    num_axis = ni,
    max_diff = cs_diff,
    clust_num1 = c1out,
    clust_num2 = c2out
  )
  return(out_df)
}

plot_max_diff <- function(max_df, iter_traits, pw = 16, ph = 4) {
  method_str <- method_str(iter_traits)
  d_str <- str_funcs::desc_str(iter_traits)
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
  print(lineplot)
  ggplot2::ggsave(filename = pnme, width = pw, height = ph)
  return()
}

plot_max_diff_both <- function(max_df1, max_df2, res_dir, pw = 4, ph = 4) {
  plotname <- paste0(res_dir, "NumAxis_Vs_MaxScoreDiff_Compare.png")
  lineplot <- ggplot2::ggplot() +
    ggplot2::geom_line(data = max_df1,
              ggplot2::aes(x = "num_axis", y = "diff", group = 1),
              color = "black") +
    ggplot2::geom_point(data = max_df1,
              ggplot2::aes(x = "num_axis", y = "diff", group = 1),
              color = "black") +
    ggplot2::geom_line(data = max_df2,
              ggplot2::aes(x = "num_axis", y = "diff"),
              color = "red", linetype = "dashed") +
    ggplot2::geom_point(data = max_df2,
              ggplot2::aes(x = "num_axis", y = "diff"),
              color = "red", linetype = "dashed")
  print(lineplot)
  ggplot2::ggsave(filename = plotname, width = pw, height = ph)
}

plot_sngl_mxdf <- function(plot_iter, max_df_list, iter_df, pw = 4, ph = 4) {
  max_df <- max_df_list[max_df_list$input_iter == plot_iter, ]
  max_df["clust_typ"] <- iter_df$clust_typ[plot_iter]
  max_df["bp_on"] <- iter_df$bp_on[plot_iter]
  max_df["cp_on"] <- iter_df$clust_prob_on[plot_iter]
  lineplot <- lineplot +
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
  return(lineplot)
}

plot_mxdf_list <- function(max_df, iter_traits, pw = 4, ph = 4) {
  #' Plot the maximum difference at each iteration for the
  #' different input variables.
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

clust_scatter <- function(clusters, b_mat,
                          se_mat,
                          iter_traits,
                          num_axis = 0,
                          pw = 8,
                          ph = 4) {
  pnme <- paste0(iter_traits$res_dir, "clusters_num_axis", num_axis, ".png")
  c1 <- colnames(b_mat)[1]
  c2 <- colnames(b_mat)[2]
  snp_list <- row.names(b_mat)
  res_df <- data.frame(
    row.names = snp_list,
    bx = b_mat[, c1],
    by = b_mat[, c2],
    bxse = se_mat[, c1],
    byse = se_mat[, c2],
    clust_num = clusters[snp_list, "clust_num"],
    clust_prob = clusters[snp_list, "clust_prob"]
  )
  ggplot2::ggplot(data = res_df,
                  ggplot2::aes(x = "bx", y = "by")) +
  ggplot2::geom_point(aes(colour = "clust_num",
                  size = "clust_prob"), shape = 21) +
  ggplot2::geom_errorbarh(
    aes(xmin = "bx" - 1.96 * "bxse", xmax = "bx" + 1.96 * "bxse",
          color = "clust_num"), linetype = "solid") +
  ggplot2::geom_errorbar(
    aes(ymin = "by" - 1.96 * "byse", ymax = "by" + 1.96 * "byse",
          color = "clust_num"), linetype = "solid") +
  ggplot2::ylab("Association with PC2") +
  ggplot2::xlab("Association with PC1") +
  ggplot2::ggtitle("Clustered by principal components") +
  ggplot2::ggsave(filename = pnme, width = pw, height = ph)
  return()
}
