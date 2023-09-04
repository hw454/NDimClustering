plot_trait_heatmap <- function(c_scores, clust_typ_str,
bp_on = TRUE, clust_prob_on = TRUE) {
#  c_scores <- test1 # When the function container is commented out use this line to rename the result df
  nm <- unique(c_scores$num_axis)
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
  if (bp_on){
    bp_str <- " and bp in score"
  } else {
    bp_str <-" and bp not in score"
  }
  if (clust_prob_on){
    clust_prob_str <- " and clust_prob in score"
  } else{
    clust_prob_str <- " and clust_prob not in score"
   }
  title_str <- paste("Association score for trait against cluster. 
  Cluster type", clust_typ_str, bp_str, clust_prob_str)
  for (i in nm){
    # Get the traits for this iteration
    trait_list <- get_col_list(c_scores, "num_axis", i, ignore_cols)
    trait_list_no_cnum <- trait_list[trait_list != "clust_num"]
    # FIXME calculate total score for each cluster an add as trait.
    # Extract the association scores for each clust trait pair.
    c_scores_term <- c_scores[c_scores$num_axis == i, ]
    c_scores_term <- c_scores_term[trait_list]
    long_form_df <- c_scores_term %>% gather("trait", "score", -"clust_num")
    # Use log scale on the association scores
    long_form_df$score <- log(abs(long_form_df$score))
    # Plot the scores against the traits.
    plotname <- paste0(res_dir, "trait_vs_ClustScores_iter", i, ".png")
    title_iter <- paste0(title_str, ". \n Iteration number ", i)
    heatplot <- ggplot(long_form_df,
                      aes(x = trait, y = clust_num, fill = score)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      geom_tile() +
      # The limits are fixed but the colour change points are moving
      #FIXME
      # Mark the two clusters with the highest difference at each step
      scale_fill_gradient2(low = "cyan", high = "blue", mid = "purple",
                          # na.value = "grey50",
                         midpoint = colmid,
                         breaks = seq(vmin, vmax, break_width),
                         limits = c(vmin, vmax)) +
      ggtitle(title_iter)
    pw <- 16
    ph <- 4
    ggsave(filename = plotname, width = pw, height = ph)
  }
}

get_col_list <- function(df, filter_col, N, ignore_cols=c()) {
  #' Filter the dataframe `df` by the column `filter_col` with value N.
  #' Return the columns names for the columns which are not all Nan once filters 
  #' and are not in `ignore_cols`
  filt_df <- df[df[filter_col] == N, ]
  c_name <- colnames(filt_df)
  keep_cols_list <-lapply(c_name, check_col_nans,
                        ignore_cols = ignore_cols,
                        df = filt_df)
  keep_cols_list <- keep_cols_list[!sapply(keep_cols_list,is.null())]
  return(keep_cols_list)
}

check_col_nans <- function(cn, df, ignore_cols) {
  if (cn %in% ignore_cols) {
    return()
  } else if (!all(is.na(filt_df[cn]))) {
    return(cn)
  } else {
    return()
  }
}

max_diff_single_axis <- function(ni, ignore_cols, c_scores, trait_list, norm_typ = "F"){
  trait_list <- get_col_list(c_scores, "num_axis", ni, ignore_cols = ignore_cols)
  trait_list_no_cnum <- trait_list[trait_list != "clust_num"]
  c_scores_term <- c_scores[c_scores$num_axis == ni, trait_list]
  clust_nums <- unique(c_scores_term$clust_num)
  cs_diff <- 0.0
  c_scores_term["total_score"] <- numeric()
  c1out <- 0
  c2out <- 0
  for (cn1 in clust_nums){
    cs1 <- as.matrix(c_scores_term[c_scores_term == cn1, trait_list_no_cnum])
    tot_score <- norm(ax.matrix(cs1), norm_typ)
    c_scores_term[ c_scores_term$clust_num==cn1, "total_score"] <- tot_score
    for (cn2 in clust_nums){
      cs2 <- as.matrix(c_scores_term[c_scores_term$clust_num == cn2,
                                    trait_list_no_cnum])
      cs_diff0 <- clust_metric(cs1, cs2, norm_typ)
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
    Max_Diff = cs_diff,
    clust_num1 = c1out,
    clust_num2 = c2out
  )
  return(out_df)
}

create_max_diff <- function(c_scores, norm_typ){
  nm <- unique(c_scores$num_axis)
  max_diff_df <- data.frame(
    num_axis = integer(),
    Max_Diff = numeric(),
    clust_num1 = integer(),
    clust_num2 = integer()
  )
  ignore_cols <- c("clust_size", "id", "num_axis", "total_score")
  max_diff_list <- lapply(nm,max_diff_single_axis,
                          ignore_cols = ignore_cols,
                          c_scores = c_scores,
                          trait_list = trait_list,
                          norm_typ = norm_typ)
  max_diff_df<- Reduce(rbind, max_diff_list)
  return(max_diff_df)
}

plot_max_diff <- function(max_df, clust_typ_str = "basic",
                          bp_on = TRUE, clust_prob_on = TRUE) {
  plotname <- paste0(res_dir, "NumAxis_Vs_MaxScoreDiff")
  if (bp_on) {
    bp_str <- " and bp on"
    bp_name_str <- "_bpON"
    } else{
    bp_str <- " and bp off"
    bp_name_str <- "_bpOFF"
  }
  if (clust_prob_on){
    clust_prob_str <- " and clust prob on"
    clust_prob_name_str <- "_clustprobON"
  } else {
    clust_prob_str <- " and clust prob off"
    clust_prob_name_str <- "_clustprobOFF"
  }
  plotname <- paste0(plotname, clust_typ_str, bp_name_str,
                     clust_prob_name_str, ".png")
  plot_title <- paste("Number of axis against max cluster difference score. 
  \n Clustering type", clust_typ_str, bp_str, clust_prob_str)
  lineplot <- ggplot() +
    geom_line(data = max_df, aes(x = num_axis, y = Max_Diff, group = 1),
              color = 'black') +
    geom_point(data = max_df, aes(x = num_axis, y = Max_Diff, group = 1),
              color = 'black') +
    ggtitle(plot_title)
  print(lineplot)
  pw <- 16
  ph <- 4
  ggsave(filename = plotname, width = pw, height = ph)
  return()
}

plot_max_diff_both <- function(max_df1, max_df2){
  plotname <- paste0(res_dir, "NumAxis_Vs_MaxScoreDiff_Compare.png")
  lineplot <- ggplot() +
    geom_line(data = max_df1,
              aes(x = num_axis, y = Max_Diff, group = 1),
              color= "black") +
    geom_point(data = max_df1,
              aes(x = num_axis, y = Max_Diff, group = 1),
              color= "black") +
    geom_line(data = max_df2, 
              aes(x = num_axis, y = Max_Diff),
              color = "red", linetype = "dashed") +
    geom_point(data = max_df2,
              aes(x = num_axis, y = Max_Diff),
              color = "red", linetype = "dashed")
  print(lineplot)
  pw <- 4
  ph <- 4
  ggsave(filename = plotname, width = pw, height = ph)
}


plot_max_diff_list <- function(max_df_list, iter_df) {
  N_sets <- unique(max_df_list$input_iter)
  plotname <- paste0(res_dir, "NumAxis_Vs_MaxScoreDiff_Compare.png")
  lineplot <- ggplot()
  print(N_sets)
  for (plot_iter in N_sets){
    max_df0 <- max_df_list[max_df_list$input_iter == plot_iter, ]
    max_df0["clust_typ"] <- iter_df$clust_typ[plot_iter]
    max_df0["bp_on"] <- iter_df$bp_on[plot_iter]
    max_df0["cp_on"] <- iter_df$clust_prob_on[plot_iter]
    lineplot<-lineplot+
    geom_line(data = max_df0,
              aes(x = num_axis, y = Max_Diff, col = clust_typ, linetype = cp_on)
              ) +
    geom_point(data = max_df0, aes(x = num_axis, y = Max_Diff,
                                 col = clust_typ,shape = bp_on)
               )
  }
  lineplot<-lineplot+ labs(x = "Number of iterations",
                          y = "Maximum difference",
                          shape= "bp_on",
                          color= "clust type",
                          linetype= "clust prob on")
  print(lineplot)
  pw <- 4
  ph <- 4
  ggsave(filename = plotname, width = pw, height = ph)
}
