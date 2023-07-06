plot_trait_heatmap <- function(c_scores){
#c_scores <- test # When the function container is commented out use this line to rename the result df
  NM <- unique(c_scores$num_axis)
  for (i in NM){
    Ni <- 4+i
    trait_list <- colnames(c_scores)[5:Ni]
    c_scores_term <- c_scores[c_scores$num_axis==i,]
    c_scores_term <- c_scores_term[,!(names(c_scores_term) %in% c('num_axis','clust_size','id'))]
    c_scores_term <- c_scores_term[ , colSums(is.na(c_scores_term))==0]
    #c_scores_term <- c_scores_term[order(c_scores_term$clust_num),]
    # Extract the association scores for each clust trait pair.
    #Scores_matrix <- data.matrix(c_scores_term[trait_list])
    #rownames(Scores_matrix) <- c_scores_term$clust_num
    #ggplot(Scores_matrix, aes(colnames(Scores_mtrix), rownames(Scores_matrix), fill= Scores_matrix)) + 
     # geom_tile()
    long_form_df <- c_scores_term %>% gather('trait','score',-'clust_num')
    #idnum<-which(colnames(c_scores_term)=='clust_num')
    #measurenum <- which(colnames(c_scores_term) %in% trait_list)
    #long_form_df2 <- melt(c_scores_term,id.vars=idnum,measure.vars=measurenum)
    plotname=paste0(res_dir,'trait_vs_ClustScores_iter',i,'.png')
    heatplot<-ggplot(long_form_df, aes(x = trait, y = clust_num, fill = score)) +
      theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
      geom_tile()
    print(heatplot)
    pw=32
    ph=4
    ggsave(filename=plotname,width=pw,height=ph)
    dev.off()
  }
}

plot_max_diff <- function(c_scores,norm_typ){
  NM <- unique(c_scores$num_axis)
  max_diff_df <- data.frame(
    num_axis = integer(),
    Max_Diff = numeric()
  )
  Max_Diff <- list(length(NM))
  for (Ni in NM){
    c_scores_term <- c_scores[c_scores$num_axis==Ni,]
    clust_nums <- unique(c_scores_term$clust_num)
    Nt=4+Ni
    trait_list <- colnames(c_scores_term)[5:Nt]
    cs_diff <- 0.0
    for (cn1 in clust_nums){
      cs1<-as.matrix(c_scores_term[c_scores_term$clust_num==cn1,trait_list])
      for (cn2 in clust_nums){
        cs2<-as.matrix(c_scores_term[c_scores_term$clust_num==cn2,trait_list])
        cs_diff0<-clust_metric(cs1,cs2,norm_typ)
        if (is.na(cs_diff0)){}
        else if (cs_diff0>cs_diff){cs_diff <- cs_diff0}
      }
    }
    max_diff_df <- max_diff_df %>% add_row(num_axis=Ni,Max_Diff=cs_diff)
  }
  #print(max_diff_df)
  plotname <- paste0(res_dir,"NumAxis_Vs_MaxScoreDiff.png")
  lineplot <- ggplot(data=max_diff_df, aes(x=num_axis, y=Max_Diff, group=1)) +
    geom_line()+
    geom_point()
  print(lineplot)
  pw=4
  ph=4
  ggsave(filename=plotname,width=pw,height=ph)
  #dev.off()
}

#test %>% plot_max_diff(thresh_norm)

