#plot_trait_heatmap <- function(c_scores,trait_info){
c_scores <- test
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
    ggsave(filename=plotname)
    dev.off()
  }
#}

#plot_trait_heatmap(c_scores,trait_info)

