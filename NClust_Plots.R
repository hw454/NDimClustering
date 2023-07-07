plot_trait_heatmap <- function(c_scores){
#c_scores <- test # When the function container is commented out use this line to rename the result df
  NM <- unique(c_scores$num_axis)
  for (i in NM){
    Ni <- 2+i
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
    pw=16
    ph=4
    ggsave(filename=plotname,width=pw,height=ph)
    dev.off()
  }
}

get_col_list <- function(df,filter_col,N,ignore_cols=c()){
  #' Filter the dataframe `df` by the column `filter_col` with value N.
  #' Return the columns names for the columns which are not all Nan once filters 
  #' and are not in `ignore_cols`
  filt_df=df[df[filter_col]==N,]
  c_name=colnames(filt_df)
  keep_cols=c()
  for (cn in c_name){
    if(cn %in% ignore_cols){}
    else if(!all(is.na(filt_df[cn]))){keep_cols<-c(keep_cols,cn)}
  }
  return(keep_cols)
}
plot_max_diff <- function(c_scores,norm_typ){
norm_typ=thresh_norm
c_scores <- test
  NM <- unique(c_scores$num_axis)
  max_diff_df <- data.frame(
    num_axis = integer(),
    Max_Diff = numeric()
  )
  ignore_cols=c('clust_size','id','num_axis','total_score')
  for (Ni in NM){
    trait_list <- get_col_list(c_scores,'num_axis',Ni,ignore_cols)
    trait_list_no_cnum<-trait_list[trait_list != 'clust_num']
    c_scores_term <- c_scores[c_scores$num_axis==Ni,trait_list]
    c_scores_term['clust_num']<-c_scores[c_scores$num_axis==Ni,'clust_num']
    clust_nums <- unique(c_scores_term$clust_num)
    cs_diff <- 0.0
    c_scores_term['total_score']=numeric()
    for (cn1 in clust_nums){
      cs1<-as.matrix(c_scores_term[c_scores_term$clust_num==cn1,trait_list_no_cnum])
      tot_score=norm(as.matrix(cs1),norm_typ)
      c_scores_term[c_scores_term$clust_num==cn1,'total_score']=tot_score
      for (cn2 in clust_nums){
        cs2<-as.matrix(c_scores_term[c_scores_term$clust_num==cn2,trait_list_no_cnum])
        cs_diff0<-clust_metric(cs1,cs2,norm_typ)
        tot_score0<-norm(as.matrix(cs2),norm_typ)
        c_scores_term[c_scores_term$clust_num==cn2,'total_score']=tot_score0
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

