#plotting functions

# barplot_pathways -----------------------------------------------------------------
#saving barplot for the glmnet models
pathways_barplot_save<-function(smodel,nam,path,cof_cutoff=0){
  #require(ggplot2)
  mfit.ctcf.2.coef = as.matrix(coef(smodel$finalModel, smodel$bestTune$lambda))
  nams<-annot_correction(rownames(mfit.ctcf.2.coef)[which(abs(mfit.ctcf.2.coef)>cof_cutoff)])
  mfit.ctcf.2.coef<-mfit.ctcf.2.coef[which(abs(mfit.ctcf.2.coef)>cof_cutoff)] 
  nams<-nams[-1]
  mfit.ctcf.2.coef<-mfit.ctcf.2.coef[-1]
  mfit.ctcf.2.sign = rep("positive", length(mfit.ctcf.2.coef)) 
  mfit.ctcf.2.sign[mfit.ctcf.2.coef < 0] = "negative"
  mfit.ctcf.2.dat = data.frame(pathways.names = nams,
                               coefficient = mfit.ctcf.2.coef, 
                               sign = mfit.ctcf.2.sign)
  pn<-factor(nams)
  plot_mfit.ctcf.2 = ggplot2::ggplot(mfit.ctcf.2.dat, ggplot2::aes(reorder(pn,-mfit.ctcf.2.coef), mfit.ctcf.2.coef, fill=sign)) + 
    ggplot2::geom_bar(position="dodge", stat="identity", width=0.75) + 
    ggplot2::scale_fill_manual(values=c("forestgreen","red")) + 
    ggplot2::scale_y_continuous( breaks = seq(min(mfit.ctcf.2.coef),max(mfit.ctcf.2.coef),by=0.1)) + 
    ggplot2::labs(x = "Pathways", y = "Coefficients") +
    ggplot2::theme_bw() + ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(),
                                         axis.line = ggplot2::element_line(colour = "black"),
                                         axis.text.x = ggplot2::element_text(angle = 65, hjust = 2, size = 4),
                                         axis.title = ggplot2::element_text(size=7,face="bold")) 
  plot_mfit.ctcf.2<-plot_mfit.ctcf.2 + ggplot2::ggtitle(nam)
  plot_mfit.ctcf.2
  png(paste(path,'/',nam,'.png',sep=''),
      width     = 8,
      height    = 8,
      units     = "in",
      res       = 1200,
      pointsize = 1
  )
  # par( mar= c(5, 5, 2, 2),  xaxs     = "i",  yaxs = "i",cex.axis = 2, cex.lab  = 2 )
  print(plot_mfit.ctcf.2)
  dev.off()
} #end function

pathways_barplot_visualize<-function(smodel,nam,cof_cutoff=0){
  #require(ggplot2)
  mfit.ctcf.2.coef = as.matrix(coef(smodel$finalModel, smodel$bestTune$lambda))
  nams<-annot_correction(rownames(mfit.ctcf.2.coef)[which(abs(mfit.ctcf.2.coef)>cof_cutoff)])
  mfit.ctcf.2.coef<-mfit.ctcf.2.coef[which(abs(mfit.ctcf.2.coef)>cof_cutoff)] 
  nams<-nams[-1]
  mfit.ctcf.2.coef<-mfit.ctcf.2.coef[-1]
  mfit.ctcf.2.sign = rep("positive", length(mfit.ctcf.2.coef)) 
  mfit.ctcf.2.sign[mfit.ctcf.2.coef < 0] = "negative"
  mfit.ctcf.2.dat = data.frame(pathways.names = nams,
                               coefficient = mfit.ctcf.2.coef, 
                               sign = mfit.ctcf.2.sign)
  pn<-factor(nams)
  plot_mfit.ctcf.2 = ggplot2::ggplot(mfit.ctcf.2.dat, ggplot2::aes(reorder(pn,-mfit.ctcf.2.coef), mfit.ctcf.2.coef, fill=sign)) + 
    ggplot2::geom_bar(position="dodge", stat="identity", width=0.75) + 
    ggplot2::scale_fill_manual(values=c("forestgreen","red")) + 
    ggplot2::scale_y_continuous( breaks = seq(min(mfit.ctcf.2.coef),max(mfit.ctcf.2.coef),by=0.1)) + 
    ggplot2::labs(x = "Pathways", y = "Coefficients") +
    ggplot2::theme_bw() + ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(),
                                         axis.line = ggplot2::element_line(colour = "black"),
                                         axis.text.x = ggplot2::element_text(angle = 65, hjust = 1, size = 6),
                                         axis.title = ggplot2::element_text(size=7,face="bold")) 
  plot_mfit.ctcf.2<-plot_mfit.ctcf.2 + ggplot2::ggtitle(nam)
  plot_mfit.ctcf.2
  print(plot_mfit.ctcf.2)
}#end function

# upset plot --------------------------------------------------------------
#visualize upset
upset_visualize<-function(sel,numb=NA){
  list_pathways <- sel
  all_pathways <- unique(unlist(list_pathways))
  int_data <- lapply(list_pathways, function(x){
    idx <- which(all_pathways %in% x)
    res<-rep(0, length(all_pathways))
    res[idx]<-1
    res})
  upset_df <- as.data.frame(int_data , as.is=T, stringsAsFactors=F, check.names=F)
  rownames(upset_df) <- all_pathways
  
  UpSetR::upset(upset_df, nintersects=numb, nsets=ncol(upset_df), 
                sets.bar.color = "#56B4E9", order.by=c('freq'), text.scale=.9, 
                point.size = 1.5, line.size = 0.7, 
                sets.x.label = "Pathways selected in each network", mainbar.y.label = "Pathway intersection size")
  
} #end function

# bubble plot functions -------------------------------------------------------------
#bubble plot based on fixed bootstap and percent methods
relevant_pathways<-function(modls,type='fixed',cutof){
  if (type=='fixed'){
    all_layers<-lapply(modls, function(x){sel=x$cof[which(abs(x$cof)>cutof),];
    names(sel)=rownames(x$cof)[which(abs(x$cof)>cutof)];
    sel})
  } #end if fixed
  
  if (type=='bootstrap'){
    for (i in 1:length(modls)){
      modls[[i]]$cof<-colSums(modls[[i]]$sel)
    }
    print('bootstrap based relavent pathways')
    all_layers<-lapply(modls, function(x){sel=x$cof[which(x$cof>cutof)];
    names(sel)=names(x$cof)[which(x$cof>cutof)];
    sel}) 
  } #end if bootstrap
  
  if (type=='percent'){
    # bubble plot based on the high percentile for each layer -----------------
    #load("/research/groups/fortino/amirhs/amir/EDCMET/customized_pipe/analysis_result/glmnet/glm_models_all_pathways.RData")
    percent_level<-cutof
    all_layers<-lapply(modls, function(x){
      cutof<-quantile(abs(x$cof[-1,]),probs = percent_level);
      print(cutof)
      sel=x$cof[which(abs(x$cof)>cutof),];
      names(sel)=rownames(x$cof)[which(abs(x$cof)>cutof)];
      sel})
  } #end if percent
  
  all_pathways<-unique(unlist(lapply(all_layers, function(x)names(x))))
  
  data<-matrix(0, nrow = length(all_pathways), ncol = length(names(all_layers)))
  rownames(data)<-all_pathways
  colnames(data)<-names(all_layers)
  for (i in 1:length(names(all_layers))){
    data[intersect(names(all_layers[[i]]),rownames(data)),i]<-all_layers[[i]][intersect(names(all_layers[[i]]),rownames(data))]
  }
  if (type=='fixed'|type=='percent'){
    data<-data[-1,] #eliminating intercept as a feature for the case of fixed and percent method
  } #end if 
  data<-as.data.frame(data)
  data$pathways<-rownames(data)
  return(data)
  
} #end function

bubble_plot<-function(data,pathways,network,coef,stability,plot_name='',fill_values=rep('#B3D9FF',length(unique(network)))){
  library(ggplot2)
  p<-ggplot(data,aes(name=plot_name,pathways, network, size = coef, alpha = stability,fill=network)) + 
    geom_point(shape=21) + theme_classic() + scale_fill_manual(values=fill_values)+
    #scale_size(limits = c(min(abs(data$coef)))) +  
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 7))
          #panel.grid.major = element_line(color='black'))
    return(p)
}



# general heat map plot -----------------------------------------------------------
heat_map<-function(data,value,Var1,Var2){
  require(ggplot2)
   p<-ggplot(data = data, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile()+scale_fill_gradient(low='black',high='white')+
  ggtitle('')+xlab('')+
  ylab('')+theme_bw()+theme(
    panel.border  = element_rect(fill=NA,color='darkred',size=0.15),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
    # axis.line = element_line(colour = 'darkblue')
    )#theme
  return(p)
}

# general barplot -----------------------------------------------------------------
bar_plot<-function(data,x,y,title,x_title,y_title,leg_title=NULL,fill_values=rep('#999999',length(unique(x)))){# a data frame where x is the name  of the vriables, y the values 
  require(ggplot2)
   p<-ggplot(data,aes(x,y,fill=x))+geom_bar(stat = 'identity')+
  theme_classic()+ ggtitle(title)+xlab(x_title)+ylab(y_title)+scale_fill_manual(values=fill_values)+
  labs(x='',fill=leg_title)+
  theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          #axis.text.x = element_text(angle = 45, hjust = 1,size = 8),    #labels in x axis
          plot.margin = unit(c(2,2,2,2),'cm'))
  return(p)
}


# general boxplot ---------------------------------------------------------
box_plot<-function(data,x,y,title,x_title,y_title,leg_title=NULL,fill_values=rep('#999999',length(unique(x))),axxis_breaks=seq(0, 1, .1),axxis_limits=c(0,1)){
  require(ggplot2)
  p<-ggplot(data, aes(x=x, y=y,fill=x)) +
  geom_boxplot() +theme_minimal()+ ggtitle(title)+xlab(x_title)+ylab(y_title)+scale_fill_manual(values=fill_values)+
  labs(x='',fill=leg_title)+
  scale_y_continuous(breaks = axxis_breaks, limits = axxis_limits) +
  #coord_cartesian(ylim = axxis)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.minor = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1,size = 8),    #labels in x axis
        plot.margin = unit(c(2,2,2,2),'cm'))
return(p)
}


# ven diagram ---------------------------------------------------------------------

ven_diag<-function(x,names_sets,fname){# x is a list x=list(set1,set2,set3)
  library(VennDiagram)
  library(RColorBrewer)
  myCol <- brewer.pal(3, "Pastel2")
  venn.diagram(
    x,
    category.names = names_sets,
    filename = fname,
    output=TRUE,
    
    # Output features
    # imagetype="png" ,
    # height = 900 , 
    # width = 900 , 
    # resolution = 300,
    # compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.4,
    # cat.fontface = "bold",
    # cat.default.pos = "outer",
    # cat.pos = c(7, 27, 100),
    cat.dist = c(0.0, 0, 0),
    # cat.fontfamily = "sans",
    rotation = 1)
}


# #heatplot ---------------------------------------------------------------


heatplot.enrichResult <- function(geneSets, foldChange=NULL) {
  d <- list2df(geneSets)
  if (!is.null(foldChange)) {
    d$foldChange <- foldChange[as.character(d[,2])]
    ## palette <- fc_palette(d$foldChange)
    p <- ggplot(d, aes_(~Gene, ~categoryID)) +
      geom_tile(aes_(fill = ~foldChange), color = "white") +
      scale_fill_continuous(low="blue", high="red", name = "fold change")
    ## scale_fill_gradientn(name = "fold change", colors = palette)
  } else {
    p <- ggplot(d, aes_(~Gene, ~categoryID)) + geom_tile(color = 'white')
  }
  p + xlab(NULL) + ylab(NULL) + theme_minimal() +
    theme(panel.grid.major = element_blank(),
          axis.text.x=element_text(angle = 60, hjust = 1))
}
list2df <- function(inputList) {
  ldf <- lapply(1:length(inputList), function(i) {
    data.frame(categoryID=rep(names(inputList[i]),
                              length(inputList[[i]])),
               Gene=inputList[[i]])
  })
  do.call('rbind', ldf)
}

# heatplot
#fc = c(G1=0.5,G2=0.1,G3=1)   # fold change for each gene and the list of the genes should be provided
#heatplot.enrichResult(list(A=c('G1','G2'),B=c('G2','G3')), fc)
heatplot.enrichResult <- function(geneSets, foldChange=NULL) {
  d <- list2df(geneSets)
  if (!is.null(foldChange)) {
    d$foldChange <- foldChange[as.character(d[,2])]
    ## palette <- fc_palette(d$foldChange)
    p <- ggplot(d, aes_(~Gene, ~categoryID)) +
      geom_tile(aes_(fill = ~foldChange), color = "white") +
      scale_fill_continuous(low="orange", high="blue", name = "Pubmed Citations")
    ## scale_fill_gradientn(name = "fold change", colors = palette)
  } else {
    p <- ggplot(d, aes_(~Gene, ~categoryID)) + geom_tile(color = 'white')
  }
  p + xlab(NULL) + ylab(NULL) + theme_minimal() +
    theme(panel.grid.major = element_blank(),
          axis.text.x=element_text(angle = 60, hjust = 1))
}
list2df <- function(inputList) {
  ldf <- lapply(1:length(inputList), function(i) {
    data.frame(categoryID=rep(names(inputList[i]),
                              length(inputList[[i]])),
               Gene=inputList[[i]])
  })
  do.call('rbind', ldf)
}