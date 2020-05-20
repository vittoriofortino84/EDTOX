rm(list=ls())
load('outputData/patway_glm_coefficients/GLM_univariate_all_scores.RData')  # For moa pathways


glm_cutoff=0
Roc_cutoff=.7
number_layers_cutoff=3



# plot function
plot_function<-function(data_x,Title=''){
  require(ggplot2)
  require(RColorBrewer)
  p<-ggplot(data_x,aes(x = reorder(annotated_pathways,-mean_of_avg_nes_for_each_pathway),
                       factor(network,level=rev(nams)),size =glm_coefs ))+
    geom_point(shape=21,aes(fill=Average_NES_EDC))  +
    scale_size(range=c(6,10),breaks=c(0.1,0.5,1,1.5),limits=c(0.01,1.5),labels = c('<= 0.1','0.5','1','1.5'))+
    #scale_size(
    scale_fill_gradientn(colours = brewer.pal(n = 9, name ='YlOrRd'),
                         breaks=seq(0,max(data_x$Average_NES_EDC),0.2),
                         limits=c(0,max(data_x$Average_NES_EDC)))+
    labs(size= 'GLM-Coef',fill='Activation Score',alpha='Pathway CV-Stability')+
    ylab('')+ xlab('')+ guides(fill= guide_colorbar(label.theme = element_text(angle = 45,size=10)  ),size=guide_legend(nrow = 1,order = 2))+
    #facet_grid(  factor(layer,level=layer_levels) ~ ., scale="free", space = 'free')+
    #facet_grid_sc(cols = vars(is_invitro),scales = 'free')+
    #facet_grid( . ~ factor(is_invitro,level=invitro_levels) , scale="free", space = 'free')+
    facet_wrap(. ~ factor(is_invitro,level=invitro_levels),nrow = 2,scales='free' )+
    theme_minimal() + 
    theme(strip.text.y = element_text(angle=0),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 16),
          legend.position  = 'bottom',
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(angle = 35, hjust = 1,size = 10))+
    ggtitle(Title)
}



name_order<-function(inp){
  nams<-factor(inp[c(grep(x=inp,pattern = 'Drug_Matrix_Rat_invitro'),
                     grep(x=inp,pattern = 'Drug_Matrix_Rat_invivo_Single_Dose_1_day'),
                     grep(x=inp,pattern = 'Drug_Matrix_Rat_invivo_Repeated_Dose_3_days'),
                     grep(x=inp,pattern = 'Drug_Matrix_Rat_invivo_Repeated_Dose_5_days'),
                     grep(x=inp,pattern = 'TG_GATEs_Human_invitro'),
                     grep(x=inp,pattern = 'TG_GATEs_Rat_invitro'),
                     grep(x=inp,pattern = 'Low'),grep(x=inp,pattern = 'Middle'),
                     grep(x=inp,pattern = 'High'),
                     grep(x=inp,pattern = '8_days'),grep(x=inp,pattern = '15_days'),
                     grep(x=inp,pattern = '29_days'),
                     grep(x=inp,pattern = 'Consensus'),
                     grep(x=inp,pattern = 'PPI'))])
  return(nams)
} #end function


# treatment
all_data_final$glm_coefs<-abs(all_data_final$glm_coefs)
all_data_final<-all_data_final[all_data_final$glm_coefs>glm_cutoff,]  
nams<-name_order(unique(all_data_final$network))

data<-all_data_final[which(all_data_final$ROC_AUC >Roc_cutoff &
                             all_data_final$number_of_most_informative_layers_for_pathway>number_layers_cutoff),]  # the pathway in at least n most informative networks


avg_avg_nes_pathway<-by(data$Average_NES_EDC,data$pathways,mean)
mean_of_avg_nes_for_each_pathway<-as.numeric(avg_avg_nes_pathway)
names(mean_of_avg_nes_for_each_pathway)<-names(avg_avg_nes_pathway)
min_avg_nes_pathway<-by(data$Average_NES_EDC,data$pathways,min)
min_of_avg_nes_for_each_pathway<-as.numeric(min_avg_nes_pathway)
names(min_of_avg_nes_for_each_pathway)<-names(min_avg_nes_pathway)
mean_of_avg_nes_for_each_pathway<-as.data.frame(list(mean_of_avg_nes_for_each_pathway,min_of_avg_nes_for_each_pathway))
names(mean_of_avg_nes_for_each_pathway)<-c('mean_of_avg_nes_for_each_pathway','min_of_avg_nes_for_each_pathway')
mean_of_avg_nes_for_each_pathway$pathways<-rownames(mean_of_avg_nes_for_each_pathway)
data<-merge(mean_of_avg_nes_for_each_pathway,data,all=T)
data<-data[which(data$min_of_avg_nes_for_each_pathway>0),] 
data$layer<-rep('PPI',nrow(data))
data$layer[grep(data$network,pattern = 'Drug_matrix',ignore.case = T)]<-'Drug_Matrix'
data$layer[grep(data$network,pattern = 'TG_GATEs',ignore.case = T)]<-'TG_GATEs'
data$layer[grep(data$network,pattern = 'Consensus',ignore.case = T)]<-'Consensus'
layer_levels<-c('Drug_Matrix','TG_GATEs','Consensus','PPI')
invitro_levels<-c('in vitro & PPI','in vivo')
data$is_invitro<-'in vitro & PPI'
data$is_invitro[grep(data$network,pattern = 'vivo',ignore.case = T)]<-'in vivo'

p<-plot_function(data)
print(p)

