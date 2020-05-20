# bubble plot for the AOPs
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
# plot function
plot_function<-function(data_x,Title='',y_color){
  require(ggplot2)
  require(RColorBrewer)
  require(facetscales)
  p<-ggplot(data_x,aes(x = reorder(annotated_pathways,-mean_of_avg_nes_for_each_pathway),
                       factor(network,level=rev(nams)),size =glm_coefs ))+
    geom_point(shape=21,aes(fill=Average_NES_positive))  +
    scale_size(range=c(6,10),breaks=c(0.1,0.5,1,1.5),limits=c(0.01,1.5),labels = c('<= 0.1','0.5','1','1.5'))+
    #scale_size(
    scale_fill_gradientn(colours = brewer.pal(n = 9, name ='YlOrRd'),
                         breaks=seq(0,max(data_x$Average_NES_positive),0.2),
                         limits=c(0,max(data_x$Average_NES_positive)))+
    labs(size= 'GLM-Coef',fill='Activation Score',alpha='Pathway CV-Stability')+
    ylab('')+ xlab('')+ guides(fill= guide_colorbar(label.theme = element_text(angle = 45,size=10)  ),size=guide_legend(nrow = 1,order = 2))+
    #facet_grid(  factor(layer,level=layer_levels) ~ ., scale="free", space = 'free')+
    #facet_grid_sc(cols = vars(is_invitro),scales = 'free')+
    #facet_grid( . ~ factor(is_invitro,level=invitro_levels) , scale="free", space = 'free')+
    #facet_wrap(. ~ factor(is_invitro,level=invitro_levels),scales='free' )+
    theme_minimal() + 
    theme(strip.text.y = element_text(angle=0),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 16),
          
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(angle = 35, hjust = 1,size = 10))+
    ggtitle(Title)
  return(p)
}

preprcess<-function(all_data,ROC_cutoff_in_vivo=0,
                    ROC_cutoff_in_vitro=0,
                    glm_invivo_cutoff=0,nes_in_vivo_cutoff=0,
                    glm_in_vitro_cutoff=0,nes_invitro_cutoff=0,
                    how_many_layer_cutoff=0){
# treatment
all_data<-all_data[!all_data$glm_coefs==0,] 
all_data$network<-gsub(x=all_data$network,pattern = 'Consensus_LINCS_HEPG2',replacement = 'Consensus_LINCS_HEPG2_1_day')

nams<<-name_order(unique(all_data$network))
all_data$glm_coefs<-abs(all_data$glm_coefs)

avg_avg_nes_pathway<-by(all_data$Average_NES_positive,all_data$pathways,mean)
mean_of_avg_nes_for_each_pathway<-as.numeric(avg_avg_nes_pathway)
names(mean_of_avg_nes_for_each_pathway)<-names(avg_avg_nes_pathway)

min_avg_nes_pathway<-by(all_data$Average_NES_positive,all_data$pathways,min)
min_of_avg_nes_for_each_pathway<-as.numeric(min_avg_nes_pathway)
names(min_of_avg_nes_for_each_pathway)<-names(min_avg_nes_pathway)

mean_of_avg_nes_for_each_pathway<-as.data.frame(list(mean_of_avg_nes_for_each_pathway,min_of_avg_nes_for_each_pathway))
names(mean_of_avg_nes_for_each_pathway)<-c('mean_of_avg_nes_for_each_pathway','min_of_avg_nes_for_each_pathway')
mean_of_avg_nes_for_each_pathway$pathways<-rownames(mean_of_avg_nes_for_each_pathway)

 data<-merge(mean_of_avg_nes_for_each_pathway,all_data,all=T)
# [1] "pathways"                         "mean_of_avg_nes_for_each_pathway" "min_of_avg_nes_for_each_pathway" 
# [4] "network"                          "Average_NES_positive"             "SD_NES_positive"                 
# [7] "Average_NES_negative"             "SD_NES_negative"                  "glm_coefs"                       
# [10] "CV_stability"                     "CV_coefficients_SD"               "annotated_pathways"              
# [13] "sum_of_voting_layers_for_pathway" "ROC_AUC"                          "average_pathway_coefficient"     
# [16] "average_pathway_ROC"              "max_pathway_ROC"                  "average_pathway_NES_positive"    
# [19] "average_pathway_NES_negative"     "pathway_length" 

data<-data[which( data$sum_of_voting_layers_for_pathway>how_many_layer_cutoff ),] 


data$layer<-rep('NA',nrow(data))
data$layer[grep(data$network,pattern = 'Drug_matrix',ignore.case = T)]<-'Drug_Matrix'
data$layer[grep(data$network,pattern = 'TG_GATEs',ignore.case = T)]<-'TG_GATEs'
data$layer[grep(data$network,pattern = 'Consensus',ignore.case = T)]<-'Consensus'
data$layer[grep(data$network,pattern = 'PPI',ignore.case = T)]<-'PPI'

layer_levels<-c('Drug_Matrix','TG_GATEs','Consensus','PPI')
invitro_levels<-c('in vitro & PPI','in vivo')
data$is_invitro<-'in vitro & PPI'
data$is_invitro[grep(data$network,pattern = 'vivo',ignore.case = T)]<-'in vivo'

# in vivo
net_colr<-c('#99CCFF','#99CCFF','#99CCFF',                     
            '#99CC66','#99CC66','#99CC66','#99CC66','#99CC66','#99CC66') 

data_invivo<-data[data$is_invitro=='in vivo', ]
data_invivo<-data_invivo[which(data_invivo$ROC_AUC >ROC_cutoff_in_vivo &
                                 data_invivo$min_of_avg_nes_for_each_pathway>nes_in_vivo_cutoff &
                                 data_invivo$glm_coefs>=glm_invivo_cutoff),]  
p1<-plot_function(data_invivo,'In Vivo',net_colr)

# in vitro
net_colr<-c('#99CCFF',                     
            '#99CC66','#99CC66',
            '#FF6640','#FF6640','#FF6640')
data_invitro<-data[data$is_invitro=='in vitro & PPI',]
data_invitro<-data_invitro[which(data_invitro$ROC_AUC >ROC_cutoff_in_vitro &
                                   data_invitro$min_of_avg_nes_for_each_pathway>nes_invitro_cutoff &
                                   data_invitro$glm_coefs>=glm_in_vitro_cutoff),]  
p2<-plot_function(data_invitro,'In Vitro & PPI',net_colr)

fp<-ggarrange(p1,p2,ncol=1,nrow = 2,common.legend = T,legend = 'bottom')
return(fp)
}




# Metabolic syndrome
load('output/metabolic_syndrome_coefficients_stabilities_two_class_NES_ROC.RData')  # For moa pathways
metabolic_s<-all_data_final;rm(all_data_final)
fp_metabolic<-preprcess(metabolic_s,ROC_cutoff_in_vivo = .55,ROC_cutoff_in_vitro = .55,how_many_layer_cutoff=5)
print(fp_metabolic)


#Diabetes
load('output/diabetes_coefficients_stabilities_two_class_NES_ROC.RData')  # For moa pathways
diabetes<-all_data_final;rm(all_data_final)
fp_diabetes<-preprcess(diabetes,ROC_cutoff_in_vivo = .69,ROC_cutoff_in_vitro = .69,how_many_layer_cutoff=5)
print(fp_diabetes)



# Atherosclerosis
load('output/atherosclerosis_coefficients_stabilities_two_class_NES_ROC.RData')  # For moa pathways
all_data_final$annotated_pathways<-gsub(all_data_final$annotated_pathways,pattern = '(NOD)',replacement = '')
atherosclerosis<-all_data_final;rm(all_data_final)
fp_atherosclerosis<-preprcess(atherosclerosis,ROC_cutoff_in_vivo = .58,ROC_cutoff_in_vitro = .58,how_many_layer_cutoff=7)
print(fp_atherosclerosis)

#20 *17 aspect ratio for all figures




