
# 1. Dose response curves (class probabilities of elastic GLM) and time response Curves -------------------------------


rm(list=ls())
library(ggpubr)
library(RColorBrewer)
load('outputData/plot_data/Dose_response.RData')
load('outputData/plot_data/time_exposure.RData')

lin_lvls_dose<-factor(c('No_EDC - No_EDC < EDC', #G1
                        'No_EDC < EDC - EDC',    #G2
                        'EDC < EDC <= EDC',      #G3
                        'EDC - EDC - EDC'))      #G4   dose levels

lin_lvls_time<-factor(c('No_EDC - No_EDC < EDC', #G1
                        'No_EDC < EDC - EDC',    #G2
                        'EDC < EDC <= EDC',      #G3
                        'EDC - EDC - EDC',       #G4
                        'EDC - EDC > No_EDC'     #G5 time levels
                        ))      


p1<-ggplot(dose_final, 
       aes(x = factor(data_layer,levels = c('Low','Middle','High')), 
           y = len, colour=factor(pattern,levels = lin_lvls_dose)))+
  geom_errorbar(aes(ymin = len-sd, ymax = len+sd),width = .2,position=position_dodge(0.04)) +
  geom_point(position=position_dodge(0.04),size = 2)+
  #scale_color_manual(values=rev(brewer.pal(n = 4, name ='RdBluGn')))+
  scale_color_manual(values=c('green3', 'forestgreen',  'darkgoldenrod2', 'darkorange'))+
  xlab('TG_GATEs_Single_Dose_1_day')+ylab('Class Probability')+ theme_minimal()+
  scale_y_continuous(breaks = seq(0, 1.01, .1), limits = c(0, 1.01)) +
  labs(colour='Groups')+
  geom_line(aes(group=pattern),show.legend = T)+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=10),
        axis.title.x = element_blank(),
        legend.position = 'right',
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        panel.spacing = unit(2,'lines'))

time_final$data_layer<-gsub(x=time_final$data_layer,pattern = '_',replacement = ' ')
p2<-ggplot(time_final, 
           aes(x = factor(data_layer,levels = c('8 days','15 days','29 days')), 
               y = len,colour=factor(pattern,levels = lin_lvls_time)))+
  geom_errorbar(aes(ymin = (len+0.01)-sd, ymax = len+sd),width = .2,position=position_dodge(0.04)) +
  geom_point(position=position_dodge(0.04),size = 2)+
  scale_color_manual(values=c('green3', 'forestgreen', 'firebrick2', 'darkgoldenrod2', 'darkorange'))+
  xlab('TG_GATEs_Repeated_Dose')+ylab('Class Probability')+ theme_minimal()+
 scale_y_continuous(breaks = seq(0, 1, .1), limits = c(0, 1.01)) +
  labs(colour='Groups')+
  geom_line(aes(group=pattern),show.legend = T)+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.text = element_text(size=8),legend.title = element_text(size=10),
        axis.title.x = element_blank(),
        legend.position = 'right',
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        panel.spacing = unit(2,'lines'))

g1<-ggarrange(p1,p2,ncol = 2,common.legend = F,legend = 'bottom')
print(g1)



# 2.  preselection of pathways based on ROC curve 
# Heatplot of Significant pathways for doses (low middle and high) and time of exposure 8,15,29 days TG-Gates-------------------------------------

heat_plot<-function(significant_pathwawys,list_groups){
  load('outputData/toxdb2gene_final.RData')
  pathsize<-as.data.frame(sapply(moa_pathways, length))
  colnames(pathsize)<-"length_pathway"
  pathsize$categoryID<-rownames(pathsize)
  source('functions/annotation_functions.R')
  library('org.Hs.eg.db')
  library(ggplot2)
  selected_list_pathways<-moa_pathways[significant_pathwawys]
  #selected_list_pathways<-selected_list_pathways[which(sapply(selected_list_pathways,length)<10)] # we select the pathways with the number of genes less than 10
  selected_pathway_list_symbol<-sapply(selected_list_pathways, function(x)mapIds(org.Hs.eg.db, x,  'SYMBOL','ENTREZID'))
  all_genes<-unique(unlist(selected_pathway_list_symbol))
  load('inputData/text_mining/endocrine_disruption_result.RData') # text mining results
  gene_freq$entrez<-as.character(gene_freq$entrez)
  gene_freq$entrez<-mapIds(org.Hs.eg.db, gene_freq$entrez,  'SYMBOL','ENTREZID')
  ind<-which(gene_freq$entrez %in% intersect(gene_freq$entrez,all_genes))
  pubgenes=gene_freq$Freq[ind]
  names(pubgenes)<-gene_freq$entrez[ind]
  DEGs<-read.csv(file='inputData/DEGs.csv',header = T,stringsAsFactors = F) # essential genes
  DEGs<-DEGs$Gene[which(DEGs$Organism=='Homo sapiens')]
  essential_genes<-intersect(all_genes,DEGs)
  Deg_genes<-essential_genes[!essential_genes %in% names(pubgenes)]
  Deg_list<-rep(0,length(Deg_genes))
  names(Deg_list)<-Deg_genes
  citation<-c(pubgenes,Deg_list)
  ldf <- lapply(1:length(selected_pathway_list_symbol), function(i) {
    data.frame(categoryID=rep(names(selected_pathway_list_symbol[i]),
                              length(selected_pathway_list_symbol[[i]])),
               Gene=selected_pathway_list_symbol[[i]])
  })
  heatplot<-do.call('rbind', ldf)
  heatplot$citation <- citation[as.character(heatplot[,2])]
  heatplot$group<-'Group'
  for (i in 1:length(list_groups)){
    for (j in 1:nrow(heatplot)){
      if (as.character(heatplot$categoryID[j]) %in% list_groups[[i]])heatplot$group[j]<-paste(heatplot$group[j],i,sep = '_')
      
    }
  }
  heatplot$essential<-heatplot$Gene %in% essential_genes
  heatplot<-merge(heatplot,pathsize,all.x=T)
  heatplot$categoryID<-pathway_annotate(as.character(heatplot$categoryID))
  heatplot$categoryID<-paste(heatplot$categoryID,heatplot$group,sep = '_')
  # print(ggplot(heatplot, aes_(~Gene, ~categoryID)) +
  #         geom_tile(aes_(fill = ~citation), color = "white") +
  #         scale_fill_continuous(low="orange", high="blue", name = "Pubmed Citations")+
  #         xlab(NULL) + ylab(NULL) + theme_minimal() +
  #         theme(panel.grid.major = element_blank(),
  #               axis.text.x=element_text(angle = 60, hjust = 1)))
  return(heatplot)
}

load('outputData/new_time_dose_sensitivity/ANOVA_Times.RData')
time_significant_pathways_uniq<-unique(do.call(c,time_significant_pathways))
time_heat_plot_data<-heat_plot(time_significant_pathways_uniq,time_significant_pathways)

load('outputData/new_time_dose_sensitivity/ANOVA_DOSEs.RData')
dose_significant_pathways_uniq<-unique(do.call(c,dose_significant_pathways))
dose_heat_plot_data<-heat_plot(dose_significant_pathways_uniq,dose_significant_pathways)

save(dose_heat_plot_data,time_heat_plot_data,file = 'outputData/new_time_dose_sensitivity/heatplot_data.RData')
write.csv(dose_heat_plot_data,file = 'outputData/excel_files/dose_sensitive_TG_gates_pathways.csv')
write.csv(time_heat_plot_data,file = 'outputData/excel_files/time_sensitive_TG_gates_pathways.csv')

# visualize the plot
load('outputData/new_time_dose_sensitivity/heatplot_data.RData')
dose_heat_plot_data<-dose_heat_plot_data[dose_heat_plot_data$length_pathway<=10,]
time_heat_plot_data<-time_heat_plot_data[time_heat_plot_data$length_pathway<=10,]
plot_heat_pl<-function(heat_plot_data){
  require(ggplot2)
  ggplot(heat_plot_data, aes_(~Gene, ~categoryID)) +
    geom_tile(aes_(fill = ~citation), color = "white") +
    scale_fill_continuous(low="orange", high="blue", name = "Pubmed Citations")+
    xlab(NULL) + ylab(NULL) + theme_minimal() +
    theme(panel.grid.major = element_blank(),
          axis.text.x=element_text(angle = 60, hjust = 1))
}
lapply(list(dose_heat_plot_data,time_heat_plot_data), function(x)plot_heat_pl(x))

