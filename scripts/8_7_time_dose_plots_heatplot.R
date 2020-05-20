
# 1. Dose response curves (class probabilities of elastic GLM) and time response Curves -------------------------------


rm(list=ls())
load('outputData/plot_data/time_exposure.RData')
load('outputData/plot_data/Dose_response.RData')

lin_lvls_dose<-factor(c('No-EDC < No-EDC < EDC ','No-EDC < EDC = EDC','EDC < EDC < EDC','EDC = EDC = EDC'))

library(RColorBrewer)
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

lin_lvls_time<-factor(c('No-EDC < No-EDC < EDC ','No-EDC < EDC = EDC','EDC > EDC > No-EDC','EDC < EDC < EDC','EDC = EDC = EDC'))
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
# #rm(list = ls())
# ## Dose Bio markers
# load('outputData/plot_data/dose_bio_marker_scenario_3.RData')
# source('functions/annotation_functions.R')
# #dose_biomarker$pathway<-pathway_annotate(dose_biomarker$pathway)
# dose_biomarker$Dose<-gsub(dose_biomarker$Dose,pattern = 'low',replacement = 'Low')
# dose_biomarker$Dose<-gsub(dose_biomarker$Dose,pattern = 'middle',replacement = 'Middle')
# dose_biomarker$Dose<-gsub(dose_biomarker$Dose,pattern = 'high',replacement = 'High')
# #library(RColorBrewer)
# p3<-ggplot(dose_biomarker, aes(x=factor(Dose,levels = c('Low','Middle','High')),
#                            y=activation_score,colour=factor(pathway))) +
#   geom_point(show.legend = F)+
#   geom_line(aes(group=pathway),show.legend = F)+
#   #guides(col=guide_legend(nrow = 2))+
#   labs(colour='Pathways')+
#   #scale_color_manual(values=rev(brewer.pal(n = 4, name ='RdYlGn')))+
#   #scale_color_gradient(low='yellow',high = 'red')+
#   xlab('TG_GATEs_Repeated_Dose')+ylab('Pathway Activation Scores')+ theme_minimal()+
#   theme(legend.text = element_text(size=10),
#         axis.title.x = element_blank(),
#         legend.position = 'right',axis.text.x = element_text(size=16),
#         axis.title = element_text(size = 16),axis.text.y = element_text(size=16),
#         panel.spacing = unit(2,'lines'))

# 
# 
# 
# ##Time Bio markers plot of the most significant pathways after exposure based on doses and time points
# load('outputData/plot_data/time_bio_marker_scenario5.RData')
# source('functions/annotation_functions.R')
# #time_biomarker$pathway<-pathway_annotate(time_biomarker$pathway)
# time_biomarker$Dose<-gsub(time_biomarker$Dose,pattern = 'eight_day',replacement = '8 day')
# time_biomarker$Dose<-gsub(time_biomarker$Dose,pattern = 'fifteen_day',replacement = '15 day')
# time_biomarker$Dose<-gsub(time_biomarker$Dose,pattern = 'twentynine_day',replacement = '29 day')
# #library(RColorBrewer)
# p4<-ggplot(time_biomarker, aes(x=factor(Dose,levels = c('8 day','15 day','29 day')),
#                            y=activation_score,colour=factor(pathway))) +
#   geom_point(show.legend = F)+
#   geom_line(aes(group=pathway),show.legend = F)+
#   #guides(col=guide_legend(nrow = 2))+
#   labs(colour='Pathways')+
#   #scale_color_manual(values=rev(brewer.pal(n = 4, name ='RdYlGn')))+
#   #scale_color_gradient(low='yellow',high = 'red')+
#   xlab('TG_GATEs_Repeated_Dose')+ylab('Pathway Activation Scores')+ theme_minimal()+
#   theme(legend.text = element_text(size=10),
#         axis.title.x = element_blank(),
#         legend.position = 'right',axis.text.x = element_text(size=16),
#         axis.title = element_text(size = 16),axis.text.y = element_text(size=16),
#         panel.spacing = unit(2,'lines'))
# ggarrange(p3,p4,ncol = 2,common.legend = F,legend = 'right')


# 2. Heatplot of Significant pathways time of exposure 8,15,29 days TG-Gates-------------------------------------
library('org.Hs.eg.db')
library(ggplot2)
source('functions/annotation_functions.R')
load('outputData/plot_data/time_bio_marker_all_increase_1_2_5_scenarios.RData') #script 9_6_2
time_biomarker<-data;rm(data)

load('outputData/toxdb2gene_final.RData')
significant_pathwawys<-unique(time_biomarker$pathway)
selected_list_pathways<-moa_pathways[significant_pathwawys]
selected_list_pathways<-selected_list_pathways[which(sapply(selected_list_pathways,length)<10)] # we select the pathways with the number of genes less than 10

names(selected_list_pathways)<-pathway_annotate(names(selected_list_pathways))
selected_pathway_list_symbol<-sapply(selected_list_pathways, function(x)mapIds(org.Hs.eg.db, x,  'SYMBOL','ENTREZID'))
all_genes<-unique(unlist(selected_pathway_list_symbol))
load('inputData/text_mining/endocrine_disruption_result.RData')
gene_freq$entrez<-as.character(gene_freq$entrez)
gene_freq$entrez<-mapIds(org.Hs.eg.db, gene_freq$entrez,  'SYMBOL','ENTREZID')
ind<-which(gene_freq$entrez %in% intersect(gene_freq$entrez,all_genes))
pubgenes=gene_freq$Freq[ind]
names(pubgenes)<-gene_freq$entrez[ind]
DEGs<-read.csv(file='inputData/DEGs.csv',header = T,stringsAsFactors = F)
DEGs<-DEGs$Gene[which(DEGs$Organism=='Homo sapiens')]
Deg_genes<-intersect(all_genes,DEGs)
Deg_genes<-Deg_genes[!Deg_genes %in% names(pubgenes)]
Deg_list<-rep(0,length(Deg_genes))
names(Deg_list)<-Deg_genes
#heatplot.enrichResult(selected_pathway_list_symbol,c(pubgenes,Deg_list))
citation<-c(pubgenes,Deg_list)
ldf <- lapply(1:length(selected_pathway_list_symbol), function(i) {
  data.frame(categoryID=rep(names(selected_pathway_list_symbol[i]),
                            length(selected_pathway_list_symbol[[i]])),
             Gene=selected_pathway_list_symbol[[i]])
})
d<-do.call('rbind', ldf)
d$citation <- citation[as.character(d[,2])]
time_heatplot<-d
save(time_heatplot,file = 'outputData/plots_tables_functions/Dose_response/time_heatplot.RData')

p3<-ggplot(d, aes_(~Gene, ~categoryID)) +
    geom_tile(aes_(fill = ~citation), color = "white") +
     scale_fill_continuous(low="orange", high="blue", name = "Pubmed Citations")+
     xlab(NULL) + ylab(NULL) + theme_minimal() +
     theme(panel.grid.major = element_blank(),
           axis.text.x=element_text(angle = 60, hjust = 1))
 

# 3. Heatplot of significant pathways dose exposure Low, Middle and High for TG-Gates ---------------------


library('org.Hs.eg.db')
library(ggplot2)
source('functions/annotation_functions.R')

load('outputData/plot_data/dose_bio_marker_all_increase_scenario_1_2_3.RData')
dose_biomarker<-data;rm(data)
load('outputData/toxdb2gene_final.RData')
significant_pathwawys<-unique(dose_biomarker$pathway)
selected_list_pathways<-moa_pathways[significant_pathwawys]
selected_list_pathways<-selected_list_pathways[which(sapply(selected_list_pathways,length)<10)]

names(selected_list_pathways)<-pathway_annotate(names(selected_list_pathways))
selected_pathway_list_symbol<-sapply(selected_list_pathways, function(x)mapIds(org.Hs.eg.db, x,  'SYMBOL','ENTREZID'))
all_genes<-unique(unlist(selected_pathway_list_symbol))
load('inputData/text_mining/endocrine_disruption_result.RData')
gene_freq$entrez<-as.character(gene_freq$entrez)
gene_freq$entrez<-mapIds(org.Hs.eg.db, gene_freq$entrez,  'SYMBOL','ENTREZID')
ind<-which(gene_freq$entrez %in% intersect(gene_freq$entrez,all_genes))
pubgenes=gene_freq$Freq[ind]
names(pubgenes)<-gene_freq$entrez[ind]
DEGs<-read.csv(file='inputData/DEGs.csv',header = T,stringsAsFactors = F)
DEGs<-DEGs$Gene[which(DEGs$Organism=='Homo sapiens')]
Deg_genes<-intersect(all_genes,DEGs)
Deg_genes<-Deg_genes[!Deg_genes %in% names(pubgenes)]
Deg_list<-rep(0,length(Deg_genes))
names(Deg_list)<-Deg_genes
#heatplot.enrichResult(selected_pathway_list_symbol,c(pubgenes,Deg_list))
citation<-c(pubgenes,Deg_list)
ldf <- lapply(1:length(selected_pathway_list_symbol), function(i) {
  data.frame(categoryID=rep(names(selected_pathway_list_symbol[i]),
                            length(selected_pathway_list_symbol[[i]])),
             Gene=selected_pathway_list_symbol[[i]])
})
d<-do.call('rbind', ldf)
d$citation <- citation[as.character(d[,2])]
dose_heatplot<-d
save(dose_heatplot,file = 'outputData/plots_tables_functions/Dose_response/dose_heatplot.RData')
p4<-ggplot(d, aes_(~Gene, ~categoryID)) +
  geom_tile(aes_(fill = ~citation), color = "white") +
  scale_fill_continuous(low="orange", high="blue", name = "Pubmed Citations")+
  xlab(NULL) + ylab(NULL) + theme_minimal() +
  theme(panel.grid.major = element_blank(),
        axis.text.x=element_text(angle = 60, hjust = 1))

 g2<-ggarrange(p4,p3,nrow = 2,common.legend = T,legend = 'bottom')
print(g2)
ggarrange(g1,g2,nrow = 2,common.legend = F,legend = 'bottom')