# Comparing the results of cross validaton across all layers and performing ANOVA test
source('functions/plot_functions.R')
source('functions/general_functions.R')
#setwd("outputData/glm/stratified_networks_CV") # CV for networks layers all pathways
setwd("outputData/glm/stratified_moa_pathways") # CV for networks layers refined pathways (2k pathways)

li<-list.files('.',pattern = '.RData')
load(li[[1]])
F1<-lapply(modls, function(x)x$f1)
for (i in 2:length(li)){ #for each CV iteration we used 4 times 5-repeated-5-fold-cv=100 data points
  load(li[[i]])
  for (j in 1:length(F1)){ #for each network 
    
    F1[[j]]<-c(F1[[j]],modls[[j]]$f1)
  }
}
setwd("../../../outputData/glm/stratified_MIE_CV/")  # CV at MIES level
li<-list.files('.',pattern = '.RData')
load(li[[1]])
F1_mie<-stability$f1
for (i in 2:length(li)){
  load(li[[i]])
  F1_mie<-c(F1_mie,stability$f1)
}
setwd('../../..')
#names(F1)<-name_correct(names(F1))
names(F1)<-gsub(x=names(F1),pattern = 'Consensus_LINCS_HEPG2',replacement = 'Consensus_LINCS_HEPG2_1_day')
nams<-factor(names(F1)[c(1:4,11,12,6,7,5,15,13,14,8,9,10)])
dat<-as.data.frame(list(F1,F1_mie))
levels(nams)<-c(levels(nams),'MIEs (Gene Level)')
nams[length(nams)+1]<-'MIEs (Gene Level)'
colnames(dat)[length(nams)]<-'MIEs (Gene Level)'
dat<-stack(dat)
res.aov<-aov(values ~ ind,data=dat)
anova_summary<-summary(res.aov)
tk_nonscaled_NES<-TukeyHSD(res.aov)
tk_nonscaled_NES<-as.data.frame(tk_nonscaled_NES$ind)

# saving the results
write.csv(tk_nonscaled_NES,file='outputData/excel_files/ANOVA_EDC_DECOY_stability_all_layers_moa.csv')
save(anova_summary,tk_nonscaled_NES,file='outputData/ANOVA_EDC_decoy_stablities_across_all_networks.RData')

# plot
color_values<-c('#E6F2FF','#CCE6FF','#B3D9FF','#99CCFF',                     # Drug matrix hep,1,3,5
                '#F2F9EC','#E6F2D9',                                         # in vitro rat and human from TGgates
                '#D9ECC6','#CCE6B3','#BFDF9F','#B3D98C','#A6D279','#99CC66', # in vivo  rat low,high,moddle,8,15,29 TGgates
                '#FFCCBF','#FF9980','#FF6640',                               # consensus hep, lINCS and PPI
                '#999999')                                                   # MIE
library(ggplot2)
library(ggpubr)
library(imager)
library(grid)
p<-ggplot(dat, aes(x=factor(dat$ind,level=nams), y=dat$values,fill=factor(dat$ind,level=nams))) +
  geom_boxplot() +theme_minimal()+ ggtitle('')+xlab('Data Layers')+ylab('F1-score')+scale_fill_manual(values=color_values)+
  labs(x='',fill='Data Layers')+
  scale_y_continuous(breaks = seq(0.7,1,0.05), limits = c(0.7,1)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.minor = element_blank(),legend.title = element_text(size = 10),
        axis.text.y   = element_text(size=16), axis.title.y   = element_text(size=16),
        legend.text = element_text(size = 16))
print(p) # 12*5

save(dat,nams,color_values,file = 'outputData/plots_tables_functions/F1_scores_EDC_vs_Decoys/F1_scores.RData')



#a_legend<-paste('a) Two levels of data were used in order to classify EDCs and negative controls', 
#                  'one is matrices compounds and  pathway activation scores from random walk with restart and gene set enrichment analysis on different gene co-expression networks ',
#                  'and the gene level which is the binary matrix of compounds and MIEs',sep='')
#b_legend<-'b) The box plot of accuracy across all data layers'
#text<-paste(a_legend,b_legend,sep = ' ')
#text.p<-ggparagraph(text = text,face='bold',size = 11,color = 'black')
#im<-load.image('images/Untitled.png')
#im2<-load.image('images/prediciotn.jpg')
#theme_set(theme_pubr())
#g1<-rasterGrob(im)
#g2<-rasterGrob(im2)

# figure 2
#ggarrange(g1,p,ncol = 2,nrow = 1,labels = c('a)','b)'),font.label = list(size=25,color='black'))

#tiff('plots/figure2.tiff',units = 'in',width = 25,height = 10,res = 300)
#ggarrange(g1,p,ncol = 2,nrow = 1,labels = c('a)','b)'),font.label = list(size=25,color='black'))
#dev.off()


#Figure 1
img2<-imager::load.image('images/figure1.png')
tiff('plots/figure1.tiff',units = 'in',width = 25,height = 10,res = 300)
plot(img2,axes = F)
dev.off()
