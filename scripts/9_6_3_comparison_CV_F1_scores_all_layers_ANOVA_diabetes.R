# Box plot of Cross Validation across all layers
#input
path_to_cv<-'output/glm_models/CV_two_class/diabetes_2/'     # path to CV folder : we performed 4 times 5 repeated 5-fold cross validation on sampo server
plot_title<-'20 Repeats of 5-fold Cross Validation for Diabetes Type 2'  #title of the box plot
# output
anova_result<-'output/excel_files/ANOVA_diabetes_Type_2.csv'      # file to save
######
setwd(path_to_cv)
li<-list.files('.',pattern = '.RData')
load(li[[1]])

F1<-lapply(CV_results, function(x)x$f1)
for (i in 2:length(li)){
  load(li[[i]])
  for (j in 1:length(F1)){
    
    F1[[j]]<-c(F1[[j]],CV_results[[j]]$f1)
  }
}
setwd('../../../..')
library(plotly)
library(dplyr)
source('functions/plot_functions.R')
source('functions/general_functions.R')
source('functions/annotation_functions.R')
names(CV_results)<-name_correct(names(CV_results))  #
names(F1)<-name_correct(names(F1))                  #
names(CV_results)<-gsub(x=names(CV_results),pattern = 'Consensus_LINCS_HEPG2',replacement = "Consensus_LINCS_HEPG2_1_day")
names(F1)<-gsub(x=names(F1),pattern = 'Consensus_LINCS_HEPG2',replacement = "Consensus_LINCS_HEPG2_1_day")
nams<-name_order(names(CV_results))

dat<-as.data.frame(sapply(F1, function(x)x))

dat<-stack(dat)
res.aov<-aov(values ~ ind,data=dat)                                                                           #ANOVA and tukey post test
disease_CV_all_layers<-TukeyHSD(res.aov)
disease_CV_all_layers<-as.data.frame(disease_CV_all_layers$ind)

write.csv(disease_CV_all_layers,file = anova_result)                   # SAVING ANOVA RESULT
dat<-arrange(dat)
dat$AO<-'Diabetes T2'
save(dat,file ='output/CV_data/T2Diabetes.RData')
#plot                                                                           # color codes
color_values<-c('#E6F2FF','#CCE6FF','#B3D9FF','#99CCFF',                        # Drug matrix hep,1,3,5
                '#F2F9EC','#E6F2D9',                                            # in vitro rat and human from TGgates
                '#D9ECC6','#CCE6B3','#BFDF9F','#B3D98C','#A6D279','#99CC66'   , # in vivo  rat low,high,moddle,8,15,29 TGgates
                '#FFCCBF','#FF9980','#FF6640')                                  # consensus hep, lINCS and PPI

# plot(box_plot(dat,factor(dat$ind,level=nams),dat$values,plot_title,'Data Layers','F1-score',leg_title='Data Layers',
#               fill_values = color_values,axxis_breaks = seq(0.91,.97,0.01),axxis_limits = c(0.91,.97))) 
# 


require(ggplot2)
ggplot(dat, aes(x=factor(dat$ind,level=nams), y=dat$values,fill=factor(dat$ind,level=nams))) +
  geom_boxplot() +theme_minimal()+ ggtitle('')+xlab('')+ylab('F1-Score')+scale_fill_manual(values=color_values)+
  labs(x='',fill='Data Layers')+
  scale_y_continuous(breaks = seq(0.92,0.95,0.01), limits = c(0.92,.95)) +
  #coord_cartesian(ylim = axxis)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.minor = element_blank(),legend.title = element_text(size = 20),
        axis.text.y   = element_text(size=20), axis.title.y   = element_text(size=24),
        plot.margin = unit(c(2,2,2,2),'cm'),legend.text = element_text(size = 24))
