# input  from  the scripts 2_3, 3_3, 4_3, 5_3
load("output/CV_data/atherosclerosis.RData"); atheroscerosis<-dat
load("output/CV_data/CAD.RData"); CAD<-dat
load("output/CV_data/metabolic_syndrome.RData"); metabolic_syndrome<-dat
load("output/CV_data/T2Diabetes.RData"); T2D<-dat;rm(dat)

all_CV_results<-do.call(rbind,list(atheroscerosis,metabolic_syndrome,T2D))
#all_CV_results<-do.call(rbind,list(CAD,atheroscerosis,metabolic_syndrome,T2D))      # Coronary artery disease is also added

data_layers<-unique(all_CV_results$ind)
disease_CV_all_layers<-list()
for (i in 1:length(data_layers)){
  dat<-all_CV_results[all_CV_results$ind==data_layers[i],]
  res.aov<-aov(values ~ AO,data=dat)
  res<-TukeyHSD(res.aov)
  disease_CV_all_layers[[i]]<-as.data.frame(res$AO)
  disease_CV_all_layers[[i]]$data_layer<-data_layers[i]
}
names(disease_CV_all_layers)<-data_layers
   
all_anova_disease_bsed<-do.call(rbind,disease_CV_all_layers)                    # ANOVA Comparing between different AOs for each data layer
all_anova_disease_bsed$pair_wise_comparison<-rownames(all_anova_disease_bsed)
write.csv(all_anova_disease_bsed,file = 'output/excel_files/ANOVA_between_disease.csv')

##
source('functions/annotation_functions.R')

nams<-name_order(unique(all_CV_results$ind))
color_values<-c('#E6F2FF','#CCE6FF','#B3D9FF','#99CCFF',                        # Drug matrix hep,1,3,5
                '#F2F9EC','#E6F2D9',                                            # in vitro rat and human from TGgates
                '#D9ECC6','#CCE6B3','#BFDF9F','#B3D98C','#A6D279','#99CC66'   , # in vivo  rat low,high,moddle,8,15,29 TGgates
                '#FFCCBF','#FF9980','#FF6640')                                  # consensus hep, lINCS and PPI
fill_values = color_values
axxis_breaks = seq(.66,.97,.01)
axxis_limits = c(.66,.97)
all_CV_results$AO<-gsub(x=all_CV_results$AO,pattern = 'Diabetes T2',replacement = 'T2 Diabetes')

require(ggplot2)
  p<-ggplot(all_CV_results, aes(x=factor(all_CV_results$ind,level=nams), y=all_CV_results$values,fill=factor(all_CV_results$ind,level=nams))) +
    geom_boxplot() +theme_minimal()+
    #ggtitle(title)+xlab(x_title)+
    ylab('F1-score')+
    scale_fill_manual(values=fill_values)+
    labs(x='',fill='Data Layers')+
    scale_y_continuous(breaks = axxis_breaks, limits = axxis_limits) +
    facet_grid(.~ factor(all_CV_results$AO,level=unique(all_CV_results$AO)) , scale="free", space = 'free')+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y = element_text(size=16),axis.title.y = element_text(size = 24,face = 'bold'),
          axis.ticks.x=element_blank(),legend.text = element_text(size = 12),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size=28,face = 'bold'),
          plot.margin = unit(c(2,2,2,2),'cm'),legend.position = 'right')
  plot(p)
  img<-imager::load.image('images/figure6.png');img<-grid::rasterGrob(img)
  library(gridExtra)
  grid.arrange(grobs=list(img,p),layout_matrix=rbind(c(1,1,NA,NA,NA,NA),c(2,2,2,2,2,2,2,2),c(2,2,2,2,2,2,2,2)))
  # figure 6
  tiff('plots/figure6.tiff',units = 'in',width = 30,height = 15,res = 300)
  grid.arrange(grobs=list(img,p),layout_matrix=rbind(c(1,1,NA,NA,NA,NA),c(2,2,2,2,2,2,2,2),c(2,2,2,2,2,2,2,2)))
  dev.off()
  
  
