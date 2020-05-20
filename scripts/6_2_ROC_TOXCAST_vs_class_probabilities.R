# 1. ROC curve analysis with hitc-values of TOXCAST 3.1 and class probilities from RWR-FGSEA-GLM pipeline  -------------------------------------
library(PRROC)
source('functions/plot_functions.R')
source('functions/general_functions.R')
load('outputData/class_probilities/class_probilities_all_compounds_final_moa.RData') #class-probabilities for all compounds in CTD
load('outputData/hitc_chemid_NR_coreg.RData') # assay endpoints related to nuclear receptors from toxcast
names(all_compounds_prob)<-name_correct(names(all_compounds_prob)) #names of the layers are corrected
assay<-chemid_hic_NR
rocmat_length<-rocmat<-ratio<-matrix(NA, nrow = length(all_compounds_prob), ncol = ncol(assay)) #ROCmat will be the matrix with rows as layers and endpoints columns
colnames(ratio)<-colnames(rocmat)<-colnames(rocmat_length)<-colnames(assay)                     #ROCmat_length is for the number of none NA 
rownames(ratio)<-rownames(rocmat)<-rownames(rocmat_length)<-names(all_compounds_prob)
for (i in 1:length(all_compounds_prob)){
  cp<-all_compounds_prob[[i]]                                          # for each data layer
  ac<-assay[intersect(rownames(assay),rownames(cp)),]
  cp<-cp[intersect(rownames(assay),rownames(cp)),]
  print(names(all_compounds_prob)[[i]])
  for (j in 1:ncol(assay)){                                           # for each assay
    m<-cbind(cp[,'edc'],as.numeric(ac[,j]))                           # class probability for edcs
    m<-na.omit(m)
    rocmat_length[i,j]<-nrow(m)
    ratio[i,j]<-sum(m[,2])/nrow(m)
    if (nrow(m)>100 ){                                                # at least 100 compounds with tested result for each assay for ROC curve
      roc<-roc.curve(scores.class0 = m[,1], weights.class0 = m[,2])
      rocmat[i,j]<-roc$auc
    }#end if
  }#enf for j assay
}#end for i   network
rocmat[is.na(rocmat)]<-0
most_reliable_assay_endpoint<-which(apply(rocmat, 2, mean) >= 0.5)  # most_reliable endpoints are those with average of 0.5 for all networks

toxcast<-read.csv('inputData/Assay_Summary_190226.csv',header = TRUE, stringsAsFactors = F)      # TOXCAST 3.1 assay_summary
descriptions<-do.call(rbind,lapply(names(most_reliable_assay_endpoint), function(x)toxcast[toxcast$assay_component_endpoint_name==x,]))
rownames(descriptions)<-names(most_reliable_assay_endpoint)
data<-t(rocmat[,most_reliable_assay_endpoint]) 
ratio<-t(ratio[,most_reliable_assay_endpoint])
k_means<-kmeans(t(data),centers = 2)                       # most informative layers are those with higher correlatio to most related endpoints using kmeans
most_informative_layers<-which(k_means$cluster==1)    

data<-as.data.frame(data);data$endpoints<-rownames(data)
descriptions$endpoints<-rownames(descriptions)
data<-merge(data,descriptions)           # integration of toxcast descriptions for assays with ROC curve analysis

# Saving
write.csv(data,file='outputData/excel_files/roc_most_related_endpoints_toxcast_moa.csv') # https://www.sciencedirect.com/science/article/pii/S1471489214001271?via%3Dihub
save(data,most_informative_layers,most_reliable_assay_endpoint,rocmat,file='outputData/most_informative_layers_final_moa.RData')


####BOX plot
# data<-as.data.frame(data)
# library(tidyverse)
# source('functions/plot_functions.R')
# source('functions/general_functions.R')
# data<-data%>% gather(network,roc, names(data))
# nams<-name_order(unique(data$network))
# color_values<-c('#E6F2FF','#CCE6FF','#B3D9FF','#99CCFF',                     # Drug matrix hep,1,3,5
#                 '#F2F9EC','#E6F2D9',                                         # in vitro rat and human from TGgates
#                 '#D9ECC6','#CCE6B3','#BFDF9F','#B3D98C','#A6D279','#99CC66', # in vivo  rat low,high,moddle,8,15,29 TGgates
#                 '#FFCCBF','#FF9980','#FF6640')                               # consensus hep, lINCS and PPI
#                                                                
# plot(box_plot(data,factor(data$network,level=nams),data$roc,'','Networks','ROC-AUC',leg_title = 'Network',
#               fill_values = color_values,axxis_breaks = seq(0.55,1,0.025),axxis_limits = c(0.55,0.725)))

###  ANOVA and tukey test between layers
  # res.aov<-aov(roc ~ network,data=data)
  # post_test_layers<-TukeyHSD(res.aov)
  # post_test_layers<-as.data.frame(post_test_layers$network)
  # write.csv(post_test_layers,file='outputData/excel_files/roc_ANOVA_toxcast.csv')

 