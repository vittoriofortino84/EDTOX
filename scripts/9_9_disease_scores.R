source('functions/glm_functions.R')
source('functions/general_functions.R')
source('functions/annotation_functions.R')

load('input/integrated_NES_scores_for_disease_biomarkers.RData') #NES scores from phase 1 pipeline for all 12000 compounds and all networks

load("output/glm_models/two_class_glm_models_metabolic_syndrome.RData")
metabolic_syndrome_scores<-prob_predict(metabolic_syndrome_two_class,fgs) #prediction of the disease scores


load('output/glm_models/two_class_glm_models_atherosclerosis.RData')
atheroscleorosis_scores<-prob_predict(atherosclerosis_two_class,fgs)

load('output/glm_models/two_class_glm_models_diabetes2.RData')
diabetes_scores<-prob_predict(diabetes_2_two_class,fgs)

load('output/chem2disease.RData')
# in the bellow function syndrome score is the list of predicted class probabilties for all data layers
disease_score_compile<-function(syndrome_scores,disease_name){
  disease_harmonic_average_score<-syndrome_scores[[1]]
  disease_harmonic_average_score<-disease_harmonic_average_score[,grep(x=colnames(disease_harmonic_average_score),
                                                                       pattern='positive')]
  disease_harmonic_average_score<-as.data.frame(disease_harmonic_average_score)
  disease_harmonic_average_score$mesh<-rownames(disease_harmonic_average_score)
  names(disease_harmonic_average_score)<-c(names(syndrome_scores)[1],'mesh')
  for (i in 2:length(syndrome_scores)){
    temp<-syndrome_scores[[i]]
    temp<-temp[,grep(x=colnames(temp),
                     pattern='positive')]
    temp<-as.data.frame(temp)
    temp$mesh<-rownames(temp)
    names(temp)<-c(names(syndrome_scores)[i],'mesh')
    disease_harmonic_average_score<-merge(disease_harmonic_average_score,temp,all = T)
    
  }
  col.ind<-which(!colnames(disease_harmonic_average_score)=='mesh')
  disease_harmonic_average_score$harmonic_disease_score<-harmonic_sum(disease_harmonic_average_score,col.ind)
  disease_harmonic_average_score$harmonic_disease_score<-scales:::rescale(disease_harmonic_average_score$harmonic_disease_score,
                                                                          to = c(0, 1))
  
  disease_harmonic_average_score$average_disease_score<-mean_function(disease_harmonic_average_score,col.ind)
  disease_harmonic_average_score$cas<-mesh2cas(as.character(disease_harmonic_average_score$mesh))
  disease_harmonic_average_score$comp_name<-mesh2name(as.character(disease_harmonic_average_score$mesh))
  
  rownames(disease_harmonic_average_score)<-disease_harmonic_average_score$mesh
  chemicals<-intersect(rownames(disease_harmonic_average_score),rownames(chem_disease))
  
  disease_labels<-chem_disease[chemicals,disease_name]
  disease_harmonic_average_score<-disease_harmonic_average_score[chemicals,]
  if (all(rownames(disease_harmonic_average_score)==rownames(disease_labels)))disease_harmonic_average_score$status_ctd<-disease_labels
  return(disease_harmonic_average_score)
}


disease_score_metabolic_syndrome<-disease_score_compile(metabolic_syndrome_scores,'metabolic syndrome')

disease_score_diabetes<-disease_score_compile(diabetes_scores,'diabetes mellitus type 2')

disease_score_atherosclerosis<-disease_score_compile(atheroscleorosis_scores,'atherosclerosis')

write.csv(disease_score_metabolic_syndrome,file = 'output/excel_files/metabolic_syndrome_scores.csv')
write.csv(disease_score_atherosclerosis,file = 'output/excel_files/atherosclerosis_scores.csv')
write.csv(disease_score_diabetes,file = 'output/excel_files/diabetes_scores.csv')


# ROC vlaidation with CTD -------------------------------------------------


#performing ROC analysis between binray data at chem disease and the predicted class probabilities
load('output/chem2disease.RData')
# x is the list of disease scores (class probabilties) realted to each data layer
# disease is the name of the disease
roc_auc<-function(x,disease){
compounds<-intersect(rownames(x),rownames(chem_disease))
x<-x[,grep(x=colnames(x),pattern='positive')]
x<-x[compounds]
labels<-chem_disease[compounds,disease]
require(PRROC)
roc<-roc.curve(scores.class0 = x, weights.class0 = labels,curve = T)
return(roc$auc)
}


Roc_metabolic_syndrome<-sapply(metabolic_syndrome_scores, roc_auc,'metabolic syndrome')
Roc_diabetes<-sapply(diabetes_scores, roc_auc,'diabetes mellitus type 2')
Roc_atherosclerosis<-sapply(atheroscleorosis_scores, roc_auc,'atherosclerosis')

roc_aucs<-do.call(rbind,list(Roc_metabolic_syndrome,Roc_diabetes,Roc_atherosclerosis))
save(roc_aucs,file = 'output/roc_disease_scores.RData')



# heatmap plot ------------------------------------------------------------

load('output/roc_disease_scores.RData')
library(gplots)
rownames(roc_aucs)<-c('metabolic_syndrome','diabetes','atherosclerosis')
heatmap.2(roc_aucs,trace = 'none', col = colorRampPalette(c('white', 'red'))(12),margins = c(20,20),
              density.info = 'none',srtCol = 45,keysize = 1,key.xlab = 'AUC-ROC',key.title = NA,cexRow = 1,cexCol = 1.2)


write.csv(roc_aucs,file = 'output/excel_files/Roc_analysis_CTD_class_prob.csv')
