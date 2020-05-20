source('functions/glm_functions.R')                            #repeated_k_fold_CV function
load("output/train_set_metabolic_syndrome_binary.RData") 
# Training the models from all layers
library(caret)
library(doParallel)

train_control <- trainControl(method = "repeatedcv",
                              number = 5,
                              repeats =2,
                              search = "random",
                              verboseIter = F,
                              classProbs=TRUE,  
                              summaryFunction=twoClassSummary)
# 1 performed as seperate experiments on sampo
cl<-makeCluster(length(train_set_metabolic_syndrome))
registerDoParallel(cl)
CV_results<-foreach(i=1:length(train_set_metabolic_syndrome),.packages=c('caret')) %dopar% caret_based_two_class_glm_repeated_k_fold_CV(train_set_metabolic_syndrome[[i]]$x,train_set_metabolic_syndrome[[i]]$y,5,5,neg_class = 'edc_negative_metabolic_syndrome',train_control = train_control)
stopCluster(cl)
names(CV_results)<-names(train_set_metabolic_syndrome)
save(CV_results,file="output/glm_models/CV_two_class/metabolic_syndrome/metabolic_syndrome_stratified_glm_models_kfold_CV_two_class1.RData") # for atherosclerosis 
#
# 2 performed as seperate experiments on sampo
cl<-makeCluster(length(train_set_metabolic_syndrome))
registerDoParallel(cl)
CV_results<-foreach(i=1:length(train_set_metabolic_syndrome),.packages=c('caret')) %dopar% caret_based_two_class_glm_repeated_k_fold_CV(train_set_metabolic_syndrome[[i]]$x,train_set_metabolic_syndrome[[i]]$y,5,5,neg_class = 'edc_negative_metabolic_syndrome',train_control = train_control)
stopCluster(cl)
names(CV_results)<-names(train_set_metabolic_syndrome)
save(CV_results,file="output/glm_models/CV_two_class/metabolic_syndrome/metabolic_syndrome_stratified_glm_models_kfold_CV_two_class2.RData") # for atherosclerosis 
#
# 3 performed as seperate experiments on sampo
cl<-makeCluster(length(train_set_metabolic_syndrome))
registerDoParallel(cl)
CV_results<-foreach(i=1:length(train_set_metabolic_syndrome),.packages=c('caret')) %dopar% caret_based_two_class_glm_repeated_k_fold_CV(train_set_metabolic_syndrome[[i]]$x,train_set_metabolic_syndrome[[i]]$y,5,5,neg_class = 'edc_negative_metabolic_syndrome',train_control = train_control)
stopCluster(cl)
names(CV_results)<-names(train_set_metabolic_syndrome)
save(CV_results,file="output/glm_models/CV_two_class/metabolic_syndrome/metabolic_syndrome_stratified_glm_models_kfold_CV_two_class3.RData") # for atherosclerosis 
#
# 4 performed as seperate experiments on sampo
cl<-makeCluster(length(train_set_metabolic_syndrome))
registerDoParallel(cl)
CV_results<-foreach(i=1:length(train_set_metabolic_syndrome),.packages=c('caret')) %dopar% caret_based_two_class_glm_repeated_k_fold_CV(train_set_metabolic_syndrome[[i]]$x,train_set_metabolic_syndrome[[i]]$y,5,5,neg_class = 'edc_negative_metabolic_syndrome',train_control = train_control)
stopCluster(cl)
names(CV_results)<-names(train_set_metabolic_syndrome)
save(CV_results,file="output/glm_models/CV_two_class/metabolic_syndrome/metabolic_syndrome_stratified_glm_models_kfold_CV_two_class4.RData") # for atherosclerosis 
