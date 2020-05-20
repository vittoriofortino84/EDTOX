source('functions/glm_functions.R')                            #multi_class_glm_repeated_k_fold_CV function
load("output/train_set_atherosclerosis_binary.RData")
# Training the models from all layers
library(caret)
library(doParallel)
library(MLmetrics)

train_control <- trainControl(method = "repeatedcv",
                              number = 5,
                              repeats =2,
                              search = "random",
                              verboseIter = F,
                              classProbs=TRUE,  
                              summaryFunction=twoClassSummary)
# 1 performed as seperate experiments on sampo
cl<-makeCluster(length(train_set_athero))
registerDoParallel(cl)
CV_results<-foreach(i=1:length(train_set_athero),.packages=c('caret','MLmetrics','ROCR')) %dopar% caret_based_two_class_glm_repeated_k_fold_CV(train_set_athero[[i]]$x,train_set_athero[[i]]$y,5,5,neg_class = "edc_negative_atherosclerosis",train_control = train_control)
stopCluster(cl)
names(CV_results)<-names(train_set_athero)
save(CV_results,file="output/glm_models/CV_two_class/atherosclerosis/stratified/stratified_glm_models_atherosclerosis_kfold_CV_two_class_1.RData") # for atherosclerosis 
# 2
cl<-makeCluster(length(train_set_athero))
registerDoParallel(cl)
CV_results<-foreach(i=1:length(train_set_athero),.packages=c('caret','MLmetrics','ROCR')) %dopar% caret_based_two_class_glm_repeated_k_fold_CV(train_set_athero[[i]]$x,train_set_athero[[i]]$y,5,5,neg_class = "edc_negative_atherosclerosis",train_control = train_control)
stopCluster(cl)
names(CV_results)<-names(train_set_athero)
save(CV_results,file="output/glm_models/CV_two_class/atherosclerosis/stratified/stratified_glm_models_atherosclerosis_kfold_CV_two_class_2.RData") # for atherosclerosis 
# 3
cl<-makeCluster(length(train_set_athero))
registerDoParallel(cl)
CV_results<-foreach(i=1:length(train_set_athero),.packages=c('caret','MLmetrics','ROCR')) %dopar% caret_based_two_class_glm_repeated_k_fold_CV(train_set_athero[[i]]$x,train_set_athero[[i]]$y,5,5,neg_class = "edc_negative_atherosclerosis",train_control = train_control)
stopCluster(cl)
names(CV_results)<-names(train_set_athero)
save(CV_results,file="output/glm_models/CV_two_class/atherosclerosis/stratified/stratified_glm_models_atherosclerosis_kfold_CV_two_class_3.RData") # for atherosclerosis 
# 4
cl<-makeCluster(length(train_set_athero))
registerDoParallel(cl)
CV_results<-foreach(i=1:length(train_set_athero),.packages=c('caret','MLmetrics','ROCR')) %dopar% caret_based_two_class_glm_repeated_k_fold_CV(train_set_athero[[i]]$x,train_set_athero[[i]]$y,5,5,neg_class = "edc_negative_atherosclerosis",train_control = train_control)
stopCluster(cl)
names(CV_results)<-names(train_set_athero)
save(CV_results,file="output/glm_models/CV_two_class/atherosclerosis/stratified/stratified_glm_models_atherosclerosis_kfold_CV_two_class_4.RData") # for atherosclerosis 


## BOOTSTRAP
# source('functions/glm_functions.R')
# load('output/train_set_atherosclerosis_binary.RData')
# library(caret)
# train_control <- trainControl(method = "repeatedcv",
#                               number = 5,
#                               repeats =2,
#                               search = "random",
#                               verboseIter = FALSE,
#                               classProbs=TRUE,  
#                               summaryFunction=twoClassSummary)
# 
# 
# cl<-makeCluster(length(train_set_athero))
# registerDoParallel(cl)
# bootstrap_results<-foreach(i=1:length(train_set_athero),.packages=c('caret','MLmetrics','ROCR')) %dopar% two_class_elastic_net_bootstrap(train_set_athero[[i]]$x,as.factor(train_set_athero[[i]]$y),train_control,percentage = 70,iter=100) #  for each layer
# stopCluster(cl)
# names(bootstrap_results)<-names(train_set_athero)
# save(bootstrap_results,file="output/glm_models/glm_models_atherosclerosis_bootstrap_two_class.RData") #models for atherosclerosis and pathways NES scores
