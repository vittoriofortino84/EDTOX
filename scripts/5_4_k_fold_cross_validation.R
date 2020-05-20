
# FOR pathway level data --------------------------------------------------


library(caret)
library(doParallel)
source("functions/glm_functions.R")
source('functions/general_functions.R')
load("outputData/new199edc_1462dec.RData")

#load("outputData/glm/integrated_fgsea_results_edc_decoy.RData")

load("outputData/glm/integrated_fgsea_results_edc_decoy_moa.RData")
train_control <- trainControl(method = "repeatedcv",
                              number = 5,
                              repeats =2,
                              search = "random",
                              verboseIter = FALSE,
                              classProbs=TRUE,  
                              summaryFunction=twoClassSummary)
# #stability test based on k-fold cross_validation for all layers
cl<-makeCluster(length(fgs))
registerDoParallel(cl)
modls<-foreach(i=1:length(fgs),.packages='caret') %dopar% caret_based_two_class_glm_repeated_k_fold_CV(fgs[[i]]$x,fgs[[i]]$y,5,5,neg_class = 'decoy',train_control = train_control) # 
names(modls)<-name_correct(names(fgs))
stopCluster(cl)
save(modls,file="outputData/glm/stratified_networks_CV/k_fold_cross_validation_networks1.RData") 
cl<-makeCluster(length(fgs))
registerDoParallel(cl)
modls<-foreach(i=1:length(fgs),.packages='caret') %dopar% caret_based_two_class_glm_repeated_k_fold_CV(fgs[[i]]$x,fgs[[i]]$y,5,5,neg_class = 'decoy',train_control = train_control) # 
names(modls)<-name_correct(names(fgs))
stopCluster(cl)
save(modls,file="outputData/glm/stratified_networks_CV/k_fold_cross_validation_networks2.RData") 
cl<-makeCluster(length(fgs))
registerDoParallel(cl)
modls<-foreach(i=1:length(fgs),.packages='caret') %dopar% caret_based_two_class_glm_repeated_k_fold_CV(fgs[[i]]$x,fgs[[i]]$y,5,5,neg_class = 'decoy',train_control = train_control) # 
names(modls)<-name_correct(names(fgs))
stopCluster(cl)
save(modls,file="outputData/glm/stratified_networks_CV/k_fold_cross_validation_networks3.RData") 
cl<-makeCluster(length(fgs))
registerDoParallel(cl)
modls<-foreach(i=1:length(fgs),.packages='caret') %dopar% caret_based_two_class_glm_repeated_k_fold_CV(fgs[[i]]$x,fgs[[i]]$y,5,5,neg_class = 'decoy',train_control = train_control) # 
names(modls)<-name_correct(names(fgs))
stopCluster(cl)
save(modls,file="outputData/glm/stratified_networks_CV/k_fold_cross_validation_networks4.RData") 

#####old unstratified splitting CV
# modls<-foreach(i=1:length(fgs),.export = 'glm_model_stability',.packages='caret') %dopar% glm_model_stability(fgs[[i]]$x,fgs[[i]]$y,5,10,neg_class = 'decoy') # 5 fold cross_validation and one repeat for each layer
# save(modls,file="outputData/glm/k_fold_cross_validation_networks.RData") 



# FOR MIE_LEVEL_DATA ------------------------------------------------------


source("functions/glm_functions.R")
library(caret)
load('outputData/chem2gene_no_out.RData')
load('outputData/new199edc_1462dec.RData')
all_genenes<-unique(unlist(sapply(chem2gene,function(x)x)))
binary_MIE<-matrix(0, nrow = length(chem2gene), ncol = length(all_genenes))
for  (i in 1:length(chem2gene)){
  binary_MIE[i,which(all_genenes %in% chem2gene[[i]])]<-1
}
colnames(binary_MIE)<-all_genenes
rownames(binary_MIE)<-names(chem2gene)
training_set<-binary_MIE[names(new_edc199_decoy_1462),]
response<-rep('decoy',nrow(training_set)) # y vector
response[1:199]<-'edc'
train_control <- trainControl(method = "repeatedcv",
                              number = 5,
                              repeats =2,
                              search = "random",
                              verboseIter = FALSE,
                              classProbs=TRUE,
                              summaryFunction=twoClassSummary)
stability<-caret_based_two_class_glm_repeated_k_fold_CV(training_set,response,5,5,neg_class = 'decoy',train_control=train_control)
save(stability,file='outputData/glm/stratified_MIE_CV/k_fold_cross_validation_MIEs1.RData')
stability<-caret_based_two_class_glm_repeated_k_fold_CV(training_set,response,5,5,neg_class = 'decoy',train_control=train_control)
save(stability,file='outputData/glm/stratified_MIE_CV/k_fold_cross_validation_MIEs2.RData')
stability<-caret_based_two_class_glm_repeated_k_fold_CV(training_set,response,5,5,neg_class = 'decoy',train_control=train_control)
save(stability,file='outputData/glm/stratified_MIE_CV/k_fold_cross_validation_MIEs3.RData')
stability<-caret_based_two_class_glm_repeated_k_fold_CV(training_set,response,5,5,neg_class = 'decoy',train_control=train_control)
save(stability,file='outputData/glm/stratified_MIE_CV/k_fold_cross_validation_MIEs4.RData')

# stability<-glm_model_stability(training_set,response,5,10,neg_class = 'decoy')  #unstratified
# save(stability,file='outputData/glm/k_fold_cross_validation_MIEs.RData')       

