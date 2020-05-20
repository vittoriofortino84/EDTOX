

# 1.Glm modeling on all toxicogenimics data layers ---------------------------------------------------------------------
source("functions/glm_functions.R") #elastic_net_model_function

#load('outputData/glm/integrated_fgsea_results_edc_decoy.RData')    # Training set of the models from all layers
load('outputData/glm/integrated_fgsea_results_edc_decoy_moa.RData')    # Training set of the models from all layers

library(caret)
library(doParallel)
# source('functions/annotation_functions.R')   #all pathways led to positive NES for at least one EDC so these lines were commented
# names(fgs)<-name_correct(names(fgs))
# for (i in 1:length(fgs)){
#   net<-fgs[[i]]
#   ind_edc<-which(net$y %in% 'edc')
#   x<-net$x[ind_edc,]
#   max_ptw<-apply(x,2,max)
#   keep_ind<-which(max_ptw>0)
#   rmv_ptw<-names(which(max_ptw<0))
#   fgs[[i]]$x<-fgs[[i]]$x[,keep_ind]
#   fgs[[i]]$removed_pathways<-rmv_ptw
# }
train_control <- trainControl(method = "repeatedcv",
                              number = 5,
                              repeats =2,
                              search = "random",
                              verboseIter = FALSE,
                              classProbs=TRUE,  
                              summaryFunction=twoClassSummary)
cl<-makeCluster(length(fgs))
registerDoParallel(cl)
modls<-foreach(i=1:length(fgs),.packages='caret') %dopar% elastic_net_model(fgs[[i]]$x,fgs[[i]]$y,train_control) # model for each layer
stopCluster(cl)
names(modls)<-names(fgs)
#save(modls,file="outputData/glm/glm_models_final.RData") #NES none_scaled
save(modls,file="outputData/glm/glm_models_final_moa.RData") #NES none_scaled with only MOA pathways
parameters<-as.data.frame(sapply(modls, function(x)x$modl$bestTune))
xlsx::write.xlsx(parameters,file = 'outputData/excel_files/glm_params.xlsx',sheetName ='hyper_parameters')


# the mies are in a binary matrix
# 2.GLM modeling on MIEs binary  data-------------------------------------------------

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
modl_mie<- elastic_net_model(training_set,response,train_control) # model for MIEs
save(modl_mie,file="outputData/glm/glm_models_mie.RData") #NES none_scaled


