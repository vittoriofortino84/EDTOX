load("output/train_set_atherosclerosis_binary.RData")
source("functions/glm_functions.R") #elastic_net_model_function
# 2.Glm modeling on all toxicogenimics data layers ---------------------------------------------------------------------
# Training the models from all layers

library(caret)
library(doParallel)
train_control <- trainControl(method = "repeatedcv",
                              number = 5,
                              repeats =2,
                              search = "random",
                              verboseIter = FALSE,
                              classProbs=TRUE,  
                              summaryFunction=twoClassSummary)
cl<-makeCluster(length(train_set_athero))
registerDoParallel(cl)
atherosclerosis_two_class<-foreach(i=1:length(train_set_athero),.packages='caret') %dopar% elastic_net_model_two_class(train_set_athero[[i]]$x,
                                                                                                                       train_set_athero[[i]]$y,
                                                                                                                       train_control) # model for each layer
stopCluster(cl)
names(atherosclerosis_two_class)<-names(train_set_athero)
save(atherosclerosis_two_class,file="output/glm_models/two_class_glm_models_atherosclerosis.RData") #models for atherosclerosis and pathways NES scores