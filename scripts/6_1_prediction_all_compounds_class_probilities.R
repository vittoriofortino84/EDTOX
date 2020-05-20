# 2.prediction of class probility for all CTD compounds based on glm models of EDC decoy for each layer --------------------------------------------------------------

source('functions/glm_functions.R')
#load("outputData/class_probilities/integrated_fgsea_allcompounds_results_final.RData") #fgs
#load("outputData/glm/glm_models_final.RData") #modls
load("outputData/class_probilities/integrated_fgsea_allcompounds_results_final_moa.RData") #fgsea NES scores
load("outputData/glm/glm_models_final_moa.RData") #modls  
all_compounds_prob<-prob_predict(modls,fgs)
#save(all_compounds_prob,file='outputData/class_probilities/class_probilities_all_compounds_final.RData')
save(all_compounds_prob,file='outputData/class_probilities/class_probilities_all_compounds_final_moa.RData')
# 
