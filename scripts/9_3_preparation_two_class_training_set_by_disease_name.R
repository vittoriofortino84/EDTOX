source('functions/general_functions.R')                         #vam2comp, comp2disease and prepare_train_set  functions
source('functions/annotation_functions.R')                      #name_correct function
#load('input/integrated_fgsea_allcompounds_results_final.RData') #NES scores from phase 1 pipeline for all 12000 compounds and all networks
load('input/integrated_NES_scores_for_disease_biomarkers.RData') #NES scores from phase 1 pipeline for all 12000 compounds and all networks

names(fgs)<-name_correct(names(fgs))

edc<-vam2comp(probmin = 0.85,probmax=1)                           #harmonic vam score between 0.85 and 1.0 for edcs

# atherosclerosis
disease_finder('athero')                                         # to find the index and exact name of the disease in the binary matrix
athero<-comp2disease(edc,decoy = NULL,"atherosclerosis")         # making the response label vector for disease
train_set_athero<-lapply(fgs,prepare_train_set,athero,type='twoclass')    # making x,y data set edc + and -
save(train_set_athero,file = 'output/train_set_atherosclerosis_binary.RData')

#diabetes mellitus type 2
disease_finder('diabet')
diabetes_2<-comp2disease(edc,decoy = NULL,"diabetes mellitus type 2",
                       disease_label = 'diabetes_mellitus_type_2')         # making the response vector for disease 
train_set_diabetes_2<-lapply(fgs,prepare_train_set,diabetes_2,type='twoclass')    # making x,y data set edc + and - for all layers
save(train_set_diabetes_2,file = 'output/train_set_diabetes_2_binary.RData')

#metabolic syndrome
disease_finder('metabolic')
metabolic<-comp2disease(edc,decoy = NULL,"metabolic syndrome",
                         disease_label = 'metabolic_syndrome')         # making the response vector for disease 
train_set_metabolic_syndrome<-lapply(fgs,prepare_train_set,metabolic,type='twoclass')    # making x,y data set edc + and - for all layers
save(train_set_metabolic_syndrome,file = 'output/train_set_metabolic_syndrome_binary.RData')

#Coronary Artery Disease
disease_finder('coronary')
cad<-comp2disease(edc,decoy = NULL,"coronary artery disease",
                        disease_label = 'coronary_artery_disease')         # making the response vector for disease 
train_set_CAD<-lapply(fgs,prepare_train_set,cad,type='twoclass')    # making x,y data set edc + and - for all layers
save(train_set_CAD,file = 'output/train_set_coronary_artery_Disease_binary.RData')

