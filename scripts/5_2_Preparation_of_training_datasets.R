
# 1.For only EDCs and Decoys preparation of data from FGSEA files from toxicogenomics pipeline -------------------------------------------------------------
#Data preparation from customized networks
source("functions/glm_functions.R")                    # data_preparation function
source('functions/general_functions.R')
load("outputData/new199edc_1462dec.RData")             # edc decoys 
load("outputData/toxdb2gene_final.RData")              # pathways
setwd('outputData/edc_decoy_fgsea/')
li<-list.files(path='.',pattern = 'RData')
edcs<-names(new_edc199_decoy_1462)[1:199]
#fgs<-data_preparation(li,edcs,toxdb_names)             #1st run with all pathways
fgs<-data_preparation(li,edcs,names(moa_pathways))      #2nd with selected pathways

used_pathways<-lapply(fgs, function(x)colnames(x$x))
names(fgs)<-name_correct(names(fgs))
setwd('../..')
#save(fgs,file='outputData/glm/integrated_fgsea_results_edc_decoy.RData')
#save(used_pathways,file='outputData/glm/used_pathways_glm.RData')
GO<-sapply(fgs,function(x)length(grep(x=colnames(x$x),pattern = 'GO',ignore.case = F)))
wiki<-sapply(fgs,function(x)length(grep(x=colnames(x$x),pattern = 'wiki',ignore.case = T)))
reactome<-sapply(fgs,function(x)length(grep(x=colnames(x$x),pattern = 'HSA',ignore.case = T)))
kegg<-sapply(fgs,function(x)length(grep(x=colnames(x$x),pattern = 'kegg',ignore.case = T)))
all_paths<-do.call(cbind,list(GO,kegg,reactome,wiki))
colnames(all_paths)<-c('GO','Kegg','reactome','wiki')
write.csv(all_paths,file = 'outputData/excel_files/pathways_per_network.csv')
save(fgs,file='outputData/glm/integrated_fgsea_results_edc_decoy_moa.RData')
save(used_pathways,file='outputData/glm/used_pathways_glm_moa.RData')


# 2. For all around 12000 compounds of CTD --------------------------------
# preparation of integrated dataframe from fgsea files performed on SAMPO cluster at UEF--------------------------------------------------------------------
source('functions/glm_functions.R')
source('functions/general_functions.R')
load('outputData/new199edc_1462dec.RData')
load('outputData/toxdb2gene_final.RData')
setwd('outputData/all_compounds_fgsea/')
li<-list.files(path='.',pattern = 'RData')
edcs<-names(new_edc199_decoy_1462)[1:199]
#fgs<-data_preparation(li,edcs,toxdb_names) #NES_scores with  no scaling
fgs<-data_preparation(li,edcs,names(moa_pathways)) #NES_scores with  no scaling the function is in glm_functions.R
names(fgs)<-name_correct(names(fgs)) #gneral_functions.R
setwd("../..")
#save(fgs,file='outputData/class_probilities/integrated_fgsea_allcompounds_results_final.RData')    #all compounds in CTD with NES scores
save(fgs,file='outputData/class_probilities/integrated_fgsea_allcompounds_results_final_moa.RData') #all compounds in CTD with NES scores
#

