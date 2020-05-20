# For the data set for disease the pathways with length of more than 50 were removed
load('outputData/class_probilities/integrated_fgsea_allcompounds_results_final_moa.RData')
source('functions/annotation_functions.R')
all_pathways<-as.data.frame(unique(unlist(lapply(fgs, function(x)colnames(x$x)))));colnames(all_pathways)<-'pathways'
all_pathways$annotations<-pathway_annotate(as.character(all_pathways$pathways))
load('outputData/toxdb2gene_final.RData')
pathway_length<-as.data.frame(sapply(moa_pathways, length));colnames(pathway_length)<-'pathways_length';pathway_length$pathways<-rownames(pathway_length)
all_pathways<-merge(all_pathways,pathway_length)
write.csv(all_pathways,file = 'outputData/excel_files/for_disease_manual_curation.csv')
####
pathway_disease<-xlsx::read.xlsx('outputData/excel_files/for_disease_manual_curation.xlsx',sheetIndex = 1)
pathway_disease$status_in_model[grep(x=pathway_disease$annotations,pattern ='cancer',ignore.case = T) ]=0
pathway_disease$status_in_model[grep(x=pathway_disease$annotations,pattern ='diabet',ignore.case = T) ]=0
pathway_disease$status_in_model[grep(x=pathway_disease$annotations,pattern ='atheros',ignore.case = T) ]=0
pathway_disease$status_in_model[grep(x=pathway_disease$annotations,pattern ='disease',ignore.case = T) ]=0
pathway_disease$status_in_model[which(pathway_disease$pathways_length>50)]=0
pathway_disease<-pathway_disease$pathways[which(pathway_disease$status_in_model==1)];pathway_disease<-as.character(pathway_disease)
source('functions/glm_functions.R')
load("outputData/new199edc_1462dec.RData")             # edc decoys list
setwd('outputData/all_compounds_fgsea/')
li<-list.files(path='.',pattern = 'RData')
edcs<-names(new_edc199_decoy_1462)[1:199]
fgs<-data_preparation(li,edcs,pathway_disease)      #2nd with selected pathways
setwd('../..')

names(fgs)<-name_correct(names(fgs))
save(fgs,file='outputData/integrated_NES_scores_for_disease_biomarkers.RData')


