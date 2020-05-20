
# harmonic and average ----------------------------------------------------------------

rm(list=ls())
source('functions/vam_functions.R')
source('functions/annotation_functions.R')
edcs<-read.csv('inputData/DEDuCT_ChemicalBasicInformation.csv',header = TRUE) #list of edcs from DEDuCT
edcs<-unique(edcs$CAS.Number)

vam<-harmonic_vam_score(as.character(edcs),type = 'cas')
vam<-vam[!is.na(vam$VAM_score),]
vam$mesh<-cas2mesh(rownames(vam))
load('outputData/new199edc_1462dec.RData')
vam$mesh<-cas2mesh(rownames(vam))
vam$is_in_training<-rep('unk',nrow(vam))
vam$is_in_training[vam$mesh %in% names(new_edc199_decoy_1462)[1:199]]<-'edc'
vam$is_in_training[vam$mesh %in% names(new_edc199_decoy_1462)[200:1661]]<-'decoy'
vam$comp_names<-mesh2name(vam$mesh)
print(paste('accuracy=',length(which(vam$VAM_score[vam$is_in_training=='unk']>0.85))/length(which(vam$is_in_training=='unk'))))
xlsx::write.xlsx2(vam,file='outputData/excel_files/DeDucts_predicted_harmonic_moa.xlsx')
colnames(vam)[which(colnames(vam) %in% 'VAM_score')]<-'harmonic_VAM_score'
vam$average_VAM_score<-average_vam_score(vam$mesh)$VAM_score
xlsx::write.xlsx2(vam,file='outputData/excel_files/DeDucts_predicted_all_scores_moa.xlsx')

# 
# # average -----------------------------------------------------------------
# rm(list=ls())
# source('functions/vam_functions.R')
# source('functions/annotation_functions.R')
# edcs<-read.csv('inputData/DEDuCT_ChemicalBasicInformation.csv',header = TRUE) #list of edcs from DEDuCT
# edcs<-unique(edcs$CAS.Number)
# vam<-average_vam_score(as.character(edcs),type = 'cas')
# 
# vam<-vam[!is.na(vam$VAM_score),]
# vam$mesh<-cas2mesh(rownames(vam))
# load('outputData/new199edc_1462dec.RData')
# vam$mesh<-cas2mesh(rownames(vam))
# vam$is_in_training<-rep('unk',nrow(vam))
# vam$is_in_training[vam$mesh %in% names(new_edc199_decoy_1462)[1:199]]<-'edc'
# vam$is_in_training[vam$mesh %in% names(new_edc199_decoy_1462)[200:1661]]<-'decoy'
# vam$comp_names<-mesh2name(vam$mesh)
# print(paste('accuracy=',length(which(vam$VAM_score[vam$is_in_training=='unk']>0.8))/length(which(vam$is_in_training=='unk'))))
# xlsx::write.xlsx2(vam,file='outputData/excel_files/DeDucts_predicted_average_moa.xlsx')
# 
