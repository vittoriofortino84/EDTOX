rm(list=ls())
library(magrittr)
source('functions/annotation_functions.R')
source('functions/vam_functions.R')
nick<-xlsx::read.xlsx2(file='inputData/Nick_Amir_Testing_30_10_2019.xlsx',sheetIndex = 1)
colnames(nick)<-c('CAS.Number','Name','SMILES..Canonical.')
all_nick<-as.character(nick$CAS.Number) %>% strsplit('&') %>% unlist() %>% 
  gsub(pattern = ' ',replacement = '') %>% unique()
all_nick<-cas2mesh(all_nick,return_back_NA = 'F')
all_nick[ all_nick %in% "3825-26-1"]<-'C023036'    #manually found in pubmed
all_nick[ all_nick %in% "307-55-1"]<-'C522391'      #manually found in pubmed
all_nick[ all_nick %in% "56558-16-8"]<-'C417207'    #manually found in pubmed
all_nick[ all_nick %in% "31508-00-6"]<-'C070055' #manually found in pubmed
all_nick[ all_nick %in% "152969-11-4"]<-'C111120' #manually found in pubmed
all_nick[ all_nick %in% "068515-48-0"]<-'C012125' #manually found in pubmed
all_nick[ all_nick %in% '028553-12-0']<-'C012125' #manually found in pubmed
all_nick[ all_nick %in% "84852-15-3"]<-'C025256' #manually found in pubmed
all_nick[ all_nick %in% "38380-01-7"]<-'C500336' #manually found in pubmed
all_nick[ all_nick %in% "111988-49-9"]<-'C417209' #manually found in pubmed
all_nick[ all_nick %in% "843-55-0" ]<-'C570106' #manually found in pubmed
all_nick[ all_nick %in% "13595-25-0"]<-'C109393' #manually found in pubmed
all_nick[ all_nick %in% "56-35-9"]<-'C005961' #manually found in pubmed
all_nick[ all_nick %in% "1461-22-9"]<-''
all_nick[ all_nick %in%  "1478-61-1"]<-''
all_nick[ all_nick %in% "5613-46-7"]<-''
all_nick[ all_nick %in% "66640-59-3"]<-''
all_nick[ all_nick %in% "375-73-5"]<-''
all_nick[ all_nick %in% "355-43-1"]<-''
all_nick[ all_nick %in% "1895-26-7" ]<-''
all_nick[ all_nick %in% "2795-39-3"]<-''
all_nick[ all_nick %in% "507-63-1"]<-''
all_nick[ all_nick %in% "57678-03-2"]<-''
all_nick[ all_nick %in% "1461-23-0"]<-''
all_nick[ all_nick %in% "620-92-8"]<-''
all_nick[ all_nick %in% "2081-08-5"]<-''
all_nick[ all_nick %in% "1571-75-1"]<-''
all_nick[ all_nick %in% "2081-08-5"]<-''
all_nick[ all_nick %in% "29426-78-6"]<-''
all_nick[ all_nick %in% "6073-11-6"]<-''
all_nick[ all_nick %in% "77-40-7"]<-''
all_nick<-all_nick[!all_nick=='']
external_validation<-as.data.frame(all_nick)
colnames(external_validation)<-c('mesh_ID')
external_validation<-cbind(external_validation,harmonic_vam_score(as.character(external_validation$mesh_ID)))
colnames(external_validation)[which(colnames(external_validation) %in% 'VAM_score')]<-'Harmonic_VAM_score'
external_validation$Average_VAM_score<-average_vam_score(as.character(external_validation$mesh_ID))$VAM_score
load('outputData/new199edc_1462dec.RData')
external_validation$status<-'unk'
external_validation$status[which(external_validation$mesh_ID %in% names(new_edc199_decoy_1462)[1:199])]<-'edc'
external_validation$status[which(external_validation$mesh_ID %in% names(new_edc199_decoy_1462)[200:1661])]<-'decoy'
external_validation$compname<-mesh2name(as.character(external_validation$mesh_ID))
external_validation$cas<-mesh2cas(as.character(external_validation$mesh_ID))
external_validation<-unique(external_validation)
external_validation<-external_validation[!is.na(external_validation$Average_VAM_score),]
colnames(external_validation)<-gsub(x=colnames(external_validation),pattern = '__','_')
colnames(external_validation)<-gsub(x=colnames(external_validation),pattern = '_HEPG2_HEPG2','_HEPG2_1_day')
library(tidyr)
external_validation<-external_validation[external_validation$status=='unk',]
Deduct<-xlsx::read.xlsx2('outputData/excel_files/DeDucts_predicted_all_scores_moa.xlsx',sheetIndex = 1)
Deduct<-Deduct[Deduct$is_in_training=='unk',]
common_Deduct_experts_compounds<-intersect(Deduct$mesh,external_validation$mesh_ID)
Deduct<-Deduct[-which(Deduct$mesh %in% common_Deduct_experts_compounds),]
colnames(Deduct)<-gsub(x=colnames(Deduct),pattern = '_HEPG2_HEPG2','_HEPG2_1_day')
colnames(Deduct)<-gsub(x=colnames(Deduct),pattern = 'Consensus_Rat_invitro_Drug.Matrix__TG_GATEs','Consensus_Rat_invitro_Drug Matrix_TG_GATEs')
colnames(Deduct)<-gsub(x=colnames(Deduct),pattern = 'comp_names','compname')
colnames(Deduct)<-gsub(x=colnames(Deduct),pattern = 'harmonic_VAM_score','Harmonic_VAM_score')
colnames(Deduct)<-gsub(x=colnames(Deduct),pattern = 'average_VAM_score','Average_VAM_score')
external_validation<-external_validation[,-which(colnames(external_validation) %in% c("mesh_ID","cas","status"))]
Deduct<-Deduct[,-which(colnames(Deduct) %in% c("mesh","X.","is_in_training"))]
Deduct$status<-'dimgrey';external_validation$status<-'blue'
Deduct<-Deduct[,intersect(colnames(external_validation),colnames(Deduct))]
external_validation<-external_validation[,intersect(colnames(external_validation),colnames(Deduct))]
merged<-rbind(external_validation,Deduct)
merged[merged=='']<-NA
merged<-na.omit(merged)
merged$Harmonic_EDC_score<-as.numeric(merged$Harmonic_VAM_score)
merged$Average_EDC_score<-as.numeric(merged$Average_VAM_score)
rownames(merged)<-name2mesh(as.character(merged$compname))
# adding disease scores
source('functions/disease_scores.R')
merged$Average_diabetes_T2_score<-mesh2disease_score(rownames(merged),'Diabetes Type 2')$average_disease_score
merged$Harmonic_diabetes_T2_score<-mesh2disease_score(rownames(merged),'Diabetes Type 2')$harmonic_disease_score
merged$Labels_diabetes_T2_CTD<-mesh2disease_score(rownames(merged),'Diabetes Type 2')$status_ctd


merged$Average_atherosclerosis_score<-mesh2disease_score(rownames(merged),'atherosclerosis')$average_disease_score
merged$Harmonic_atherosclerosis_score<-mesh2disease_score(rownames(merged),'atherosclerosis')$harmonic_disease_score
merged$Labels_atherosclerosis_CTD<-mesh2disease_score(rownames(merged),'atherosclerosis')$status_ctd


merged$Average_metabolic_syndrome_score<-mesh2disease_score(rownames(merged),'metabolic syndrome')$average_disease_score
merged$Harmonic_metabolic_syndrome_score<-mesh2disease_score(rownames(merged),'metabolic syndrome')$harmonic_disease_score
merged$Labels_metabolic_syndrome_CTD<-mesh2disease_score(rownames(merged),'metabolic syndrome')$status_ctd

heat_map_plot<-function(merged,nams,legend_title){
  columns_to_keep<-c('compname','status',nams)
  merged<-merged[,columns_to_keep]
  
  merged<-merged[order(merged$Harmonic_EDC_score,decreasing = T),]
  merged$turn<-1;merged$turn[75:nrow(merged)]<-2
  
  st1<-merged$status[merged$turn==1]
  st2<-merged$status[merged$turn==2]
 
  
  comp_lvls<-merged$compname
  subdata<-gather(merged,key='data_layer',
                  value = 'class_prob',-c(compname,status,turn))
  subdata$class_prob<-as.numeric(subdata$class_prob)

  library(reshape2)
  library(RColorBrewer)
  library(ggplot2)
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
  ggplot(subdata,aes(y = factor(data_layer,levels = rev(nams)), x = factor(compname,levels = comp_lvls),
                     fill = class_prob))+
    geom_tile()+
    scale_fill_gradientn(colours = myPalette(100))+xlab('')+ylab('')+labs(fill=legend_title)+
    theme_minimal()+ 
    guides(fill = guide_colourbar(barwidth = 2, barheight = 10))+
    facet_wrap( . ~ subdata$turn,nrow = 2,scales = 'free')+
    theme(axis.title = element_text(size = 16),strip.text = element_text(size = 16),
          legend.title  = element_text(size = 16),legend.text = element_text(size = 14),legend.position = 'right',
          axis.text.y = element_text(size=16),
          axis.text.x = element_text(size=16,angle = 45,hjust = 1))}




nams_edc_scores<-c(
               'Drug_Matrix_Rat_invitro_Single_Dose_1_day','Drug_Matrix_Rat_invivo_Single_Dose_1_day',    
               'Drug_Matrix_Rat_invivo_Repeated_Dose_3_days',  'Drug_Matrix_Rat_invivo_Repeated_Dose_5_days', 
               'TG_GATEs_Human_invitro_Single_Dose_1_day','TG_GATEs_Rat_invitro_Single_Dose_1_day',      
               'TG_GATEs_Rat_invivo_Single_Dose_Low_1_day','TG_GATEs_Rat_invivo_Single_Dose_Middle_1_day',
               'TG_GATEs_Rat_invivo_Single_Dose_High_1_day','TG_GATEs_Rat_invivo_Repeated_Dose_8_days',    
               'TG_GATEs_Rat_invivo_Repeated_Dose_15_days', 'TG_GATEs_Rat_invivo_Repeated_Dose_29_days',  
               'Consensus_Rat_invitro_Drug Matrix_TG_GATEs','Consensus_LINCS_HEPG2_1_day',                 
               'PPI_STRINGdb','Average_EDC_score','Harmonic_EDC_score' )




nams_disease_scores<-c(
               'Average_EDC_score','Harmonic_EDC_score',
               "Average_diabetes_T2_score", "Harmonic_diabetes_T2_score",'Labels_diabetes_T2_CTD',
               "Average_atherosclerosis_score","Harmonic_atherosclerosis_score", 'Labels_atherosclerosis_CTD',             
               "Average_metabolic_syndrome_score","Harmonic_metabolic_syndrome_score",'Labels_metabolic_syndrome_CTD'   )


save(heat_map_plot,merged,nams_disease_scores,nams_edc_scores,file = 'outputData/plots_tables_functions/Heatmap_validation/heatmap_edc_scores_disease.RData')
  
############### plot visualization

  rm(list=ls())
  load('outputData/plots_tables_functions/Heatmap_validation/heatmap_edc_scores_disease.RData')
  heat_map_plot(merged,nams_edc_scores,'class probability')
  heat_map_plot(merged,nams_disease_scores,'score')
