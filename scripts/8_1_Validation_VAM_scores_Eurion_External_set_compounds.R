# 1.PLOTTING merging test compounds from Nick Plant and DEDuCT--------------------------------------------------------------
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
merged$Harmonic_VAM_score<-as.numeric(merged$Harmonic_VAM_score)
merged$Average_VAM_score<-as.numeric(merged$Average_VAM_score)
merged<-merged[order(merged$Harmonic_VAM_score,decreasing = T),]
merged$turn<-1;merged$turn[75:nrow(merged)]<-2
print(length(which(merged$Harmonic_VAM_score > 0.6))/nrow(merged))

st1<-merged$status[merged$turn==1]
st2<-merged$status[merged$turn==2]
data_1<-merged[merged$turn==1,]
data_2<-merged[merged$turn==2,]

pplott<-function(data){
comp_lvls<<-merged$compname
subdata<-gather(data,key='data_layer',
                 value = 'class_prob',-c(compname,status,turn))
subdata$class_prob<-as.numeric(subdata$class_prob)
subdata$data_layer<-gsub(x=subdata$data_layer,pattern = "Harmonic_VAM_score",
                          replacement  ="Harmonic_EDC_score" )

subdata$data_layer<-gsub(x=subdata$data_layer,pattern = "Average_VAM_score",
                          replacement  ="Average_EDC_score" )

library(reshape2)
library(RColorBrewer)
library(ggplot2)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
source('functions/general_functions.R')
nams<<-factor(c('Drug_Matrix_Rat_invitro_Single_Dose_1_day','Drug_Matrix_Rat_invivo_Single_Dose_1_day',    
 'Drug_Matrix_Rat_invivo_Repeated_Dose_3_days',  'Drug_Matrix_Rat_invivo_Repeated_Dose_5_days', 
 'TG_GATEs_Human_invitro_Single_Dose_1_day','TG_GATEs_Rat_invitro_Single_Dose_1_day',      
 'TG_GATEs_Rat_invivo_Single_Dose_Low_1_day','TG_GATEs_Rat_invivo_Single_Dose_Middle_1_day',
 'TG_GATEs_Rat_invivo_Single_Dose_High_1_day','TG_GATEs_Rat_invivo_Repeated_Dose_8_days',    
 'TG_GATEs_Rat_invivo_Repeated_Dose_15_days', 'TG_GATEs_Rat_invivo_Repeated_Dose_29_days',  
 'Consensus_Rat_invitro_Drug Matrix_TG_GATEs','Consensus_LINCS_HEPG2_1_day',                 
 'PPI_STRINGdb','Average_EDC_score','Harmonic_EDC_score'))
net_colr<<-c('#99CCFF','#99CCFF','#99CCFF','#99CCFF',                     # Drug matrix hep,1,3,5
            '#99CC66','#99CC66',                                         # in vitro rat and human from TGgates
            '#99CC66','#99CC66','#99CC66','#99CC66','#99CC66','#99CC66', # in vivo  rat low,high,moddle,8,15,29 TGgates
            '#FF6640','#FF6640','#FF6640','purple','purple')
subdata<<-subdata
p1 <- ggplot(subdata,aes(y = factor(data_layer,levels = rev(nams)), x = factor(compname,levels = comp_lvls),
                         fill = class_prob))+
  geom_tile()+
  scale_fill_gradientn(colours = myPalette(100))+xlab('')+ylab('')+labs(fill='class probability')+
  #coord_equal()+ 
  #coord_fixed(ratio = .5)+
  #guides(shape=guide_legend())+
  theme_minimal()+ 
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10))+
  facet_wrap( . ~ merged$turn,nrow = 2,scales = 'free')+
  theme(axis.title = element_text(size = 16),strip.text = element_text(size = 16),
        legend.title  = element_text(size = 16),legend.text = element_text(size = 14),legend.position = 'right',
        axis.text.y = element_text(size=16,colour=rev(net_colr)),
        #axis.text.x = element_text(size=16,angle = 45,hjust = 1,colour=st))
        axis.text.x = element_text(size=16,angle = 45,hjust = 1))
return(p1)


}
#g1<-pplott(data_1,st1);g2<-pplott(data_2,st2)
#ggarrange(g1,g2,nrow = 2,common.legend = T,legend = 'bottom')
p<-pplott(merged)
print(p)
validation_data<-subdata
save(validation_data,nams,comp_lvls,file = 'outputData/plots_tables_functions/Heatmap_validation/validation_set.RData')

# 2.Saving ----------------------------------------------------------------


rm(list=ls())
library(magrittr)
source('functions/annotation_functions.R')
load('outputData/VAM_scores_hs_most_informative_layers_moa.RData')

nick<-xlsx::read.xlsx2(file='inputData/Nick_Amir_Testing_30_10_2019.xlsx',sheetIndex = 1)
colnames(nick)<-c('CAS.Number','Name','SMILES..Canonical.')
all_nick<-as.character(nick$CAS.Number) %>% strsplit('&') %>% unlist() %>% gsub(pattern = ' ',replacement = '') %>% unique()
all_nick<-cas2mesh(all_nick,return_back_NA = 'F')
all_nick[ all_nick %in% "3825-26-1"]<-'C023036'    #manually found in pubmed
all_nick[ all_nick %in% "307-55-1"]<-'C522391'
all_nick[ all_nick %in% "56558-16-8"]<-'C417207'
all_nick[ all_nick %in% "31508-00-6"]<-'C070055'
all_nick[ all_nick %in% "152969-11-4"]<-'C111120'
all_nick[ all_nick %in% "068515-48-0"]<-'C012125'
all_nick[ all_nick %in% '028553-12-0']<-'C012125'
all_nick[ all_nick %in% "84852-15-3"]<-'C025256'
all_nick[ all_nick %in% "38380-01-7"]<-'C500336'
all_nick[ all_nick %in% "111988-49-9"]<-'C417209'
all_nick[ all_nick %in% "843-55-0" ]<-'C570106'
all_nick[ all_nick %in% "13595-25-0"]<-'C109393'
all_nick[ all_nick %in% "56-35-9"]<-'C005961'
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

test_train_prdictions<-as.data.frame(all_VAM_score[which(names(all_VAM_score) %in% all_nick)])
test_train_prdictions<-merge(all_prob,test_train_prdictions,by=0,all=F)
test_train_prdictions$cas<-mesh2cas(test_train_prdictions$Row.names)
test_train_prdictions$comname<-mesh2name(test_train_prdictions$Row.names)
test_train_prdictions$is_in_training_set<-test_train_prdictions$Row.names %in% names(edc_VAM_scores)
all_test_length<-length(test_train_prdictions$`all_VAM_score[which(names(all_VAM_score) %in% all_nick)]`[test_train_prdictions$is_in_training_set=='FALSE'])
pred<-length(which(test_train_prdictions$`all_VAM_score[which(names(all_VAM_score) %in% all_nick)]`[test_train_prdictions$is_in_training_set=='FALSE']>.85))
pred/all_test_length
xlsx::write.xlsx2(test_train_prdictions,file = 'outputData/excel_files/harmonic_nicks_predicted.xlsx')

# 

rm(list=ls())
library(magrittr)
source('functions/annotation_functions.R')

load('outputData/VAM_scores_average_most_informative_layers_moa.RData')
nick<-xlsx::read.xlsx2(file='inputData/Nick_Amir_Testing_30_10_2019.xlsx',sheetIndex = 1)
colnames(nick)<-c('CAS.Number','Name','SMILES..Canonical.')
all_nick<-as.character(nick$CAS.Number) %>% strsplit('&') %>% unlist() %>% gsub(pattern = ' ',replacement = '') %>% unique()
all_nick<-cas2mesh(all_nick,return_back_NA = 'F')
all_nick[ all_nick %in% "3825-26-1"]<-'C023036'    #manually found in pubmed
all_nick[ all_nick %in% "307-55-1"]<-'C522391'
all_nick[ all_nick %in% "56558-16-8"]<-'C417207'
all_nick[ all_nick %in% "31508-00-6"]<-'C070055'
all_nick[ all_nick %in% "152969-11-4"]<-'C111120'
all_nick[ all_nick %in% "068515-48-0"]<-'C012125'
all_nick[ all_nick %in% '028553-12-0']<-'C012125'
all_nick[ all_nick %in% "84852-15-3"]<-'C025256'
all_nick[ all_nick %in% "38380-01-7"]<-'C500336'
all_nick[ all_nick %in% "111988-49-9"]<-'C417209'
all_nick[ all_nick %in% "843-55-0" ]<-'C570106'
all_nick[ all_nick %in% "13595-25-0"]<-'C109393'
all_nick[ all_nick %in% "56-35-9"]<-'C005961'
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

test_train_prdictions<-as.data.frame(all_VAM_score[which(names(all_VAM_score) %in% all_nick)])
test_train_prdictions<-merge(all_prob,test_train_prdictions,by=0,all=F)
test_train_prdictions$cas<-mesh2cas(test_train_prdictions$Row.names)
test_train_prdictions$comname<-mesh2name(test_train_prdictions$Row.names)
test_train_prdictions$is_in_training_set<-test_train_prdictions$Row.names %in% names(edc_VAM_scores)
all_test_length<-length(test_train_prdictions$`all_VAM_score[which(names(all_VAM_score) %in% all_nick)]`[test_train_prdictions$is_in_training_set=='FALSE'])
pred<-length(which(test_train_prdictions$`all_VAM_score[which(names(all_VAM_score) %in% all_nick)]`[test_train_prdictions$is_in_training_set=='FALSE']>.85))
pred/all_test_length
xlsx::write.xlsx2(test_train_prdictions,file = 'outputData/excel_files/average_nicks_predicted.xlsx')


# harmonic and average----------------------------------------------------------------

rm(list=ls())
library(magrittr)
source('functions/annotation_functions.R')
source('functions/vam_functions.R')
nick<-xlsx::read.xlsx2(file='inputData/Nick_Amir_Testing_30_10_2019.xlsx',sheetIndex = 1)
colnames(nick)<-c('CAS.Number','Name','SMILES..Canonical.')
all_nick<-as.character(nick$CAS.Number) %>% strsplit('&') %>% unlist() %>% gsub(pattern = ' ',replacement = '') %>% unique()
all_nick<-cas2mesh(all_nick,return_back_NA = 'F')
all_nick[ all_nick %in% "3825-26-1"]<-'C023036'    #manually found in pubmed
all_nick[ all_nick %in% "307-55-1"]<-'C522391'
all_nick[ all_nick %in% "56558-16-8"]<-'C417207'
all_nick[ all_nick %in% "31508-00-6"]<-'C070055'
all_nick[ all_nick %in% "152969-11-4"]<-'C111120'
all_nick[ all_nick %in% "068515-48-0"]<-'C012125'
all_nick[ all_nick %in% '028553-12-0']<-'C012125'
all_nick[ all_nick %in% "84852-15-3"]<-'C025256'
all_nick[ all_nick %in% "38380-01-7"]<-'C500336'
all_nick[ all_nick %in% "111988-49-9"]<-'C417209'
all_nick[ all_nick %in% "843-55-0" ]<-'C570106'
all_nick[ all_nick %in% "13595-25-0"]<-'C109393'
all_nick[ all_nick %in% "56-35-9"]<-'C005961'
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
xlsx::write.xlsx2(external_validation,file = 'outputData/excel_files/nicks_predicted_moa.xlsx')