
# 1. Plot of harmonic EDC score VS toxpi for DeDuCT compounds ----------------------------------------------------------------

rm(list=ls())
library(ggplot2)
library(ggrepel)
source('functions/annotation_functions.R')
load('outputData/VAM_scores_hs_most_informative_layers_moa.RData')
names(edc_VAM_scores)<-mesh2cas(names(edc_VAM_scores))  # converting mesh to cas numbers
toxpi_edc<-read.csv('inputData/toxpi_scores.csv')
toxpi_edc$CASRN<-gsub(toxpi_edc$CASRN,pattern = '\\.',replacement = '')
rownames(toxpi_edc)<-toxpi_edc$CASRN
vam_score<-edc_VAM_scores[intersect(names(edc_VAM_scores),toxpi_edc$CASRN)]
toxpi_score<-toxpi_edc[intersect(names(edc_VAM_scores),toxpi_edc$CASRN),"score"]
names(toxpi_score)<-intersect(names(edc_VAM_scores),toxpi_edc$CASRN)
df_edc<-as.data.frame(cbind(vam_score,toxpi_score))
df_edc$compound_cas<-names(toxpi_score)
rownames(df_edc)<-cas2name(rownames(df_edc))
deduct_smiles<-read.csv('inputData/DEDuCT_ChemicalBasicInformation.csv')
deduct_smiles<-deduct_smiles[which(deduct_smiles$CAS.Number %in% df_edc$compound_cas),c(2,6)]
colnames(deduct_smiles)<-c('compound_cas','SMILES')
all_data_h<-merge(deduct_smiles,df_edc)
all_data_h$comp_name<-cas2name(as.character(all_data_h$compound_cas))
ind<-which((all_data_h$toxpi_score >0.15 & all_data_h$vam_score>0.75)|(all_data_h$toxpi_score<0.055 & all_data_h$vam_score>0.95))
all_data_h$comp_name[-ind]<-NA
all_data_h$colr<-ifelse(all_data_h$toxpi_score >0.15 & all_data_h$vam_score>0.75, 'r', 'g')
all_data_h$type<-'Harmonic Sum Score'
load('outputData/VAM_scores_average_most_informative_layers_moa.RData')
names(edc_VAM_scores)<-mesh2cas(names(edc_VAM_scores))  # converting mesh to cas numbers
toxpi_edc<-read.csv('inputData/toxpi_scores.csv')
toxpi_edc$CASRN<-gsub(toxpi_edc$CASRN,pattern = '\\.',replacement = '')
rownames(toxpi_edc)<-toxpi_edc$CASRN
vam_score<-edc_VAM_scores[intersect(names(edc_VAM_scores),toxpi_edc$CASRN)]
toxpi_score<-toxpi_edc[intersect(names(edc_VAM_scores),toxpi_edc$CASRN),"score"]
names(toxpi_score)<-intersect(names(edc_VAM_scores),toxpi_edc$CASRN)
df_edc<-as.data.frame(cbind(vam_score,toxpi_score))
df_edc$compound_cas<-names(toxpi_score)
rownames(df_edc)<-cas2name(rownames(df_edc))
deduct_smiles<-read.csv('inputData/DEDuCT_ChemicalBasicInformation.csv')
deduct_smiles<-deduct_smiles[which(deduct_smiles$CAS.Number %in% df_edc$compound_cas),c(2,6)]
colnames(deduct_smiles)<-c('compound_cas','SMILES')
all_data<-merge(deduct_smiles,df_edc)
all_data$comp_name<-cas2name(as.character(all_data$compound_cas))
txpi_edc_score<-all_data
ind<-which((txpi_edc_score$toxpi_score >0.15 & txpi_edc_score$vam_score>0.75)|(txpi_edc_score$toxpi_score<0.06 & txpi_edc_score$vam_score>0.95))
txpi_edc_score$comp_name[-ind]<-NA
txpi_edc_score$colr<-ifelse(all_data_h$toxpi_score >0.15 & all_data_h$vam_score>0.75, 'r', 'g')
txpi_edc_score$type<-'Average Score'
all_data_h<-rbind(all_data_h,txpi_edc_score)
# plot
p<-ggplot(all_data_h,aes(x=vam_score,y=toxpi_score,colour=factor(all_data_h$colr,levels = c('r','g'))))+
  geom_point(size=ifelse(all_data_h$toxpi_score >0.15 & all_data_h$vam_score>0.75, 4, 2),
             show.legend = T)+
  geom_text_repel(label=all_data_h$comp_name,size=5,segment.size = 0.2, segment.color = "grey",show.legend = F) +
  scale_color_manual(name='DEDuCT EDCs',values = c('red','grey'),labels=c('High ToxPi score','Low  ToxPi score'))+
  ggtitle('')+ xlab('EDC scores')+ylab('ToxPi scores')+
  facet_grid( . ~  type , scale="free", space = 'free')+
  theme_minimal() +
  theme(axis.text = element_text(size = 20),axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
        strip.text = element_text(size = 20),legend.text = element_text(size = 16))
plot(p)
toxpi_DeDuct<-all_data_h
save(toxpi_DeDuct,file='outputData/plots_tables_functions/ToxPi/Toxpi_reesults.RData')


# 2. Plot of harmonic sum edc score VS TOXPI for unknown compounds -----------------------------------------------------------
rm(list=ls())
source('functions/annotation_functions.R')
toxpi<-read.csv('inputData/toxpi_scores.csv')
toxpi$CASRN<-gsub(toxpi$CASRN,pattern = '\\.',replacement = '')
load('outputData/VAM_scores_hs_most_informative_layers_moa.RData') # VAM scores harmonic sum of the most informative layers
tox<-toxpi[,c(1,4)]
colnames(tox)<-c('cas','toxpi')
vam<-as.data.frame(unk_VAM_scores)   # only compounds with unknown labels not edc and decoy

#vam<-as.data.frame(all_VAM_score)     # all compounds 


colnames(vam)<-'unk_VAM_scores'




vam$cas<-mesh2cas(rownames(vam))
vam$mesh<-rownames(vam)
all_mat<-merge(tox,vam,all = F)
all_mat$compound_name<-cas2name(all_mat$cas)
# plot
library(ggplot2)
library(ggrepel)
#all_mat$compound_name[which(all_mat$toxpi<0.1)]<-''
ggplot(all_mat,aes(x=unk_VAM_scores,y=toxpi))+
  geom_point(size=1,color=ifelse(all_mat$toxpi >0.1 & all_mat$unk_VAM_scores>0.85, "red", "grey50"))+
  geom_text_repel(label=all_mat$compound_name,size=3,segment.size = 0.2, segment.color = "grey")+
  ylab('TOXPI score')+ xlab('Harmonic EDC score')+
  theme_classic()+  ggtitle('Harmonic EDC score VS TOXPI scores For Unknown compounds ')+xlab('Harmonic EDC score')+ylab('TOXPI score')

# 3. Plot of Average VS TOXPI for unknown compounds --------------------------------------------------------------
rm(list=ls())
source('functions/annotation_functions.R')
toxpi<-read.csv('inputData/toxpi_scores.csv')
toxpi$CASRN<-gsub(toxpi$CASRN,pattern = '\\.',replacement = '')
load('outputData/VAM_scores_average_most_informative_layers_moa.RData')  #average of the most informative layeres
tox<-toxpi[,c(1,4)]
colnames(tox)<-c('cas','toxpi')
vam<-as.data.frame(unk_VAM_scores)
colnames(vam)<-'unk_VAM_scores'
vam$cas<-mesh2cas(rownames(vam))
vam$mesh<-rownames(vam)
all_mat<-merge(tox,vam,all = F)
all_mat$compound_name<-cas2name(all_mat$cas)
# plot
library(ggplot2)
library(ggrepel)
all_mat$compound_name[which(all_mat$toxpi<0.1)]<-''
ggplot(all_mat,aes(x=unk_VAM_scores,y=toxpi))+
  geom_point(size=1,color=ifelse(all_mat$toxpi >0.1 & all_mat$unk_VAM_scores>0.85, "red", "grey50"))+
  geom_text_repel(label=all_mat$compound_name,size=3,segment.size = 0.2, segment.color = "grey")+ylab('TOXPI score')+ xlab('Average EDC score')+
  theme_classic() +  ggtitle('Average EDC score VS TOXPI scores For Unknown compounds ')+ xlab('Average EDC score')+ylab('TOXPI score')


# 4. excel file of harmonic sum edc score VS TOXPI for all compounds -----------------------------------------------------------
rm(list=ls())
source('functions/annotation_functions.R')
toxpi<-read.csv('inputData/toxpi_scores.csv')
toxpi$CASRN<-gsub(toxpi$CASRN,pattern = '\\.',replacement = '')
load('outputData/VAM_scores_hs_most_informative_layers_moa.RData') # VAM scores harmonic sum of the most informative layers
tox<-toxpi[,c(1,4)]
colnames(tox)<-c('cas','toxpi')
vam<-as.data.frame(all_VAM_score)   # only compounds with unknown labels not edc and decoy
colnames(vam)<-'harmonic_sum_VAM_scores'
vam$cas<-mesh2cas(as.character(rownames(vam)))
vam$mesh<-rownames(vam)
load('outputData/new199edc_1462dec.RData')
vam$status<-'unk'
vam$status[which(rownames(vam) %in% names(new_edc199_decoy_1462)[1:199])]<-'edc'
vam$status[which(rownames(vam) %in% names(new_edc199_decoy_1462)[200:1661])]<-'decoy'
all_mat<-merge(tox,vam,all = F)
all_mat$compound_name<-cas2name(as.character(all_mat$cas))
write.csv(all_mat,file='outputData/excel_files/toxpi_harmonic_vam_moa.csv')



# 5. Excel file of Average VS TOXPI for all compounds --------------------------------------------------------------
rm(list=ls())
source('functions/annotation_functions.R')
toxpi<-read.csv('inputData/toxpi_scores.csv')
toxpi$CASRN<-gsub(toxpi$CASRN,pattern = '\\.',replacement = '')
load('outputData/VAM_scores_average_most_informative_layers_moa.RData')  #average of the most informative layeres
tox<-toxpi[,c(1,4)]
colnames(tox)<-c('cas','toxpi')
vam<-as.data.frame(all_VAM_score)
colnames(vam)<-'average_VAM_scores'
vam$cas<-mesh2cas(as.character(rownames(vam)))
vam$mesh<-rownames(vam)
load('outputData/new199edc_1462dec.RData')
vam$status<-'unk'
vam$status[which(rownames(vam) %in% names(new_edc199_decoy_1462)[1:199])]<-'edc'
vam$status[which(rownames(vam) %in% names(new_edc199_decoy_1462)[200:1661])]<-'decoy'
all_mat<-merge(tox,vam,all = F)
all_mat$compound_name<-cas2name(as.character(all_mat$cas))
write.csv(all_mat,file='outputData/excel_files/toxpi_average_vam_moa.csv')



