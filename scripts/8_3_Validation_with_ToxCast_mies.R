
# 1. Retreiving ToxCast MIES for all compounds in TOXCAST ------------------------------------------------------------

# MIES for all  compounds in ToxCast based on toxcast assays: if a toxcast assay is active in hitcall for a compound its target gene is considered as
# MIE  for that compound

load('outputData/ToxCast/hitcall.RData')           # ToxCast HITCALL matrix
load('outputData/ToxCast/assay_target_gene.RData') # Toxcast target gene ids
tst_comp<-rownames(toxcast_3_1)
mies<-list()
for (i in 1:length(tst_comp)){
  print(tst_comp[i])
  assays<-colnames(toxcast_3_1)[which(toxcast_3_1[rownames(toxcast_3_1) %in% tst_comp[i],]==1)]
  if(length(assays)>0){
    mie<-c()   
    mie<-unique(unlist(lapply(assays, function(x)toxcast$intended_target_entrez_gene_id[toxcast$assay_component_endpoint_name==x])))
    
  }else{
    mie<-NA
  }
  mies[[i]]<-mie
}
names(mies)<-tst_comp

mies[which(lapply(mies, length) %in% 0)]<-c()
mies[which(lapply(mies, is.na) %in% TRUE)]<-c()
cas_mies<-mies
mesh_mies<-mies
source('functions/annotation_functions.R')
names(mesh_mies)<-cas2mesh(names(mesh_mies))
save(mesh_mies,cas_mies,file='outputData/all_toxcast_mies.RData')



# 2. calculation of EDC score for for all componuds in the intersection of CTD and TOXCAST ----------------------------------------------------------------------

load('outputData/all_toxcast_mies.RData')   #MIES from TOXcast
load('outputData/chem2gene_no_out.RData')
int_names<-intersect(names(chem2gene),names(mesh_mies))
mesh_mies_tox_ctd<-mesh_mies[int_names]
source('functions/vam_functions.R')                 # mie2vam function
predicted_tox_based_vam<-mie2vam(mesh_mies_tox_ctd) # all pipeline components
predicted_from_toxcast<-predicted_tox_based_vam[,c('harmonic_VAM_score','average_VAM_score')]
colnames(predicted_from_toxcast)<-c('harmonic_from_toxcast','average_from_toxcast')
predicted_from_toxcast$average_from_CTD<-average_vam_score(rownames(predicted_from_toxcast))$VAM_score
integrated_tox_with_ctd<-na.omit(predicted_from_toxcast[,2:3])
load('outputData/new199edc_1462dec.RData')
integrated_tox_with_ctd$status<-'unk'
integrated_tox_with_ctd$status[which(rownames(integrated_tox_with_ctd) %in% names(new_edc199_decoy_1462)[1:199])]<-'edc'
integrated_tox_with_ctd$status[which(rownames(integrated_tox_with_ctd) %in% names(new_edc199_decoy_1462)[200:1661])]<-'decoy'
source('functions/annotation_functions.R')
integrated_tox_with_ctd$comp_names<-mesh2name(rownames(integrated_tox_with_ctd))
integrated_tox_with_ctd$cas<-mesh2cas(rownames(integrated_tox_with_ctd))

save(predicted_tox_based_vam,integrated_tox_with_ctd,file='outputData/predicted_toxcast_mie2vam.RData')
write.csv(integrated_tox_with_ctd,file = 'outputData/excel_files/predicted_toxcast_mie2vam.csv')
ctd_driven<-xlsx::read.xlsx('outputData/excel_files/all_compounds_average_vam_moa.xlsx',sheetIndex = 1) # predicted EDC scores from CTD MIES
save(ctd_driven,file='outputData/ctd_driven_probability_scores_all_layers.RData')



# 3. plot Visualization of the results of ToxCast MIEs for EDCs and Decoys --------------------------------------------------------
rm(list=ls());gc()
library(ggplot2)
library(ggrepel)
load('outputData/predicted_toxcast_mie2vam.RData') # predicted scores from ToxCast MIEs

Toxcast_driven<-predicted_tox_based_vam[,c("Drug_Matrix_Rat_invivo_Single_Dose_1_day",
                                           "Drug_Matrix_Rat_invivo_Repeated_Dose_5_days"
                                            )] #selection of two in vivo models from ToxCast predictions
colnames(Toxcast_driven)<-paste('using_tox_cast_mies',colnames(Toxcast_driven),sep = '')
Toxcast_driven$mesh<-rownames(Toxcast_driven)
load('outputData/ctd_driven_probability_scores_all_layers.RData') #prediceted from CTD MIEs
ctd_driven<-ctd_driven[,c('Drug_Matrix_Rat_invivo_Single_Dose_1_day',
                          'Drug_Matrix_Rat_invivo_Repeated_Dose_5_days',
                          'mesh','comp_names',
                          'is_in_training')]
merged_tox_ctd_prdictions<-merge(ctd_driven,Toxcast_driven,all=F)
merged_tox_ctd<-merged_tox_ctd_prdictions[which(merged_tox_ctd_prdictions$is_in_training=='edc' | merged_tox_ctd_prdictions$is_in_training=='decoy'),]
write.csv(merged_tox_ctd,file = 'outputData/excel_files/from_toxcast_mies.csv')
# for DM1 and DM5 invivo
merged_tox_ctd_prdictions_dm1<-merged_tox_ctd[,c('comp_names','Drug_Matrix_Rat_invivo_Single_Dose_1_day',
                                                 'using_tox_cast_miesDrug_Matrix_Rat_invivo_Single_Dose_1_day',
                                                 'is_in_training')]
merged_tox_ctd_prdictions_dm5<-merged_tox_ctd[,c('comp_names',"Drug_Matrix_Rat_invivo_Repeated_Dose_5_days",
                                                 "using_tox_cast_miesDrug_Matrix_Rat_invivo_Repeated_Dose_5_days",
                                                 'is_in_training')]
colnames(merged_tox_ctd_prdictions_dm1)[2:3]<-colnames(merged_tox_ctd_prdictions_dm5)[2:3]<-c('CTD_based',
                                                                                              'ToxCast_based')

merged_tox_ctd_prdictions_dm1$data_layer<-'Drug_Matrix_Rat_invivo_Single_Dose_1_day'
merged_tox_ctd_prdictions_dm5$data_layer<-'Drug_Matrix_Rat_invivo_Repeated_Dose_5_days'
all_data<-rbind(merged_tox_ctd_prdictions_dm1,merged_tox_ctd_prdictions_dm5)
all_data_invivo<-all_data
all_data_invivo$in_vit<-'In Vivo'


load('outputData/predicted_toxcast_mie2vam.RData') # predicted scores from ToxCast MIEs
colnames(predicted_tox_based_vam)<-gsub(x=colnames(predicted_tox_based_vam),pattern = '__',replacement = '_')
colnames(predicted_tox_based_vam)<-gsub(x=colnames(predicted_tox_based_vam),pattern = '_HEPG2_HEPG2',replacement = '_HEPG2')
Toxcast_driven<-predicted_tox_based_vam[,c("Drug_Matrix_Rat_invitro_Single_Dose_1_day",
                                           "Consensus_Rat_invitro_Drug Matrix_TG_GATEs",
                                           "Consensus_LINCS_HEPG2",
                                           "PPI_STRINGdb"
)] #selection of  in vitro models from ToxCast predictions
layer_names<-colnames(Toxcast_driven)
colnames(Toxcast_driven)<-paste('ToxCast',colnames(Toxcast_driven),sep = '_')
Toxcast_driven$mesh<-rownames(Toxcast_driven)
# CTD
load('outputData/ctd_driven_probability_scores_all_layers.RData') #prediceted from CTD MIEs
colnames(ctd_driven)<-gsub(x=colnames(ctd_driven),pattern = 'Drug.Matrix__TG_GATEs',replacement = 'Drug Matrix_TG_GATEs')
colnames(ctd_driven)<-gsub(x=colnames(ctd_driven),pattern = '_HEPG2_HEPG2',replacement = '_HEPG2')
ctd_driven<-ctd_driven[,c("Drug_Matrix_Rat_invitro_Single_Dose_1_day",
                          "Consensus_Rat_invitro_Drug Matrix_TG_GATEs",
                          "Consensus_LINCS_HEPG2",
                          "PPI_STRINGdb",
                          'mesh','comp_names',
                          'is_in_training')]

merged_tox_ctd_prdictions<-merge(ctd_driven,Toxcast_driven,all=F)
merged_tox_ctd<-merged_tox_ctd_prdictions[which(merged_tox_ctd_prdictions$is_in_training=='edc' | merged_tox_ctd_prdictions$is_in_training=='decoy'),]

write.csv(merged_tox_ctd,file = 'outputData/excel_files/from_toxcast_mies_invitro.csv')

# for all in vitro layers
mixed<-list()
for (i in layer_names){
  
  tox_ind<-which(colnames(merged_tox_ctd) %in% paste('ToxCast',i,sep='_'))
  ctd_ind<-which(colnames(merged_tox_ctd) %in% i)
  inds<-c(6,ctd_ind,tox_ind,7)
  dat<-merged_tox_ctd[,inds]
  colnames(dat)<-c('comp_names','CTD_based','ToxCast_based','is_in_training')
  dat$data_layer<-i
  mixed[[i]]<-dat
}

all_data_invitro<-do.call(rbind,mixed)
all_data_invitro$in_vit<-'In Vitro'
all_data<-rbind(all_data_invitro,all_data_invivo)
all_data$is_in_training<-gsub(all_data$is_in_training,pattern = 'decoy',replacement = 'Negative Control')
all_data$is_in_training<-gsub(all_data$is_in_training,pattern = 'edc',replacement = 'EDCs')
all_data$comp_names[all_data$is_in_training=='EDCs' & all_data$ToxCast_based>0.5]<-''
all_data$comp_names[all_data$is_in_training=='Negative Control' & all_data$ToxCast_based<0.5]<-''


p<-ggplot(all_data,aes(x=CTD_based ,y=ToxCast_based,
                       shape=factor(all_data$is_in_training),
                       colour=factor(all_data$is_in_training)))+
  geom_point(size=ifelse(is.na(all_data$comp_names), 5,2))+
  ylab('Class probability (MIES from ToxCast)')+ xlab('Class probability (MIES from CTD)')+
  scale_shape_manual(name='Compounds',values = c(19,20),labels=c('EDCs','Negative Control'))+
  scale_color_manual(name='Compounds',values = c('red','grey'),labels=c('EDCs','Negative Control'))+
  theme_minimal() +   
  geom_text_repel(show.legend = F,label=all_data$comp_names,size=6,segment.size = 0.2, segment.color = "grey")+
  #ggtitle('In vivo Drug Matrix ')+
  facet_wrap( in_vit ~  data_layer,nrow = 4 )+
  theme(axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 16))
plot(p)

Mie_2_class_prob<-all_data
save(Mie_2_class_prob,file='outputData/plots_tables_functions/MIEs_TOXCast_to_Class_prob/toxcast_mies_class_prob.RData')











# 4. Accuracy of the 6 models based on F1-scores----------------------------------------------------------------
rm(list = ls())
#ToxCast
load('outputData/predicted_toxcast_mie2vam.RData') # predicted scores from ToxCast MIEs
colnames(predicted_tox_based_vam)<-gsub(x=colnames(predicted_tox_based_vam),pattern = '__',replacement = '_')
colnames(predicted_tox_based_vam)<-gsub(x=colnames(predicted_tox_based_vam),pattern = '_HEPG2_HEPG2',replacement = '_HEPG2')
Toxcast_driven<-predicted_tox_based_vam[,1:6]
colnames(Toxcast_driven)<-paste('ToxCast',colnames(Toxcast_driven),sep = '_')
Toxcast_driven$mesh<-rownames(Toxcast_driven)
# CTD
load('outputData/ctd_driven_probability_scores_all_layers.RData') #prediceted from CTD MIEs
colnames(ctd_driven)<-gsub(x=colnames(ctd_driven),pattern = 'Drug.Matrix__TG_GATEs',replacement = 'Drug Matrix_TG_GATEs')
colnames(ctd_driven)<-gsub(x=colnames(ctd_driven),pattern = '_HEPG2_HEPG2',replacement = '_HEPG2')
ctd_driven<-ctd_driven[,c(2,3,4,9,10,11,19,18,20)]
merged_tox_ctd_prdictions<-merge(ctd_driven,Toxcast_driven,all=F)                      
merged_tox_ctd<-merged_tox_ctd_prdictions[which(merged_tox_ctd_prdictions$is_in_training=='edc' | merged_tox_ctd_prdictions$is_in_training=='decoy'),]
# calculation of the F1 accuracy scores from confusion matrix
source('functions/glm_functions.R')
acuracy_mat<-merged_tox_ctd
acuracy_real<-as.character(acuracy_mat[,"is_in_training"])
acuracy_mat[,1:9]<-c()
evaluation<-cm<-list()
for (i in 1:ncol(acuracy_mat)){
  test_mat<-cbind(acuracy_mat[,i],acuracy_real)
  test_mat<-na.omit(test_mat)
  test_mat[test_mat[,1]>0.5]<-'edc'
  test_mat[,1][test_mat[,1]<=0.5]<-'decoy'
  cm[[i]]<-confusion_matrix(test_mat[,2],test_mat[,1])
  evaluation[[i]]<-cm_based_evaluation(cm[[i]])
  }
names(evaluation)<-names(cm)<-colnames(acuracy_mat)
knitr::kable(sapply(evaluation, function(x)mean(x$f1)),col.names = 'F1-score')
tbl<-as.data.frame(sapply(evaluation, function(x)mean(x$f1)))
colnames(tbl)<-'F1-score-From-ToxCast-MIEs'
tbl$data_layer<-rownames(tbl)
write.csv(tbl,file = 'outputData/excel_files/F1_From_ToxCast_MIES.csv')


# tbl$`F1-score`<-round(tbl$`F1-score`,digits = 2)
# tbl$data_layers<-c('Drug_Matrix_Rat_invivo_Single_Dose_1_day','Drug_Matrix_Rat_invivo_Repeated_Dose_5_days')
# tbl<-tbl[,c(2,1)]
# rownames(tbl)<-c(' ','  ')
# tbl<-ggpubr::ggtexttable(tbl)
# img<-imager::load.image('images/figure4.png');img<-grid::rasterGrob(img)
# #ggarrange(img,p,tbl,ncol = 1,nrow = 3,font.label = list(size=25,color='black'),heights = c(1,2,1))
# library(gridExtra)
# grid.arrange(grobs=list(img,p,tbl),layout_matrix=rbind(c(1,1,NA,NA,NA,NA),c(2,2,2,2,2,2,2,2),c(2,2,2,2,2,2,2,2),3))
# # figure 4
# tiff('plots/figure4.tiff',units = 'in',width = 18,height = 12,res = 300)
# grid.arrange(grobs=list(img,p,tbl),layout_matrix=rbind(c(1,1,NA,NA,NA,NA),c(2,2,2,2,2,2,2,2),c(2,2,2,2,2,2,2,2),c(3,3,NA,NA,NA,NA)))
# dev.off()
# 


