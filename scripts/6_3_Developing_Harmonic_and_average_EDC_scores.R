
# 1_harmonic --------------------------------------------------------------


library(ggplot2)
library(ggrepel)
source('functions/plot_functions.R')
source('functions/annotation_functions.R')
source('functions/general_functions.R')

#load('outputData/class_probilities/class_probilities_all_compounds_final.RData')       # all pathways
load('outputData/class_probilities/class_probilities_all_compounds_final_moa.RData')    # moa pathways

#load('outputData/most_informative_layers_final.RData')                                # all pathways
load('outputData/most_informative_layers_final_moa.RData')                             # moa pathways

load('outputData/new199edc_1462dec.RData')
all_prob<-as.data.frame(all_compounds_prob[[1]][,'edc'])
all_prob$compname<-rownames(all_prob)
names(all_prob)<-c(names(all_compounds_prob)[[1]],'compname')

for (i in 1:length(all_compounds_prob)){
  print(names(all_compounds_prob)[i])
  p<-as.data.frame(all_compounds_prob[[i]][,'edc'])
  p$compname<-rownames(p)
  names(p)<-c(names(all_compounds_prob)[[i]],'compname')
  all_prob<-merge(all_prob,p,all = T)
}
colnames(all_prob)<-name_correct(colnames(all_prob))
rownames(all_prob)<-all_prob[,'compname']
all_prob<-all_prob[,which(!colnames(all_prob) %in% 'compname')]
ind<-which(colnames(all_prob) %in% names(most_informative_layers))  # index for the most informative layers


hs<-harmonic_sum(all_prob,ind)            # #calculation of harmonc sums as VAM scores for the most informative layers
hs<-scales:::rescale(hs,to = c(0, 1))


all_VAM_score<-hs       #as harmonic sum of most informative layers

# VAM for edcs
edc_VAM_scores<-edc_probs<-sort(all_VAM_score[which(names(all_VAM_score) %in% names(new_edc199_decoy_1462)[1:199])],decreasing = TRUE)
# VAM for unknown compounds edcs and decoys are excluded
unk_VAM_scores<-ctd_probs<-sort(all_VAM_score[which(!names(all_VAM_score) %in% names(new_edc199_decoy_1462))],decreasing = TRUE)
# saving
#save(all_prob,all_VAM_score,edc_VAM_scores,unk_VAM_scores,file='outputData/VAM_scores_hs_most_informative_layers.RData')
#save(all_prob,all_VAM_score,edc_VAM_scores,unk_VAM_scores,file='outputData/VAM_scores_average_most_informative_layers.RData')
save(all_prob,all_VAM_score,edc_VAM_scores,unk_VAM_scores,file='outputData/VAM_scores_hs_most_informative_layers_moa.RData')
par( mar= c(30, 5, 2, 2))
#barplot(edc_probs,las=2,cex.names = 0.4,ylab = 'TOXICITY SCORE',main = 'VAM toxicity scores for EDCs' )
barplot(ctd_probs,las=2,cex.names = 0.6,ylab = 'TOXICITY SCORE',main='VAM score for the CTD compounds' )
edc_plot<-as.data.frame(edc_VAM_scores)
rownames(edc_plot)<-mesh2cas(rownames(edc_plot))
plot(bar_plot(edc_plot,rownames(edc_plot),edc_plot$edc_VAM_scores,'VAM scores for EDCs ','Compounds','VAM_score',leg_title='Compounds'))


# 2_average ---------------------------------------------------------------

rm(list=ls())
library(ggplot2)
library(ggrepel)
source('functions/plot_functions.R')
source('functions/annotation_functions.R')
source('functions/general_functions.R')

#load('outputData/class_probilities/class_probilities_all_compounds_final.RData')       # all pathways
load('outputData/class_probilities/class_probilities_all_compounds_final_moa.RData')    # moa pathways

#load('outputData/most_informative_layers_final.RData')                                # all pathways
load('outputData/most_informative_layers_final_moa.RData')                             # moa pathways

load('outputData/new199edc_1462dec.RData')
all_prob<-as.data.frame(all_compounds_prob[[1]][,'edc'])
all_prob$compname<-rownames(all_prob)
names(all_prob)<-c(names(all_compounds_prob)[[1]],'compname')

for (i in 1:length(all_compounds_prob)){
  print(names(all_compounds_prob)[i])
  p<-as.data.frame(all_compounds_prob[[i]][,'edc'])
  p$compname<-rownames(p)
  names(p)<-c(names(all_compounds_prob)[[i]],'compname')
  all_prob<-merge(all_prob,p,all = T)
}
colnames(all_prob)<-name_correct(colnames(all_prob))
rownames(all_prob)<-all_prob[,'compname']
all_prob<-all_prob[,which(!colnames(all_prob) %in% 'compname')]
ind<-which(colnames(all_prob) %in% names(most_informative_layers))  # index for the most informative layers

means<-mean_function(all_prob,ind)

all_VAM_score<-means   #as mean         of most informative layers

# VAM for edcs
edc_VAM_scores<-edc_probs<-sort(all_VAM_score[which(names(all_VAM_score) %in% names(new_edc199_decoy_1462)[1:199])],decreasing = TRUE)

# VAM for unknown compounds edcs and decoys are excluded
unk_VAM_scores<-ctd_probs<-sort(all_VAM_score[which(!names(all_VAM_score) %in% names(new_edc199_decoy_1462))],decreasing = TRUE)
save(all_prob,all_VAM_score,edc_VAM_scores,unk_VAM_scores,file='outputData/VAM_scores_average_most_informative_layers_moa.RData')


par( mar= c(30, 5, 2, 2))
#barplot(edc_probs,las=2,cex.names = 0.4,ylab = 'TOXICITY SCORE',main = 'VAM toxicity scores for EDCs' )
barplot(ctd_probs,las=2,cex.names = 0.6,ylab = 'TOXICITY SCORE',main='VAM score for the CTD compounds' )
edc_plot<-as.data.frame(edc_VAM_scores)
rownames(edc_plot)<-mesh2cas(rownames(edc_plot))
plot(bar_plot(edc_plot,rownames(edc_plot),edc_plot$edc_VAM_scores,'VAM scores for EDCs ','Compounds','VAM_score',leg_title='Compounds'))





# excel files of EDC scores for all compounds in CTD-------------------------------------------------------------
#harmonic and average

source('functions/annotation_functions.R')
source('functions/vam_functions.R')  #harmonic_vam_score and average_vam_score functions
load('outputData/chem2gene_no_out.RData')
all_comp<-names(chem2gene)
all_vam<-harmonic_vam_score(all_comp)
all_vam<-all_vam[!is.na(all_vam$VAM_score),]
all_vam$comp_names<-mesh2name(rownames(all_vam))
all_vam$mesh<-rownames(all_vam)
all_vam$is_in_training<-rep('unk',nrow(all_vam))
load('outputData/new199edc_1462dec.RData')
all_vam$is_in_training[all_vam$mesh %in% names(new_edc199_decoy_1462)[1:199]]='edc'
all_vam$is_in_training[all_vam$mesh %in% names(new_edc199_decoy_1462)[200:1661]]='decoy'
colnames(all_vam)[which(colnames(all_vam)=="VAM_score")]='harmonic_sum_edc_score'
all_vam$average_edc_score<-average_vam_score(all_vam$mesh)$VAM_score
all_vam$cas<-mesh2cas(all_vam$mesh)
colnames(all_vam)[colnames(all_vam)=="Consensus_LINCS_HEPG2_HEPG2"]="Consensus_LINCS_HEPG2"
colnames(all_vam)[colnames(all_vam)=="Consensus_Rat_invitro_Drug Matrix__TG_GATEs"]="Consensus_Rat_invitro_Drug Matrix_TG_GATEs"
save(all_vam,file='outputData/CTD_edc_scores_dictionary.RData')
xlsx::write.xlsx2(all_vam,file='outputData/excel_files/all_compounds_vam_moa.xlsx')

