
# 1. Integration of stabilties from K-FOLD-CV with GLM coefficients -----------------------------------------------------------


#setwd("outputData/glm/stratified_networks_CV/")     # CV results for all models all pathways (the first run with 7000pathways)
setwd("outputData/glm/stratified_moa_pathways/")     # CV results for all models with moa pathways (around 2000)

li<-list.files('.',pattern = '.RData')
load(li[[1]])
cof<-lapply(modls, function(x)x$coefs)
for (i in 2:length(li)){
  load(li[[i]])
  for (j in 1:length(cof)){
    
    cof[[j]]<-c(cof[[j]],modls[[j]]$coefs)
  }
}
setwd('../../..')
source('functions/general_functions.R')
#load('outputData/glm/glm_models_final.RData')         # final glm models for all layers all pathways
load('outputData/glm/glm_models_final_moa.RData')      # final glm models for all layers only pathways related to MOAs

#names(modls)<-name_correct(names(modls))
#names(cof)<-name_correct(names(cof))
if (all(names(modls)==names(cof))){
stability_cof<-list()                                  # making a list of stabilites and coefiicients
  for (i in 1:length(cof)){
   coffs<-as.data.frame(cof[[i]])                      
   stability<-apply((abs(coffs)>0)*1, 1, sum)          # sum of the coefficients more than zero for all CV folds
   standard_dev<-apply(coffs, 1, sd)                   # standard deviation of CV results
   main_coefs<-modls[[i]]$cof[-1]                      # coefficients  from the final model the first one is intercept and is removed
   stability_cof[[i]]<-as.data.frame(list(main_coefs,stability,standard_dev),col.names = c('coef','stability','SD'))
   stability_cof[[i]]$network<-names(cof)[[i]]
   stability_cof[[i]]$pathways<-rownames(stability_cof[[i]])
  } #end for
names(stability_cof)<-names(cof)
data<-do.call(rbind,stability_cof)
} # end if


#save(data,file='outputData/patway_glm_coefficients/cv_based_glm_coefficients_stabilities_two_class.RData')
#write.csv2(data,file='outputData/excel_files/glm_coefficients_stabilitis.csv') # glm coefiicients and their stabilities are restored
save(data,file='outputData/patway_glm_coefficients/cv_based_glm_coefficients_stabilities_two_class_moa.RData')
write.csv(data,file='outputData/excel_files/glm_coefficients_stabilitis_moa.csv') # glm coefiicients and their stabilities are restored




# 2. Integration of NES scores and certatinty levels for pathways with GLM coefficients-----------------------------


source('functions/annotation_functions.R')
#load('outputData/patway_glm_coefficients/cv_based_glm_coefficients_stabilities_two_class.RData')  #all 7k pathways approach
#load('outputData/class_probilities/integrated_fgsea_allcompounds_results_final.RData')
#names(fgs)<-name_correct(names(fgs))
load('outputData/most_informative_layers_final_moa.RData') #most informative pathways
load('outputData/patway_glm_coefficients/cv_based_glm_coefficients_stabilities_two_class_moa.RData') #moa based pathways
load('outputData/glm/integrated_fgsea_results_edc_decoy_moa.RData')

avg_nes<-list()
for (i in 1:length(fgs)){               # for each network in fgs 
  edc_ind<-which(fgs[[i]]$y=='edc')     # ONLY AVERAGE of NES scores AND their SD FOR EDCS
  x<-fgs[[i]]$x[edc_ind,]               # ONLY AVERAGE of NES scores AND their SD FOR EDC
  avg_nes_edc<-apply(x,2,mean)          # average 
  sd_nes_edc<-apply(x,2,sd)             # standard deviation 
  
  dec_ind<-which(fgs[[i]]$y=='decoy')   # ONLY AVERAGE of NES scores AND their SD FOR decoys
  xx<-fgs[[i]]$x[dec_ind,]              # ONLY AVERAGE of NES scores AND their SD FOR decoys
  avg_nes_decoy<-apply(xx,2,mean)       # average 
  sd_nes_decoy<-apply(xx,2,sd)          # standard deviation 
  
  dt<-as.data.frame(list(avg_nes_edc,sd_nes_edc,avg_nes_decoy,sd_nes_decoy))
  names(dt)<-c('Average_NES_EDC','SD_NES_EDC','Average_NES_decoy','SD_NES_decoy')
  dt$network<-names(fgs)[[i]]
  dt$pathways<-rownames(dt)
  dt$n_edc<-length(edc_ind)
  dt$n_decoy<-length(dec_ind)
  avg_nes[[i]]<-dt
}
names(avg_nes)<-names(fgs)
all_nes_avg_sd<-do.call(rbind,avg_nes)
all_data<-merge(all_nes_avg_sd,data,by=c('network','pathways'))
colnames(all_data)<-c('network','pathways','Average_NES_EDC','SD_NES_EDC','Average_NES_decoy','SD_NES_decoy','n_edc','n_decoy',
'glm_coefs','CV_stability','CV_coefficients_SD')   
source('functions/annotation_functions.R')
all_data$annotated_pathways<-pathway_annotate(all_data$pathways)
#load("outputData/glm/glm_models_final.RData")
load("outputData/glm/glm_models_final_moa.RData")
list_pathways  <- lapply(modls, function(x){y=x$sel
y=y[-1] #the first element is intercept so we remove it for upset plot
y})
all_pathways <- unique(unlist(list_pathways))
int_data <- lapply(list_pathways, function(x){
  idx <- which(all_pathways %in% x)
  res<-rep(0, length(all_pathways))
  res[idx]<-1
  res})
upset_df <- as.data.frame(int_data , as.is=T, stringsAsFactors=F, check.names=F)
rownames(upset_df) <- all_pathways
number_of_networks_for_pathway<-apply(upset_df,1 ,sum)
number_of_most_informative_layers_for_pathway<-apply(upset_df[,most_informative_layers],1 ,sum)
dd<-as.data.frame(list(number_of_networks_for_pathway,number_of_most_informative_layers_for_pathway))
colnames(dd)<-c('number_of_networks_for_pathway','number_of_most_informative_layers_for_pathway')
dd$pathways<-rownames(dd)
final_data<-merge(all_data,dd,all=T)
#final_data<-final_data[!final_data$glm_coefs==0,]                                              #the pathways with GLM coefficients of zero are removed
mean_stability_each<-by(final_data$CV_stability,list(final_data$pathways),mean)                # avaerage of stability for each pathway across all networks
mean_stability_each_pathway<-as.data.frame(as.numeric(mean_stability_each),row.names=names(mean_stability_each))
colnames(mean_stability_each_pathway)<-'average_pathway_stability'
mean_stability_each_pathway$pathways<-rownames(mean_stability_each_pathway)
final_data<-merge(final_data,mean_stability_each_pathway)
final_data$pathway_certainty_percentage<-
(final_data$number_of_networks_for_pathway/length(unique(final_data$network)))*final_data$average_pathway_stability

# saving
save(final_data,all_data,file='outputData/patway_glm_coefficients/NESavg_sd_EDCs_stability_coefficients_pathways_moa.RData') #new run 2k pathways related to moa
write.csv(final_data,file='outputData/excel_files/NESavg_sd_EDCs_stability_coefficients_pathways_moa.csv')


# 3. t_test normalization test---------------------------------------------------------------
# shapiro test
# load('outputData/glm/integrated_fgsea_results_edc_decoy_moa.RData')
# shapiro_test_all<-lapply(fgs, function(y)sapply(apply(y$x,2,shapiro.test),function(x)x$p.value))
# #shapiro_test_edc<-lapply(fgs, function(y)sapply(apply(y$x[y$y=='edc',],2,shapiro.test),function(x)x$p.value))
# shapiro_test_decoy<-lapply(fgs, function(y)sapply(apply(y$x[y$y=='decoy',],2,shapiro.test),function(x)x$p.value))
# # since there was no normal distribution we decidede to integrate the ROC-AUC to the table of coefficients
# #jarque bra test
# load('outputData/glm/integrated_fgsea_results_edc_decoy_moa.RData')
# library(tseries)
# jarque_bra_test_all<-lapply(fgs, function(y)sapply(apply(y$x,2,jarque.bera.test),function(x)x$p.value))
# jarque_bra_test_edc<-lapply(fgs, function(y)sapply(apply(y$x[y$y=='edc',],2,jarque.bera.test),function(x)x$p.value))
# jarque_bra_test_decoy<-lapply(fgs, function(y)sapply(apply(y$x[y$y=='decoy',],2,jarque.bera.test),function(x)x$p.value))
# 



# 4. integration with ROC_curve analysis : to avoid the problem of  lack of  normality in  data based on shapiro t_test normalization test -------------------------------------------
load('outputData/glm/integrated_fgsea_results_edc_decoy_moa.RData')
library(PRROC)
roc_auc<-list()
for (i in 1:length(fgs)){               # for each network in fgs  
  edc_ind<-which(fgs[[i]]$y=='edc')     # ONLY AVERAGE of NES scores AND their SD FOR EDCS
  x<-fgs[[i]]$x[edc_ind,]  
  y<-rep(1,length(edc_ind))
  
  dec_ind<-which(fgs[[i]]$y=='decoy')   # ONLY AVERAGE of NES scores AND their SD FOR decoys
  xx<-fgs[[i]]$x[dec_ind,]              
  yy<-rep(0,length(dec_ind))
  mat<-rbind(x,xx)
  response<-c(y,yy)
  
     rocs<-rep(NA,ncol(mat))           # calculation of ROC curve using EDC and decoy as label VS the vector of NES scores
     
     for (j in 1:ncol(mat)){
       roc<-roc.curve(scores.class0 = mat[,j], weights.class0 = response,curve = F)
      # print(roc$auc)
       rocs[j]<-roc$auc
       }
  
  dt<-as.data.frame(rocs)
  names(dt)<-'ROC_AUC'
  dt$network<-names(fgs)[[i]]
  dt$pathways<-colnames(mat)

  roc_auc[[i]]<-dt
}
names(roc_auc)<-names(fgs)
roc_data<-do.call(rbind,roc_auc)
load('outputData/patway_glm_coefficients/NESavg_sd_EDCs_stability_coefficients_pathways_moa.RData')
rm(all_data_roc,all_data_final)
all_data_roc<-merge(final_data,roc_data,by=c('network','pathways'))

mean_cof<-by(all_data_roc$glm_coefs,list(all_data_roc$pathways),mean)
mean_nes_edc_across_net<-by(all_data_roc$Average_NES_EDC,list(all_data_roc$pathways),mean)
mean_nes_dec_across_net<-by(all_data_roc$Average_NES_decoy,list(all_data_roc$pathways),mean)
mean_roc<-by(all_data_roc$ROC_AUC,list(all_data_roc$pathways),mean)
maxx_roc<-by(all_data_roc$ROC_AUC,list(all_data_roc$pathways),max)

if(all(names(mean_cof)==names(mean_nes_edc_across_net)) &                            # Checking the order of names before integration to be unique
   all(names(mean_cof)==names(mean_nes_dec_across_net)) &
   all(names(mean_cof)==names(mean_roc)) &
   all(names(mean_cof)==names(maxx_roc))){
   mean_cof_each_pathway<-as.data.frame(list(as.numeric(mean_cof),                    # mean of GLM coefficients across all networks for each pathway
                                             as.numeric(mean_roc),                    # mean of ROC across all networks for each pathway
                                             as.numeric(maxx_roc),                    # max  of network with maximum network ROC
                                             as.numeric(mean_nes_edc_across_net),     # mean of NES for EDCs across all layers for each pathway
                                             as.numeric(mean_nes_dec_across_net),     # mean of NES DECOYS   across all layers for each pathway
                                             names(mean_cof)))}
colnames(mean_cof_each_pathway)<-c('average_pathway_coefficient',
                                   'average_pathway_ROC',
                                   'max_pathway_ROC',
                                   'average_pathway_NES_edc',
                                   'average_pathway_NES_decoy',
                                   'pathways')
all_data_final<-merge(all_data_roc,mean_cof_each_pathway)

# integration with the pathway size
load('outputData/toxdb2gene_final.RData')                                             # all pathways 
pathways_length<-as.data.frame(sapply(moa_pathways, length));colnames(pathways_length)='pathway_length' #length of pathway
pathways_length$pathways<-rownames(pathways_length)
all_data_final<-merge(all_data_final,pathways_length)

# save
save(all_data_final,file='outputData/patway_glm_coefficients/GLM_univariate_all_scores.RData')
write.csv(all_data_final,file='outputData/excel_files/NESavg_sd_EDCs_stability_coefficients_pathways_moa.csv')


# 5. For presentation to Janni the pathways were seleceted based on the ROC curve of more than .7--------------------------------------------
load('outputData/patway_glm_coefficients/GLM_univariate_all_scores.RData')
all_data_final<-all_data_final[!all_data_final$glm_coefs==0,]  
library(magrittr)
all_data_final$pathway_type<-all_data_final$annotated_pathways %>% strsplit('_') %>% lapply('[',1) %>% unlist() # pathway category wiki, kegg,etc
all_data_final$pathway_type[grep(x=all_data_final$pathway_type,pattern = 'aop:',ignore.case = T)]<-'GO'
all_data_final<-all_data_final[all_data_final$ROC_AUC>.7,] #pathways with ROC_auc of more than 0.7
for (i in unique(all_data_final$pathway_type)){
  to_save<-all_data_final[all_data_final$pathway_type==i,c("annotated_pathways",
                                                           "average_pathway_NES_edc",
                                                           "average_pathway_NES_decoy",
                                                           "pathway_length")]
  to_save<-unique(to_save)
  xlsx::write.xlsx(to_save,file = 'outputData/excel_files/for_janni.xlsx',sheetName = i,append = T,row.names = F)  # appending all in one excel file 
}





# 6. Spearman correlation comparison between the  profile of pathway activation score for EDCs across in vitro and in vivo data layers --------

load('outputData/patway_glm_coefficients/GLM_univariate_all_scores.RData')
in_vivo_networks<-unique(all_data_final$network)[grep(x=unique(all_data_final$network),pattern = 'invivo')]   # in vivo
in_vitro_networks<-unique(all_data_final$network)[-grep(x=unique(all_data_final$network),pattern = 'invivo')] # in vitro
spearman_pathway_activation_profile<-matrix(NA, nrow = length(in_vivo_networks), ncol = length(in_vitro_networks),dimnames = list(in_vivo_networks,in_vitro_networks))

# rows are in vivo and columns are in vitro 
for (i in 1:length(in_vivo_networks)){
  for (j in 1:length(in_vitro_networks)){
    vv<-all_data_final[all_data_final$network==in_vivo_networks[i],c('pathways','network','Average_NES_EDC')];colnames(vv)<-c('pathways','network_vivo','Average_NES_EDC_vivo')
    vt<-all_data_final[all_data_final$network==in_vitro_networks[j],c('pathways','network','Average_NES_EDC')];colnames(vt)<-c('pathways','network_vitro','Average_NES_EDC_vitro')
    vv_vt<-merge(vv,vt)
    spearman_pathway_activation_profile[i,j]<-cor(vv_vt$Average_NES_EDC_vivo,vv_vt$Average_NES_EDC_vitro,method = 'spearman')
  }
}
colnames(spearman_pathway_activation_profile)<-gsub(x=colnames(spearman_pathway_activation_profile),pattern = "Consensus_LINCS_HEPG2" ,replacement = "Consensus_LINCS_HEPG2_1_day" )
# gplots::heatmap.2(spearman_pathway_activation_profile,trace = 'none', col = colorRampPalette(c('blue', 'red'))(12),margins = c(20,25),
#                   density.info = 'none',srtCol = 45,keysize = 1,key.xlab = 'Spearman correlation',key.title = NA,cexRow = 1,cexCol = 1.2)
# 

library(ComplexHeatmap)
library(circlize)
mat=spearman_pathway_activation_profile
col_fun<-colorRamp2(c(0,1),c('yellow','red'))
cn<-colnames(mat)
colnames(mat)<-c()
column_ha<-HeatmapAnnotation(spearman=anno_boxplot(mat,height = unit(4,'cm')),show_annotation_name =  F)
botom_ha<-HeatmapAnnotation(text=anno_text(cn,rot=45),show_annotation_name = F)
row_ha<-rowAnnotation(spearman=anno_boxplot(mat,width  = unit(6,'cm'),axis_param = list(side='bottom',labels_rot=45)),show_annotation_name =  F)
ht<-Heatmap(mat,name='Spearman Correlation',col = col_fun,
           # top_annotation = column_ha,
            bottom_annotation = botom_ha,
            #right_annotation = row_ha,
            height = unit(12,'cm'),border = )
draw(ht,padding=unit(c(0,4,0,10),'cm'),heatmap_legend_side='left')














