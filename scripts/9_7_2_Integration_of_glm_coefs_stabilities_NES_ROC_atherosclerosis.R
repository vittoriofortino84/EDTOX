
# 1. preparation of coefficients and stabilites for GLM  -------------------------
# cross validation was done four times using 5-repeated 5-fold
# needed path and files for each disease to be changed 
# input
load('output/glm_models/two_class_glm_models_atherosclerosis.RData')                       # glm model for disease
model_name<-atherosclerosis_two_class                                                      # name of models
path_to_CV_results_for_disease<-'output/glm_models/CV_two_class/atherosclerosis/'          # path to CV results
# output
RData_file<-'output/atherosclerosis_coefficients_stabilities_two_class.RData'              # OUTPUT save file name RData
## Just the above lines should be changed for each diesease

source('functions/annotation_functions.R')
setwd(path_to_CV_results_for_disease)                                 # Cross validation results path
li<-list.files('.',pattern = '.RData')
load(li[[1]])
cof<-lapply(CV_results, function(x)x$coefs)                           # 
for (i in 2:length(li)){                                              # concatanating stabilities from different experiments OF cv (4 RESULTS)
  load(li[[i]])
  for (j in 1:length(cof)){
    cof[[j]]<-c(cof[[j]],CV_results[[j]]$coefs)
  } # end for j networks
} #end for i number of 5-fold-CV experiments

setwd('../../../..')

if (all(names(model_name)==names(cof))) {
  
  stability_cof<-list()                                      # making a list of stabilites and coefiicients
  for (i in 1:length(cof)){
    coffs<-as.data.frame(cof[[i]])                      
    stability<-apply((abs(coffs)>0)*1, 1, sum)            # sum of the coefficients more than zero for all CV folds
    standard_dev<-apply(coffs, 1, sd)                     # standard deviation of CV results for each pathway
    main_coefs<-model_name[[i]]$cof[-1]                # coefficients  from the final model   inercept was removed
    stability_cof[[i]]<-as.data.frame(list(main_coefs,stability,standard_dev),col.names = c('coef','stability','SD'))
    stability_cof[[i]]$network<-names(cof)[[i]]
    stability_cof[[i]]$pathways<-rownames(stability_cof[[i]])
  } #end for
  names(stability_cof)<-names(cof)
  data<-do.call(rbind,stability_cof)
}
#data$pathways<-pathway_annotate(data$pathways)
save(data,file=RData_file)


# 2. merging of coefficients with average of NES scores for positive disease labeled edcs-------------------------------------------------
rm(list=ls())
# input 
load('output/glm_models/two_class_glm_models_atherosclerosis.RData')    #GLM coefficients and stabilities
model_name<-atherosclerosis_two_class
load('output/train_set_atherosclerosis_binary.RData') #TRAINING SET (nes SCORES)
fgs<-train_set_athero
disease_positive_label<-"edc_positive_atherosclerosis"                  # disease positive label 
disease_negative_label<-"edc_negative_atherosclerosis"                  # disease negative label
glm_stability_file_path<-'output/atherosclerosis_coefficients_stabilities_two_class.RData' # output of previuos step
# output
RData_file<-'output/atherosclerosis_coefficients_stabilities_two_class_NES.RData'            # OUTPUT
## Just the above lines should be changed for each diesease

source('functions/annotation_functions.R')
avg_nes<-list()
for (i in 1:length(fgs)){                                                # for each network in train_set_athero 
  
  edc_ind<-which(fgs[[i]]$y==disease_positive_label)                    # ONLY AVERAGE AND SD FOR EDCS with positive disease labels
  x<-fgs[[i]]$x[edc_ind,]                                                # ONLY AVERAGE AND SD FOR EDCs  with positive disease labels
  avg_pos<-apply(x,2,mean)                                                         # average of all pathway scores EDCs +
  sd_nes_pos<-apply(x,2,sd)                                                        # standard deviation of all pathway scores EDCs +
  
  neg_ind<-which(fgs[[i]]$y==disease_negative_label)                     # ONLY AVERAGE AND SD FOR EDCS with NEGATIVE disease labels
  xx<-fgs[[i]]$x[neg_ind,]                                               # ONLY AVERAGE AND SD FOR EDCs  with NEGATIVE disease labels
  avg_neg<-apply(xx,2,mean)                                                        # average of all pathway scores EDCs -
  sd_nes_neg<-apply(xx,2,sd)                                                       # standard deviation of all pathway scores EDCs -
  
  dt<-as.data.frame(list(avg_pos,sd_nes_pos,avg_neg,sd_nes_neg))
  names(dt)<-c('Average_NES_positive_labels',
               'SD_NES_positive_labels',
               'Average_NES_negative_labels',
               'SD_NES_negative_labels')
  dt$network<-names(fgs)[[i]]
  dt$pathways<-rownames(dt)
  avg_nes[[i]]<-dt
}
names(avg_nes)<-names(fgs)
all_nes_avg_sd<-do.call(rbind,avg_nes)
load(glm_stability_file_path)




all_data<-merge(all_nes_avg_sd,data,by=c('network','pathways'))
colnames(all_data)<-c("network",                             # columns 
                      "pathways",
                      "Average_NES_positive",
                      "SD_NES_positive",
                      "Average_NES_negative",
                      "SD_NES_negative",
                      "glm_coefs",
                      "CV_stability",
                      "CV_coefficients_SD")   

source('functions/annotation_functions.R')
all_data$annotated_pathways<-pathway_annotate(all_data$pathways)
list_pathways  <- lapply(model_name, function(x){y=x$sel
y=y[-1] #the first element is intercept so we remove it for upset plot
y})
all_pathways <- unique(unlist(list_pathways))
int_data <- lapply(list_pathways, function(x){
  idx <- which(all_pathways %in% x)
  res<-rep(0, length(all_pathways))
  res[idx]<-1
  res})
upset_df<-as.data.frame(int_data , as.is=T, stringsAsFactors=F, check.names=F)
rownames(upset_df) <- all_pathways
sum_of_voting_layers_for_pathway<-apply(upset_df,1 ,sum)
dd<-as.data.frame(sum_of_voting_layers_for_pathway)
dd$pathways<-rownames(dd)
final_data<-merge(all_data,dd,all=T)
final_data$network<-name_correct(final_data$network)   #
final_data<-final_data[!final_data$glm_coefs==0,]

save(final_data,file=RData_file)                # SAVING the file coeffieicnts and NES averages for both EDC + and - labels


# 3. Performing ROC analysis on CAS as a univariate test and integ --------
rm(list=ls())
# input
load('output/train_set_atherosclerosis_binary.RData')                   # training set of the  disease
fgs<-train_set_athero                                                            # training set of the  disease
disease_positive_label<-"edc_positive_atherosclerosis"                  # disease positive label 
disease_negative_label<-"edc_negative_atherosclerosis"                  # disease negative label
path_to_final_data<-'output/atherosclerosis_coefficients_stabilities_two_class_NES.RData'   # path to final data of glm coefs, stabilites and NES scores
# output
final_Rdata<-'output/atherosclerosis_coefficients_stabilities_two_class_NES_ROC.RData'          # RDATA file
final_excel<-'output/excel_files/final_data_Atherosclerosis.csv'                                        # Excel file
## Just the above lines should be changed for each diesease 

library(PRROC)
roc_auc<-list()
for (i in 1:length(fgs)){                                # for each network in fgs  
  edc_ind<-which(fgs[[i]]$y==disease_positive_label)     # ONLY AVERAGE of NES scores AND their SD FOR EDCS
  x<-fgs[[i]]$x[edc_ind,]  
  y<-rep(1,length(edc_ind))
  
  dec_ind<-which(fgs[[i]]$y==disease_negative_label)   # ONLY AVERAGE of NES scores AND their SD FOR decoys
  xx<-fgs[[i]]$x[dec_ind,]              
  yy<-rep(0,length(dec_ind))
  mat<-rbind(x,xx)
  response<-c(y,yy)
  
  rocs<-rep(NA,ncol(mat))                               # calculatio of ROC curve using EDC + and EDC - as label VS the vector of NES scores
  
  for (j in 1:ncol(mat)){
    roc<-roc.curve(scores.class0 = mat[,j], weights.class0 = response,curve = F)
    # print(roc$auc)
    rocs[j]<-roc$auc
  } # for j 
  
  dt<-as.data.frame(rocs)
  names(dt)<-'ROC_AUC'
  dt$network<-names(fgs)[[i]]
  dt$pathways<-colnames(mat)
  
  roc_auc[[i]]<-dt
} # for i
names(roc_auc)<-names(fgs)
roc_data<-do.call(rbind,roc_auc)
roc_data$network<-gsub(x=roc_data$network,pattern = "Consensus_LINCS_HEPG2_HEPG2",replacement = "Consensus_LINCS_HEPG2")
roc_data$network<-gsub(x=roc_data$network,pattern = "Consensus_Rat_invitro_Drug Matrix__TG_GATEs",replacement = "Consensus_Rat_invitro_Drug Matrix_TG_GATEs" )
load(path_to_final_data)

all_data_roc<-merge(final_data,roc_data,by=c('network','pathways'))

mean_cof<-by(all_data_roc$glm_coefs,list(all_data_roc$pathways),mean)
mean_nes_edc_pos_across_net<-by(all_data_roc$Average_NES_positive,list(all_data_roc$pathways),mean)
mean_nes_edc_neg_across_net<-by(all_data_roc$Average_NES_negative,list(all_data_roc$pathways),mean)
mean_roc<-by(all_data_roc$ROC_AUC,list(all_data_roc$pathways),mean)
maxx_roc<-by(all_data_roc$ROC_AUC,list(all_data_roc$pathways),max)

if(all(names(mean_cof)==names(mean_nes_edc_pos_across_net)) &                            # Checking the order of names before integration to be unique
   all(names(mean_cof)==names(mean_nes_edc_neg_across_net)) &
   all(names(mean_cof)==names(mean_roc)) &
   all(names(mean_cof)==names(maxx_roc))){
  mean_cof_each_pathway<-as.data.frame(list(as.numeric(mean_cof),                    # mean of GLM coefficients across all networks for each pathway
                                            as.numeric(mean_roc),                    # mean of ROC across all networks for each pathway
                                            as.numeric(maxx_roc),                    # max  of network with maximum network ROC
                                            as.numeric(mean_nes_edc_pos_across_net),     # mean of NES for EDCs across all layers for each pathway
                                            as.numeric(mean_nes_edc_neg_across_net),     # mean of NES DECOYS   across all layers for each pathway
                                            names(mean_cof)))}
colnames(mean_cof_each_pathway)<-c('average_pathway_coefficient',
                                   'average_pathway_ROC',
                                   'max_pathway_ROC',
                                   'average_pathway_NES_positive',
                                   'average_pathway_NES_negative',
                                   'pathways')
all_data_final<-merge(all_data_roc,mean_cof_each_pathway)

# integration with the pathway size
load('output/toxdb2gene_final.RData')                 # all pathways from phase 1 
pathways_length<-as.data.frame(sapply(moa_pathways, length));colnames(pathways_length)='pathway_length' # length of pathway
pathways_length$pathways<-rownames(pathways_length)
all_data_final<-merge(all_data_final,pathways_length)

# save
save(all_data_final,file=final_Rdata)
write.csv(all_data_final,file=final_excel)



# 4. Data presentation to Janni for interpreatation -----------------------


rm(list=ls())
#input
all_data_final<-read.csv('output/excel_files/final_data_Atherosclerosis.csv',stringsAsFactors = F)
# output
output_file<-'output/excel_files/for_janni_Atherosclerosis.xlsx'

all_data_final<-all_data_final[!all_data_final$glm_coefs==0,]  
library(magrittr)
all_data_final$pathway_type<-all_data_final$annotated_pathways %>% strsplit('_') %>% lapply('[',1) %>% unlist() # pathway category wiki, kegg,etc
all_data_final$pathway_type[grep(x=all_data_final$pathway_type,pattern = 'aop:',ignore.case = T)]<-'GO'
all_data_final<-all_data_final[all_data_final$ROC_AUC>.7,] #pathways with ROC_auc of more than 0.7
for (i in unique(all_data_final$pathway_type)){
  to_save<-all_data_final[all_data_final$pathway_type==i,c("annotated_pathways",
                                                           "average_pathway_NES_positive",
                                                           "average_pathway_NES_negative",
                                                           "pathway_length")]
  to_save<-unique(to_save)
  xlsx::write.xlsx(to_save,file = output_file,sheetName = i,append = T,row.names = F)  # appending all in one excel file 
}































