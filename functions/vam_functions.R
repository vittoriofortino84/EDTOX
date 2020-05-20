# retireving all AVAERAGE scores for 6 MOST INFORMATIVE NETWORKS FOR the compounds in our library
average_vam_score<-function(inp,type='mesh'){
  #load('outputData/VAM_scores_average_most_informative_layers.RData')
  load('outputData/VAM_scores_average_most_informative_layers_moa.RData')  
  vam=matrix(NA, nrow = length(inp), ncol = ncol(all_prob)+1)
  colnames(vam)=c(colnames(all_prob),'VAM_score')
  rownames(vam)=inp
  if (type=='cas'|type=='name'){
    source('functions/annotation_functions.R')
    if (type=='cas') inp=cas2mesh(inp)
    if (type=='name') inp=name2mesh(inp)
      } #end if
  for (i in (1:nrow(vam))){
    tmp=as.numeric(all_prob[inp[i],c(1:ncol(all_prob))])
    vam[i,1:ncol(all_prob)]=tmp
    vam[i,'VAM_score']=all_VAM_score[inp[i]]
  } # end for
  return(as.data.frame(vam))
} # end function

#retireving all HARMONIC SUM scores for 6 MOST INFORMATIVE NETWORKS FOR the compounds in our library
harmonic_vam_score<-function(inp,type='mesh'){
  #load('outputData/VAM_scores_hs_most_informative_layers.RData')
  load('outputData/VAM_scores_hs_most_informative_layers_moa.RData')
  vam=matrix(NA, nrow = length(inp), ncol = ncol(all_prob)+1)
  colnames(vam)=c(colnames(all_prob),'VAM_score')
  rownames(vam)=inp
  if (type=='cas'|type=='name'){
    source('functions/annotation_functions.R')
    if (type=='cas') inp=cas2mesh(inp)
    if (type=='name') inp=name2mesh(inp)
  } #end if
  for (i in (1:nrow(vam))){
    tmp=as.numeric(all_prob[inp[i],c(1:ncol(all_prob))])
    vam[i,1:ncol(all_prob)]=tmp
    vam[i,'VAM_score']=all_VAM_score[inp[i]]
  } # end for
  return(as.data.frame(vam))
} # end function

# this function repeats the whole RWR-FGSEA-GLM_prediction pipeline using a list of chemicals with MIES 
mie2vam<-function(mie_list){
  source("functions/pipeline.R")                                                 # pipeline functions
  load('outputData/toxdb2gene_final.RData')                                      # pathways file
  #load('outputData/glm/used_pathways_glm.RData')                                # pathways used in GLM modeling (script 5_2)
  load('outputData/glm/used_pathways_glm_moa.RData')                             # pathways used in GLM modeling (script 5_2) MOA pathways
  tox2db<-c(go2gene,kegg2gene,reactome2gene,wikiptw2gene,msigdb2gene)            # tox2db is all pathways
  library(BiocParallel)
  ##
   load("outputData/network/wTO_DM.RData")
   #1. DM_hepatocytes.1
   ptw<-tox2db[which(names(tox2db) %in% used_pathways[[1]])]
   p<-bpstart(MulticoreParam(40))
   fgsea_res<-pipeline(wTO_DM[[1]],ptw,0.1,1000,mie_list,p)
   bpstop(p)
   save(fgsea_res,file = 'outputData/mie2vam_pred/DM_hep.RData')
   gc()
   ## 
   #2. DM_liver.1
   ptw<-tox2db[which(names(tox2db) %in% used_pathways[[2]])]
   p<-bpstart(MulticoreParam(40))
   fgsea_res<-pipeline(wTO_DM[[2]],ptw,0.1,1000,mie_list,p)
   bpstop(p)
   save(fgsea_res,file = 'outputData/mie2vam_pred/DM1.RData')
   gc()
   ## 
   #3. DM_liver.5
   ptw<-tox2db[which(names(tox2db) %in% used_pathways[[4]])]
   p<-bpstart(MulticoreParam(40))
   fgsea_res<-pipeline(wTO_DM[[4]],ptw,0.05,1000,mie_list,p)
   bpstop(p)
   save(fgsea_res,file = 'outputData/mie2vam_pred/DM5.RData')
   gc()
   ## 
   
   load('outputData/network/wTO_hep_cons.RData')
   #4. Consensus hepatocytes
   ptw<-tox2db[which(names(tox2db) %in% used_pathways[[8]])]
   p<-bpstart(MulticoreParam(40))
   fgsea_res<-pipeline_consensus(wTO_hep_cons,ptw,0.05,1000,mie_list,p)
   bpstop(p)
   save(fgsea_res,file = 'outputData/mie2vam_pred/hep_cons.RData')
   gc()
   #
   
   load('outputData/network/wTO_LINCS_cons.RData')
   #5. Consensus lincs
   ptw<-tox2db[which(names(tox2db) %in% used_pathways[[9]])]
   p<-bpstart(MulticoreParam(40))
   fgsea_res<-pipeline_consensus(wTO_LINCS_cons,ptw,0.1,700,mie_list,p)
   bpstop(p)
   save(fgsea_res,file = 'outputData/mie2vam_pred/lincs_cons.RData')
   gc()
   
   # 6. PPI
   ptw<-tox2db[which(names(tox2db) %in% used_pathways[[10]])]
   p<-bpstart(MulticoreParam(40))
   fgsea_res<-pipeline_ppi(ptw,0.85,1000,mie_list,p)
   bpstop(p)
   save(fgsea_res,file = 'outputData/mie2vam_pred/ppi.RData')
   gc()
  # after RWR-FGSEA we need NES treatment
  source('functions/glm_functions.R')
  source('functions/general_functions.R')
  load('outputData/toxdb2gene_final.RData')
  setwd('outputData/mie2vam_pred/')
  li<-list.files(path='.',pattern = 'RData')
  fgs<-data_preparation(li,'n',toxdb_names) #NES_scores with  no scaling
  names(fgs)<-name_correct(names(fgs))
  setwd("../..")
  # prediction of class probilities
  source('functions/general_functions.R')
  # load("outputData/glm/glm_models_final.RData") #modls from edcs and decoy
  # names(modls)<-name_correct(names(modls))
   load("outputData/glm/glm_models_final_moa.RData") #new modls from edcs and decoys based on moa pathways
   
   modls<-modls[names(fgs)]
   for (i in 1:length(modls)){        # checkig the missing columns in the results with coefficents names and corrects missing ones
   compound_number<-nrow(fgs[[i]]$x)  
   nc<-length(which(!rownames(modls[[i]]$cof) %in% colnames(fgs[[i]]$x)))
   new_mat<-matrix(0,nrow = compound_number,ncol = nc-1)
   new_features_names<-rownames(modls[[i]]$cof)[which(!rownames(modls[[i]]$cof) %in% colnames(fgs[[i]]$x))]
   new_features_names<-new_features_names[!new_features_names %in% "(Intercept)"]
   colnames(new_mat)<-new_features_names
   fgs[[i]]$x<-cbind(fgs[[i]]$x,new_mat)
   fgs[[i]]$x<-as.matrix(fgs[[i]]$x[,rownames(modls[[i]]$cof)[!rownames(modls[[i]]$cof) %in% "(Intercept)"]])
      } #end for
   
   all_compounds_prob<-prob_predict(modls,fgs)  #prdiction of the obtained class probilities 
    # VAM score as harmonic sum of the predicted values
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
    #hs<-scales:::rescale(hs,to = c(0, 1))
    all_prob$harmonic_VAM_score<-harmonic_sum(all_prob,1:6)
    all_prob$average_VAM_score<-mean_function(all_prob,1:6)
    return(all_prob)
}# end function
