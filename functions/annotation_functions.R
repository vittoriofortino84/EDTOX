

#Corrects and uniforms layers 
name_correct<-function(x){
  x<-gsub('.RData','',x)
  x<-gsub('ppi','PPI_STRINGdb',x)
  x<-gsub('DM','Drug Matrix_',x)
  x<-gsub('TGG|Fgsea_tg','TG_Gate',x)
  x<-gsub('_hep','Hepatocytes',x)
  x<-gsub('hep_cons','Consensus_Rat_invitro_DM_TG_gates',x)
  x<-gsub('TGG|Fgsea_tg','TG_Gate',x)
  x<-gsub('lincs_cons','Consensus_LINCS',x)
  x<-gsub('Consensus_LINCS','Consensus_LINCS_HEPG2',x)
  x<-gsub('hep_cons','Consensus_Hepatocytes(DM,TG-gates)',x)
  x<-gsub("Drug Matrix_Hepatocytes",'Drug_Matrix_Rat_invitro_Single_Dose_1_day',x)
  x<-gsub('Drug Matrix_1','Drug_Matrix_Rat_invivo_Single_Dose_1_day',x)
  x<-gsub('Drug Matrix_3','Drug_Matrix_Rat_invivo_Repeated_Dose_3_days',x)
  x<-gsub('Drug Matrix_5','Drug_Matrix_Rat_invivo_Repeated_Dose_5_days',x)
  x<-gsub('TG_Gate_single.high','TG_GATEs_Rat_invivo_Single_Dose_High_1_day',x)
  x<-gsub('TG_Gate_single.low','TG_GATEs_Rat_invivo_Single_Dose_Low_1_day',x)
  x<-gsub('TG_Gate_single.middle','TG_GATEs_Rat_invivo_Single_Dose_Middle_1_day',x)
  x<-gsub('Consensus_Hepatocytes','Consensus_Rat_invitro',x)
  x<-gsub('TG_Gate_human_in_vitro','TG_GATEs_Human_invitro_Single_Dose_1_day',x)
  x<-gsub('TG_Gate_rat_in_vitro','TG_GATEs_Rat_invitro_Single_Dose_1_day',x)
  x<-gsub('TG_Gate_rat_in_vivo_rep_15_day','TG_GATEs_Rat_invivo_Repeated_Dose_15_days',x)
  x<-gsub('TG_Gate_rat_in_vivo_rep_29_day','TG_GATEs_Rat_invivo_Repeated_Dose_29_days',x)
  x<-gsub('TG_Gate_rat_in_vivo_rep_8_day','TG_GATEs_Rat_invivo_Repeated_Dose_8_days',x)
  x<-gsub('PPI (STRINGdb)','PPI_STRINGdb',x,fixed = T)
  x<-gsub('Consensus_Hepatocytes (DM,TG-gata)','Consensus_Rat_invitro_DM_TG_gates',x,fixed = T)
  x<-gsub('Consensus_Rat_invitro (Drug Matrix_,TG-gata)','Consensus_Rat_invitro_DM_TG_GATEs',x,fixed = T)
  x<-gsub("Consensus_LINCS_HEPG2(I,II)",'Consensus_LINCS_HEPG2',x,fixed = T)
  x<-gsub("Consensus_Rat_invitro_DM_TG_gates","Consensus_Rat_invitro_DM_TG_GATEs",x,fixed = T)
  return(x) 
}#end function

# 1. Compound conversions -------------------------------------------------

mesh2name<-function(inp){
  ixns=read.csv('inputData/annotaion/chem_mesh_cas_dictionary.csv',stringsAsFactors = F)
  nn<-inp
  for (i in 1:length(inp)){
    if (length(ixns$name[which(ixns$mesh==inp[i])])>0){
    nn[i]=ixns$name[which(ixns$mesh==inp[i])]
    }
    else {
      nn[i]=''
    }
  }
return(nn)
} #end function

name2mesh<-function(inp){
  ixns=read.csv('inputData/annotaion/chem_mesh_cas_dictionary.csv',stringsAsFactors = F)
  nn<-inp
  for (i in 1:length(inp)){
    if (length(ixns$name[which(tolower(ixns$name)==tolower(inp[i]))])>0){
      nn[i]=ixns$mesh[which(tolower(ixns$name)==tolower(inp[i]))]
    }
    else {
      nn[i]=''
    }
  }
  return(nn)
} #end function

name2cas<-function(inp){
  ixns=read.csv('inputData/annotaion/chem_mesh_cas_dictionary.csv',stringsAsFactors = F)
  nn<-inp
  for (i in 1:length(inp)){
    if (length(ixns$name[which(tolower(ixns$name)==tolower(inp[i]))])>0){
      nn[i]=ixns$cas[which(tolower(ixns$name)==tolower(inp[i]))]
    }
    else {
      nn[i]=''
    }
  }
  return(nn)
} #end function

cas2name<-function(inp,return_back_NA='TRUE'){
  ixns=read.csv('inputData/annotaion/chem_mesh_cas_dictionary.csv',stringsAsFactors = F)
  nn<-inp
  for (i in 1:length(inp)){
    if (length(ixns$name[which(ixns$cas==inp[i])])>0){
    nn[i]=ixns$name[which(ixns$cas==inp[i])]
    }
    else {
      if(return_back_NA=='TRUE')nn[i]=''
    }
  }
  return(nn)
} #end function

cas2mesh<-function(inp,return_back_NA='TRUE'){
  ixns=read.csv('inputData/annotaion/chem_mesh_cas_dictionary.csv',stringsAsFactors = F)
  nn<-inp
  for (i in 1:length(inp)){
    if (length(ixns$mesh[which(ixns$cas==inp[i])])>0){
    nn[i]=ixns$mesh[which(ixns$cas==inp[i])]
    }
    else {
      if(return_back_NA=='TRUE')nn[i]=''
    }
  }
  return(nn)
} #end function

mesh2cas<-function(inp){
  ixns=read.csv('inputData/annotaion/chem_mesh_cas_dictionary.csv',stringsAsFactors = F)
  nn<-inp
  for (i in 1:length(inp)){
    if (length(ixns$cas[which(ixns$mesh==inp[i])])>0){
    nn[i]=ixns$cas[which(ixns$mesh==inp[i])]
    }
    else {
      nn[i]=''
    }
  }
  return(nn)
} #end function

# 2.pathway conversions ---------------------------------------------------

pathway_annotate<-function(inp){
  library(GO.db)
  Go_annotations<-sapply(GOTERM,function(x)x@Term) #Dictionary of GO terms
  keggs<-read.delim(file='inputData/kegg.txt')     #Dictionary of KEGG
  keggs<-keggs[,1:2]
  reacts<-read.csv(file='inputData/reactome.csv')  #Dictionary of reactome
  aops<-read.csv(file='outputData/aop_mat.csv')    #Dictionary of AOPs made by script 1_2
  
  nn<-inp
  
  for (i in 1:length(inp)){
    # GO annotaions
    if (any(grep(pattern = 'GO:',x = inp[i] ))==TRUE){
      nam<-Go_annotations[inp[i]]
      nam<-gsub(' ','_',nam)
      nn[i]<-paste('GO_',nam,sep = '')
    } #end if
    
    # wiki annotations
    if (any(grep(pattern = 'WikiPathways',x=inp[i]))==TRUE){
      nam<-inp[i]
      nam<-gsub('%WikiPathways_20190510%','',nam)
      nam<-gsub('%Homo sapiens','',nam)
      nam<-gsub(' ','_',nam)
      nam<-gsub('WP','_WP',nam)
      nn[i]<-paste('WIKI_',nam,sep = '')
    } #end if
    
    # kegg pathways
    if (any(grep(pattern = 'KEGG:',x=inp[i]))==TRUE){
      nam<-gsub('KEGG:','',inp[i])
      if (length(which(keggs$TermID==nam))>0) {
        nam<-keggs$Term[which(keggs$TermID==nam)]
      }else {
        warning(paste('No name for the identifier',nam))
      } #end if warning
      nn[i]<-paste('KEGG_',nam,sep = '')
    } #end if
    
    # Reactome pathways
    if (any(grep(pattern = 'R-HSA-',x=inp[i]))==TRUE){
      nam<-inp[i]
      nam<-reacts$V1[which(reacts$V2==nam)]
      nn[i]<-paste('REACTOME_',nam,sep = '')
    } #end if
    
    nn[i]<-toupper(nn[i])
    
    # AOP annotaions
    if (any(colnames(aops) %in% nn[i])==T){
      aopind<-which(colnames(aops) == nn[i])
      aopc<-aops$X[aops[,aopind]==1]
      aopc<-gsub(' ',',',do.call(paste,as.list(aopc)))
      nn[i]<-paste(aopc,'_',nn[i],sep = '')
    } #end if
  }#end for
  
  return(nn)
} #end function

kegg_reactome_2_mainclass<-function(list_of_pathways,type='KEGG'){
  require(qdapRegex)
  if(type=='GO'){
    # number of AOPs
    Go_p<-list_of_pathways
    Go_p<-gsub(x=Go_p,pattern = '_',',')
    p_name<-qdapRegex::rm_between(Go_p, 'Aop:', ',', extract=TRUE)
    aops<-unique(unlist(p_name[!is.na(p_name)]))
    n_aops<-length(aops)
    return(list(Number_aops=n_aops,AOPS=aops))
  }
  
  if(type=='REACTOME'){
    load('inputData/annotaion/reactome.RData')           # reactome dictionary
    res<-(rep(NA,length(list_of_pathways)))      
    for (i  in 1:length(list_of_pathways)){
      inp<-list_of_pathways[i]
      find_rel = table(react_hie[which(react_hie[,2] %in% inp),1])
      reactome[names(find_rel),]
      test_val = table(react_top[which(react_top$pathway %in% c(names(find_rel),inp)),'top_level_pathway'])
      if (length(reactome[names(test_val),'term'])==1)res[i]<-reactome[names(test_val),'term']
      res[is.na(res)]<-'Other Category'
    }
    wh<-length(list_of_pathways)-length(which(res=='Other Category'))
    t<-table(res)
    t<-t[!names(t) %in% 'Other Category']
    t<-as.data.frame(t/wh*100)    # percemtage of the pathways
    #t<-as.data.frame(t)          # number of the pathways
    if(ncol(t)>1)colnames(t)<-c('class','prop')
    return(list(names=res,prop=t))
  }
  
  
  if(type=='KEGG'){
    load('inputData/annotaion/kegg_category_lib.RData')  # kegg dictionary       
    res<-(rep(NA,length(list_of_pathways)))      
    for (i  in 1:length(list_of_pathways)){
      inp<-list_of_pathways[i]
      inp<-gsub(x=inp,pattern ='KEGG:|_|hsa',replacement = '')
      hh<-names(unlist(sapply(kegg_category_lib, function(x) which(inp %in% x))))
      if(length(hh)==1){res[i]=hh}
      res[is.na(res)]<-'Other Category'
    }
    wh<-length(list_of_pathways)-length(which(res=='Other Category'))
    t<-table(res)
    t<-t[!names(t) %in% 'Other Category']
    t<-as.data.frame(t/wh*100) #percentage of the pathways
    #t<-as.data.frame(t)       # number of the pathways
    if(ncol(t)>1)colnames(t)<-c('class','prop')
    return(list(names=res,prop=t)) #returning the percentage for all  and the category for each pathway
  }
  
}


  #  res<-gsub(x=res,
  #            pattern = 'Cell Cycle|Autophagy|Protein localization|Chromatin organization|Programmed Cell Death|Cell-Cell communication|Transport of small molecules|Cellular responses to external stimuli|Vesicle-mediated transport',
  #            replacement ='Cellular Processes' )
  #  res<-gsub(x=res,pattern = '[()]',replacement ='')
  #  res<-gsub(x=res,pattern = 'DNA Repair|DNA Replication|Gene expression Transcription',replacement ='Genetic Information Processing')
  #  res<-gsub(x=res,pattern = 'Signal Transduction',replacement ='Environmental Information Processing')
  # # #res<-gsub(x=res,pattern = '',replacement ='Diseases')
  #  res<-gsub(x=res,pattern = 'Developmental Biology|Immune System|Reproduction|Neuronal System|Muscle contraction|Circadian Clock',replacement ='Organismal Systems')
  #  res<-gsub(x=res,pattern = 'Metabolism of RNA|Metabolism of proteins',replacement ='Metabolism')

# 
# 
# unify_kegg_reactome_categories<-function(inp){
#   
#   #Dictionary
#   pathway_dictionary<-list(
#     Cell_growth_and_death = list(               KEGG = c("Cellular Processes, Cell Growth and Death"),#ok
#                                                 REACTOME = c("Cell Cycle","Programmed Cell Death")),
#   
#     Cell_transport_and_catabolism = list(       KEGG = c("Cellular Processes, Transport and catabolism"), # we do not use
#                                                 REACTOME = c("Autophagy","Protein localization")),
# 
#     Organelle_biogenesis_and_maintenance = list(KEGG =c("Genetic Information Processing, Translation", # we do not
#                                                         "Genetic Information Processing, Folding, sorting and degradation"),
#                                                 REACTOME = c("Organelle biogenesis and maintenance")),
#     
#     Replication_and_repair = list(              KEGG =c("Genetic Information Processing, Replication and repair"), 
#                                                 REACTOME = c("DNA Repair",
#                                                              "DNA Replication")),
#     
#     Transcription = list(                       KEGG =c("Genetic Information Processing, Transcription"), 
#                                                 REACTOME = c("Gene expression (Transcription)")),
#     
#     Metabolism = list(                          KEGG =c("Metabolism"), 
#                                                 REACTOME = c("Metabolism",
#                                                              "Metabolism of proteins","Metabolism of RNA")),
#     
#     Environmental_information_processing = list(KEGG =c("Environmental Information Processing"),
#                                                 REACTOME = c("Signal Transduction",
#                                                              "Extracellular matrix organization",
#                                                              "Transport of small molecules")),
#     
#     Chromatin_organization = list(              REACTOME = c("Chromatin organization")),
#     
#     Cellular_responses_stimuli = list(          REACTOME = c("Cellular responses to external stimuli")),
#     # Endocrine_system = list(                    KEGG =c("Organismal Systems, Endocrine system"), 
#     #                                             REACTOME = c("Reproduction")),
#     # Immune_system = list(                       KEGG =c("Organismal Systems, Immune system"), 
#     #                                             REACTOME = c("Immune systems")),
#     # Nervous_system = list(                      KEGG = c("Organismal Systems, Nervous system"), 
#     #                                             REACTOME = c("Neuronal System")),
#     # Environ_adaptation = list(                  KEGG = c("Organismal Systems, Environmental adaptation"), 
#     #                                             REACTOME = c("Circadian Clock")),
#     #Digestive_system = list(                    REACTOME = c("Digestion and absorption")),
#     Organismal_systems = list(                  KEGG = c("Organismal Systems"),
#                                                 REACTOME = c("Immune System",
#                                                              "Neuronal System",
#                                                              "Digestion and absorption",
#                                                              "Circadian Clock",
#                                                              "Reproduction")),
#     
#     Human_disease =   list(                     KEGG = c("Human Diseases"),
#                                                 REACTOME = c("Disease")),
#     
#     Cell_Communication =  list(                 KEGG = c("Cell Communication"),
#                                                 REACTOME = c("Cell-Cell communication"))
#     
#      )
#   res<-rep(NA,length(inp))
#   for (i in 1:length(inp)){
#     hh<-names(unlist(sapply(pathway_dictionary, function(x) which(inp[i] %in% unlist(x)))))
#     if(length(hh)==1){res[i]=hh}else{res[i]='Other_Category'}
#   }
#   return(res)
#   
# }
# 
# pathway2kegg_translation<-function(list_pathways,fam_inp='Wiki'){
#   all_mappings<-read.table(file = 'inputData/annotaion/all_mappings.tsv',stringsAsFactors = F,sep='\t',header = F)
#   colnames(all_mappings)<-c('pathway_name','ID','family','mapp_type','category','mapped_ID','mapped_family')
#   all_mappings<-all_mappings[which(all_mappings$family=='kegg' | all_mappings$mapped_family=='kegg'),]
#   POOL_LIST_PATHWAYS<-unique(c(as.character(all_mappings$ID),as.character(all_mappings$mapped_ID)))
#   
#   equivalent_in_kegg<-rep(NA,length(list_pathways))
#   for (i in 1: length(list_pathways)){
#     # Given a pathway p
#     res<-NA
#     p<-list_pathways[i]
#     # name and type preprocessing
#     if(fam_inp=='REACTOME'){p_name<-p}
#     if(fam_inp=='Wiki'){
#     p_name<-paste('WP',qdapRegex::rm_between(p, '%WP', '%', extract=TRUE)[[1]],sep='')}
#     # checking for mapping
#  
#       # if the pathways in wiki and reactome is in the pool list of mappings with kegg
#       if (p_name %in% POOL_LIST_PATHWAYS){
#         
#         row_ind<-which(all_mappings$mapped_ID==p_name | all_mappings$ID==p_name)[1]
#         col_ind<-colnames(all_mappings)[which(all_mappings[row_ind,]==p_name)]
#         if (col_ind =="mapped_ID")res<-all_mappings[row_ind,'ID']
#         if (col_ind =="ID")res<-all_mappings[row_ind,"mapped_ID"]
#         
#         
#       } # end if in the pool list of pathways
#  
#     equivalent_in_kegg[i]<-res
#   } # end for
#   equivalent_in_kegg<-gsub(x=equivalent_in_kegg,pattern = 'path:',replacement = '')
#   return(equivalent_in_kegg)
# } # end function
# 
# path2mainclass_old<-function(list_of_pathways){
#   res<-(rep(NA,length(list_of_pathways)))
#   
#   for (i  in 1:length(list_of_pathways)){
#     inp<-list_of_pathways[i]
#     if (any(grep(x=inp,pattern = 'GO')))res[i]<-'GO'
#     if (any(grep(x=inp,pattern = 'Wiki')))res[i]<-keg2category(pathway2kegg_translation(inp,fam_inp='Wiki'))
#     if (any(grep(x=inp,pattern = 'KEGG')))res[i]<-keg2category(inp)
#     # if (any(grep(x=inp,pattern = 'R-HSA')))res[i]<-unify_kegg_reactome_categories(reactom2category(inp))
#     if (any(grep(x=inp,pattern = 'R-HSA')))res[i]<-keg2category(pathway2kegg_translation(inp,fam_inp='REACTOME'))
#   }
#   
#   return(res)
# }










