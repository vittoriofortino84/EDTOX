
# #chembl to cas conversion -----------------------------------------------

# this function was modified from the package webchem to retrieve the cas numbers from drug names using webchem package

library(webchem)
newgetid<-function (query, language = "en", match = c("best", "first", 
                                                      "all", "ask", "na"), verbose = TRUE) 
{
  require(httr)
  require(jsonlite)
  match <- match.arg(match)
  foo <- function(query, language, match, verbose) {
    query <- URLencode(query)
    limit <- 50
    qurl <- paste0("wikidata.org/w/api.php?action=wbsearchentities&format=json&type=item")
    qurl <- paste0(qurl, "&language=", language, "&limit=", 
                   limit, "&search=", query)
    if (verbose) 
      message("Querying ", qurl)
    Sys.sleep(0.3)
    cont <- fromJSON(content(GET(qurl, user_agent("webchem (https://github.com/ropensci/webchem)")), 
                             "text"))
    search <- cont$search
    if (length(search) == 0) {
      if (verbose) 
        message("Substance not found! Returing NA. \n")
      id <- NA
      attr(id, "matched") <- NA
      attr(id, "distance") <- NA
      return(id)
    }
    search <- search[search$match$type %in% c("label", "alias"), 
                     ]
    search <- search[tolower(iconv(search$match$text, "latin1", 
                                   "ASCII", sub = "")) == tolower(URLdecode(query)), 
                     ]
    if (nrow(search) > 1) {
      if (verbose) 
        message("More then one Link found. \n")
      if (match == "na") {
        if (verbose) 
          message("Returning NA. \n")
        id <- NA
        matched_sub <- NA
        d <- NA
      }
      if (match == "all") {
        if (verbose) 
          message("Returning all matches. \n")
        id <- search$id
        matched_sub <- search$label
        d <- "all"
      }
      if (match == "first") {
        if (verbose) 
          message("Returning first match. \n")
        id <- search$id[1]
        matched_sub <- search$label[1]
        d <- "first"
      }
      if (match == "best") {
        if (verbose) 
          message("Returning best match. \n")
        dd <- adist(URLdecode(query), search$label)/nchar(search$label)
        id <- search$id[which.min(dd)]
        d <- round(dd[which.min(dd)], 2)
        matched_sub <- search$label[which.min(dd)]
      }
      if (match == "ask") {
        tochoose <- data.frame(match = search$label, 
                               url = search$url, description = search$description)
        print(tochoose)
        message("\nEnter rownumber of compounds (other inputs will return 'NA'):\n")
        take <- as.numeric(scan(n = 1, quiet = TRUE))
        if (length(take) == 0) {
          id <- NA
          matched_sub <- NA
          d <- NA
        }
        if (take %in% seq_len(nrow(tochoose))) {
          id <- search$id[take]
          matched_sub <- search$label[take]
          d <- "interactive"
        }
      }
    }
    else {
      id <- search$id
      d <- 0
      matched_sub <- search$label
    }
    names(id) <- NULL
    attr(id, "matched") <- matched_sub
    attr(id, "distance") <- d
    return(id)
  }
  if (length(query)>0){   #added to the source code
    out <- lapply(query, foo, language = language, match = match, 
                  verbose = verbose)
    if (match != "all") {
      out <- data.frame(t(sapply(out, function(y) {
        c(y, attr(y, "matched"), attr(y, "distance"))
      })), stringsAsFactors = FALSE)
      # names(out) <- c("id", "match", "distance") #commented in the main code
      out[["query"]] <- query
    }
    return(out)
  }    # added to the main code from webchem package
}

toxdb<-read.csv('inputData/tox2db.csv') # Toxdb:  γ-linolenate biosynthesis pathway
comps<-unique(toxdb$drug)
ids<-rep(NA,length(comps))
for(i in 1:length(comps)){  #getting the id for each compound and saving them in ids 
  ids[[i]]<-newgetid(as.character(comps[i]), match = 'best')[[1]]
  print(comps[i])            
}
ids<-na.omit(ids)
ids<-ids[!ids==0]
res<-wd_ident(ids,verbose = F)
toxdb<-merge(res,toxdb)
save(toxdb,file='outputData/ids_toxdb.RData')
save(toxdb,file='inputData/ToxDB/annotation_toxdb_dictionary.RData')




# 2. Intergating with toxdb -----------------------------------------------


load("outputData/class_probilities/integrated_fgsea_allcompounds_results_final_moa.RData");fgs<-fgs[c(1:7,11:15)] #our predicted pathway scores RWR-FGSEA
tox_db2_class_pro<-function(pathw,toxdb_pat){# pathw is the annottaion for the pathway, toxdb_pat is the file retireved from ToxDB
  source('functions/annotation_functions.R') # mesh2cas function 
  
  load('inputData/ToxDB/annotation_toxdb_dictionary.RData')# chembl2cas dictionary 
  dict<-unique(toxdb[,c('chembl','cas')])
  
  toxdb_pat$combinedID<-paste(toxdb_pat$study,toxdb_pat$celltype,toxdb_pat$time,sep = '_')
  colnames(toxdb_pat)[colnames(toxdb_pat)=='score']<-'activation_score_ToxDB'
  layers<- c( "DrugMatrix_hepatocyte_1 d","DrugMatrix_liver_1 d", "DrugMatrix_liver_3 d","DrugMatrix_liver_5 d",
              rep("TG-GATEs_inVivoRatLiverSingleDose_1 d",3),"TG-GATEs_inVitroHumanLiver_1 d","TG-GATEs_inVitroRatLiver_1 d",
              "TG-GATEs_inVivoRatLiverRepeatDose_15 d","TG-GATEs_inVivoRatLiverRepeatDose_29 d","TG-GATEs_inVivoRatLiverRepeatDose_8 d") 
  
  integrated_data<-list()
  for (i in 1:length(fgs)){
    print(names(fgs)[i])
    layer<-layers[i]
    temp_tox<-toxdb_pat[toxdb_pat$combinedID==layer,]
    temp_tox<-merge(temp_tox,dict,all = F)
    temp_our<-as.data.frame(fgs[[i]]$x[,pathw])
    colnames(temp_our)<-'activation_score_RWR_FGSEA'
    temp_our$mesh<-rownames(temp_our)
    
    temp_our$cas<-mesh2cas(temp_our$mesh)
    res<-merge(temp_our,temp_tox,all = F)
    res$nam<-names(fgs)[i]
    integrated_data[[i]]<-res
    
  }
  all_layers<-do.call(rbind,integrated_data) #combining 12 layers
  all_layers$pathway<-pathw
  return(all_layers)
}

# wiki estrogen receptor
estr_name<-'Estrogen Receptor Pathway%WikiPathways_20190510%WP2881%Homo sapiens'
estr_toxdb<-read.csv2('inputData/ToxDB/ToxDB_wiki_Estrogen_Receptor_Pathway.csv',stringsAsFactors = F)
estr_integrated<-tox_db2_class_pro(estr_name,estr_toxdb)

# wiki aryl hydrocarbon receptor
aryl_toxdb<-read.csv2('inputData/ToxDB/ToxDB_wiki_Aryl_Hydrocarbon_Receptor.csv',stringsAsFactors = F)
aryl_name<-'Aryl Hydrocarbon Receptor%WikiPathways_20190510%WP2586%Homo sapiens'
aryl_integrated<-tox_db2_class_pro(aryl_name,aryl_toxdb)

# wiki integrated breast cancer
breast_toxdb<-read.csv2('inputData/ToxDB/ToxDB_wiki_Integrated_Breast_Cancer_Pathway.csv',stringsAsFactors = F)
breast_name<-'Integrated Breast Cancer Pathway%WikiPathways_20190510%WP1984%Homo sapiens'
breast_integrated<-tox_db2_class_pro(breast_name,breast_toxdb)

# wiki Tryptophan metabolism
tryptophan_toxdb<-read.csv2('inputData/ToxDB/ToxDB_wiki_Tryptophan_metabolism.csv',stringsAsFactors = F)
tryptophan_name<-'Tryptophan metabolism%WikiPathways_20190510%WP465%Homo sapiens'
tryptophan_integrated<-tox_db2_class_pro(tryptophan_name,tryptophan_toxdb)

final_res<-do.call(rbind,list(aryl_integrated,breast_integrated,estr_integrated,tryptophan_integrated))
write.csv(final_res,file = 'outputData/excel_files/toxdb_patway_sores.csv')


# 3. plot analysis -------------------------------------------------------------

pat_scores<-read.csv('outputData/excel_files/toxdb_patway_sores.csv',stringsAsFactors = F)
pat_scores$is_invitro<-'in vivo'
pat_scores$is_invitro[grep(pat_scores$nam,pattern = 'vitro')]<-'in vitro'
pat_scores$study_type<-'DrugMatrix'
pat_scores$study_type[grep(pat_scores$nam,pattern = 'TG_')]<-'TG-GATES_Single'
pat_scores$study_type[grep(pat_scores$nam,pattern = 'TG_GATEs_Rat_invivo_Repeated')]<-'TG-GATES_Repeated'

pat_scores<-pat_scores[pat_scores$activation_score_RWR_FGSEA>0 & pat_scores$activation_score_ToxDB>0 ,]
pat_scores$pathway<-gsub(pat_scores$pathway,pattern = '%WikiPathways_20190510%WP2586%Homo sapiens',replacement = '')
pat_scores$pathway<-gsub(pat_scores$pathway,pattern = '%WikiPathways_20190510%WP1984%Homo sapiens',replacement = '')
pat_scores$pathway<-gsub(pat_scores$pathway,pattern = '%WikiPathways_20190510%WP2881%Homo sapiens',replacement = '')
pat_scores$pathway<-gsub(pat_scores$pathway,pattern = '%WikiPathways_20190510%WP465%Homo sapiens',replacement = '')
pat_scores$time<-gsub(pat_scores$time,pattern = ' d',replacement = '')
pat_scores$time<-as.numeric(as.character(pat_scores$time))
load('outputData/new199edc_1462dec.RData')
pat_scores$status<-''
pat_scores$status[which(pat_scores$mesh %in% names(new_edc199_decoy_1462)[1:199])]<-'EDC'

#pat_scores<-pat_scores[pat_scores$status=='EDC',]
library(ggrepel);library(ggplot2);require(RColorBrewer)
#lvls<-unique(pat_scores$combinedID)[c(1,7,6,2:4,5,10,8,9)]
lvls<-unique(pat_scores$nam)[c(1:4,8,9,6,7,5,12,10,11)]
color_values<-c('Red','#CCE6FF','#B3D9FF','#99CCFF',                     # Drug matrix hep,1,3,5
                'Orange','Pink',                                         # in vitro rat and human from TGgates
                '#D9ECC6','#CCE6B3','#BFDF9F','#B3D98C','#A6D279','#99CC66' # in vivo  rat low,high,moddle,8,15,29 TGgates
)                                               
ggplot(pat_scores, aes(x=activation_score_ToxDB, y=activation_score_RWR_FGSEA,colour=factor(nam,levels = lvls)),size=2) +
  #geom_point(size=ifelse(pat_scores$status=='edc',0,0.5))+
  geom_point()+
  labs(color='Data Layer')+
  #scale_size_manual(name='Time (Day)',values = seq(2,3.5,by=0.25)) +
  #scale_shape_manual(name='Time (Day)',values=c(17,25,18,24,16,15))+
  geom_text_repel(label=pat_scores$status,size=3,segment.size = 0.1, segment.color = "grey") +
  facet_wrap(pathway ~ study, nrow  = 2)+
  ggtitle('')+xlab('Pathway activation Scores ToxDB')+ylab('Pathway activation Scores RWR-FGSEA')+
  scale_color_manual(values = color_values)+
  theme_minimal()+ theme(axis.title.x = element_text(size = 20),axis.title = element_text(size = 20),axis.text.x = element_text(size=15),
                         axis.text.y = element_text(size=15),strip.text = element_text(size = 20),legend.position = 'bottom',
                         legend.text = element_text(size = 15))
# 

# 
# 
save(pat_scores,lvls,color_values,file = 'outputData/plots_tables_functions/ToxDB/patways_scores.RData')




# 4.  ---------------------------------------------------------------------
pat_scores<-read.csv('outputData/excel_files/toxdb_patway_sores.csv',stringsAsFactors = F)
#pat_scores<-pat_scores[!pat_scores$activation_score_RWR_FGSEA==0,]
pat_scores$pathway<-gsub(pat_scores$pathway,pattern = '%WikiPathways_20190510%WP2586%Homo sapiens',replacement = '')
pat_scores$pathway<-gsub(pat_scores$pathway,pattern = '%WikiPathways_20190510%WP1984%Homo sapiens',replacement = '')
pat_scores$pathway<-gsub(pat_scores$pathway,pattern = '%WikiPathways_20190510%WP2881%Homo sapiens',replacement = '')
pat_scores$pathway<-gsub(pat_scores$pathway,pattern = '%WikiPathways_20190510%WP465%Homo sapiens',replacement = '')
pat_scores$time<-gsub(pat_scores$time,pattern = ' d',replacement = '')
pat_scores$time<-as.numeric(as.character(pat_scores$time))
cmp<-list()
layers_name<-unique(tst$nam)
for (i in 1:length(layers_name)){
  print(layers_name[i])
  cmp[[i]]<-unique(tst$chembl[tst$nam==layers_name[i]])
}
Reduce('intersect',cmp)



##### OLD wrong approach for ToxDB

# 1.harmonic edc_score VS toxdb scores -------------------------------------------------------------

# 
# source('functions/annotation_functions.R')
# load("outputData/ids_toxdb.RData")     #output of the above code
# load("outputData/VAM_scores_hs_most_informative_layers_moa.RData") 
# vam_ctd<-as.data.frame(all_VAM_score)
# vam_ctd$cas<-mesh2cas(names(all_VAM_score))
# colnames(vam_ctd)<-c('VAM_score','cas')
# all_scores<-merge(vam_ctd,toxdb)
# library(ggplot2)
# library(ggrepel)
# 
# all_scores$dose=gsub(x = all_scores$dose,pattern = '[?]',replacement = '')
# all_scores = all_scores[-which(all_scores$celltype == "heart" | all_scores$celltype == "kidney" | all_scores$celltype == "muscle"), 
#           c("drug","VAM_score","score","study","celltype","time","dose")]
# all_scores = all_scores[grep("d",all_scores$time),]
# all_scores = all_scores[-grep("4 d",all_scores$time),];all_scores = all_scores[-grep("7 d",all_scores$time),]
# all_scores$time = droplevels(all_scores$time)
# all_scores$drug = paste(all_scores$drug, all_scores$study, all_scores$celltype, all_scores$time, all_scores$dose,sep="_")
# all_scores$time = as.numeric(gsub(all_scores$time,pattern=" d",replacement=""))
# write.csv(all_scores,file = 'outputData/excel_files/toxdb_harmonic_vam_moa.csv')
# vam_gp = all_scores; vam_gp$drug[which(vam_gp$VAM_score < 0.75 | vam_gp$score < 3)] = ""
# vam_gp$lbl<-rep("grey50",nrow(vam_gp))
# vam_gp$lbl[which(vam_gp$VAM_score > 0.75 & vam_gp$score > 3)] = "blue"
# names(vam_gp)[3]<-'TOXDB_score'
# vam_gp$drug=gsub(x = vam_gp$drug,pattern = '_TG-GATEs_|_DrugMatrix_',replacement = '_')
# vam_gp$drug=gsub(x = vam_gp$drug,pattern = '_TG-GATEs_|_DrugMatrix_',replacement = '_')
# vam_gp$drug=gsub(x = vam_gp$drug,pattern = 'inVivo',replacement = 'VV_')
# vam_gp$drug=gsub(x = vam_gp$drug,pattern = 'inVitro',replacement = 'VT_')
# vam_gp$drug=gsub(x = vam_gp$drug,pattern = 'hepatocyte',replacement = 'HP')
# vam_gp$drug=gsub(x = vam_gp$drug,pattern = 'RepeatDose',replacement = 'RD')
# vam_gp$drug=gsub(x = vam_gp$drug,pattern = 'Rat',replacement = 'Ra')
# vam_gp$drug=gsub(x = vam_gp$drug,pattern = 'Liver',replacement = 'Li')
# ggplot(vam_gp, aes(x=VAM_score, y=TOXDB_score)) +
#   geom_point(aes(shape=factor(time),size=factor(time)), color=vam_gp$lbl,alpha = 0.7) + 
#   scale_size_manual(name='Time (Day)',values = seq(2,3.5,by=0.25)) +
#   scale_shape_manual(name='Time (Day)',values=c(17,25,18,24,16,15))+
#   geom_text_repel(aes(label = drug), size = 2, segment.size = 0.2, segment.color = "grey") +
#   facet_wrap( ~study, nrow  = 2)+
#   ggtitle('Harmonic Sum EDC score VS TOXDB scores')+xlab('Harmonic EDC score')+ylab('TOXDB score')+
#   theme_minimal()
# 
# 
# 
# # 2.average EDC_scores VS toxdb scores---------------------------------------------------------------
# 
# 
# 
# source('functions/annotation_functions.R')
# load("outputData/ids_toxdb.RData")     #output of the script 1_5
# load("outputData/VAM_scores_average_most_informative_layers_moa.RData") 
# vam_ctd<-as.data.frame(all_VAM_score)
# vam_ctd$cas<-mesh2cas(names(all_VAM_score))
# colnames(vam_ctd)<-c('VAM_score','cas')
# all_scores<-merge(vam_ctd,toxdb)
# library(ggplot2)
# library(ggrepel)
# 
# all_scores$dose=gsub(x = all_scores$dose,pattern = '[?]',replacement = '')
# all_scores = all_scores[-which(all_scores$celltype == "heart" | all_scores$celltype == "kidney" | all_scores$celltype == "muscle"), 
#                         c("drug","VAM_score","score","study","celltype","time","dose")]
# all_scores = all_scores[grep("d",all_scores$time),]
# all_scores = all_scores[-grep("4 d",all_scores$time),];all_scores = all_scores[-grep("7 d",all_scores$time),]
# all_scores$time = droplevels(all_scores$time)
# all_scores$drug = paste(all_scores$drug, all_scores$study, all_scores$celltype, all_scores$time, all_scores$dose,sep="_")
# all_scores$time = as.numeric(gsub(all_scores$time,pattern=" d",replacement=""))
# write.csv(all_scores,file = 'outputData/excel_files/toxdb_average_vam_moa.csv')
# vam_gp = all_scores; vam_gp$drug[which(vam_gp$VAM_score < 0.75 | vam_gp$score < 3.5)] = ""
# vam_gp$lbl<-rep("grey50",nrow(vam_gp))
# vam_gp$lbl[which(vam_gp$VAM_score > 0.75 & vam_gp$score > 3.5)] = "blue"
# names(vam_gp)[3]<-'TOXDB_score'
# vam_gp$drug=gsub(x = vam_gp$drug,pattern = '_TG-GATEs_|_DrugMatrix_',replacement = '_')
# vam_gp$drug=gsub(x = vam_gp$drug,pattern = '_TG-GATEs_|_DrugMatrix_',replacement = '_')
# #vam_gp$drug=gsub(x = vam_gp$drug,pattern = 'inVivo',replacement = '')
# #vam_gp$drug=gsub(x = vam_gp$drug,pattern = 'inVitro',replacement = '')
# #vam_gp$drug=gsub(x = vam_gp$drug,pattern = 'hepatocyte',replacement = '')
# vam_gp$drug=gsub(x = vam_gp$drug,pattern = 'RepeatDose',replacement = '')
# vam_gp$drug=gsub(x = vam_gp$drug,pattern = 'Rat',replacement = '')
# #vam_gp$drug=gsub(x = vam_gp$drug,pattern = 'Liver',replacement = '')
# ggplot(vam_gp, aes(x=VAM_score, y=TOXDB_score),size=1) +
#   #geom_point(aes(shape=factor(time)), color=vam_gp$lbl,alpha = 0.7) + 
#   geom_point(aes(color=factor(time)))+
#   labs(color='Time (day)')+
#   scale_size_manual(name='Time (Day)',values = seq(2,3.5,by=0.25)) +
#   scale_shape_manual(name='Time (Day)',values=c(17,25,18,24,16,15))+
#   geom_text_repel(aes(label = drug), size = 4, segment.size = 0.2, segment.color = "grey") +
#   facet_wrap( ~study, nrow  = 2)+
#   ggtitle('Average EDC score VS TOXDB scores')+xlab('Average EDC score')+ylab('TOXDB score')+
#   theme_minimal()+ theme(axis.title.x = element_text(size = 20),axis.title = element_text(size = 20),
#                          strip.text = element_text(size = 20))
# 
