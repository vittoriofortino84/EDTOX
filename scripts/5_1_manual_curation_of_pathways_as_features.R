#Preparation of the excel file and pretreatment of the pathways 
# load("outputData/class_probilities/integrated_fgsea_allcompounds_results_final.RData")
# pathways_list<-Reduce('union',lapply(fgs, function(x)colnames(x$x)))
# source('functions/annotation_functions.R')
# intg<-as.data.frame(list(pathways_list,annotated_pathways))
# colnames(intg)<-c('pathway','annotation')
### The fle was manually curated in excels and those related to diseases were signed as 0 and those realted to MOA as 1 then
# 
load('outputData/toxdb2gene_final.RData')
moa_pathways<-xlsx::read.xlsx('outputData/excel_files/pathways_annotations_manually_curated.xlsx',sheetIndex = 1)
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='GO_',ignore.case = F) ]=0         # GO_terms 
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='cancer',ignore.case = T) ]=1      # Evidence for the cancer for EDCs
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='metabolism',ignore.case = T) ]=1  # Emphasis on metabolism
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='aop:',ignore.case = T) ]=1        # Emphasis on AOPs related
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='development',ignore.case = T)]=0
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='behavior',ignore.case = T) ]=0
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='disease',ignore.case = T) ]=0
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='contraction',ignore.case = T) ]=0
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='differentiation',ignore.case = T) ]=0
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='infect',ignore.case = T) ]=0
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='bacter',ignore.case = T) ]=0
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='atrophy',ignore.case = T) ]=0
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='malaria',ignore.case = T) ]=0
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='viral',ignore.case = T) ]=0
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='hypertrophy',ignore.case = T) ]=0
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='morphogenesis',ignore.case = T) ]=0
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='apopt',ignore.case = T) ]=0
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='ribosom',ignore.case = T) ]=0
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='polymerase',ignore.case = T) ]=0
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='radiation',ignore.case = T) ]=0
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='heat_shock',ignore.case = T) ]=0
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern ='_rna',ignore.case = T) ]=0
moa_pathways$MOA0_disease1[grep(x=moa_pathways$annotation,pattern =' rna',ignore.case = T) ]=0

# Considering only those pathways with moa=1 
moa_pathways<-moa_pathways$pathway[moa_pathways$MOA0_disease1==1]
moa_pathways<-as.character(moa_pathways)
all_pathways<-c(go2gene,go2gene_noDesc,kegg2gene,msigdb2gene,reactome2gene,wikiptw2gene)
moa_pathways<-all_pathways[moa_pathways]

# calculation of jaccard index between the pathways to removoe duplicated pathways based on similrity of genes
source('functions/general_functions.R')
jacard_sim_moa<-list_jaccrd_similarity(moa_pathways) #pairwise jaccard similarity
save(jacard_sim_moa,file='outputData/jaccrad_between_moa_based_pathways.RData')
#

load('outputData/jaccrad_between_moa_based_pathways.RData')
diag(jacard_sim_moa)<-0
jacard_sim_moa[lower.tri(jacard_sim_moa)]<-0
melted_jac<-reshape::melt(jacard_sim_moa)
duplicates<-melted_jac[which(melted_jac$value==1),]
colnames(duplicates)<-c('first_pathways','second_pathways','jaccard_sim')
source('functions/annotation_functions.R')
duplicates$first_annotations<-pathway_annotate(as.character(duplicates$first_pathways))
duplicates$second_annotations<-pathway_annotate(as.character(duplicates$second_pathways))
duplicates$first_category<-sapply(strsplit(duplicates$first_annotations,'_'),'[',1)
duplicates$second_category<-sapply(strsplit(duplicates$second_annotations,'_'),'[',1)
to_remove_pathways<-duplicates$second_pathways
for (i in 1:length(to_remove_pathways)){
  if (duplicates$first_category[i]=='REACTOME' & duplicates$second_category[i]=='KEGG'){
    to_remove_pathways[i]=duplicates$first_pathways[i]
  }
}
to_remove_pathways<-unique(to_remove_pathways)
to_remove_pathways<-as.character(to_remove_pathways)
moa_pathways[to_remove_pathways]<-c() # Removing pathways with completely similar genes
#selection of pathways based on genes expressed in the liver 
pr<-xlsx::read.xlsx2('inputData/journal.pcbi.1004847.s026.XLS',1,stringsAsFactors=F)   #list of probe ids and their corresponding gene ids in microarray
human_orthologs_expressed_in_rat_liver<-unique(pr$HUMAN_ORTHOLOG_ENTREZ[pr$in_9071_rat_liver_expressed_set=='Yes'])
percents<-sapply(moa_pathways, function(x)length(intersect(x,human_orthologs_expressed_in_rat_liver))/length(x))
percent_tbl<-as.data.frame(percents)
percent_tbl$annotation<-pathway_annotate(rownames(percent_tbl))
to_remove_pathways<-rownames(percent_tbl)[which(percent_tbl$percents==0)]
moa_pathways[to_remove_pathways]<-c()
save(moa_pathways,toxdb_names,go2gene,go2gene_noDesc,msigdb2gene,reactome2gene,kegg2gene,wikiptw2gene,duplicates,
     file = 'outputData/toxdb2gene_final.RData')

