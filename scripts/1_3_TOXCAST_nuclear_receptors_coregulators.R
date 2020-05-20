# 3.    Retrieving Nuclear_receptor_ENDPOINTS from TOXCAST --------------------------------------
# Subseting nuclear receptors and their coregulators from the toxcast hitc binary file based on two data sources
# one from NRs and coregs (NURSA) and one from our experts
# The genes entrez ids related to each TOXCAST endpoint was retrieved from TOXCAST Assay_Summary
# input :  inputData/Assay_Summary_190226.csv' from ToxCast 3.1
#          https://nursa.org/nursa/molecules/index.jsf

toxcast<-read.csv('inputData/Assay_Summary_190226.csv',header = TRUE, stringsAsFactors = F)      # TOXCAST 3.1 assay_summary
TOXBASED_ASSAYS<-toxcast$assay_component_endpoint_name[which(toxcast$intended_target_family=='nuclear receptor')] #TOXcast data of nuclear receptor 1
toxcast<-subset(toxcast,select=c(8,39,53))
library(tidyr)
library(dplyr)
toxcast<- toxcast %>% mutate(intended_target_entrez_gene_id=strsplit(as.character(intended_target_entrez_gene_id),'|',fixed = TRUE)) %>% unnest(intended_target_entrez_gene_id)
toxcast<-na.omit(toxcast)
toxcast<-distinct(toxcast)                                                                     
mie_assay_names<-unique(toxcast$assay_component_endpoint_name)
nr<-read.csv('inputData/nr_and_coreg.csv',header = TRUE,stringsAsFactors = F) # https://nursa.org/nursa/molecules/index.jsf Data of nuclear receptors and coregulators
nr_genes<-unique(nr$entrez)
# Howmany MIES in  CTD are Coregulators
coreg_genes<-unique(nr$entrez[nr$Type=="Coregulator"])
load('outputData/chem2gene_no_out.RData')
all_possible_mies<-unique(unlist(lapply(chem2gene, function(x)x)))
print(paste(length(intersect(all_possible_mies,coreg_genes)),'genes are Co-regulators from',length(all_possible_mies),'MIES in CTD'))
expert_list<-read.csv('inputData/NR_EDCs.csv',header = TRUE,stringsAsFactors = F)  # Data of nuclear receptors from our experts 
library('org.Hs.eg.db')
exp_gids<-unique(expert_list$Gene)
exp_gids<-exp_gids[-which(exp_gids=='')]
exp_gids<-mapIds(org.Hs.eg.db, exp_gids, 'ENTREZID', 'SYMBOL')
nr_genes<-union(nr_genes,exp_gids)
nr_genes<-as.character(na.omit(nr_genes))
toxcast_nr<-toxcast[which(as.character(toxcast$intended_target_entrez_gene_id) %in% as.character(sort(nr_genes))),]
nr_assays_names<-unique(toxcast_nr$assay_component_endpoint_name)
nr_assays_names<-union(nr_assays_names,TOXBASED_ASSAYS)
toxcast_3_1<-read.csv('inputData/hitc_Matrix_190226.csv',header = TRUE, stringsAsFactors = F)      # TOXCAST 3.1 assay_summary
library(readxl)
all_cas<-read_excel('inputData/DSSTox_Identifiers_and_CASRN.xlsx') # all compounds tox_ids and cas dictionary
all_cas<-all_cas$casrn
all_cas2<-gsub('[-]','',all_cas) #transforming tox ids into CAS number according to page 6 of: 
all_cas2<-gsub('^','C',all_cas2) #https://lri.americanchemistry.com/Users-Guide-for-Accessing-and-Interpreting-ToxCast-Data.pdf
# manual merging takes a bit time
toxcast_3_1<-toxcast_3_1[which(toxcast_3_1$X %in% all_cas2),]
for (i in 1:dim(toxcast_3_1)[1]){
  print(i)
  nam<-toxcast_3_1$X[i]
  rownames(toxcast_3_1)[i]<-all_cas[which(all_cas2==nam)]
}
hitc_mat_NR_Coreg<-toxcast_3_1[,which(colnames(toxcast_3_1) %in% nr_assays_names)]
hitc_mat_MIE<-toxcast_3_1[,which(colnames(toxcast_3_1) %in% mie_assay_names)]
save(toxcast_3_1,hitc_mat_NR_Coreg,hitc_mat_MIE,file='outputData/toxcast_direct_based_NR_coreg_from_TOXCAST_binary.RData')# Toxcast chemicals,NR endpoints
ixns<-read.csv('inputData/CTD_chem_gene_ixns.csv.gz',comment.char = c("#"),stringsAsFactors = F) #all compounds from CTD
library(data.table)
ixns<-as.data.table(ixns)
ixns<-ixns %>% distinct(ixns$CasRN,ixns$ChemicalID)
colnames(ixns)<-c('cas','chem')
hitc_mat_NR_Coreg<-as.data.frame(hitc_mat_NR_Coreg)
hitc_mat_NR_Coreg$cas<-rownames(hitc_mat_NR_Coreg)
chemid_hic_NR<-merge(as.data.frame(ixns),hitc_mat_NR_Coreg)
chemid_hic_NR<-as.matrix(chemid_hic_NR)
rownames(chemid_hic_NR)<-chemid_hic_NR[,2]
save(chemid_hic_NR,file='outputData/hitc_chemid_NR_coreg.RData') #all compounds in CTD which were tested for nuclear assay endpoints in ToxCast
write.csv(chemid_hic_NR,file = 'outputData/excel_files/nuclear_recptors_coregulators.csv')
