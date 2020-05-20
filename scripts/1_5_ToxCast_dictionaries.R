# preparation of TOXCAST assay2target_gene and Toxcast hitcall
# 1. Assay-target_gene
toxcast<-read.csv('inputData/Assay_Summary_190226.csv',header = TRUE, stringsAsFactors = F)      # TOXCAST 3.1 assay_summary
toxcast<-subset(toxcast,select=c(8,39,53))
library(tidyr)
library(dplyr)
toxcast<- toxcast %>%
  mutate(intended_target_entrez_gene_id=strsplit(as.character(intended_target_entrez_gene_id),'|',fixed = TRUE)) %>% 
  unnest(intended_target_entrez_gene_id)
toxcast<-na.omit(toxcast)
toxcast<-distinct(toxcast)
save(toxcast,file = 'outputData/ToxCast/assay_target_gene.RData')

# 2. Hitcall of ToxCast 3.1
toxcast_3_1<-read.csv('inputData/hitc_Matrix_190226.csv',header = TRUE, stringsAsFactors = F)      # TOXCAST 3.1 assay_summary
library(readxl)
all_cas<-read_excel('inputData/DSSTox_Identifiers_and_CASRN.xlsx') # all compounds tox_ids
all_cas<-all_cas$casrn
all_cas2<-gsub('[-]','',all_cas) #transforming tox ids into CAS number according to page 6 of: 
all_cas2<-gsub('^','C',all_cas2) #https://lri.americanchemistry.com/Users-Guide-for-Accessing-and-Interpreting-ToxCast-Data.pdf
toxcast_3_1<-toxcast_3_1[which(toxcast_3_1$X %in% all_cas2),]
for (i in 1:dim(toxcast_3_1)[1]){
  print(i)
  nam<-toxcast_3_1$X[i]
  rownames(toxcast_3_1)[i]<-all_cas[which(all_cas2==nam)]
}
toxcast_3_1[,1]<-c()
save(toxcast_3_1,file='outputData/ToxCast/hitcall.RData')

