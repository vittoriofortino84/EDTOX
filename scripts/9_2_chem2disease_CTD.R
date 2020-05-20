
# Binary matrix of chemicals and disease ----------------------------------
#setwd("/research/groups/fortino/amirhs/EDCmet_phase2")
load("input/new199edc_1462dec.RData")
ixns<-read.csv('input/CTD_chemicals_diseases.csv.gz',comment.char = c("#"),header = FALSE,stringsAsFactors = F) 
ixns<-ixns[,c(2,4)]
ixns<-unique(ixns)
colnames(ixns)<-c('chemical','disease')
ixns$disease<-gsub('[()]+','',ixns$disease)    # Omitting paranthese characters in  disease column to make uniform annoations
ixns$disease<-gsub('[,]+','',ixns$disease)
ixns$disease<-tolower(ixns$disease)            # Lower case converting of all words in the disease column
chemicals<-unique(ixns$chemical)               # all chemicala
all_diseases<-unique(ixns$disease)             # all diseaes types
chem_disease<-matrix(0, nrow = length(chemicals), ncol = length(all_diseases),dimnames = list(chemicals,all_diseases))  #binary matrix chemicals disease types
for (i in 1:length(chemicals)){
  print(i)
  chem_disease[i,which(all_diseases %in% ixns[ixns$chemical %in% chemicals[i],'disease'])]<-1
}
chem_disease_edc_dec<-chem_disease[rownames(chem_disease) %in% names(new_edc199_decoy_1462),]
save(chem_disease_edc_dec,chem_disease,file='output/chem2disease.RData')

# Ratio for edcs and decoys used in classifier of phase 1-----------------------------------------------

load('output/chem2disease.RData')
load("input/new199edc_1462dec.RData")
disease_sum_edc<-disease_sum_dec<-rep(NA,ncol(chem_disease_edc_dec))
for (j in 1:ncol(chem_disease_edc_dec)){
  disease_sum_edc[j]<-sum(chem_disease_edc_dec[rownames(chem_disease_edc_dec) %in% names(new_edc199_decoy_1462)[1:199],j])   #for edcs
  disease_sum_dec[j]<-sum(chem_disease_edc_dec[rownames(chem_disease_edc_dec) %in% names(new_edc199_decoy_1462)[200:1661],j])#for decoys
}

print(paste(length(which(rownames(chem_disease_edc_dec) %in% names(new_edc199_decoy_1462)[1:199])),'edcs'))
print(paste(length(which(rownames(chem_disease_edc_dec) %in% names(new_edc199_decoy_1462)[200:1661])),'decoys'))


names(disease_sum_edc)<-names(disease_sum_edc)<-colnames(chem_disease_edc_dec)
barplot(disease_sum_edc,las=2)
lines(disease_sum_dec,col='red')
ratio<-disease_sum_edc/disease_sum_dec
#ratio<-ratio[!is.na(ratio)]
sum_edc<-as.data.frame(do.call(cbind,list(disease_sum_edc,disease_sum_dec,ratio)))
names(sum_edc)<-c('sum_edc_disease','sum_dec_disease','ratio_disease')
xlsx::write.xlsx2(sum_edc,file = 'output/sum_edcs_decs_diseases.xlsx')



# validation of predicted compounds from phase 1 based on disease approach -------------------------
load('output/chem2disease.RData')
source('functions/general_functions.R')
actives<-vam2comp(probmin = 0.85,probmax=1,type='harmonic')
inactives<-vam2comp(probmin = 0,probmax = 0.4,type = 'harmonic')
disease_sum_edc<-disease_sum_dec<-rep(NA,ncol(chem_disease_edc_dec))
for (j in 1:ncol(chem_disease_edc_dec)){
  disease_sum_edc[j]<-sum(chem_disease[rownames(chem_disease) %in% actives,j])   #for compounds EDC_score > 0.85
  disease_sum_dec[j]<-sum(chem_disease[rownames(chem_disease) %in% inactives,j]) #for compounds EDC_score < 0.4
}
names(disease_sum_edc)<-names(disease_sum_edc)<-colnames(chem_disease_edc_dec)
barplot(disease_sum_edc,las=2)
lines(disease_sum_dec,col='red')
ratio<-disease_sum_edc/disease_sum_dec
sum_edc<-as.data.frame(do.call(cbind,list(disease_sum_edc,disease_sum_dec,ratio)))
names(sum_edc)<-c('sum_edc_disease','sum_dec_disease','ratio_disease')
xlsx::write.xlsx2(sum_edc,file = 'output/sum_edcs_decs_diseases_vam_0.85_0.4.xlsx')


