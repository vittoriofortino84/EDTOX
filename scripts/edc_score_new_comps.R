# edc score for new compounds from their MIEs
source('functions/vam_functions.R')
newcomps<-list(comp1=c('196','2101' ,'3065'),comp2=c('196','2101')) #entrez ids of genes as MIEs for each compound 
edc_score=mie2vam(newcomps)