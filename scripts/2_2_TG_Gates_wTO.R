#  wTO TG-GATEs ------------------------------------------------------------
require(xlsx)

file_nam<-read.xlsx2('inputData/journal.pcbi.1004847.s026.XLS',1,stringsAsFactors=F) #For assigning micro-array probeids to genes
get_TGG_final_mtxs<-function(x,pheno,taxon,in_vitro,prd,split_by, liver_expr=file_nam){
  pheno<-pheno[!pheno$DOSE_LEVEL=='Control',]                                 # The controls are removed becaue they are log2fold ratios
  pheno<-pheno[pheno$SACRIFICE_PERIOD%in%prd,]
  if(split_by=='dose')pheno<-split(x=pheno,f=as.factor(pheno$DOSE_LEVEL))
  else if (split_by=='period')pheno<-split(x=pheno,f=as.factor(pheno$SACRIFICE_PERIOD))
  in_9071<-liver_expr$in_9071_rat_liver_expressed_set=='Yes'
  if(taxon=='human' & in_vitro==T) probe_id<-liver_expr$HPH.HGU133.2          #probe set HPH.HGU133.2 : human in vitro
  else if(taxon=='rat' & in_vitro==T)probe_id<-liver_expr$RPH.RG230.2         #probe set RPH.RG230.2  : rat in vitro
  else probe_id<-liver_expr$RAT.LIVER.RG230.2                                 #probe set LIVER.RG230.2: rat in vivo
  has_orth<-liver_expr$HUMAN_ORTHOLOG_ENTREZ!=''                              
  logi<-has_orth&in_9071&(probe_id!='')
  probe_set<-probe_id[logi]
  human_orth<-liver_expr$HUMAN_ORTHOLOG_ENTREZ[logi]
  rownames(x)<-x[,1]
  if (class(pheno)=='list'){
    x<-sapply(pheno,function(y){
      mtx<-x[probe_set,colnames(x)%in%paste('00',y$BARCODE,sep='')];    #performing the data of annotations on the main matrix
      rownames(mtx)<-human_orth;
      mtx})
    return(x)
  }
  else{
    x<-x[probe_set,colnames(x)%in%paste('00',pheno$BARCODE,sep='')]
    rownames(x)<-human_orth
    return(list(x))
  }
}
path_l2r<-paste(getwd(),'/inputData/log2ratio/TG/',sep='')
files<-list.files(path_l2r)
TGG<-lapply(files,function(x)read.csv(paste(path_l2r,x,sep=''),sep='\t',stringsAsFactors=F,check.names=F))
names(TGG)<-c('TGG_human_in_vitro','TGG_rat_in_vitro','TGG_rat_in_vivo_rep','TGG_rat_in_vivo_single')
path_cell<-paste(getwd(),'/inputData/log2ratio/cell/',sep='')    #TG-gate annotations results of normalized microarray data
files<-list.files(path_cell)
TGG_cell<-lapply(files,function(x)read.csv(paste(path_cell,x,sep=''),stringsAsFactors=F,check.names=F))
names(TGG_cell)<-c('TGG_human_in_vitro','TGG_rat_in_vitro','TGG_rat_in_vivo_rep','TGG_rat_in_vivo_single')
taxon<-c('human',rep('rat',3))                           # The first network for human and the three last for rat
in_vitro<-c(T,T,F,F)                                     #Two net works are in vitro and the two last are in vivo
period<-list('24 hr','24 hr',c('8 day','15 day','29 day'),'24 hr')  

split_by<-list('NULL','NULL','period','dose')           #No splitting in the first and second splitting based on period and dose for the last two networks 
TGG_final_mtx<-mapply(get_TGG_final_mtxs,TGG,TGG_cell,taxon,in_vitro,period,split_by)
# x=TGG list of csv log2 ratios
# pheno=TGG_cell csv files (annotations for tg-gate log2ratios)
#taxon =taxon 'human' 'rat' 'rat' 'rat'
#invitro =invitro T,T,F,F
#period= prd '24 hr','24 hr',c('8 day','15 day','29 day'),'24 hr'
#split_by= 'NULL','NULL','period','dose')
TGG_final_mtx<-unlist(TGG_final_mtx,recursive = F)
library(wTO)
library(doParallel)
cl<-makeCluster(8)
registerDoParallel(cl)
wTO_TGG<-foreach(i=1:length(TGG_final_mtx),.export = 'wTO.fast') %dopar% wTO.fast(Data=TGG_final_mtx[[i]],
                                                                                  method = 'p', sign = 'sign',method_resampling = "Bootstrap",delta = 0.2)
stopCluster(cl)
names(wTO_TGG)<-names(TGG_final_mtx)
save(wTO_TGG,file='outputData/network/wTO_TGG.RData')
