
#  wTO Drug Matrix ---------------------------------------------------------

path_l2r<-paste(getwd(),'/inputData/log2ratio/DM/',sep='')
files<-list.files(path_l2r)
DM<-lapply(files,function(x)read.csv(paste(path_l2r,x,sep=''),sep='\t',stringsAsFactors=F,check.names=F))
names(DM)<-c('DM_hepatocytes','DM_liver')
DM_annotation<-read.csv('inputData/DrugMatrix Exp annotations.txt',sep='\t',stringsAsFactors=F)
require(xlsx)
pr<-read.xlsx2('inputData/journal.pcbi.1004847.s026.XLS',1,stringsAsFactors=F)   #list of probe ids and their corresponding gene ids in microarray

#in function x= DMs two lists from csv file, pheno= DM annotations and in_vitro is a logical vector T,F
get_DM_final_mtxs<-function(x,pheno,in_vitro,liver_expr=pr){
  rownames(x)<-x$Name
  pheno<-pheno[-grep('control',pheno$Source.Name.s_),]                 # removes rows containing "control" in source name because data are log2ratio
  all_durations<-unique(pheno$Dose.Duration..C.)                       # i.e: 3 1 5 4 0.25 7 14  duration times for DM
  all_durations<-all_durations[all_durations>=1 & all_durations<10]    # Durations between 1 and 10 are selected :i.e 1 3  4 5 7
  pheno<-pheno[pheno$Dose.Duration..C.%in%all_durations,]              # elimination of annotaions with i.e 0.25 and 14                 
  pheno<-split(pheno,f=as.factor(pheno$Dose.Duration..C.))             # transfroming reamined dose durations as different lists
  pheno<-pheno[sapply(pheno, nrow)>100]                                # i.e 4 and 7 were removed because nrows<100;to prevent wrong correlations 
  in_9071<-liver_expr$in_9071_rat_liver_expressed_set=='Yes'
  if(in_vitro==T) probe_id<-liver_expr$RPH.RG230.2                     # Assigning probeids to related gene for in vitro and in vivo
  else probe_id<-liver_expr$RAT.LIVER.RG230.2                          # First item in DM list is human hepatocytes (in vitro) and the second item is Rat liver (in vivo)
  has_orth<-liver_expr$HUMAN_ORTHOLOG_ENTREZ!=''
  logi<-has_orth&in_9071&(probe_id!='')
  probe_set<-probe_id[logi]
  human_orth<-liver_expr$HUMAN_ORTHOLOG_ENTREZ[logi]
  #x= DMs two lists from csv file, pheno= modified DM annotations and in_vitro is a logical vector T,F
  x<-lapply(pheno,function(y){
    mtx<-x[probe_set,colnames(x)%in%y[,1]];                            #Row names are gene ids
    rownames(mtx)<-human_orth;
    mtx})
}

in_vitro<-c(T,F)         #the first item in DM list is human hepatocytes (in vitro) and the second item is Rat liver (in vivo)
DM_final_mtx<-mapply(function(x,y)get_DM_final_mtxs(x,DM_annotation,y),DM,in_vitro,SIMPLIFY = F) #The final list with rows as geneids and columns as drugs
DM_final_mtx<-unlist(DM_final_mtx,recursive = F)
DM_final_mtx<-DM_final_mtx[lapply(DM_final_mtx,ncol)>100]                      
set.seed(123)
library(doParallel)               #For paralellization 
library(wTO)                      # For doing weghted gene co-expression with topology overlap (wTO) 
cl<-makeCluster(4)                #use to parallel 4 nodes because there are 4 lists inside the DM (including DM_hepatocyte, DM_liver.1, DM_liver.3, DM_liver.5)
registerDoParallel(cl)            #used to register the parallel backend with the foreach package using 4 nodes
wTO_DM<-foreach(i=1:length(DM_final_mtx),.export = 'wTO.fast') %dopar% wTO.fast(Data=DM_final_mtx[[i]],
                                                                                method = 'p', sign = 'sign',method_resampling = "Bootstrap",delta = 0.2)
stopCluster(cl)
names(wTO_DM)<-names(DM_final_mtx)
save(wTO_DM,file='outputData/network/wTO_DM.RData')

