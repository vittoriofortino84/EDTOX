
# 3. wTO LINCS ---------------------------------------------------------------
library(cmapR)

ds_path<-c(paste(getwd(),"/inputData/log2ratio/LINCS/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx",sep=''),
           paste(getwd(),  "/inputData/log2ratio/LINCS/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx",sep=''))

sign_met_path<-c(paste(getwd(),'/inputData/log2ratio/LINCS/GSE70138/GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz',sep=''),
                 paste(getwd(),  "/inputData/log2ratio/LINCS/GSE92742/GSE92742_Broad_LINCS_sig_metrics.txt.gz",sep=''))
sign_met<-lapply(sign_met_path,function(x)read.csv(x,sep='\t',stringsAsFactors=F))
sign_met<-lapply(sign_met,function(x)x[grep('HEPG2',x$sig_id),])            # Selecting HEPG2
sign_met<-lapply(sign_met,function(x)x[grep('24H',x$sig_id),])              # selecting 24 Hrs
sign_met<-lapply(sign_met,function(x)x[!grepl("ctl_vehicle",x$pert_type),]) # Removing control samples
my_ds<-mapply(function(x,i)parse.gctx(i,cid=x$sig_id),sign_met,ds_path,SIMPLIFY = F)
liver_expr<-read.csv('inputData/journal.pcbi.1004847.s026.csv',stringsAsFactors = F)
in_9071<-liver_expr$in_9071_rat_liver_expressed_set=='Yes'
has_orth<-!is.na(liver_expr$HUMAN_ORTHOLOG_ENTREZ)
human_orth<-as.character(liver_expr$HUMAN_ORTHOLOG_ENTREZ[in_9071&has_orth])
my_ds<-lapply(my_ds,function(x){y=x@mat;y[rownames(y)%in%human_orth,]})
library(doParallel)
library(wTO)
cl<-makeCluster(2)
registerDoParallel(cl)
wTO_LINCS<-foreach(i=1:length(my_ds),.export = 'wTO.fast') %dopar% wTO.fast(Data=as.data.frame(my_ds[[i]]),
                                                                            method = 'p', sign = 'sign',method_resampling = "Bootstrap",delta = 0.2)
stopCluster(cl)
# consensus network
get_final_wTO<-function(x){
  
  x$Node.1<-as.character(x$Node.1)
  x$Node.2<-as.character(x$Node.2)
  return(x)
}
wTO_LINCS<-lapply(wTO_LINCS,get_final_wTO)
wTO_LINCS_cons<-wTO.Consensus(data=list(wTO_Data1=data.frame(Node.1=wTO_LINCS[[1]]$Node.1,
                                                             Node.2=wTO_LINCS[[1]]$Node.2,
                                                             wTO=wTO_LINCS[[1]]$wTO,
                                                             pval=wTO_LINCS[[1]]$pval),
                                        wTO_Data2=data.frame(Node.1=wTO_LINCS[[2]]$Node.1,
                                                             Node.2=wTO_LINCS[[2]]$Node.2,
                                                             wTO=wTO_LINCS[[2]]$wTO,
                                                             pval=wTO_LINCS[[2]]$pval)))

save(wTO_LINCS,file = 'outputData/network/wTO_LINCS_cons.RData')

