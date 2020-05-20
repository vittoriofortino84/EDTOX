
# 4. wTO Consensus -----------------------------------------------------------

load('outputData/network/wTO_DM.RData')
load('outputData/network/wTO_TGG.RData')
wTO_DM_hep<-wTO_DM$DM_hepatocytes.1
wTO_TGG_hep<-wTO_TGG$TGG_rat_in_vitro
wTO_hep<-list(wTO_DM_hep=wTO_DM_hep,wTO_TGG_hep=wTO_TGG_hep)
get_final_wTO<-function(x){
  x$Node.1<-as.character(x$Node.1)
  x$Node.2<-as.character(x$Node.2)
  return(x)
}
wTO_hep<-lapply(wTO_hep,get_final_wTO)
library(wTO)
wTO_hep_cons<-wTO.Consensus(data=list(wTO_Data1=data.frame(Node.1=wTO_hep[[1]]$Node.1,
                                                           Node.2=wTO_hep[[1]]$Node.2,
                                                           wTO=wTO_hep[[1]]$wTO,
                                                           pval=wTO_hep[[1]]$pval),
                                      wTO_Data2=data.frame(Node.1=wTO_hep[[2]]$Node.1,
                                                           Node.2=wTO_hep[[2]]$Node.2,
                                                           wTO=wTO_hep[[2]]$wTO,
                                                           pval=wTO_hep[[2]]$pval)))
save(wTO_hep_cons,file='outputData/network/wTO_hep_cons.RData')



