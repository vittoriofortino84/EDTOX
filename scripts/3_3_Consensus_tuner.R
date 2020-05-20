
source("functions/general_functions.R") #functions  jaccard dissimilarity
load("outputData/new199edc_1462dec.RData")
load("outputData/network/wTO_hep_cons.RData")
load("outputData/network/wTO_LINCS_cons.RData")
all_cons<-list(wTO_hep_cons,wTO_LINCS_cons)
names(all_cons)<-c('hep_cons','lincs_cons') 
new_edc199_decoy_1462<-as.array(new_edc199_decoy_1462)
library(dnet)
library(igraph)
library(cluster)
perclist<-c(0.02,0.03,0.05,0.1)        #which percent of the edges to extract for randomwalking
nglist<-c(200,500,700,1000)            #The number of most top visited genes to consider for jaccard dissimilarity


#network preparation and edge cutoffs
for (l in 1:length(all_cons)){
  x<-all_cons[[l]]
  x$Node.1<-as.character(x$Node.1)               #pairwise calculation so length(node1)=length(node2)=length(wTO)=length(pval)=length(pval.adjusted)      
  x$Node.2<-as.character(x$Node.2)
  x$CN[x$pval.fisher>5*10^-3]<-0
  x$CN<-abs(x$CN)
  silm<-matrix(NA, nrow = length(perclist), ncol = length(nglist))
  rownames(silm)<-perclist
  colnames(silm)<-nglist
  ind_per<-0
  for (perc in perclist ){
    ind_ng<-0
    ind_per<-ind_per+1
    gr<-graph_from_data_frame(x[order(x$CN,decreasing = T),1:2][1:round(nrow(x)*perc),],directed=F)  #gettting the   percent of weighted edge values
    d_sign<-lapply(new_edc199_decoy_1462,function(x)intersect(x,V(gr)$name))    #intersection of the genes in the graph and the genes related to the compunds
    d_sign<-d_sign[sapply(d_sign,length)>0]                     #omitting chemicals with no intersected genes inside the graph 
    mat<-sapply(d_sign,function(x)(V(gr)$name%in%x)*1)           #matrix of the genes in the row and chemicals in the columns
    rownames(mat)<-V(gr)$name 
    #random walk
    probm<-dRWR(gr,setSeeds=mat,normalise.affinity.matrix='quantile',parallel=T,verbose = F) 
    probm<-as.matrix(probm)
    colnames(probm)<-names(d_sign)                        #Columns are chemicals
    rownames(probm)<-V(gr)$name 
    for (ng in nglist){
      ind_ng<-ind_ng+1
      print(paste(l,'-',perc,'_',ng))
      pro<-probm
      # cutoff based on  top visited genes
      for (i in 1:dim(pro)[2]){
        cutof_q<-sort(pro[,i],decreasing = T)[ng]
        pro[pro[,i]<cutof_q,i]<-0
        pro[pro[,i]>=cutof_q,i]<-1
      }
      #jaccard and heatmap matrix
      pro<-t(pro)
      jac<-jaccrd_dissimilarity(pro)
      #silhouette plot
      ndecoy<-length(which(names(new_edc199_decoy_1462)[200:1661] %in% rownames(pro)))
      nedc<-length(which(names(new_edc199_decoy_1462)[1:199] %in% rownames(pro)))
      resp<-c(rep(1,nedc),rep(0,ndecoy)) 
      km <- kmeans(jac, centers = 2, nstart=25)
      km$cluster<-resp
      ss<-silhouette(km$cluster, jac)
      silm[ind_per,ind_ng]<-mean(ss[1:nedc,3])
    } # end for ng
  } # end for percent
  save(silm,file =paste("outputData/tuning/Sil_all_cons_",names(all_cons)[l],'.RData',sep = ''))
  #assign(paste('Sil_all_cons',names(all_cons)[l]),silm)
} #end for l 
