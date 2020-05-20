
source("functions/general_functions.R") #functions  jaccard dissimilarity
library(STRINGdb)
library(igraph)
library(org.Hs.eg.db)
library(dnet)
library(cluster)
load("outputData/new199edc_1462dec.RData")
new_edc199_decoy_1462<-as.array(new_edc199_decoy_1462)


perclist<-c(0.6,0.65,0.7,0.75,0.8,0.85)        #which percent of the edges to extract for randomwalking
nglist<-c(200,500,700,1000)      #The number of most top visited genes to consider for jaccard dissimilarity


string_db<-STRINGdb$new(version="10",species=9606)
entz_id<-data.frame(entz_id=keys(org.Hs.egGENENAME),stringsAsFactors = F)
string_ids<-string_db$map(my_data_frame=entz_id,my_data_frame_id_col_names='entz_id',removeUnmappedRows=T)
string_inter<-string_db$get_interactions(string_ids$STRING_id)
string_inter<-string_inter[,-grep('coexpression',names(string_inter))]
sel_inter<-1-(string_inter[,3:(ncol(string_inter)-1)]/1000)
string_inter$combined_score<-1-Reduce('*',sel_inter)

silm<-matrix(NA, nrow = length(perclist), ncol = length(nglist))
rownames(silm)<-perclist
colnames(silm)<-nglist
ind_per<-0


for (perc in perclist ){
  ind_ng<-0
  ind_per<-ind_per+1
  
  edge_list<-string_inter[string_inter$combined_score>=perc,1:2]
  edge_list<-merge(edge_list,string_ids,by.x='from',by.y='STRING_id')
  edge_list<-merge(edge_list,string_ids,by.x='to',by.y='STRING_id')[,3:4]
  colnames(edge_list)<-c('Node.1','Node.2')
  gr<-graph_from_edgelist(as.matrix(edge_list),directed=F)
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
    print(paste(perc,'_',ng))
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

    #silhouette score
    ndecoy<-length(which(rownames(new_edc199_decoy_1462)[200:1661] %in% rownames(pro)))
    nedc<-length(which(rownames(new_edc199_decoy_1462)[1:199] %in% rownames(pro)))
    resp<-c(rep(1,nedc),rep(0,ndecoy)) 
    km <- kmeans(jac, centers = 2, nstart=25)
    km$cluster<-resp
    ss<-silhouette(km$cluster, jac)
    silm[ind_per,ind_ng]<-mean(ss[1:nedc,3])
  } # end for ng
  
} # end for percent
save(silm,file =paste("outputData/tuning/ppi/Sil_PPI.RData",sep = ''))


