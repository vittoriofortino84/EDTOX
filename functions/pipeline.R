
# RWR-FGSEA pipeline_DM_TG for drug matrix and tg-gates networks ----------------------------------------------------------


# x is the network, 
# ptw is pathways
# perc is the percentage of the edges from the network to be used and 
# ng is the number of sorted genes after random walk to be used for fgsea
# com_seed is the seeds or MIEs for the compounds (list of Entrez genes one set of genes per compound)
# p is parralellization  option
pipeline<-function(x,ptw,perc,ng,comp_seeds,p){
  
  
  library(fgsea)
  library(dnet)
  library(igraph)
  
  #graph optimization
  x$Node.1<-as.character(x$Node.1)               #pairwise calculation so length(node1)=length(node2)=length(wTO)=length(pval)=length(pval.adjusted)      
  x$Node.2<-as.character(x$Node.2)
  x$wTO[x$pval.adj>5*10^-3]<-0                   #non-significant weighted topological overlap values with p-adjusted values of more than 0.005 are taken to zero 
  x$wTO<-abs(x$wTO)                              #absolute signs for wto values are taken
  gr<-graph_from_data_frame(x[order(x$wTO,decreasing = T),1:2][1:round(nrow(x)*perc),],directed=F)  #gettting the   percent of weighted edge values
  d_sign<-lapply(comp_seeds,function(x)intersect(x,V(gr)$name))    #intersection of the genes in the graph and the genes related to the compunds
  d_sign<-d_sign[sapply(d_sign,length)>0]                     #omitting chemicals with no intersected genes inside the graph 
  mat<-sapply(d_sign,function(x)(V(gr)$name%in%x)*1)           #matrix of the genes in the row and chemicals in the columns
  rownames(mat)<-V(gr)$name 
  
  #random walk with restrart
  probm<-dRWR(gr,setSeeds=mat,normalise.affinity.matrix='quantile',parallel=T,verbose = F) 
  probm<-as.matrix(probm)
  colnames(probm)<-names(d_sign)                        #Columns are chemicals
  rownames(probm)<-V(gr)$name
  for (i in 1:dim(probm)[2]){
    cutof_q<-sort(probm[,i],decreasing = T)[ng]
    probm[probm[,i]<cutof_q,i]<-0
  }
  
  
  #Fgsea
  q<-apply(probm,2,function(y)sum(y!=0)) # Counts the number of genes for each chemical which are at least once seen in random walk (probility >0)
  probm<-probm[,q>=ng]                     # We take the chemicals with at least more than ng observed genes after random walk on the GCN
  probm<-as.matrix(probm)
  colnames(probm)<-names(q)[which(q>=ng)]
  newscore<-ES<-NES<-pval<-matrix(NA,nrow=length(ptw),ncol=ncol(probm),dimnames=list(names(ptw),colnames(probm))) #making a matrix with null values  the rows will be the pathways and the columns will be the compounds
 
   for(i in 1:ncol(probm)){
    #for each chemical and its sets of genes we do this procedure
    res<-fgsea(pathways=ptw,stats=probm[probm[,i]!=0,i],nperm=10000,minSize=1,maxSize=200,BPPARAM=p)
    #the pathways more than 500 genes are not considered and those with at least one gene will be considered
    # in the case of nuclear receptors there are sometimes pathways with only one gene
    ES[res$pathway,i]<-res$ES
    NES[res$pathway,i]<-res$NES
    pval[res$pathway,i]<-res$pval
    newscore[res$pathway,i]<-(res$NES)*(-log10(res$pval))
    if(i%%100==0) gc()           # at the end of the loop clear the memory
  }
  
  return(list(ES=ES,NES=NES,pval=pval,newscore=newscore))
}# end function




# pipeline_PPI for PPI network-----------------------------------------------------------

pipeline_ppi<-function(ptw,combscore,ng,comp_seeds,p){
  
  library(fgsea)
  library(dnet)
  library(igraph)
 
  load('outputData/network/ppi_network.RData')                   # PPI nwtwork
  d_sign<-lapply(comp_seeds,function(x)intersect(x,V(gr)$name))  #intersection of the genes in the graph and the genes related to the compunds
  d_sign<-d_sign[sapply(d_sign,length)>0]                        #omitting chemicals with no intersected genes inside the graph 
  mat<-sapply(d_sign,function(x)(V(gr)$name%in%x)*1)             #matrix of the genes in the row and chemicals in the columns
  rownames(mat)<-V(gr)$name 
  
  #random walk
  probm<-dRWR(gr,setSeeds=mat,normalise.affinity.matrix='quantile',parallel=T,verbose = F) 
  probm<-as.matrix(probm)
  colnames(probm)<-names(d_sign)                        #Columns are chemicals
  rownames(probm)<-V(gr)$name 
  for (i in 1:dim(probm)[2]){
    cutof_q<-sort(probm[,i],decreasing = T)[ng]
    probm[probm[,i]<cutof_q,i]<-0
  }
  
  
  #Fgsea
  q<-apply(probm,2,function(y)sum(y!=0)) # Counts the number of genes for each chemical which are at least once seen in random walk (probility >0)
  probm<-probm[,q>=ng]                     # We take the chemicals with at least more than 2000 observed genes after random walk on the GCN
  probm<-as.matrix(probm)
  colnames(probm)<-names(q)[which(q>=ng)]
  newscore<-ES<-NES<-pval<-matrix(NA,nrow=length(ptw),ncol=ncol(probm),dimnames=list(names(ptw),colnames(probm))) #making a matrix with null values  the rows will be the pathways and the columns will be the compounds
  for(i in 1:ncol(probm)){
    #for each chemical and its sets of genes we do this procedure
    res<-fgsea(pathways=ptw,stats=probm[probm[,i]!=0,i],nperm=10000,minSize=1,maxSize=200,BPPARAM=p)
    #the pathways more than 500 genes are not considered and those with at least one gene will be considered
    # in the case of nuclear receptors there are sometimes pathways with only one gene
    ES[res$pathway,i]<-res$ES
    NES[res$pathway,i]<-res$NES
    pval[res$pathway,i]<-res$pval
    newscore[res$pathway,i]<-(res$NES)*(-log10(res$pval))
    if(i%%100==0) gc()           # at the end of the loop clear the memory
  }
  
  return(list(ES=ES,NES=NES,pval=pval,newscore=newscore))
}# end function



# pipeline_consensus for LINCS AND DM-TG CONSENSUS NETWORK------------------------------------------------------


# x is the consensus network, 
#ptw is pathways
# perc is the percentage of the edges from the network to be used and 
# ng is the number of genes in the network to be used
#com_seed is the seeds for the compounds
#p is parralellization  option
pipeline_consensus<-function(x,ptw,perc,ng,comp_seeds,p){
  
  
  library(fgsea)
  library(dnet)
  library(igraph)
  
  #graph optimization
  x$Node.1<-as.character(x$Node.1)               #pairwise calculation so length(node1)=length(node2)=length(wTO)=length(pval)=length(pval.adjusted)      
  x$Node.2<-as.character(x$Node.2)
  x$CN[x$pval.fisher>5*10^-3]<-0                   #non-significant weighted topological overlap values with p-adjusted values of more than 0.005 are taken to zero 
  x$CN<-abs(x$CN)                              #absolute signs for wto values are taken
  gr<-graph_from_data_frame(x[order(x$CN,decreasing = T),1:2][1:round(nrow(x)*perc),],directed=F)  #gettting the   percent of weighted edge values
  d_sign<-lapply(comp_seeds,function(x)intersect(x,V(gr)$name))    #intersection of the genes in the graph and the genes related to the compunds
  d_sign<-d_sign[sapply(d_sign,length)>0]                     #omitting chemicals with no intersected genes inside the graph 
  mat<-sapply(d_sign,function(x)(V(gr)$name%in%x)*1)           #matrix of the genes in the row and chemicals in the columns
  rownames(mat)<-V(gr)$name 
  
  #random walk with restrart
  probm<-dRWR(gr,setSeeds=mat,normalise.affinity.matrix='quantile',parallel=T,verbose = F) 
  probm<-as.matrix(probm)
  colnames(probm)<-names(d_sign)                        #Columns are chemicals
  rownames(probm)<-V(gr)$name
  for (i in 1:dim(probm)[2]){
    cutof_q<-sort(probm[,i],decreasing = T)[ng]
    probm[probm[,i]<cutof_q,i]<-0
  }
  
  
  #Fgsea
  q<-apply(probm,2,function(y)sum(y!=0)) # Counts the number of genes for each chemical which are at least once seen in random walk (probility >0)
  probm<-probm[,q>=ng]                     # We take the chemicals with at least more than 2000 observed genes after random walk on the GCN
  probm<-as.matrix(probm)
  colnames(probm)<-names(q)[which(q>=ng)]
  newscore<-ES<-NES<-pval<-matrix(NA,nrow=length(ptw),ncol=ncol(probm),dimnames=list(names(ptw),colnames(probm))) #making a matrix with null values  the rows will be the pathways and the columns will be the compounds
  
  for(i in 1:ncol(probm)){
    #for each chemical and its sets of genes we do this procedure
    res<-fgsea(pathways=ptw,stats=probm[probm[,i]!=0,i],nperm=10000,minSize=1,maxSize=200,BPPARAM=p)
    #the pathways more than 500 genes are not considered and those with at least one gene will be considered
    # in the case of nuclear receptors there are sometimes pathways with only one gene
    ES[res$pathway,i]<-res$ES
    NES[res$pathway,i]<-res$NES
    pval[res$pathway,i]<-res$pval
    newscore[res$pathway,i]<-(res$NES)*(-log10(res$pval))
    #if(i%%100==0) gc()           # at the end of the loop clear the memory
  }
  
  return(list(ES=ES,NES=NES,pval=pval,newscore=newscore))
}# end function
