# 1.    Extracting MIEs from CTD for all chemicals to use as seeds in RWR --------------
#
# Getting the compounds and their related genes to use as seeds in the random-walk procedure
# input:  http://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz  in inputData/CTD_chem_gene_ixns.csv.gz
# output: list of chemical and their related MIEs (genes) outputData/chem2gene_no_out.RData 
#
# compound gene interarctions release of june 2020
ixns<-read.csv('inputData/CTD_june_2020/CTD_chem_gene_ixns.csv.gz',comment.char = c("#"),
               stringsAsFactors = F) 
colnames(ixns)<-c("ChemicalName","ChemicalID", "CasRN","GeneSymbol", "GeneID","GeneForms",
                  "Organism","OrganismID","Interaction" ,"InteractionActions", "PubMedIDs")
library(data.table)
ixns<-as.data.table(ixns)
ixns<-ixns[,paste(InteractionActions,collapse = '|'),by=list(chemID=ixns$ChemicalID,cas=ixns$CasRN,geneID=ixns$GeneID,geneSymbol=ixns$GeneSymbol)]
InterAction<-lapply(ixns$V1,function(x)strsplit(x,'|',fixed = T))
InterAction<-lapply(InterAction,function(x)sapply(x,strsplit,split='^',fixed=T))
type<-sort(unique(unlist(sapply(InterAction,function(x)sapply(x,function(y)y[2])))))
#making a matrix (m*n) with m=1044988 the number of gene-compound and n=51 interaction types 
InterAction_mat<-matrix(FALSE,nrow=nrow(ixns),ncol=length(type),dimnames=list(NULL,type))
#Interaction_mat: a logical matrix with rows: chemical-gene interatction and columns as Interaction types
for(i in 1:length(InterAction)){
  for(j in 1:length(InterAction[[i]])){
    InterAction_mat[i,InterAction[[i]][[j]][2]]<-TRUE
  }
}
# METABOLIC INTERACTIONS:
# "acetylation"            "acylation"              "ADP-ribosylation"       "alkylation"            
# "amination"              "carbamoylation"         "carboxylation"          "chemical synthesis"    
# "cleavage"               "degradation"            "ethylation"             "farnesylation"         
# "geranoylation"          "glucuronidation"        "glutathionylation"      "glycation"             
# "glycosylation"          "hydrolysis"             "hydroxylation"          "lipidation"            
# "methylation"            "N-linked glycosylation" "nitrosation"            "O-linked glycosylation"
# "oxidation"              "palmitoylation"         "phosphorylation"        "prenylation"           
# "reduction"              "ribosylation"           "sulfation"              "sumoylation"           
# "ubiquitination"
idx_met_proc<-c(2,4:7,9:12,14,15,18,20:26,28,31,33:39,41,43,47,48,50)
#TRANSPORT INTERACTIONS:
# "export"    "import"    "secretion" "uptake" 
idx_transport<-c(16,27,44,51)
InterAction_mat[,"metabolic processing"]<-InterAction_mat[,"metabolic processing"]|
  apply(InterAction_mat[,idx_met_proc],1,any) #IF THE CHEMICAL-GENE PERTURBATION WAS 1 IN ANY OF THE METABOLIC INTERAACTIONS GETS 1
InterAction_mat[,"transport"]<-InterAction_mat[,"transport"]|apply(InterAction_mat[,idx_transport],1,any) ##IF THE CHEMICAL-GENE PERTURBATION WAS 1 IN ANY OF THE TRANSPORT GETS 1
InterAction_mat<-InterAction_mat[,-c(idx_met_proc,idx_transport)]  #SIMPLIFYING THE columns of the LOGICAL MATRIX into 14 interaction types
#Performing multiple correspondence analysis to get the variables which are more informative
### The plot is saved in plots folder
#Performing multiple correspondence analysis to get the variables which are more informative
### The plot is saved in plots folder
library(FactoMineR)
library(factoextra)
res.mca<-MCA(InterAction_mat,ncp = 14,graph = F)
fviz_eig(res.mca,ncp = 14,addlabels = T)
fviz_mca_var(res.mca,choice = 'var',axes = 1:2,repel = T,shape.var = 19,labelsize=8)+theme(axis.title = element_text(size=30))+labs(title='')
#reaction','binding','activity','expression','metabolic processing' were the more distant types
#

type.count<-function(x){
  if (class(x)=='logical') x*1
  else colSums(x)
}
sel_type<-c('reaction','binding','activity','expression','metabolic processing')
idx<-apply(InterAction_mat[,sel_type,drop=F],1,any)
ixns<-ixns[idx,]
InterAction_mat <- InterAction_mat[idx,sel_type]
##################################################### NEW CODE
# compile for each chemical a data.frame indicating for each gene the data sources supporting the interaction
provolone <- list()
selt <- c('R','B', 'A', 'E', 'M')
for(c in unique(ixns$chemID)){
  test = ixns[ixns$chemID==c,]
  mat = InterAction_mat[ixns$chemID==c,]
  if(nrow(test)==1) {
    provolone[[length(provolone)+1]] = test$geneID
  }
  else {
    gid = unique(test$geneID)
    sum_gid = c()
    cat_gid = c()
    for(g in gid) {
      #print(which(g == test$geneID))
      sum_gid = c(sum_gid, sum(mat[which(g == test$geneID),]))
      cat_gid = c(cat_gid,  paste(selt[which(mat[which(g == test$geneID),] == TRUE)], collapse="."))
    }
    provolone[[length(provolone)+1]] = data.frame(gid, sum_gid, cat_gid, row.names = 1)
  }
}

# compile for each chemical the chemical-gene interactions that need to be removed
mies_to_remove = lapply(provolone, function(x) {
  if(!is.numeric(x)) {
    if(nrow(x) > 50) {
      id.EM = which(x$sum_gid == 1 & (x$cat_gid == 'E' | x$cat_gid == 'M'))
      id.EM = c(id.EM, which(x$sum_gid == 2 & x$cat_gid == 'E.M'))
      if((nrow(x) - length(id.EM)) > 200) {
        id.REM = which(x$sum_gid == 1 & x$cat_gid == 'R')
        id.REM = c(id.REM, which(x$sum_gid == 2 & x$cat_gid == 'R.E'))
        id.REM = c(id.REM, which(x$sum_gid == 2 & x$cat_gid == 'R.M'))
        id.REM = c(id.REM, which(x$sum_gid == 3 & x$cat_gid == 'R.E.M'))
        if((nrow(x) - (length(id.EM) + length(id.REM))) > 200) {
          id.AB = which(x$sum_gid == 1 & (x$cat_gid == 'A' | x$cat_gid == 'B'))
          if((nrow(x) - (length(id.EM) + length(id.REM) + length(id.AB))) > 300) {
            tab = table(x$cat_gid[-c(id.EM, id.REM, id.AB)])
            print(tab[which(tab > 0)])
            print(length(x$cat_gid[-c(id.EM, id.REM, id.AB)]))
            print("-------------------------------------------------------------")
            #return(-1)
            return(c(rownames(x)[id.EM], rownames(x)[id.REM], rownames(x)[id.AB]))
          }
          return(c(rownames(x)[id.EM], rownames(x)[id.REM], rownames(x)[id.AB]))
        }
        return(c(rownames(x)[id.EM], rownames(x)[id.REM]))
      }
      return(rownames(x)[id.EM])
    }
    else NULL
  }
  else NULL
})
# determine the final set of entrez gene ids associated to each chemical
final_set_mies = list()
for(i in 1:length(provolone)) {
  if(is.data.frame(provolone[[i]])) {
    if(nrow(provolone[[i]]) > 50 & length(mies_to_remove[[i]]) > 0) {
      final_set_mies[[i]] = setdiff(rownames(provolone[[i]]), mies_to_remove[[i]])
    }
    else  final_set_mies[[i]] = rownames(provolone[[i]])
  }
  else
    final_set_mies[[i]] = provolone[[i]]
}
names(final_set_mies) = names(mies_to_remove) = names(provolone) = unique(ixns$chemID)
final_set_mies<-final_set_mies[-which(sapply(final_set_mies,length)>100)]
final_set_mies<-final_set_mies[-which(sapply(final_set_mies, length)==0)]
chem2gene<-final_set_mies
save(list = c('final_set_mies', 'mies_to_remove', 'provolone'), file = "outputData/new_MIES_set.RData")
save(chem2gene,file = 'outputData/chem2gene_no_out.RData')

                                                              
                                                              
                                                              
##makes a dictionary with three columns cas name mesh
 ixns<-read.csv('inputData/CTD_chemicals.csv.gz',comment.char = c("#"),stringsAsFactors = F,header = F) #all compounds from CTD
 ixns$V2<-gsub(pattern = 'MESH:',replacement = '',ixns$V2)
 library(data.table)
 ixns<-as.data.table(ixns)
 ixns<-ixns[,c(1,2,3)]
 colnames(ixns)<-c('name','mesh','cas')
 chem_gene<-read.csv('inputData/CTD_chem_gene_ixns.csv.gz',comment.char = c("#"),stringsAsFactors = F) #around 2000000 gene chmical interaction with 11 variables
 chem_gene<-unique(chem_gene[,c("ChemicalName","ChemicalID","CasRN")])
 colnames(chem_gene)<-c('name','mesh','cas')
 ixns<-merge(ixns,chem_gene,by = c('name','mesh','cas'),all=T)
 ixns<-unique(ixns)
 write.csv(ixns,file = 'inputData/annotaion/chem_mesh_cas_dictionary.csv')
