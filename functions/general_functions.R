
#FUNCTION jaccard_dissimilarity returns  n*n matrix of pairwise jaccard dissimilarity from n*m binary matrix 
#         the columns must have names
jaccrd_dissimilarity<-function(pro){ #row wise calculation of jaccard index
  jac<-matrix(NA, nrow = dim(pro)[1], ncol = dim(pro)[1])
  nu<-dim(pro)[1]
  nu0<-nu-1
  for (i in 1:nu0){
    er<-i+1
    for (j in er:nu){
      g1<-colnames(pro)[which(pro[i,] %in% 1)]
      g2<-colnames(pro)[which(pro[j,] %in% 1)]
      tmp<-1-(length(intersect(g1,g2))/length(union(g1,g2)))
      jac[i,j]<-tmp
      jac[j,i]<-tmp
    } #end for j
  }#end for i
  diag(jac)<-0
  rownames(jac)<-colnames(jac)<-rownames(pro)
  return(jac)
}#end function   

#retunrs n*n from a list with gene sets
list_jaccrd_similarity<-function(pro){ #row wise calculation of jaccard index
  jac<-matrix(NA, nrow = length(pro), ncol = length(pro))
  nu<-length(pro)
  nu0<-nu-1
  for (i in 1:nu0){
    er<-i+1
    for (j in er:nu){
      tmp<-length(intersect(pro[[i]],pro[[j]]))/length(union(pro[[i]],pro[[j]]))
      jac[i,j]<-tmp
      jac[j,i]<-tmp
    } #end for j
  }#end for i
  diag(jac)<-1
  rownames(jac)<-colnames(jac)<-names(pro)
  return(jac)
}#end function

# Column wise zero and one scaling of data
zero_one_scale<-function(x){
  datt<-x
  for (d in 1:ncol(datt)){
    dn<-na.omit(datt[,d])
    if (length(dn)>0){
      maxx<-max(dn)
      minn<-min(dn)
      delta<-maxx-minn
      datt[,d]<-(datt[,d]-minn)/delta
    } # end if
  } #end for
  return(datt)
} # end function




#Corrects and uniforms layers 
name_correct<-function(x){
  x<-gsub('.RData','',x)
  x<-gsub('ppi','PPI_STRINGdb',x)
  x<-gsub('DM','Drug Matrix_',x)
  x<-gsub('TGG|Fgsea_tg','TG_Gate',x)
  x<-gsub('_hep','Hepatocytes',x)
  x<-gsub('hep_cons','Consensus_Rat_invitro_DM_TG_gates',x)
  x<-gsub('TGG|Fgsea_tg','TG_Gate',x)
  x<-gsub('lincs_cons','Consensus_LINCS',x)
  x<-gsub('Consensus_LINCS','Consensus_LINCS_HEPG2',x)
  x<-gsub('hep_cons','Consensus_Hepatocytes(DM,TG-gates)',x)
  x<-gsub("Drug Matrix_Hepatocytes",'Drug_Matrix_Rat_invitro_Single_Dose_1_day',x)
  x<-gsub('Drug Matrix_1','Drug_Matrix_Rat_invivo_Single_Dose_1_day',x)
  x<-gsub('Drug Matrix_3','Drug_Matrix_Rat_invivo_Repeated_Dose_3_days',x)
  x<-gsub('Drug Matrix_5','Drug_Matrix_Rat_invivo_Repeated_Dose_5_days',x)
  x<-gsub('TG_Gate_single.high','TG_GATEs_Rat_invivo_Single_Dose_High_1_day',x)
  x<-gsub('TG_Gate_single.low','TG_GATEs_Rat_invivo_Single_Dose_Low_1_day',x)
  x<-gsub('TG_Gate_single.middle','TG_GATEs_Rat_invivo_Single_Dose_Middle_1_day',x)
  x<-gsub('Consensus_Hepatocytes','Consensus_Rat_invitro',x)
  x<-gsub('TG_Gate_human_in_vitro','TG_GATEs_Human_invitro_Single_Dose_1_day',x)
  x<-gsub('TG_Gate_rat_in_vitro','TG_GATEs_Rat_invitro_Single_Dose_1_day',x)
  x<-gsub('TG_Gate_rat_in_vivo_rep_15_day','TG_GATEs_Rat_invivo_Repeated_Dose_15_days',x)
  x<-gsub('TG_Gate_rat_in_vivo_rep_29_day','TG_GATEs_Rat_invivo_Repeated_Dose_29_days',x)
  x<-gsub('TG_Gate_rat_in_vivo_rep_8_day','TG_GATEs_Rat_invivo_Repeated_Dose_8_days',x)
  x<-gsub('PPI (STRINGdb)','PPI_STRINGdb',x,fixed = T)
  x<-gsub('Consensus_Hepatocytes (DM,TG-gata)','Consensus_Rat_invitro_DM_TG_gates',x,fixed = T)
  x<-gsub('Consensus_Rat_invitro (Drug Matrix_,TG-gata)','Consensus_Rat_invitro_DM_TG_GATEs',x,fixed = T)
  x<-gsub("Consensus_LINCS_HEPG2(I,II)",'Consensus_LINCS_HEPG2',x,fixed = T)
  x<-gsub("Consensus_Rat_invitro_DM_TG_gates","Consensus_Rat_invitro_DM_TG_GATEs",x,fixed = T)
  return(x) 
}#end function

name_order<-function(inp){
  nams<-factor(inp[c(grep(x=inp,pattern = 'Drug_Matrix_Rat_invitro'),
                                grep(x=inp,pattern = 'Drug_Matrix_Rat_invivo_Single_Dose_1_day'),
                                grep(x=inp,pattern = 'Drug_Matrix_Rat_invivo_Repeated_Dose_3_days'),
                                grep(x=inp,pattern = 'Drug_Matrix_Rat_invivo_Repeated_Dose_5_days'),
                                grep(x=inp,pattern = 'TG_GATEs_Human_invitro'),
                                grep(x=inp,pattern = 'TG_GATEs_Rat_invitro'),
                                grep(x=inp,pattern = 'Low'),grep(x=inp,pattern = 'Middle'),
                                grep(x=inp,pattern = 'High'),
                                grep(x=inp,pattern = '8_days'),grep(x=inp,pattern = '15_days'),
                                grep(x=inp,pattern = '29_days'),
                                grep(x=inp,pattern = 'Consensus'),
                                grep(x=inp,pattern = 'PPI'))])
  return(nams)
} #end function

#calculates harmonic sum data can be a matrix or a data frame, col.scores are the selected columns
#the rows must have names
harmonic_sum<-function(data,col.scores){
  HS <- rep(NA, nrow(data))
  for(i in 1:nrow(data)){
    dtemp<-as.numeric(data[i,col.scores])
    dtemp<-dtemp[!is.na(dtemp)]
    HS[i] <- sum(sort(dtemp,decreasing=T)/(1:length(dtemp))^2)
  } # end for
  names(HS)<-rownames(data)
  return(HS)
} # end function


#calculates mean data can be a matrix or a data frame, col.scores are the selected columns
#the rows must have names
mean_function<-function(data,col.scores){
 means<-rep(NA,nrow(data))
 for (i in 1:length(means)){
   dtemp<-as.numeric(data[i,col.scores]) 
   means[i]<-  mean(dtemp,na.rm = T)
 }
names(means)<-rownames(data)
return(means)
}

# ROC permutation teset https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2687965/pdf/btp211.pdf
roc_permutation<-function(tst_score,gold_standard,npermutation=1000,seed='123'){#tst_scsore should ne numeric, gold_standard should be one for positive and 0 for negative response
  set.seed(seed)
  require(PRROC)
  x_zero<-roc.curve(scores.class0 = tst_score, weights.class0 = gold_standard,curve = F)
  x_zero<-x_zero$auc
  res<-rep(NA,length(npermutation))
  for (i in 1:npermutation){
    ind<-sample(1:length(tst_score),length(tst_score),replace = F)
    gold_standard<-gold_standard[ind]
    roc<-roc.curve(scores.class0 = tst_score, weights.class0 = gold_standard,curve = F)
    res[i]<-roc$auc
  } #end for
  p_perm<-length(which(res>=x_zero))/npermutation
  return(p_perm)
} #end function

