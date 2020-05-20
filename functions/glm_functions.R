###################
#fucntions  are: 

#1.data_preparation(x) x is the list of fgsea RData files from the pipeline and edc_names is the name of compounds to be in training set 
#                      returns an integrated list for all layers and the y vectors for each layer with 'edc' and 'decoy' labels
#                      Example:  li<-list.files(path='.',pattern = 'RData')
#                              edcs<-c('cafeine','letrozole',,,,)
#                               fgs<-data_preparation(li,edcs)
#
#2.elastic_net_model(data,lab,trcntrl,seed) data is a  matrix of training set where each column is a feature,lab should be factor as response,trcntrl is model tranin_control
#                      returns model of glm
#
#3.glm_model_stability(xx,yy,kk,reppeat,neg_class='no',seed=1) performs k_fold_cross_validation for a binary glm model xx is training data matrix, y is 
#                     response vecotor as binary factor, reppeat number of cross-validation repeating, neg_class is label for negative label
#                     returns accuracy, precision, f1, mcc
#
#4.elastic_net_bootstrap(data,lab,trcntrl,percentage=70,iter=100,seed=1) performs bootstrapng test; data is a training matrix, lab is factor of 
#                                                                        labels, trcntrl is training control options,percentage is the 
#                                                                        percent of data picked for bootstrap, iter is the number of bootstrapping
#                                                                        returns a binray matrix for each feature, each row is one bootstrap
#                                                                        iteration.
# 
#5.prob_predict(x,y) predicts class probilites for y based on models from x, y is the list of all layers as output of function data_preparation 
#                                                                            x is the list of all layers as output of function elastic_net_model
#                                                                            returns a matrix for each layer of data length of x should be length y

# preprocess --------------------------------------------------------------


#preparing data,  x is the list of  files from fgsea  customized pipeline
data_preparation<-function(x,edc_names,pathway_names=NULL){
  li<-x
  inf<-c()
  cnt<-1
  for (i in li){
    nam<-tools::file_path_sans_ext(i)
    print(nam)
    #rm(fgsea_res)
    load(i)
    if (length(pathway_names)>0){
      print('with selected pathways')
      data<-fgsea_res$NES[which(as.character(rownames(fgsea_res$NES)) %in% pathway_names),]
    }
    data<-t(data)                                    # after this transpose the chemicals are rows and the pathways are columns
   
    data[is.na(data)]<-0                             #  putting NA NES scores equal to zero
    if (length(which(colSums(abs(data))==0))>=1){    #  removing the pathways with zero values for all chemicals
    data<-as.matrix(data[,!colSums(abs(data))==0])
    }else{
      data<-as.matrix(data)
      }
    
    lab<-rep('decoy',nrow(data)) # y vector
    lab[which(rownames(data) %in% edc_names)]<-'edc'
    inform<-c()
    inform$x<-data
    inform$y<-lab
    inform$n_edc<-length(which(lab=='edc'))
    inform$n_decoy<-length(which(lab=='decoy'))
    inf[[cnt]]<-inform
    cnt<-cnt+1
  }
  names(inf)<-li 
  return(inf)
}






# glm_modeling ------------------------------------------------------------

#making model from glm
elastic_net_model_two_class<-elastic_net_model <-function(data,lab,trcntrl,seed=1){ 
  set.seed(seed)
  modl<-train(data, lab,
              method = "glmnet",
              family = "binomial",
              metric = "ROC",
              tuneLength = 20,
              trControl = trcntrl)
  mfit.ctcf.2.coef = as.matrix(coef(modl$finalModel, modl$bestTune$lambda))
  sel<-rownames(mfit.ctcf.2.coef)[which(abs(mfit.ctcf.2.coef)>0.0)] 
  return(list(modl=modl,sel=sel,cof=mfit.ctcf.2.coef))
}   #two_class


#making model from glm for multi label classifications
elastic_net_model_multi_class <-function(data,lab,trcntrl,seed=1){ 
  set.seed(seed)
  modl<-train(data, lab,
              method = "glmnet",
              family = "multinomial",
              metric = "logLoss",
              tuneLength = 20,
              trControl = trcntrl)
  return(modl)
} #multi class


#stabiltiy of the models from glmnet using k-fold cross validation 
glm_model_stability<-function(xx,yy,kk,reppeat,neg_class='no',seed=1){
  n<-nrow(xx)
  set.seed(seed)
  f <- ceiling(n/kk)
  accuracy<-precision<-recall<-mcc<-f1<-sensittivity<-specifficity<-tps<-tns<-fps<-fns<-matrix(0, nrow = kk, ncol = reppeat)
  
  for (j in 1:reppeat) {                  # repeating  CV j times
    s <- sample(rep(1:kk, f), n)
    for (i in 1:kk) {                      # i=1 to k th fold
      train.index <- seq_len(n)[(s != i)]   # training data index
      train.x<-xx[train.index,]
      train.y<-yy[train.index]
      elastic_net_model <- train(train.x,train.y,
                                 method = "glmnet",
                                 family = "binomial",
                                 metric = "ROC",
                                 tuneLength = 20,
                                 trControl = train_control)
      
      test.index <- seq_len(n)[(s == i)]    #test data index
      test.x<-xx[test.index,]
      test.y<-yy[test.index]
      predtest <- predict(elastic_net_model, test.x)
      pred<-ifelse(predtest %in% neg_class,0,1)
      real<-ifelse(test.y %in% neg_class,3,5)
      sumrp<-real+pred
      tp<-length(which(sumrp==6))
      tn<-length(which(sumrp==3))
      fn<-length(which(sumrp==5))
      fp<-length(which(sumrp==4))
      accuracy[i,j]=(tp+tn)/(tp+tn+fp+fn) 
      precision[i,j]=tp/(tp+fp)
      recall[i,j]=tp/(tp+fn)
      f1[i,j]=(2*precision[i,j]*recall[i,j])/(precision[i,j]+recall[i,j])
      #mcc[i,j]=((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
      sensittivity[i,j]=tp/(tp+fn)
      specifficity[i,j]=tn/(tn+fp)
     
      tps[i,j]<-tp
      tns[i,j]<-tn
      fps[i,j]<-fp
      fns[i,j]<-fn
    }#end for k fold
  } #end for repeat_cv
  return(list(models=elastic_net_model,accuracy=accuracy,f1=f1,precision=precision,recall=recall,sensittivity=sensittivity,specifficity=specifficity,tp=tps,fp=fps,fn=fns,tn=tns))
}

#stabiltiy of the features from glmnet using bootstrap
elastic_net_bootstrap <-function(data,lab,trcntrl,percentage=70,iter=100,seed=1){ 
  set.seed(seed)
  nk<-ceiling((percentage/100)*nrow(data))
  sel<-matrix(0, nrow = iter, ncol = ncol(data))
  colnames(sel)<-colnames(data)
  for (i in 1:iter){
    
    ind<-sample(1:nrow(data),nk,replace = TRUE)
    dx<-data[ind,]
    dy<-lab[ind]
    modl<-train(dx, dy,
                method = "glmnet",
                family = "binomial",
                metric = "ROC",
                tuneLength = 20,
                trControl = trcntrl)
    mfit.ctcf.2.coef = as.matrix(coef(modl$finalModel, modl$bestTune$lambda))
    relevant<-mfit.ctcf.2.coef[-1]
    sel[i,which(abs(relevant)>0.0)]<-1
  } #end for
  return(list(sel=sel))
} #end function


# stratified kfold CV using caret
caret_based_two_class_glm_repeated_k_fold_CV<-function(xx,yy,kk,reppeat,neg_class='no',train_control=train_control){
  all_index<-caret::createMultiFolds(yy,k=kk,times=reppeat)
  accuracy<-precision<-recall<-f1<-sensittivity<-specifficity<-tps<-tns<-fps<-fns<-rep(0, length(all_index))
  coefs<-list()
  for (i in 1:length(all_index)) {                  # repeating  CV j times
    train.index <- all_index[[i]]  # training data index
    train.x<-xx[train.index,]
    train.y<-yy[train.index]
    elastic_net_model <- elastic_net_model_two_class(train.x,train.y,train_control)
    #collecting the coefficients of glm models
    mfit.ctcf.2.coef = as.vector(coef(elastic_net_model$modl$finalModel, elastic_net_model$modl$bestTune$lambda))
    relevant<-mfit.ctcf.2.coef[-1]
    names(relevant)<-colnames(xx)
    coefs[[i]]<-relevant
    # evaluations
    test.index <-     #test data index
      test.x<-xx[-train.index,]
    test.y<-yy[-train.index]
    predtest <- predict(elastic_net_model$modl, test.x)
    pred<-ifelse(predtest %in% neg_class,0,1)
    real<-ifelse(test.y %in% neg_class,3,5)
    sumrp<-real+pred
    tp<-length(which(sumrp==6))
    tn<-length(which(sumrp==3))
    fn<-length(which(sumrp==5))
    fp<-length(which(sumrp==4))
    accuracy[i]=(tp+tn)/(tp+tn+fp+fn) 
    precision[i]=tp/(tp+fp)
    recall[i]=tp/(tp+fn)
    f1[i]=(2*precision[i]*recall[i])/(precision[i]+recall[i])
    sensittivity[i]=tp/(tp+fn)
    specifficity[i]=tn/(tn+fp)
    tps[i]<-tp
    tns[i]<-tn
    fps[i]<-fp
    fns[i]<-fn
  } #end for repeat_k_cv
  return(list(indexes=all_index,coefs=coefs,accuracy=accuracy,f1=f1,precision=precision,recall=recall,sensittivity=sensittivity,specifficity=specifficity,tp=tps,fp=fps,fn=fns,tn=tns))
}




#prediction of class probilities based on the  models from the function elastic_net_model
prob_predict<-function(x,y){
  result=c()
  for (i in 1:length(x)){
    res<-predict(x[[i]]$modl,y[[i]]$x,type='prob')
    result[[i]]=as.matrix(res)
    rownames(result[[i]])=rownames(y[[i]]$x)
  }
  names(result)=names(x)
  return(result)
}


# Confusion matrix --------------------------------------------------------
confusion_matrix<-function(real,pred){#two vectors of real observations and predicted observations it can be two or multiclass labeled
  all_conditions<-unique(union(real,pred))
  cm<-matrix(NA,ncol = length(all_conditions),nrow = length(all_conditions))
  rownames(cm)<-colnames(cm)<-all_conditions
  for (i in all_conditions){
    for (j in all_conditions){
      cm[i,j]<-length(which(real %in% i ==pred %in% j & real %in% i==TRUE & pred %in% j==TRUE))
    } #end for pred
  }   # end for real
  return(cm)
}# end function

cm_based_evaluation <- function(cm, type = "basic"){ ## http://blog.revolutionanalytics.com/2016/03/com_class_eval_metrics_r.html
  n = sum(cm) # number of instances
  nc = nrow(cm) # number of classes
  diag = diag(cm) # number of correctly classified instances per class 
  rowsums = apply(cm, 1, sum) # number of instances per class
  colsums = apply(cm, 2, sum) # number of predictions per class
  p = rowsums / n # distribution of instances over the actual classes
  q = colsums / n # distribution of instances over the predicted classes
  # metrics
  accuracy = sum(diag) / n 
  precision = diag / colsums 
  recall = diag / rowsums 
  # macro-averaged Metrics
  f1 = 2 * precision * recall / (precision + recall) 
  macroPrecision = mean(precision)
  macroRecall = mean(recall)
  macroF1 = mean(f1,na.rm = T)
  # one vs all
  oneVsAll = lapply(1 : nc,function(i){
    v = c(cm[i,i],
          rowsums[i] - cm[i,i],
          colsums[i] - cm[i,i],
          n-rowsums[i] - colsums[i] + cm[i,i]);
    return(matrix(v, nrow = 2, byrow = T))})
  s = matrix(0, nrow = 2, ncol = 2)
  for(i in 1 : nc){s = s + oneVsAll[[i]]}
  #print(s)
  # - avearge accuracy
  avgAccuracy = sum(diag(s)) / sum(s)
  # - micro average metrics
  micro_prf = (diag(s) / apply(s,1, sum))[1]
  # Majority-class Metrics
  mcIndex = which(rowsums==max(rowsums))[1] # majority-class index
  mcAccuracy = as.numeric(p[mcIndex]) 
  mcRecall = 0*p;  mcRecall[mcIndex] = 1
  mcPrecision = 0*p; mcPrecision[mcIndex] = p[mcIndex]
  mcF1 = 0*p; mcF1[mcIndex] = 2 * mcPrecision[mcIndex] / (mcPrecision[mcIndex] + 1)
  if(type == "basic")
    return(list(accuracy=accuracy, precision=precision, recall=recall, f1=f1))
  else if(type == "macro")
    return(list(totAcc = accuracy, allPrec = precision, allRec = recall, allF1 = f1,
                macAcc = avgAccuracy, macPrec = macroPrecision, macRec = macroRecall, macF1 = macroF1,
                micAcc = micro_prf))
  else if(type == "majority")
    return(list(mcAccuracy=mcAccuracy, mcPrecision=mcPrecision, mcRecall=mcRecall, mcF1=mcF1))
}


