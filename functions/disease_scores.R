# it takes the argument meshs as one or multiple compounds and the name of the disease as string and returns the scores 
# related to the disease including the class probabilities as well as the harmonic and average scores
mesh2disease_score<-function(meshs,disease_name){
  if (disease_name=='Diabetes Type 2')path2score='inputData/disease_scores/diabetes_scores.csv'
  if (disease_name=='atherosclerosis')path2score='inputData/disease_scores/atherosclerosis_scores.csv'
  if(disease_name=='metabolic syndrome')path2score='inputData/disease_scores/metabolic_syndrome_scores.csv'
  disease_scores<-read.csv(path2score,stringsAsFactors = F)
  res<-matrix(NA,nrow = length(meshs),ncol = 21)
  for (i in 1:length(meshs)){
    sel<-  disease_scores[disease_scores$mesh==meshs[i],2:22]
    res[i,]<-as.character(sel[1,])
    }
  
  res<-as.data.frame(res)
  colnames(res)<-colnames(disease_scores)[2:22]
  return(res)
}