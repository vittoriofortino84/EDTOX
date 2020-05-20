# Distributio of EDC score on pie chart
edc_decoy<-readRDS('outputData/excel_files/all_edc_scores.rds')   # all edc scores
DiabetesT2=read.csv('inputData/disease_scores/diabetes_scores.csv',stringsAsFactors = F)
Atherosclerosis=read.csv('inputData/disease_scores/atherosclerosis_scores.csv',stringsAsFactors = F)
Metabolic_syndrome=read.csv('inputData/disease_scores/metabolic_syndrome_scores.csv',stringsAsFactors = F)

pie_distribution_function<-function(all_scores,harmonic_s,average_s){
dist_res_function<-function(all_scores,harmonic_s,average_s){
  harmonic_sum<-all_scores[,colnames(all_scores) %in% harmonic_s]
  average_score<-all_scores[,colnames(all_scores) %in% average_s]
  
  intervals<-seq(.2,1,.2)
  dist_res<-matrix(NA,length(intervals),3)
  j=0;k=1
  for (i in intervals){
  dist_res[k,3]<-length(which(harmonic_sum>=j & harmonic_sum<i))
  dist_res[k,2]<-length(which(average_score>=j & average_score<i))

  dist_res[k,1]<-paste(j,i,sep = '-')
  j=i;k=k+1}
dist_res<-data.frame(dist_res,stringsAsFactors = F)
colnames(dist_res)<-c('Category','average_number','harmonic_number')
dist_res$colr<-c('green3', 'pink', 'red', 'yellow', 'darkorange')

dist_res$average_number<-as.numeric(dist_res$average_number)
dist_res$harmonic_number<-as.numeric(dist_res$harmonic_number)

dist_res$prop_average<-dist_res$average_number/sum(dist_res$average_number)
dist_res$prop_harmonic<-dist_res$harmonic_number/sum(dist_res$harmonic_number)
data_average<-dist_res[,c("Category" ,"average_number","prop_average",'colr')]
data_harmonic<-dist_res[,c("Category" ,"harmonic_number","prop_harmonic",'colr')]
colnames(data_harmonic)<-colnames(data_average)<-c('Category','Number','prop','colr')
return(list(data_harmonic=data_harmonic,data_average=data_average))
}
all_vam_scores<-dist_res_function(all_scores,harmonic_s,average_s)
donut_plot<-function(data,title){
  require(dplyr)
  
  count.data<-data
  count.data <- count.data %>% 
    arrange(desc(Number)) %>%
    mutate(lab.ypos = cumsum(prop) - 0.5*prop)
  
  return(count.data)
}
pie_chart<-function(pie_data,title){
  require(ggplot2)
  ggplot(pie_data, aes(x = 2, y = prop, fill=factor(Category,levels=rev(Category)))) +
    geom_bar(stat = "identity", color = "white") +
    coord_polar(theta = "y")+labs(fill=title)+
    
    geom_text(aes(y = lab.ypos, label = Number), color = "Black",size=4)+
    scale_fill_manual(values = rev(pie_data$colr))+
    facet_wrap(.~pie_data$type)+
    theme_void()+
    xlim(0.5, 2.5)
}
average_edc<-donut_plot(all_vam_scores$data_average,'Average EDC score');average_edc$type='Average'
harmonic_edc<-donut_plot(all_vam_scores$data_harmonic,'Harmonic EDC score');harmonic_edc$type='harmonic'
harmonic_average_combined<-rbind(average_edc,harmonic_edc)
p1<-pie_chart(average_edc,'')
p2<-pie_chart(harmonic_edc,'')
ggpubr::ggarrange(p1,p2,nrow = 1,common.legend = T,legend = 'bottom')}

pie_distribution_function(edc_decoy,'harmonic_sum_edc_score','average_edc_score')  #edc score classification
pie_distribution_function(Atherosclerosis,'harmonic_disease_score','average_disease_score') #Atherosclerosis
pie_distribution_function(DiabetesT2,'harmonic_disease_score','average_disease_score') #Diabetes T2
pie_distribution_function(Metabolic_syndrome,'harmonic_disease_score','average_disease_score') #Metabolic syndrome

save(pie_distribution_function,edc_decoy,Atherosclerosis,DiabetesT2,Metabolic_syndrome,
     file = 'outputData/plots_tables_functions/disease_score_distribution/disease_scores_pie_chart.RData')