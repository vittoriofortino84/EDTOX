# Distributio of EDC score on pie chart
all_vam<-readRDS('outputData/excel_files/all_edc_scores.rds')   # all edc scores
intervals<-seq(.2,1,.2)
dist_res<-matrix(NA,length(intervals),3)
j=0;k=1

for (i in intervals){
  dist_res[k,3]<-length(which(all_vam$harmonic_sum_edc_score>=j & all_vam$harmonic_sum_edc_score<i))
  dist_res[k,2]<-length(which(all_vam$average_edc_score>=j & all_vam$average_edc_score<i))
  dist_res[k,1]<-paste(j,i,sep = '-')
  j=i;k=k+1
}

dist_res<-as.data.frame(dist_res)

colnames(dist_res)<-c('Category','average_number','harmonic_number')

dist_res$colr<-c('green3', 'pink', 'red', 'yellow', 'darkorange')


dist_res$average_number<-as.numeric(dist_res$average_number)
dist_res$harmonic_number<-as.numeric(dist_res$harmonic_number)

dist_res$prop_average<-dist_res$average_number/sum(dist_res$average_number)
dist_res$prop_harmonic<-dist_res$harmonic_number/sum(dist_res$harmonic_number)

data_average<-dist_res[,c("Category" ,"average_number","prop_average",'colr')]
data_harmonic<-dist_res[,c("Category" ,"harmonic_number","prop_harmonic",'colr')]
colnames(data_harmonic)<-colnames(data_average)<-c('Category','Number','prop','colr')
donut_plot<-function(data,title){
  require(dplyr)
  
  count.data<-data
  count.data <- count.data %>% 
    arrange(desc(Number)) %>%
    mutate(lab.ypos = cumsum(prop) - 0.5*prop)
  
  return(count.data)
}

average_edc<-donut_plot(data_average,'Average EDC score')
harmonic_edc<-donut_plot(data_harmonic,'Harmonic EDC score')
pie_chart<-function(pie_data,title){
  require(ggplot2)
  ggplot(pie_data, aes(x = 2, y = prop, fill = factor(Category,levels = rev(Category)))) +
  geom_bar(stat = "identity", color = "white") +
    coord_polar(theta = "y")+labs(fill=title)+
    
    geom_text(aes(y = lab.ypos, label = Number), color = "Black",size=4)+
    scale_fill_manual(values = rev(pie_data$colr))+
    
    theme_void()+
    xlim(0.5, 2.5)
}


pie_chart(average_edc,'Average EDC score')
pie_chart(harmonic_edc,'Harmonic EDC score')

#saveRDS(all_res,file = 'inputData/statistics/distribution_compounds_EDCs.rds')

