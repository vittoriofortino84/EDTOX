
# Kegg and Reactome -------------------------------------------------------
#save picture 17*7
rm(list=ls())
library(RColorBrewer)
source('functions/general_functions.R')
source('functions/annotation_functions.R')
load('outputData/patway_glm_coefficients/GLM_univariate_all_scores.RData')  # For moa pathways
# treatment
all_data_final$glm_coefs<-abs(all_data_final$glm_coefs)
nams<-name_order(unique(all_data_final$network))               # all networks
all_data_final<-all_data_final[which(all_data_final$network %in% nams),]

data<-all_data_final[which(all_data_final$glm_coefs>0 &
                           all_data_final$Average_NES_EDC>0),
                           c('pathways','network')]  # the pathway in at least n most informative networks

data<-data[-grep(x=data$pathways,pattern ='Wiki'),] #removing wiki
ntype<-unique(data$network)

#kegg
data_KEGG<-data[grep(x=data$pathways,pattern ='KEGG'),]
keggs_prop<-lapply(ntype, function(x)kegg_reactome_2_mainclass(data_KEGG$pathways[which(data_KEGG$network %in% x)],
type='KEGG')$prop)
names(keggs_prop)<-ntype
for (i in 1:length(keggs_prop))keggs_prop[[i]]$network=ntype[i]
keggs_prop<-do.call(rbind,keggs_prop)

keggs_prop$class<-as.character(keggs_prop$class)

kegg_lvls<- c(
  "Cellular Processes"                  , 
  "Environmental Information Processing",          
 "Genetic Information Processing"       ,
 "Organismal Systems"                   ,                            
 "Human Diseases"                       ,
 "Metabolism-amino acid "               ,                        
 "Metabolism-carbohydrate "             ,
 "Metabolism-glycan "                   ,                           
 "Metabolism-lipid "                    ,
 "Metabolism-nucleotide "               ,                        
 "Metabolism-energy "                   ,
 "Metabolism-cofactors and vitamins"    ,             
 "Metabolism-biosynthesis secondary metabolites",
 "Metabolism-xenobiotics biodegradation",         
 "Metabolism-biosynthesis polyketides,terpenoids")
correct_omitted<-function(df,class_lvls){
f<-lapply(ntype, function(x){
  dt<-NA
  temp<-setdiff(class_lvls,unique(df$class[which(df$network== x)]))
  if (length(temp)>0){dt<-as.data.frame(temp);colnames(dt)<-'class'
  
  dt$network<-x}
  return(dt)
})
dt<-do.call(rbind,f)
dt<-na.omit(dt)
dt$prop<-.1
return(dt)
}
zero_class<-correct_omitted(keggs_prop,kegg_lvls)
keggs_prop<-rbind(keggs_prop,zero_class)


g2_names<-c("TG_GATEs_Rat_invivo_Repeated_Dose_29_days","TG_GATEs_Rat_invivo_Repeated_Dose_15_days",
            "TG_GATEs_Rat_invivo_Repeated_Dose_8_days", "TG_GATEs_Rat_invivo_Single_Dose_High_1_day",
            "TG_GATEs_Rat_invivo_Single_Dose_Middle_1_day")

g3_names<-c("TG_GATEs_Rat_invitro_Single_Dose_1_day", "TG_GATEs_Human_invitro_Single_Dose_1_day",
            "Consensus_LINCS_HEPG2","Consensus_Rat_invitro_DM_TG_GATEs","PPI_STRINGdb" )

keggs_prop$gr<-'g1'
keggs_prop$gr[keggs_prop$network %in% g2_names]<-'g2'
keggs_prop$gr[keggs_prop$network %in% g3_names]<-'g3'
library(ggplot2)
 p1<-ggplot(keggs_prop, aes(x=factor(network,levels = nams),y=prop,fill=factor(class,levels = kegg_lvls))) +
  geom_bar(stat = "identity",width = .8,position = position_dodge2(preserve = "single",width = .9))+

  scale_fill_manual(values=c(brewer.pal(9, 'Blues')[c(2,4,6,7,9)],
                             brewer.pal(9, 'Reds')[c(3,4,5,6,7)],
                             brewer.pal(9, 'Oranges')[c(3,4,5,6,7)]))+
  theme_minimal()+
  xlab('')+ylab('% pathways')+ 
 scale_y_continuous(breaks = seq(0, max(keggs_prop$prop), 5), limits = c(0, max(keggs_prop$prop))) +
  labs(fill='Category')+
  theme(legend.text = element_text(size=10),
        strip.text = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'bottom',axis.text.x = element_text(size=10,angle = 0),
        axis.title = element_text(size = 12),axis.text.y = element_text(size=10))+
        facet_wrap(factor(keggs_prop$gr,levels = c('g1','g2','g3')),nrow = 3,scales = 'free')
print(p1)
save(keggs_prop,kegg_lvls,nams,file = 'outputData/plots_tables_functions/Enrichment_Result_categorization/kegg.RData')


## Reactome
data_REACTOME<-data[grep(x=data$pathways,pattern ='R-HSA'),]
react_prop<-lapply(ntype, function(x)kegg_reactome_2_mainclass(data_REACTOME$pathways[which(data_REACTOME$network %in% x)],
type='REACTOME')$prop)
names(react_prop)<-ntype
for (i in 1:length(react_prop))react_prop[[i]]$network=ntype[i]
react_prop<-do.call(rbind,react_prop)
react_prop$class<-as.character(react_prop$class)

#
 react_lvls<-c(
 "Vesicle-mediated transport"             ,
 "Cell-Cell communication"                ,
 "Cellular responses to external stimuli" ,
 "Gene expression (Transcription)"        ,
 "Immune System"                          ,
 "Transport of small molecules"           ,
 "Autophagy"                              ,
 "DNA Repair"                             ,
 "Muscle contraction"                     ,
 "Neuronal System"                        ,
 "Developmental Biology"                  ,
 "Circadian Clock"                        ,
 "DNA Replication"                        ,
 "Cell Cycle"                             ,
 "Digestion and absorption"               ,
 "Extracellular matrix organization"      ,
 "Chromatin organization"                 ,
 "Hemostasis"                             ,
 "Organelle biogenesis and maintenance"   ,
 "Protein localization"                   ,
 "Signal Transduction"                    ,
 "Programmed Cell Death"                  ,
 "Reproduction"                           ,
 "Disease"                                ,
 "Metabolism"                             ,
 "Metabolism of proteins"                 ,
 "Metabolism of RNA"
   )

zero_class<-correct_omitted(react_prop,react_lvls)
react_prop<-rbind(react_prop,zero_class)

react_prop$gr<-'g1'
react_prop$gr[react_prop$network %in% g2_names]<-'g2'
react_prop$gr[react_prop$network %in% g3_names]<-'g3'

## Reactome plot
p2<-ggplot(react_prop, aes(x=factor(network,levels = nams),y=prop,fill=factor(class,levels = react_lvls))) +
  geom_bar(stat='identity',width = .8,position = position_dodge2(preserve = "single",width = .9))+
   scale_fill_manual(values=c(brewer.pal(9, 'YlOrRd')[c(4,5,6,7)],
                              brewer.pal(9, 'YlGnBu')[c(4,5,6,7,8)],
                              brewer.pal(9, 'Greens'),
                              brewer.pal(9, 'Oranges')[c(4,5,6,7,8,9)],
                              brewer.pal(9, 'Reds')[c(7,8,9)]))+
  
  theme_minimal()+
  xlab('')+ylab('% pathways')+ 
  scale_y_continuous(breaks = seq(0, max(react_prop$prop), 5), limits = c(0, max(react_prop$prop))) +
  labs(fill='Category')+ 
  theme(legend.text = element_text(size=10),
                                strip.text = element_blank(),
                                axis.title.x = element_blank(),
                                legend.position = 'bottom',axis.text.x = element_text(size=10,angle = 0),
                                axis.title = element_text(size = 16),axis.text.y = element_text(size=14))+
  facet_wrap(factor(react_prop$gr,levels = c('g1','g2','g3')),nrow = 3,scales = 'free')
print(p2)

save(react_prop,react_lvls,nams,file = 'outputData/plots_tables_functions/Enrichment_Result_categorization/reactome.RData')

# #AOPs -------------------------------------------------------------------


data_GO<-data[grep(x=data$pathways,pattern ='GO'),]
data_GO$pathways<-pathway_annotate(data_GO$pathways)
AOPs_list<-sapply(ntype, function(x)kegg_reactome_2_mainclass(unique(data_GO$pathways[which(data_GO$network %in% x)]),
type='GO')$AOPS)



all_aops <- unique(unlist(AOPs_list))
int_data <- lapply(AOPs_list, function(x){
  idx <- which(all_aops %in% x)
  res<-rep(0, length(all_aops))
  res[idx]<-1
  res})
upset_df <- as.data.frame(int_data , as.is=T, stringsAsFactors=F, check.names=F)
rownames(upset_df) <- all_aops
#Aops upset
p3<-UpSetR::upset(upset_df, nintersects=12, nsets=ncol(upset_df), 
              sets.bar.color = "#56B4E9", order.by=c('freq'), text.scale=1.4, 
              point.size = 2.5, line.size = 0.7, 
              sets.x.label = "AOPs  in each network", mainbar.y.label = "AOPs intersection size")

print(p3)


AOPs_number<-sapply(ntype, function(x)kegg_reactome_2_mainclass(unique(data_GO$pathways[which(data_GO$network %in% x)]),
                                                               type='GO')$Number_aops)
aop_data<-as.data.frame(AOPs_number)
aop_data$network<-rownames(aop_data)

color_values<-c('#E6F2FF','#CCE6FF','#B3D9FF','#99CCFF',                     # Drug matrix hep,1,3,5
                '#F2F9EC','#E6F2D9',                                         # in vitro rat and human from TGgates
                '#D9ECC6','#CCE6B3','#BFDF9F','#B3D98C','#A6D279','#99CC66', # in vivo  rat low,high,moddle,8,15,29 TGgates
                '#FFCCBF','#FF9980','#FF6640')   

# Aops simple plot
p4<-ggplot(aop_data, aes(x = factor(network,levels = nams),
                           y = AOPs_number, fill = factor(network,levels = nams))) +
  geom_bar(stat = "identity") +ylab('')+
  #scale_fill_viridis_d() +
  scale_fill_manual(values = color_values)+
  theme_minimal()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+coord_flip()

print(p4)



 aop_data$prop<-aop_data$AOPs_number/sum(aop_data$AOPs_number)
 color_df<-as.data.frame(list(color_values,nams));colnames(color_df)<-c('color','network')
 aop_data<-merge(aop_data,color_df)
 library(dplyr)
 count.data<-aop_data
 count.data <- count.data %>%
   arrange(desc(AOPs_number)) %>%
   mutate(lab.ypos = cumsum(prop) - 0.5*prop)
 count.data

 
 # Aops donut plot
 p3<-ggplot(count.data, aes(x = 2, y = prop, fill = factor(network,levels = rev(network)))) +
   geom_bar(stat = "identity", color = "white") +
   coord_polar(theta = "y", start = 0)+labs(fill='Data Layer')+
   geom_text(aes(y = lab.ypos, label = AOPs_number), color = "Black",size=4)+
   scale_fill_manual(values=as.character(rev(count.data$color)))+
   theme_void()+
   xlim(0.5, 2.5)

 save(count.data,file = 'outputData/plots_tables_functions/Enrichment_Result_categorization/AOPS_count.RData')
print(p3) 