library(ggplot2);library(ggpubr);library(grid);library(facetscales)

# 1.  time response pattern recognition TG-GATEs 8 15 29 --------
rm(list=ls())
load('outputData/CTD_edc_scores_dictionary.RData')
all_vam<-all_vam[,grep(colnames(all_vam),pattern = 'invivo|is_in|comp_names')]
all_vam<-na.omit(all_vam)

delta<-all_vam[,c("TG_GATEs_Rat_invivo_Repeated_Dose_8_days",
                  "TG_GATEs_Rat_invivo_Repeated_Dose_15_days",
                  "TG_GATEs_Rat_invivo_Repeated_Dose_29_days"
)] # Repeated dose
groups<-list()
g1<-delta[apply(delta, 1,function(x) x[1] <.4  & x[2] < .4  & x[3] > .6),1:3]
groups[[1]]<-rownames(g1)

g2<-delta[apply(delta, 1,function(x) x[1] <.4  & x[2] > .6  & x[3] > .6),1:3]
groups[[2]]<-rownames(g2)

g3<-delta[apply (delta, 1,function(x) (x[1] >= .6 & x[1] < .8) & x[2] >=  (x[1] + .1) & 
                   x[3] >= x[2] ),1:3]
groups[[3]]<-rownames(g3)

g4<-delta[apply (delta, 1,function(x) x[1] >= 0.9 & x[2] >=  0.9 &  x[3] >= 0.9),1:3]
groups[[4]]<-rownames(g4)

g5<-delta[apply(delta, 1,function(x) x[1] > .6  & x[2] > .6  & x[3] < .4),1:3]
groups[[5]]<-rownames(g5)

gr_list<-list()
for (i in 1:length(groups)){
  ind<-groups[[i]]
  lvls<-c( "8_days", "15_days","29_days")
  Tg_rep<-all_vam[ind,grep(colnames(all_vam),pattern = 'TG_GATEs_Rat_invivo_Repeated_Dose|comp_names')]
  colnames(Tg_rep)<-gsub(x=colnames(Tg_rep),pattern ='TG_GATEs_Rat_invivo_Repeated_Dose_',replacement = '' )
  # colnames(Tg_rep)<-gsub(x=colnames(Tg_rep),pattern ='_days',replacement = '' )
  Tg_rep$mesh<-rownames(Tg_rep)
  Tg_rep<-tidyr::gather(data = Tg_rep,key='data_layer',value='class_prob',-c(comp_names,mesh))
  Tg_rep$cluster_gr<-paste('Group',i,sep = '_')
  gr_list[[i]]<-Tg_rep
}


data_prep<-function(inp){
  require(dplyr)
  df <- inp %>%
    group_by(data_layer) %>%
    summarise(
      sd = sd(class_prob),
      len = mean(class_prob)
    )
  return(df)
}
data<-lapply(gr_list, function(x)data_prep(x))
data[[1]]$pattern<-'No_EDC - No_EDC < EDC'
data[[2]]$pattern<-'No_EDC < EDC - EDC'
data[[3]]$pattern<-'EDC < EDC <= EDC'
data[[4]]$pattern<-'EDC - EDC - EDC'
data[[5]]$pattern<-'EDC - EDC > No_EDC'
lin_lvls<-factor(c('No_EDC - No_EDC < EDC','No_EDC < EDC - EDC','EDC < EDC <= EDC','EDC - EDC - EDC','EDC - EDC > No_EDC'))

final<-do.call(rbind,data)
time_final<-final;save(time_final,file='outputData/plot_data/time_exposure.RData')

library(RColorBrewer)
ggplot(final, 
       aes(x = factor(data_layer,levels = c('8_days','15_days','29_days')), 
           y = len, ymin = len-sd, ymax = len+sd,colour=factor(pattern,levels = lin_lvls)))+
  geom_pointrange()+geom_errorbar(width = .2) +geom_point(size = 1.5)+
  scale_color_manual(values=rev(brewer.pal(n = 6, name ='RdYlGn')))+
  xlab('TG_GATEs_Single_Repeated_Dose')+ylab('Class Probability')+ theme_minimal()+
  labs(colour='Dose Response Groups')+
  geom_line(aes(group=pattern),show.legend = T)+
  theme(legend.text = element_text(size=10),
        legend.position = 'right',axis.text.x = element_text(size=16),
        axis.title = element_text(size = 16),axis.text.y = element_text(size=16),
        panel.spacing = unit(2,'lines'))


Tg_rep_single<-do.call(rbind,gr_list)
Tg_rep_single$desc<-'EDC'
Tg_rep_single$desc[Tg_rep_single$cluster_gr=='Group_1']<-'No_EDC - No_EDC < EDC'
Tg_rep_single$desc[Tg_rep_single$cluster_gr=='Group_2']<-'No_EDC < EDC - EDC'
Tg_rep_single$desc[Tg_rep_single$cluster_gr=='Group_3']<-'EDC < EDC <= EDC'
Tg_rep_single$desc[Tg_rep_single$cluster_gr=='Group_4']<-'EDC - EDC - EDC'
Tg_rep_single$desc[Tg_rep_single$cluster_gr=='Group_5']<-'EDC - EDC > No_EDC'

Tg_rep_single$sec_gr<-'Group'
Tg_rep_single$sec_gr[Tg_rep_single$cluster_gr=='Group_1']<-'g1'
Tg_rep_single$sec_gr[Tg_rep_single$cluster_gr=='Group_2']<-'g1'
Tg_rep_single$sec_gr[Tg_rep_single$cluster_gr=='Group_3']<-'g2'
Tg_rep_single$sec_gr[Tg_rep_single$cluster_gr=='Group_4']<-'g2'
Tg_rep_single$sec_gr[Tg_rep_single$cluster_gr=='Group_5']<-'g3'

Tg_rep_single$f_gr<-'Group'
Tg_rep_single$f_gr[Tg_rep_single$cluster_gr=='Group_1']<-'f1'
Tg_rep_single$f_gr[Tg_rep_single$cluster_gr=='Group_2']<-'f2'
Tg_rep_single$f_gr[Tg_rep_single$cluster_gr=='Group_3']<-'f1'
Tg_rep_single$f_gr[Tg_rep_single$cluster_gr=='Group_4']<-'f2'
Tg_rep_single$f_gr[Tg_rep_single$cluster_gr=='Group_5']<-'f1'


scales_y <- list(
  `g1` = scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .1)),
  `g2` = scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, .05)),
  `g3` = scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .1))
  
)
library(RColorBrewer)
ggplot(Tg_rep_single, aes(x=factor(data_layer,levels = lvls), y=class_prob,colour=factor(desc,levels = lin_lvls))) +
  geom_point(show.legend = F)+
  geom_line(aes(group=comp_names),show.legend = T)+
  #guides(col=guide_legend(nrow = 2))+
  labs(colour='Exposure Time Groups')+
  scale_color_manual(values=rev(brewer.pal(n = 6, name ='RdYlGn')))+
  #scale_color_gradient(low='yellow',high = 'red')+
  xlab('TG_GATEs_Repeated_Dose')+ylab('Class Probability')+ theme_minimal()+
  theme(strip.text = element_blank(),legend.text = element_text(size=10),
        legend.position = 'right',axis.text.x = element_text(size=16),
        axis.title = element_text(size = 20),axis.text.y = element_text(size=16),panel.spacing = unit(2,'lines')) +
  facet_grid_sc(rows=vars(sec_gr),cols = vars(f_gr),scales = list(y=scales_y))







# 2. preselection and anova for the pathways based on time ---------------------------
load('outputData/patway_glm_coefficients/GLM_univariate_all_scores.RData')

all_data_final<-all_data_final[which(all_data_final$ROC_AUC >0.7),]
Roc_based_preselected_pathways<-all_data_final$pathways[grep(x=all_data_final$network,
                                                             pattern = 'TG_GATEs_Rat_invivo_Repeated_Dose',ignore.case = F)]


load("outputData/class_probilities/integrated_fgsea_allcompounds_results_final_moa.RData")


fgs<-fgs[c("TG_GATEs_Rat_invivo_Repeated_Dose_8_days",
           "TG_GATEs_Rat_invivo_Repeated_Dose_15_days",
           "TG_GATEs_Rat_invivo_Repeated_Dose_29_days")] #All predicted pathway scores using RWR-FGSEA strategy
names(fgs)<-c('eight_day','fifteen_day','twentynine_day')
common_pathways<-Reduce('intersect',lapply(fgs, function(x)colnames(x$x))) # COMMON PATHWAYS IN LOW MIDDLE AND HIGH
common_pathways<-intersect(common_pathways,Roc_based_preselected_pathways)

load('outputData/CTD_edc_scores_dictionary.RData')
all_vam<-all_vam[,grep(colnames(all_vam),pattern = 'invivo|is_in|comp_names')]
all_vam<-na.omit(all_vam)
delta<-all_vam[,c("TG_GATEs_Rat_invivo_Repeated_Dose_8_days",
                  "TG_GATEs_Rat_invivo_Repeated_Dose_15_days",
                  "TG_GATEs_Rat_invivo_Repeated_Dose_29_days"
)] # Single dose class probabilities from elastic net


pas_anova<-function(g_edc){
  edcs<-list();for (i in 1:length(fgs)){edcs[[i]]<-fgs[[i]]$x[g_edc,common_pathways]}
  names(edcs)<-names(fgs)
  
  # for each pathway we do ANOVA with Tukey test
  res<-matrix(NA, nrow = length(common_pathways), ncol = 6)
  for (i in 1:length(common_pathways)){
    res[i,1]<-mean(edcs$eight_day[,i])
    res[i,2]<-mean(edcs$fifteen_day[,i])
    res[i,3]<-mean(edcs$twentynine_day[,i])
    anov_inp<-as.data.frame(list(edcs$eight_day[,i],edcs$fifteen_day[,i],edcs$twentynine_day[,i]))
    colnames(anov_inp)<-c('edc_eight_day','edc_fifteen_day','edc_twentynine_day')
    anov_inp<-stack(anov_inp)
    res.aov<-aov(values ~ ind,data=anov_inp)
    tk<-TukeyHSD(res.aov)
    tk<-as.data.frame(tk$ind)
    res[i,4:6]<-tk$`p adj`
  }
  colnames(res)<-c('average_edc_eight_day','average_edc_fifteen_day','average_edc_twentynine_day',
                   paste('p_value',rownames(tk),sep='_')
  )
  rownames(res)<-common_pathways
  res<-as.data.frame(res)
  return(res)
}
# selection of the compounds as edcs 
g1_edc<-delta[apply(delta, 1,function(x) x[1] < .4 & x[2] < .4 & x[3] > .6),1:3];g1_edc<-rownames(g1_edc)
g2_edc<-delta[apply(delta, 1,function(x) x[1] < .4 & x[2] > .6 & x[3] > .6),1:3];g2_edc<-rownames(g2_edc)

g3_edc<-delta[apply(delta, 1,function(x) (x[1] >= .6 & x[1] < .8) & x[2] >=  (x[1] + .1) & x[3] >= x[2] ),1:3]
g3_edc<-rownames(g3_edc)
g4_edc<-delta[apply(delta, 1,function(x) x[1] >= .9 & x[2] >=  .9 &  x[3] >= .9),1:3];g4_edc<-rownames(g4_edc)
g5_edc<-delta[apply(delta, 1,function(x) x[1] > .6  & x[2] > .6  & x[3] < .4),1:3];g5_edc<-rownames(g5_edc)
# anova on pathway activation scores
anova_time_1_5<-lapply(list(g1_edc,g2_edc,g3_edc,g4_edc,g5_edc), function(x)pas_anova(x))
names(anova_time_1_5)<-c('g1','g2','g3','g4','g5')

# group1 8d<.4, 15d<.4 and 29d>.6
res<-anova_time_1_5$g1
g1_pathways<-rownames(res)[which( 
    res$`p_value_edc_twentynine_day-edc_eight_day`                   <= 0.05 &
    res$`p_value_edc_twentynine_day-edc_fifteen_day`                 <= 0.05 &
    res$`p_value_edc_fifteen_day-edc_eight_day`                       > 0.05 & 
    res$average_edc_twentynine_day > res$average_edc_fifteen_day             &
    res$average_edc_twentynine_day > res$average_edc_eight_day
)]

# group2 8d<.4 15d>.6 and 29d>.6
res<-anova_time_1_5$g2
g2_pathways<-rownames(res)[which( 
    res$`p_value_edc_twentynine_day-edc_eight_day`                  <= 0.05 &
    res$`p_value_edc_fifteen_day-edc_eight_day`                     <= 0.05 &
    res$`p_value_edc_twentynine_day-edc_fifteen_day`                 > 0.05 &
    res$average_edc_fifteen_day > res$average_edc_eight_day                 &
    res$average_edc_twentynine_day >   res$average_edc_eight_day        
)]



# group 3:  8d (0.6-0.8), 15d=>8d+.1 and 29d => 15d
res<-anova_time_1_5$g3
g3_pathways<-rownames(res)[which( 
    res$`p_value_edc_twentynine_day-edc_eight_day`                     <= 0.05 &
    res$`p_value_edc_fifteen_day-edc_eight_day`                        <= 0.05 &
    #res$`p_value_edc_twentynine_day-edc_fifteen_day`                   <= 0.05 &
    res$average_edc_fifteen_day > res$average_edc_eight_day                    & 
    res$average_edc_twentynine_day > res$average_edc_eight_day                 &
    res$average_edc_twentynine_day >= res$average_edc_fifteen_day     
)]

# group 4: 8d 15d and 29d > 0.9
res<-anova_time_1_5$g4
g4_pathways<-rownames(res)[which( 
    res$`p_value_edc_fifteen_day-edc_eight_day`        > 0.05 &
    res$`p_value_edc_twentynine_day-edc_eight_day`     > 0.05 &
    res$`p_value_edc_twentynine_day-edc_fifteen_day`   > 0.05 
)]


# group5 8d>0.6 15d>.6 and 29d<.4
res<-anova_time_1_5$g5
g5_pathways<-rownames(res)[which( 
  res$`p_value_edc_twentynine_day-edc_eight_day`                  <= 0.05 &
    res$`p_value_edc_twentynine_day-edc_fifteen_day`                <= 0.05 &
    res$average_edc_twentynine_day < res$average_edc_fifteen_day            &
    res$average_edc_twentynine_day <   res$average_edc_eight_day        
)]

time_significant_pathways<-list(g1_pathways,g2_pathways,g3_pathways,g4_pathways,g5_pathways)
time_edc_groups<-list(g1_edc,g2_edc,g3_edc,g4_edc,g5_edc)

save(anova_time_1_5,time_significant_pathways,time_edc_groups,file = 'outputData/new_time_dose_sensitivity/ANOVA_Times.RData')



