library(ggplot2);library(ggpubr);library(grid);library(facetscales)
# 1. Dose response pattern recognition TG-GATEs Low Middle and --------
rm(list=ls())
load('outputData/CTD_edc_scores_dictionary.RData')
all_vam<-all_vam[,grep(colnames(all_vam),pattern = 'invivo|is_in|comp_names')]
all_vam<-na.omit(all_vam)

delta<-all_vam[,c("TG_GATEs_Rat_invivo_Single_Dose_Low_1_day",
                  "TG_GATEs_Rat_invivo_Single_Dose_Middle_1_day",
                  "TG_GATEs_Rat_invivo_Single_Dose_High_1_day"
)] # Single dose
groups<-list()
# the compounds are divided in to 4 groups based on their proifles
g1<-delta[apply(delta, 1,function(x) x[1] < .4 & x[2] < .4 & x[3] > .6),1:3]
groups[[1]]<-rownames(g1)

g2<-delta[apply(delta, 1,function(x) x[1] < .4 & x[2] > .6 & x[3] > .6),1:3]
groups[[2]]<-rownames(g2)

g3<-delta[apply (delta, 1,function(x) (x[1] >= .6 & x[1] < .8) & x[2] >=  (x[1] + .1) & 
                   x[3] >= x[2] ),1:3]
groups[[3]]<-rownames(g3)

g4<-delta[apply (delta, 1,function(x) x[1] >= .9 & x[2] >=  .9 &  x[3] >= .9),1:3]
groups[[4]]<-rownames(g4)


gr_list<-list()
for (i in 1:length(groups)){
  ind<-groups[[i]]
  lvls<-c( "Low", "Middle","High")
  Tg_single<-all_vam[ind,grep(colnames(all_vam),pattern = 'TG_GATEs_Rat_invivo_Single_Dose|comp_names')]
  colnames(Tg_single)<-gsub(x=colnames(Tg_single),pattern ='TG_GATEs_Rat_invivo_Single_Dose_',replacement = '' )
  colnames(Tg_single)<-gsub(x=colnames(Tg_single),pattern ='_1_day',replacement = '' )
  Tg_single$mesh<-rownames(Tg_single)
  Tg_single<-tidyr::gather(data = Tg_single,key='data_layer',value='class_prob',-c(comp_names,mesh))
  Tg_single$cluster_gr<-paste('Group',i,sep = '_')
  gr_list[[i]]<-Tg_single
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
lin_lvls<-factor(c('No_EDC - No_EDC < EDC','No_EDC < EDC - EDC','EDC < EDC <= EDC','EDC - EDC - EDC'))
final<-do.call(rbind,data)
dose_final<-final;save(dose_final,file='outputData/plot_data/Dose_response.RData')
library(RColorBrewer)
ggplot(final, 
       aes(x = factor(data_layer,levels = c('Low','Middle','High')), 
           y = len, ymin = len-sd, ymax = len+sd,colour=factor(pattern,levels = lin_lvls)))+
  geom_pointrange()+geom_errorbar(width = .2) +geom_point(size = 1.5)+
  scale_color_manual(values=rev(brewer.pal(n = 4, name ='RdYlGn')))+
  xlab('TG_GATEs_Single_Dose_1_day')+ylab('Class Probability')+ theme_minimal()+
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

Tg_rep_single$sec_gr<-'Group'
Tg_rep_single$sec_gr[Tg_rep_single$cluster_gr=='Group_1']<-'g1'
Tg_rep_single$sec_gr[Tg_rep_single$cluster_gr=='Group_2']<-'g1'
Tg_rep_single$sec_gr[Tg_rep_single$cluster_gr=='Group_3']<-'g2'
Tg_rep_single$sec_gr[Tg_rep_single$cluster_gr=='Group_4']<-'g2'

Tg_rep_single$f_gr<-'Group'
Tg_rep_single$f_gr[Tg_rep_single$cluster_gr=='Group_1']<-'f1'
Tg_rep_single$f_gr[Tg_rep_single$cluster_gr=='Group_2']<-'f2'
Tg_rep_single$f_gr[Tg_rep_single$cluster_gr=='Group_3']<-'f1'
Tg_rep_single$f_gr[Tg_rep_single$cluster_gr=='Group_4']<-'f2'

scales_y <- list(
  `g1` = scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .1)),
  `g2` = scale_y_continuous(limits = c(0.55, 1), breaks = seq(0.4, 1, .1))
)
library(RColorBrewer)
ggplot(Tg_rep_single, aes(x=factor(data_layer,levels = lvls), y=class_prob,colour=factor(desc,levels = lin_lvls))) +
  geom_point(show.legend = F)+
  geom_line(aes(group=comp_names),show.legend = T)+
  #guides(col=guide_legend(nrow = 2))+
  labs(colour='Dose Response Groups')+
  scale_color_manual(values=rev(brewer.pal(n = 4, name ='RdYlGn')))+
  #scale_color_gradient(low='yellow',high = 'red')+
  xlab('TG_GATEs_Single_Dose_1_day')+ylab('Class Probability')+ theme_minimal()+
  theme(strip.text = element_blank(),
        legend.text = element_text(size=10),
        legend.position = 'right',axis.text.x = element_text(size=16),
        axis.title = element_text(size = 20),axis.text.y = element_text(size=11),
        panel.spacing = unit(2,'lines')) +
  facet_grid_sc(rows=vars(sec_gr),cols = vars(f_gr),scales = list(y=scales_y))




# 2. preselection and anova for the pathways based on doses ---------------------------
# we preselect the pathways with the GLM more than .7 based on ROC analysis of 
# pathway activation scores and class labels EDCs and Decoys for the training set
load('outputData/patway_glm_coefficients/GLM_univariate_all_scores.RData')
all_data_final<-all_data_final[which(all_data_final$ROC_AUC >0.7),]
Roc_based_preselected_pathways<-all_data_final$pathways[grep(x=all_data_final$network,
                                                      pattern = 'TG_GATEs_Rat_invivo_Single_Dose',ignore.case = F)]
# The common pathways in the list of Low middle and high will be used for ANOVA
load("outputData/class_probilities/integrated_fgsea_allcompounds_results_final_moa.RData")
fgs<-fgs[c("TG_GATEs_Rat_invivo_Single_Dose_Low_1_day",
           "TG_GATEs_Rat_invivo_Single_Dose_Middle_1_day",
           "TG_GATEs_Rat_invivo_Single_Dose_High_1_day")]      #All predicted pathway scores using RWR-FGSEA strategy
names(fgs)<-c('low','middle','high')
common_pathways<-Reduce('intersect',lapply(fgs, function(x)colnames(x$x))) # COMMON PATHWAYS IN LOW MIDDLE AND HIGH
common_pathways<-intersect(common_pathways,Roc_based_preselected_pathways) # Those in the list of preselected list will be used

load('outputData/CTD_edc_scores_dictionary.RData')
all_vam<-all_vam[,grep(colnames(all_vam),pattern = 'invivo|is_in|comp_names')]
all_vam<-na.omit(all_vam)
delta<-all_vam[,c("TG_GATEs_Rat_invivo_Single_Dose_Low_1_day",
                  "TG_GATEs_Rat_invivo_Single_Dose_Middle_1_day",
                  "TG_GATEs_Rat_invivo_Single_Dose_High_1_day"
)] # Single dose class probabilities from elastic net

pas_anova<-function(g_edc){
  edcs<-list();for (i in 1:length(fgs)){edcs[[i]]<-fgs[[i]]$x[g_edc,common_pathways]}
  names(edcs)<-names(fgs)
  
  # for each pathway we do ANOVA with Tukey post test
  res<-matrix(NA, nrow = length(common_pathways), ncol = 6)
  for (i in 1:length(common_pathways)){
    res[i,1]<-mean(edcs$low[,i])
    res[i,2]<-mean(edcs$middle[,i])
    res[i,3]<-mean(edcs$high[,i])
    
    anov_inp<-as.data.frame(list(edcs$low[,i],edcs$middle[,i],edcs$high[,i]))
    colnames(anov_inp)<-c('edc_low','edc_middle','edc_high')
    anov_inp<-stack(anov_inp)
    res.aov<-aov(values ~ ind,data=anov_inp)
    tk<-TukeyHSD(res.aov)
    tk<-as.data.frame(tk$ind)
    res[i,4:6]<-tk$`p adj`
  }
  colnames(res)<-c('average_edc_low','average_edc_middle','average_edc_high',
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

# anova on pathway activation scores
anova_doses_1_4<-lapply(list(g1_edc,g2_edc,g3_edc,g4_edc), function(x)pas_anova(x))
names(anova_doses_1_4)<-c('g1','g2','g3','g4')

# group1 low<.4, middle<.4 and high>.6
res<-anova_doses_1_4$g1
  g1_pathways<-rownames(res)[which( 
      res$`p_value_edc_high-edc_low`                  <= 0.05 &
      res$`p_value_edc_high-edc_middle`               <= 0.05 &
      res$`p_value_edc_middle-edc_low`                 > 0.05 &  
      res$average_edc_high > res$average_edc_middle           &
      res$average_edc_high > res$average_edc_low
)]

# group2 low<.4 middle>.6 and high>.6
res<-anova_doses_1_4$g2
g2_pathways<-rownames(res)[which( 
    res$`p_value_edc_high-edc_low`                    <= 0.05 &
    res$`p_value_edc_middle-edc_low`                  <= 0.05 &
    res$`p_value_edc_high-edc_middle`                  > 0.05 &
    res$average_edc_middle > res$average_edc_low              &
    res$average_edc_high >   res$average_edc_low      
    )]
  
# group 3  low (0.6-0.8), midle=>low+.1 and high => middle
res<-anova_doses_1_4$g3
g3_pathways<-rownames(res)[which( 
    res$`p_value_edc_high-edc_low`                    <= 0.05 &
    res$`p_value_edc_middle-edc_low`                  <= 0.05 &
 
    res$average_edc_high >= res$average_edc_middle            & 
    res$average_edc_high > res$average_edc_low                &
    res$average_edc_middle > res$average_edc_low         
)]

# group 4 low,middle and high >0.9
res<-anova_doses_1_4$g4
g4_pathways<-rownames(res)[which( 
res$`p_value_edc_middle-edc_low`   > 0.05 &
res$`p_value_edc_high-edc_low`     > 0.05 &
res$`p_value_edc_high-edc_middle`  > 0.05 
)]

dose_significant_pathways<-list(g1_pathways,g2_pathways,g3_pathways,g4_pathways)
dose_edc_groups<-list(g1_edc,g2_edc,g3_edc,g4_edc)

save(anova_doses_1_4,dose_significant_pathways,dose_edc_groups,file = 'outputData/new_time_dose_sensitivity/ANOVA_DOSEs.RData')


