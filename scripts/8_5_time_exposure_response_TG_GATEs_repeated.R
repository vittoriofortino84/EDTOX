library(ggplot2);library(ggpubr);library(grid);library(facetscales)
# 1. Manual Dose response pattern recognition TG-GATEs Low Middle and --------
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

 g3<-delta[apply(delta, 1,function(x) x[1] > .6  & x[2] > .6  & x[3] < .4),1:3]
 groups[[3]]<-rownames(g3)
# 
# 
# g4<-delta[apply (delta, 1,function(x) x[2] > 0.6 & x[2] > (x[1] + 0.5) & (x[3] + 0.5) < x[2]),1:3]
# groups[[4]]<-rownames(g4)


g4<-delta[apply (delta, 1,function(x) (x[1] >= .6 & x[1] < .8) & x[2] >=  (x[1] + .1) & 
                   x[3] >= x[2] ),1:3]
groups[[4]]<-rownames(g4)

g5<-delta[apply (delta, 1,function(x) x[1] >= 0.9 & x[2] >=  0.9 &  x[3] >= 0.9),1:3]
groups[[5]]<-rownames(g5)

#excl<-unlist(groups);groups[[8]]<-rownames(delta)[-which(rownames(delta) %in% excl)][1:100]
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
data[[1]]$pattern<-'No-EDC < No-EDC < EDC '
data[[2]]$pattern<-'No-EDC < EDC = EDC'

data[[3]]$pattern<-'EDC > EDC > No-EDC'
#data[[4]]$pattern<-'No-EDC < EDC > no-EDC'


data[[4]]$pattern<-'EDC < EDC < EDC'
data[[5]]$pattern<-'EDC = EDC = EDC'

lin_lvls<-factor(c('No-EDC < No-EDC < EDC ','No-EDC < EDC = EDC','EDC > EDC > No-EDC','EDC < EDC < EDC','EDC = EDC = EDC'))

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
Tg_rep_single$desc[Tg_rep_single$cluster_gr=='Group_1']<-'No-EDC < No-EDC < EDC '
Tg_rep_single$desc[Tg_rep_single$cluster_gr=='Group_2']<-'No-EDC < EDC = EDC'
Tg_rep_single$desc[Tg_rep_single$cluster_gr=='Group_3']<-'EDC > EDC > no-EDC '
# Tg_rep_single$desc[Tg_rep_single$cluster_gr=='Group_6']<-'EDC > EDC > EDC'
# Tg_rep_single$desc[Tg_rep_single$cluster_gr=='Group_7']<-'EDC > no-EDC < EDC'
Tg_rep_single$desc[Tg_rep_single$cluster_gr=='Group_4']<-'no-EDC < EDC > no-EDC'
Tg_rep_single$desc[Tg_rep_single$cluster_gr=='Group_5']<-'EDC < EDC < EDC'
Tg_rep_single$desc[Tg_rep_single$cluster_gr=='Group_6']<-'EDC = EDC = EDC'



lin_lvls<-factor(c('No-EDC < No-EDC < EDC ','No-EDC < EDC = EDC',
                   'EDC > EDC > no-EDC ','no-EDC < EDC > no-EDC','EDC < EDC < EDC','EDC = EDC = EDC'
                   ))

Tg_rep_single$sec_gr<-'Group'
Tg_rep_single$sec_gr[Tg_rep_single$cluster_gr=='Group_1']<-'g1'
Tg_rep_single$sec_gr[Tg_rep_single$cluster_gr=='Group_2']<-'g1'
Tg_rep_single$sec_gr[Tg_rep_single$cluster_gr=='Group_3']<-'g2'
Tg_rep_single$sec_gr[Tg_rep_single$cluster_gr=='Group_4']<-'g2'
Tg_rep_single$sec_gr[Tg_rep_single$cluster_gr=='Group_5']<-'g3'
# Tg_rep_single$sec_gr[Tg_rep_single$cluster_gr=='Group_6']<-'g3'
# Tg_rep_single$sec_gr[Tg_rep_single$cluster_gr=='Group_7']<-'g4'
Tg_rep_single$sec_gr[Tg_rep_single$cluster_gr=='Group_6']<-'g3'

Tg_rep_single$f_gr<-'Group'
Tg_rep_single$f_gr[Tg_rep_single$cluster_gr=='Group_1']<-'f1'
Tg_rep_single$f_gr[Tg_rep_single$cluster_gr=='Group_2']<-'f2'
Tg_rep_single$f_gr[Tg_rep_single$cluster_gr=='Group_3']<-'f1'
Tg_rep_single$f_gr[Tg_rep_single$cluster_gr=='Group_4']<-'f2'
Tg_rep_single$f_gr[Tg_rep_single$cluster_gr=='Group_5']<-'f1'
# Tg_rep_single$f_gr[Tg_rep_single$cluster_gr=='Group_6']<-'f2'
# Tg_rep_single$f_gr[Tg_rep_single$cluster_gr=='Group_7']<-'f1'
Tg_rep_single$f_gr[Tg_rep_single$cluster_gr=='Group_6']<-'f2'



scales_y <- list(
  `g1` = scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .1)),
  `g2` = scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .1)),
  `g3` = scale_y_continuous(limits = c(0.55, 1), breaks = seq(0.55, 1, .05))

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






# 2_1. ANOVA of pathway scores scenario 1-----------------------------------

rm(list=ls())
load("outputData/class_probilities/integrated_fgsea_allcompounds_results_final_moa.RData")

fgs<-fgs[c("TG_GATEs_Rat_invivo_Repeated_Dose_8_days",
           "TG_GATEs_Rat_invivo_Repeated_Dose_15_days",
           "TG_GATEs_Rat_invivo_Repeated_Dose_29_days")] #All predicted pathway scores using RWR-FGSEA strategy
names(fgs)<-c('eight_day','fifteen_day','twentynine_day')
common_pathways<-Reduce('intersect',lapply(fgs, function(x)colnames(x$x))) # COMMON PATHWAYS IN LOW MIDDLE AND HIGH

load('outputData/CTD_edc_scores_dictionary.RData')
all_vam<-all_vam[,grep(colnames(all_vam),pattern = 'invivo|is_in|comp_names')]
all_vam<-na.omit(all_vam)
delta<-all_vam[,c("TG_GATEs_Rat_invivo_Repeated_Dose_8_days",
                  "TG_GATEs_Rat_invivo_Repeated_Dose_15_days",
                  "TG_GATEs_Rat_invivo_Repeated_Dose_29_days"
)] # Single dose class probabilities from elastic net

# selection of the compounds as edcs and decoys
g_edc<-delta[apply(delta, 1,function(x) x[1] < .4 & x[2] < .4 & x[3] > .6),1:3];g_edc<-rownames(g_edc)

g_decoys<-delta[apply (delta, 1,function(x) x[1] <  0.4 & x[2] <   0.4 &  x[3] <  0.4),1:3];g_decoys<-rownames(g_decoys)

set.seed(123)
boot_strap<-decoys_ind<-list()
for (j in 1:100){
  print(j)
  
  
  # Decoys
  g_dec<-g_decoys[sample(1:length(g_decoys),length(g_edc))]
  decoys_ind[[j]]<-g_dec
  #intersetion of the pathway scores for low middle and high TG-GATEs single dose
  
  edcs<-list();for (i in 1:length(fgs)){edcs[[i]]<-fgs[[i]]$x[g_edc,common_pathways]}
  decs<-list();for (i in 1:length(fgs)){decs[[i]]<-fgs[[i]]$x[g_dec,common_pathways]}
  names(edcs)<-names(decs)<-names(fgs)
  
  
  # for each pathway we do ANOVA with Tukey test
  res<-matrix(NA, nrow = length(common_pathways), ncol = 21)
  for (i in 1:length(common_pathways)){
    res[i,1]<-mean(edcs$eight_day[,i])
    res[i,2]<-mean(edcs$fifteen_day[,i])
    res[i,3]<-mean(edcs$twentynine_day[,i])
    
    res[i,4]<-mean(decs$eight_day[,i])
    res[i,5]<-mean(decs$fifteen_day[,i])
    res[i,6]<-mean(decs$twentynine_day[,i])
    
    anov_inp<-as.data.frame(list(edcs$eight_day[,i],edcs$fifteen_day[,i],edcs$twentynine_day[,i],
                                 decs$eight_day[,i],decs$fifteen_day[,i],decs$twentynine_day[,i]))
    colnames(anov_inp)<-c('edc_eight_day','edc_fifteen_day','edc_twentynine_day',
                          'dec_eight_day','dec_fifteen_day','dec_twentynine_day')
    anov_inp<-stack(anov_inp)
    res.aov<-aov(values ~ ind,data=anov_inp)
    #anova_summary<-summary(res.aov)
    tk<-TukeyHSD(res.aov)
    tk<-as.data.frame(tk$ind)
    res[i,7:21]<-tk$`p adj`
  }
  colnames(res)<-c('average_edc_eight_day','average_edc_fifteen_day','average_edc_twentynine_day',
                   'average_dec_eight_day','average_dec_fifteen_day','average_dec_twentynine_day',
                   paste('p_value',rownames(tk),sep='_')
  )
  rownames(res)<-common_pathways
  res<-as.data.frame(res)
  boot_strap[[j]]<-res
}

save(edcs,decoys_ind,decs,boot_strap,file = 'outputData/bootstrap_ANOVA_Time_1.RData')




# 2_2. ANOVA of pathway scores scenario 2-----------------------------------

rm(list=ls())
load("outputData/class_probilities/integrated_fgsea_allcompounds_results_final_moa.RData")

fgs<-fgs[c("TG_GATEs_Rat_invivo_Repeated_Dose_8_days",
           "TG_GATEs_Rat_invivo_Repeated_Dose_15_days",
           "TG_GATEs_Rat_invivo_Repeated_Dose_29_days")] #All predicted pathway scores using RWR-FGSEA strategy
names(fgs)<-c('eight_day','fifteen_day','twentynine_day')
common_pathways<-Reduce('intersect',lapply(fgs, function(x)colnames(x$x))) # COMMON PATHWAYS IN LOW MIDDLE AND HIGH

load('outputData/CTD_edc_scores_dictionary.RData')
all_vam<-all_vam[,grep(colnames(all_vam),pattern = 'invivo|is_in|comp_names')]
all_vam<-na.omit(all_vam)
delta<-all_vam[,c("TG_GATEs_Rat_invivo_Repeated_Dose_8_days",
                  "TG_GATEs_Rat_invivo_Repeated_Dose_15_days",
                  "TG_GATEs_Rat_invivo_Repeated_Dose_29_days"
)] # Single dose class probabilities from elastic net

# selection of the compounds as edcs and decoys
g_edc<-delta[apply(delta, 1,function(x) x[1] <.4  & x[2] > .6  & x[3] > .6),1:3];g_edc<-rownames(g_edc)

g_decoys<-delta[apply (delta, 1,function(x) x[1] <  0.4 & x[2] <   0.4 &  x[3] <  0.4),1:3];g_decoys<-rownames(g_decoys)

set.seed(123)
boot_strap<-decoys_ind<-list()
for (j in 1:100){
  print(j)
  
  
  # Decoys
  g_dec<-g_decoys[sample(1:length(g_decoys),length(g_edc))]
  decoys_ind[[j]]<-g_dec
  #intersetion of the pathway scores for low middle and high TG-GATEs single dose
  
  edcs<-list();for (i in 1:length(fgs)){edcs[[i]]<-fgs[[i]]$x[g_edc,common_pathways]}
  decs<-list();for (i in 1:length(fgs)){decs[[i]]<-fgs[[i]]$x[g_dec,common_pathways]}
  names(edcs)<-names(decs)<-names(fgs)
  
  
  # for each pathway we do ANOVA with Tukey test
  res<-matrix(NA, nrow = length(common_pathways), ncol = 21)
  for (i in 1:length(common_pathways)){
    res[i,1]<-mean(edcs$eight_day[,i])
    res[i,2]<-mean(edcs$fifteen_day[,i])
    res[i,3]<-mean(edcs$twentynine_day[,i])
    
    res[i,4]<-mean(decs$eight_day[,i])
    res[i,5]<-mean(decs$fifteen_day[,i])
    res[i,6]<-mean(decs$twentynine_day[,i])
    
    anov_inp<-as.data.frame(list(edcs$eight_day[,i],edcs$fifteen_day[,i],edcs$twentynine_day[,i],
                                 decs$eight_day[,i],decs$fifteen_day[,i],decs$twentynine_day[,i]))
    colnames(anov_inp)<-c('edc_eight_day','edc_fifteen_day','edc_twentynine_day',
                          'dec_eight_day','dec_fifteen_day','dec_twentynine_day')
    anov_inp<-stack(anov_inp)
    res.aov<-aov(values ~ ind,data=anov_inp)
    #anova_summary<-summary(res.aov)
    tk<-TukeyHSD(res.aov)
    tk<-as.data.frame(tk$ind)
    res[i,7:21]<-tk$`p adj`
  }
  colnames(res)<-c('average_edc_eight_day','average_edc_fifteen_day','average_edc_twentynine_day',
                   'average_dec_eight_day','average_dec_fifteen_day','average_dec_twentynine_day',
                   paste('p_value',rownames(tk),sep='_')
  )
  rownames(res)<-common_pathways
  res<-as.data.frame(res)
  boot_strap[[j]]<-res
}

save(edcs,decoys_ind,decs,boot_strap,file = 'outputData/bootstrap_ANOVA_Time_2.RData')




# 2_3. ANOVA of pathway scores scenario 3-----------------------------------

rm(list=ls())
load("outputData/class_probilities/integrated_fgsea_allcompounds_results_final_moa.RData")

fgs<-fgs[c("TG_GATEs_Rat_invivo_Repeated_Dose_8_days",
           "TG_GATEs_Rat_invivo_Repeated_Dose_15_days",
           "TG_GATEs_Rat_invivo_Repeated_Dose_29_days")] #All predicted pathway scores using RWR-FGSEA strategy
names(fgs)<-c('eight_day','fifteen_day','twentynine_day')
common_pathways<-Reduce('intersect',lapply(fgs, function(x)colnames(x$x))) # COMMON PATHWAYS IN LOW MIDDLE AND HIGH

load('outputData/CTD_edc_scores_dictionary.RData')
all_vam<-all_vam[,grep(colnames(all_vam),pattern = 'invivo|is_in|comp_names')]
all_vam<-na.omit(all_vam)
delta<-all_vam[,c("TG_GATEs_Rat_invivo_Repeated_Dose_8_days",
                  "TG_GATEs_Rat_invivo_Repeated_Dose_15_days",
                  "TG_GATEs_Rat_invivo_Repeated_Dose_29_days"
)] # Single dose class probabilities from elastic net

# selection of the compounds as edcs and decoys
g_edc<-delta[apply(delta, 1,function(x) x[1] > .6  & x[2] > .6  & x[3] < .4),1:3];g_edc<-rownames(g_edc)

g_decoys<-delta[apply (delta, 1,function(x) x[1] <  0.4 & x[2] <   0.4 &  x[3] <  0.4),1:3];g_decoys<-rownames(g_decoys)

set.seed(123)
boot_strap<-decoys_ind<-list()
for (j in 1:100){
  print(j)
  
  
  # Decoys
  g_dec<-g_decoys[sample(1:length(g_decoys),length(g_edc))]
  decoys_ind[[j]]<-g_dec
  #intersetion of the pathway scores for low middle and high TG-GATEs single dose
  
  edcs<-list();for (i in 1:length(fgs)){edcs[[i]]<-fgs[[i]]$x[g_edc,common_pathways]}
  decs<-list();for (i in 1:length(fgs)){decs[[i]]<-fgs[[i]]$x[g_dec,common_pathways]}
  names(edcs)<-names(decs)<-names(fgs)
  
  
  # for each pathway we do ANOVA with Tukey test
  res<-matrix(NA, nrow = length(common_pathways), ncol = 21)
  for (i in 1:length(common_pathways)){
    res[i,1]<-mean(edcs$eight_day[,i])
    res[i,2]<-mean(edcs$fifteen_day[,i])
    res[i,3]<-mean(edcs$twentynine_day[,i])
    
    res[i,4]<-mean(decs$eight_day[,i])
    res[i,5]<-mean(decs$fifteen_day[,i])
    res[i,6]<-mean(decs$twentynine_day[,i])
    
    anov_inp<-as.data.frame(list(edcs$eight_day[,i],edcs$fifteen_day[,i],edcs$twentynine_day[,i],
                                 decs$eight_day[,i],decs$fifteen_day[,i],decs$twentynine_day[,i]))
    colnames(anov_inp)<-c('edc_eight_day','edc_fifteen_day','edc_twentynine_day',
                          'dec_eight_day','dec_fifteen_day','dec_twentynine_day')
    anov_inp<-stack(anov_inp)
    res.aov<-aov(values ~ ind,data=anov_inp)
    #anova_summary<-summary(res.aov)
    tk<-TukeyHSD(res.aov)
    tk<-as.data.frame(tk$ind)
    res[i,7:21]<-tk$`p adj`
  }
  colnames(res)<-c('average_edc_eight_day','average_edc_fifteen_day','average_edc_twentynine_day',
                   'average_dec_eight_day','average_dec_fifteen_day','average_dec_twentynine_day',
                   paste('p_value',rownames(tk),sep='_')
  )
  rownames(res)<-common_pathways
  res<-as.data.frame(res)
  boot_strap[[j]]<-res
}

save(edcs,decoys_ind,decs,boot_strap,file = 'outputData/bootstrap_ANOVA_Time_3.RData')




# 2_4. ANOVA of pathway scores scenario 4-----------------------------------

rm(list=ls())
load("outputData/class_probilities/integrated_fgsea_allcompounds_results_final_moa.RData")

fgs<-fgs[c("TG_GATEs_Rat_invivo_Repeated_Dose_8_days",
           "TG_GATEs_Rat_invivo_Repeated_Dose_15_days",
           "TG_GATEs_Rat_invivo_Repeated_Dose_29_days")] #All predicted pathway scores using RWR-FGSEA strategy
names(fgs)<-c('eight_day','fifteen_day','twentynine_day')
common_pathways<-Reduce('intersect',lapply(fgs, function(x)colnames(x$x))) # COMMON PATHWAYS IN LOW MIDDLE AND HIGH

load('outputData/CTD_edc_scores_dictionary.RData')
all_vam<-all_vam[,grep(colnames(all_vam),pattern = 'invivo|is_in|comp_names')]
all_vam<-na.omit(all_vam)
delta<-all_vam[,c("TG_GATEs_Rat_invivo_Repeated_Dose_8_days",
                  "TG_GATEs_Rat_invivo_Repeated_Dose_15_days",
                  "TG_GATEs_Rat_invivo_Repeated_Dose_29_days"
)] # Single dose class probabilities from elastic net

# selection of the compounds as edcs and decoys
g_edc<-delta[apply(delta, 1,function(x) x[2] > 0.6 & x[2] > (x[1] + 0.5) & (x[3] + 0.5) < x[2]),1:3];g_edc<-rownames(g_edc)

g_decoys<-delta[apply (delta, 1,function(x) x[1] <  0.4 & x[2] <   0.4 &  x[3] <  0.4),1:3];g_decoys<-rownames(g_decoys)

set.seed(123)
boot_strap<-decoys_ind<-list()
for (j in 1:100){
  print(j)
  
  
  # Decoys
  g_dec<-g_decoys[sample(1:length(g_decoys),length(g_edc))]
  decoys_ind[[j]]<-g_dec
  #intersetion of the pathway scores for low middle and high TG-GATEs single dose
  
  edcs<-list();for (i in 1:length(fgs)){edcs[[i]]<-fgs[[i]]$x[g_edc,common_pathways]}
  decs<-list();for (i in 1:length(fgs)){decs[[i]]<-fgs[[i]]$x[g_dec,common_pathways]}
  names(edcs)<-names(decs)<-names(fgs)
  
  
  # for each pathway we do ANOVA with Tukey test
  res<-matrix(NA, nrow = length(common_pathways), ncol = 21)
  for (i in 1:length(common_pathways)){
    res[i,1]<-mean(edcs$eight_day[,i])
    res[i,2]<-mean(edcs$fifteen_day[,i])
    res[i,3]<-mean(edcs$twentynine_day[,i])
    
    res[i,4]<-mean(decs$eight_day[,i])
    res[i,5]<-mean(decs$fifteen_day[,i])
    res[i,6]<-mean(decs$twentynine_day[,i])
    
    anov_inp<-as.data.frame(list(edcs$eight_day[,i],edcs$fifteen_day[,i],edcs$twentynine_day[,i],
                                 decs$eight_day[,i],decs$fifteen_day[,i],decs$twentynine_day[,i]))
    colnames(anov_inp)<-c('edc_eight_day','edc_fifteen_day','edc_twentynine_day',
                          'dec_eight_day','dec_fifteen_day','dec_twentynine_day')
    anov_inp<-stack(anov_inp)
    res.aov<-aov(values ~ ind,data=anov_inp)
    #anova_summary<-summary(res.aov)
    tk<-TukeyHSD(res.aov)
    tk<-as.data.frame(tk$ind)
    res[i,7:21]<-tk$`p adj`
  }
  colnames(res)<-c('average_edc_eight_day','average_edc_fifteen_day','average_edc_twentynine_day',
                   'average_dec_eight_day','average_dec_fifteen_day','average_dec_twentynine_day',
                   paste('p_value',rownames(tk),sep='_')
  )
  rownames(res)<-common_pathways
  res<-as.data.frame(res)
  boot_strap[[j]]<-res
}

save(edcs,decoys_ind,decs,boot_strap,file = 'outputData/bootstrap_ANOVA_Time_4.RData')




# 2_5. ANOVA of pathway scores scenario 5-----------------------------------

rm(list=ls())
load("outputData/class_probilities/integrated_fgsea_allcompounds_results_final_moa.RData")

fgs<-fgs[c("TG_GATEs_Rat_invivo_Repeated_Dose_8_days",
           "TG_GATEs_Rat_invivo_Repeated_Dose_15_days",
           "TG_GATEs_Rat_invivo_Repeated_Dose_29_days")] #All predicted pathway scores using RWR-FGSEA strategy
names(fgs)<-c('eight_day','fifteen_day','twentynine_day')
common_pathways<-Reduce('intersect',lapply(fgs, function(x)colnames(x$x))) # COMMON PATHWAYS IN LOW MIDDLE AND HIGH

load('outputData/CTD_edc_scores_dictionary.RData')
all_vam<-all_vam[,grep(colnames(all_vam),pattern = 'invivo|is_in|comp_names')]
all_vam<-na.omit(all_vam)
delta<-all_vam[,c("TG_GATEs_Rat_invivo_Repeated_Dose_8_days",
                  "TG_GATEs_Rat_invivo_Repeated_Dose_15_days",
                  "TG_GATEs_Rat_invivo_Repeated_Dose_29_days"
)] # Single dose class probabilities from elastic net

# selection of the compounds as edcs and decoys
g_edc<-delta[apply(delta, 1,function(x) (x[1] >= .6 & x[1] < .8) & x[2] >=  (x[1] + .1) & 
                     x[3] >= x[2] ),1:3];g_edc<-rownames(g_edc)

g_decoys<-delta[apply (delta, 1,function(x) x[1] <  0.4 & x[2] <   0.4 &  x[3] <  0.4),1:3];g_decoys<-rownames(g_decoys)

set.seed(123)
boot_strap<-decoys_ind<-list()
for (j in 1:100){
  print(j)
  
  
  # Decoys
  g_dec<-g_decoys[sample(1:length(g_decoys),length(g_edc))]
  decoys_ind[[j]]<-g_dec
  #intersetion of the pathway scores for low middle and high TG-GATEs single dose
  
  edcs<-list();for (i in 1:length(fgs)){edcs[[i]]<-fgs[[i]]$x[g_edc,common_pathways]}
  decs<-list();for (i in 1:length(fgs)){decs[[i]]<-fgs[[i]]$x[g_dec,common_pathways]}
  names(edcs)<-names(decs)<-names(fgs)
  
  
  # for each pathway we do ANOVA with Tukey test
  res<-matrix(NA, nrow = length(common_pathways), ncol = 21)
  for (i in 1:length(common_pathways)){
    res[i,1]<-mean(edcs$eight_day[,i])
    res[i,2]<-mean(edcs$fifteen_day[,i])
    res[i,3]<-mean(edcs$twentynine_day[,i])
    
    res[i,4]<-mean(decs$eight_day[,i])
    res[i,5]<-mean(decs$fifteen_day[,i])
    res[i,6]<-mean(decs$twentynine_day[,i])
    
    anov_inp<-as.data.frame(list(edcs$eight_day[,i],edcs$fifteen_day[,i],edcs$twentynine_day[,i],
                                 decs$eight_day[,i],decs$fifteen_day[,i],decs$twentynine_day[,i]))
    colnames(anov_inp)<-c('edc_eight_day','edc_fifteen_day','edc_twentynine_day',
                          'dec_eight_day','dec_fifteen_day','dec_twentynine_day')
    anov_inp<-stack(anov_inp)
    res.aov<-aov(values ~ ind,data=anov_inp)
    #anova_summary<-summary(res.aov)
    tk<-TukeyHSD(res.aov)
    tk<-as.data.frame(tk$ind)
    res[i,7:21]<-tk$`p adj`
  }
  colnames(res)<-c('average_edc_eight_day','average_edc_fifteen_day','average_edc_twentynine_day',
                   'average_dec_eight_day','average_dec_fifteen_day','average_dec_twentynine_day',
                   paste('p_value',rownames(tk),sep='_')
  )
  rownames(res)<-common_pathways
  res<-as.data.frame(res)
  boot_strap[[j]]<-res
}

save(edcs,decoys_ind,decs,boot_strap,file = 'outputData/bootstrap_ANOVA_Time_5.RData')




# 2_6. ANOVA of pathway scores scenario 6-----------------------------------
#THe compound is EDC at HIGH

rm(list=ls())
load("outputData/class_probilities/integrated_fgsea_allcompounds_results_final_moa.RData")

fgs<-fgs[c("TG_GATEs_Rat_invivo_Repeated_Dose_8_days",
           "TG_GATEs_Rat_invivo_Repeated_Dose_15_days",
           "TG_GATEs_Rat_invivo_Repeated_Dose_29_days")] #All predicted pathway scores using RWR-FGSEA strategy
names(fgs)<-c('eight_day','fifteen_day','twentynine_day')
common_pathways<-Reduce('intersect',lapply(fgs, function(x)colnames(x$x))) # COMMON PATHWAYS IN LOW MIDDLE AND HIGH

load('outputData/CTD_edc_scores_dictionary.RData')
all_vam<-all_vam[,grep(colnames(all_vam),pattern = 'invivo|is_in|comp_names')]
all_vam<-na.omit(all_vam)
delta<-all_vam[,c("TG_GATEs_Rat_invivo_Repeated_Dose_8_days",
                  "TG_GATEs_Rat_invivo_Repeated_Dose_15_days",
                  "TG_GATEs_Rat_invivo_Repeated_Dose_29_days"
)] # Single dose class probabilities from elastic net

# selection of the compounds as edcs and decoys
g_edc<-delta[apply(delta, 1,function(x) x[1] >= 0.9 & x[2] >=  0.9 &  x[3] >= 0.9),1:3];g_edc<-rownames(g_edc)

g_decoys<-delta[apply (delta, 1,function(x) x[1] <  0.4 & x[2] <   0.4 &  x[3] <  0.4),1:3];g_decoys<-rownames(g_decoys)

set.seed(123)
boot_strap<-decoys_ind<-list()
for (j in 1:100){
  print(j)
  
  
  # Decoys
  g_dec<-g_decoys[sample(1:length(g_decoys),length(g_edc))]
  decoys_ind[[j]]<-g_dec
  #intersetion of the pathway scores for low middle and high TG-GATEs single dose
  
  edcs<-list();for (i in 1:length(fgs)){edcs[[i]]<-fgs[[i]]$x[g_edc,common_pathways]}
  decs<-list();for (i in 1:length(fgs)){decs[[i]]<-fgs[[i]]$x[g_dec,common_pathways]}
  names(edcs)<-names(decs)<-names(fgs)
  
  
  # for each pathway we do ANOVA with Tukey test
  res<-matrix(NA, nrow = length(common_pathways), ncol = 21)
  for (i in 1:length(common_pathways)){
    res[i,1]<-mean(edcs$eight_day[,i])
    res[i,2]<-mean(edcs$fifteen_day[,i])
    res[i,3]<-mean(edcs$twentynine_day[,i])
    
    res[i,4]<-mean(decs$eight_day[,i])
    res[i,5]<-mean(decs$fifteen_day[,i])
    res[i,6]<-mean(decs$twentynine_day[,i])
    
    anov_inp<-as.data.frame(list(edcs$eight_day[,i],edcs$fifteen_day[,i],edcs$twentynine_day[,i],
                                 decs$eight_day[,i],decs$fifteen_day[,i],decs$twentynine_day[,i]))
    colnames(anov_inp)<-c('edc_eight_day','edc_fifteen_day','edc_twentynine_day',
                          'dec_eight_day','dec_fifteen_day','dec_twentynine_day')
    anov_inp<-stack(anov_inp)
    res.aov<-aov(values ~ ind,data=anov_inp)
    #anova_summary<-summary(res.aov)
    tk<-TukeyHSD(res.aov)
    tk<-as.data.frame(tk$ind)
    res[i,7:21]<-tk$`p adj`
  }
  colnames(res)<-c('average_edc_eight_day','average_edc_fifteen_day','average_edc_twentynine_day',
                   'average_dec_eight_day','average_dec_fifteen_day','average_dec_twentynine_day',
                   paste('p_value',rownames(tk),sep='_')
  )
  rownames(res)<-common_pathways
  res<-as.data.frame(res)
  boot_strap[[j]]<-res
}

save(edcs,decoys_ind,decs,boot_strap,file = 'outputData/bootstrap_ANOVA_Time_6.RData')





# 3. Analysis of ANOVA pathways scores across dose in TG-GATEs ---------------------------------------------------
# increase scenarios in pathway scores
#### scenario 1
library(tidyr)
rm(list = ls())
load('outputData/bootstrap_ANOVA_Time_1.RData')
final_res<-matrix(0, nrow = length(boot_strap), ncol = nrow(boot_strap[[1]]))
colnames(final_res)<-rownames(boot_strap[[1]])
for (i in 1:nrow(final_res)){
  res<-boot_strap[[i]]
  ind<-which( 
    res$`p_value_dec_twentynine_day-edc_twentynine_day`                   <= 0.05 & 
      res$`p_value_dec_fifteen_day-edc_twentynine_day`               <= 0.05 &
      res$`p_value_dec_eight_day-edc_twentynine_day`                  <= 0.05 &  
      res$`p_value_edc_twentynine_day-edc_eight_day`          <= 0.05 &
      res$`p_value_edc_twentynine_day-edc_fifteen_day`               <= 0.05 &
      res$average_edc_twentynine_day > res$average_dec_twentynine_day             & 
      res$average_edc_twentynine_day > res$average_dec_fifteen_day           &
      res$average_edc_twentynine_day > res$average_dec_eight_day              &
      res$average_edc_twentynine_day > res$average_edc_fifteen_day           &
      res$average_edc_twentynine_day > res$average_edc_eight_day
  )
  final_res[i,ind]<-1
}

sel<-colnames(final_res)[which(colSums(final_res)>=100)]
mean_edcs<-as.data.frame(lapply(edcs, function(x)apply(x[,sel],2,mean)))
#matplot(t(mean_edcs),type = 'l',lty = 4)

mean_edcs$type<-'none-EDC <= none-EDC < EDC'
mean_edcs$pathway<-rownames(mean_edcs)
data_1<-gather(mean_edcs,key='Dose',value='activation_score',-c(type,pathway))

#### scenario 2
load('outputData/bootstrap_ANOVA_Time_2.RData')
final_res<-matrix(0, nrow = length(boot_strap), ncol = nrow(boot_strap[[1]]))
colnames(final_res)<-rownames(boot_strap[[1]])
for (i in 1:nrow(final_res)){
  res<-boot_strap[[i]]
  ind<-which( 
      res$`p_value_dec_twentynine_day-edc_twentynine_day`                 <= 0.05 &
      res$`p_value_dec_fifteen_day-edc_fifteen_day`             <= 0.05 &
      res$`p_value_edc_twentynine_day-edc_eight_day`                  <= 0.05 &
      res$`p_value_edc_fifteen_day-edc_eight_day`                <= 0.05 &
      res$average_edc_twentynine_day > res$average_dec_twentynine_day             & 
      res$average_edc_twentynine_day > res$average_dec_fifteen_day           &
      res$average_edc_twentynine_day > res$average_dec_eight_day              &
      res$average_edc_fifteen_day > res$average_dec_twentynine_day           & 
      res$average_edc_fifteen_day > res$average_dec_fifteen_day         &
      res$average_edc_fifteen_day > res$average_dec_eight_day            &
      res$average_edc_fifteen_day > res$average_edc_eight_day            &
      res$average_edc_twentynine_day >   res$average_edc_eight_day      
  )
  final_res[i,ind]<-1
}

sel<-colnames(final_res)[which(colSums(final_res)>=100)]
mean_edcs<-as.data.frame(lapply(edcs, function(x)apply(x[,sel],2,mean)))
#matplot(t(mean_edcs),type = 'l',lty = 4)


mean_edcs$type<-'non-EDC < EDC <= EDC'
mean_edcs$pathway<-rownames(mean_edcs)
data_2<-gather(mean_edcs,key='Dose',value='activation_score',-c(type,pathway))

#### scenario 5
#ind<-which( apply(res[,7:21], 1,function(p_values)all(p_values<=0.05)) %in% TRUE) #all p-values shoule be significant
load('outputData/bootstrap_ANOVA_Time_5.RData')
final_res<-matrix(0, nrow = length(boot_strap), ncol = nrow(boot_strap[[1]]))
colnames(final_res)<-rownames(boot_strap[[1]])
for (i in 1:nrow(final_res)){
  res<-boot_strap[[i]]
  ind<-which(
      res$`p_value_dec_twentynine_day-edc_twentynine_day`                <= 0.05 &
      res$`p_value_dec_fifteen_day-edc_fifteen_day`                      <= 0.05 &
     # res$`p_value_dec_eight_day-edc_eight_day`                          <= 0.05 &
        
      res$`p_value_edc_twentynine_day-edc_eight_day`                     <= 0.05 &
      res$`p_value_edc_fifteen_day-edc_eight_day`                        <= 0.05 &
     # res$`p_value_edc_twentynine_day-edc_fifteen_day`                   <= 0.05 &
 
      res$average_edc_fifteen_day > res$average_edc_eight_day                &         
      res$average_edc_twentynine_day >= res$average_edc_fifteen_day 
  ) #
  final_res[i,ind]<-1
}

sel<-colnames(final_res)[which(colSums(final_res)>=100)]
mean_edcs<-as.data.frame(lapply(edcs, function(x)apply(x[,sel],2,mean)))
matplot(t(mean_edcs),type = 'l',lty = 4)
mean_edcs$type<-'EDC < EDC <= EDC'
mean_edcs$pathway<-rownames(mean_edcs)
data_3<-tidyr::gather(mean_edcs,key='Dose',value='activation_score',-c(type,pathway))
time_biomarker<-data_3;save(time_biomarker,file='outputData/plot_data/time_bio_marker_scenario5.RData')

data<-do.call(rbind,list(data_1,data_2,data_3))

#library(RColorBrewer)
ggplot(data, aes(x=factor(Dose), y=activation_score,colour=factor(pathway))) +
  geom_point(show.legend = F)+
  geom_line(aes(group=pathway),show.legend = F)+
  #guides(col=guide_legend(nrow = 2))+
  labs(colour='Pathways')+
  #scale_color_manual(values=rev(brewer.pal(n = 4, name ='RdYlGn')))+
  #scale_color_gradient(low='yellow',high = 'red')+
  xlab('TG_GATEs_Repeated_Dose')+ylab('Pathway Activation Scores')+ theme_minimal()+
  theme(legend.text = element_text(size=10),
        # strip.text = element_blank(),
        legend.position = 'bottom',
        axis.text.x = element_text(size=16,angle = 45,hjust = 1),
        axis.title = element_text(size = 11),
        axis.text.y = element_text(size=11),panel.spacing = unit(2,'lines')) +
  facet_wrap(. ~ type)

save(data,file='outputData/plot_data/time_bio_marker_all_increase_1_2_5_scenarios.RData')
#mean_decs<-as.data.frame(lapply(decs, function(x)apply(x[,sel],2,mean)))
#mean_edcs$type<-'decoy'
#boxplot(mean_edcs)
#plot(apply(mean_edcs, 2, mean))
#lines(apply(mean_decs, 2, mean))
#matplot(t(mean_decs))

