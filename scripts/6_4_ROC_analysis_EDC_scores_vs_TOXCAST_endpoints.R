
# 1. harmonic -------------------------------------------------------------

# ROC curve analysis and permutation test between nuclear receptor endpoints in toxcast and harmonic_VAM score
library(PRROC)
source('functions/general_functions.R')     # ROC permutation function
load('outputData/VAM_scores_hs_most_informative_layers_moa.RData')
load('outputData/hitc_chemid_NR_coreg.RData') # hitc binary for nuclear receptor toxcast endpoints
load("outputData/most_informative_layers_final_moa.RData")
assay<-chemid_hic_NR[,most_reliable_assay_endpoint]
ac<-assay[intersect(rownames(assay),names(all_VAM_score)),]
cp<-all_VAM_score[intersect(rownames(assay),names(all_VAM_score))]
rocmat_length<-rocmat<-perm_res<-matrix(NA, nrow = 1, ncol = ncol(assay))
colnames(rocmat)<-colnames(rocmat_length)<-colnames(assay)
for (j in 1:ncol(assay)){
  m<-cbind(cp,as.numeric(ac[,j]))
  m<-na.omit(m)
  rocmat_length[j]<-nrow(m)
  ratio<-sum(m[,2])/nrow(m)
  roc<-roc.curve(scores.class0 = m[,1], weights.class0 = m[,2],curve = T)
  rocmat[j]<-roc$auc
  perm_res[j]<-roc_permutation(m[,1],m[,2],npermutation = 10000)
}#enf for j
ROC_res<-as.data.frame(t(rbind(rocmat,perm_res)))
names(ROC_res)<-c('AUC_ROC','P_value')
save(ROC_res,file = 'outputData/permutation_ROC_hs_VAM_score_moa.RData')
write.csv(ROC_res,file = 'outputData/excel_files/permutation_ROC_hs_VAM_score_moa.csv')



# 2. average ---------------------------------------------------------------

# ROC curve analysis and permutation test between nuclear receptor endpoints in toxcast and average_VAM score
library(PRROC)
source('functions/general_functions.R')     # ROC permutation function
load('outputData/VAM_scores_average_most_informative_layers_moa.RData')
load('outputData/hitc_chemid_NR_coreg.RData') # hitc binary for nuclear receptor toxcast endpoints
#load("outputData/most_informative_layers_final.RData")
load("outputData/most_informative_layers_final_moa.RData")
assay<-chemid_hic_NR[,most_reliable_assay_endpoint]
ac<-assay[intersect(rownames(assay),names(all_VAM_score)),]
cp<-all_VAM_score[intersect(rownames(assay),names(all_VAM_score))]
rocmat_length<-rocmat<-perm_res<-matrix(NA, nrow = 1, ncol = ncol(assay))
colnames(rocmat)<-colnames(rocmat_length)<-colnames(assay)
for (j in 1:ncol(assay)){
  m<-cbind(cp,as.numeric(ac[,j]))
  m<-na.omit(m)
  rocmat_length[j]<-nrow(m)
  ratio<-sum(m[,2])/nrow(m)
  roc<-roc.curve(scores.class0 = m[,1], weights.class0 = m[,2],curve = T)
  rocmat[j]<-roc$auc
  perm_res[j]<-roc_permutation(m[,1],m[,2],npermutation = 10000)
}#enf for j
ROC_res<-as.data.frame(t(rbind(rocmat,perm_res)))
names(ROC_res)<-c('AUC_ROC','P_value')
# save
save(ROC_res,file = 'outputData/permutation_ROC_average_VAM_score_moa.RData')
write.csv(ROC_res,file = 'outputData/excel_files/permutation_ROC_average_VAM_score_moa.csv')


# 3. Heatmap of the class probabilites for different layers VS the most significant assay endpoints  --------------------------------------

Roc_harmonic_sum<-read.csv('outputData/excel_files/permutation_ROC_hs_VAM_score_moa.csv') # the result of permutation for harmonic score 
Roc_average<-read.csv('outputData/excel_files/permutation_ROC_average_VAM_score_moa.csv') # the result of permutation for average score
significant_endpoints<-intersect(Roc_average$X[Roc_average$P_value==0],Roc_harmonic_sum$X[Roc_harmonic_sum$P_value==0])
load('outputData/most_informative_layers_final_moa.RData')


#data<-t(rocmat[,significant_endpoints]) # endpoints with p-permutations=0
data<-t(rocmat);data<-data[-which(rowSums(data)==0),]
data<-data[which(apply(data, 1, mean) >= .6),]


colnames(data)<-gsub(x=colnames(data),pattern = 'HEPG2_HEPG2', replacement = 'HEPG2')
colnames(data)<-gsub(x=colnames(data),pattern = 'HEPG2', replacement = 'HEPG2_1_day')
colnames(data)<-gsub(x=colnames(data),pattern = '__TG_GATEs', replacement = '_TG_GATEs')
roc_aucs<-data
save(roc_aucs,file = 'outputData/plots_tables_functions/ROC_class_prob_Toxcast_hitc/ROC_aucs_class_prob_hitc.RData')

## using the library gplots
library(gplots)
par(mar=c(7,4,4,2)+2)
ht<-heatmap.2(roc_aucs,trace = 'none', col = colorRampPalette(c('white', 'red'))(12),margins = c(20,20),
          density.info = 'none',srtCol = 45,keysize = 1,key.xlab = 'AUC-ROC',key.title = NA,cexRow = 1,cexCol = 1.2)
print(ht)
## using complex heatmap
library(ComplexHeatmap)
library(circlize)
mat=data
col_fun<-colorRamp2(c(min(data),max(data)),c('white','red'))
cn<-colnames(mat)
colnames(mat)<-c()
column_ha<-HeatmapAnnotation(roc=anno_boxplot(mat,height = unit(4,'cm')),show_annotation_name =  F)
botom_ha<-HeatmapAnnotation(text=anno_text(cn,rot=45),show_annotation_name = F)
row_ha<-rowAnnotation(roc=anno_boxplot(mat,width  = unit(6,'cm'),axis_param = list(side='bottom',labels_rot=45)),show_annotation_name =  F)
ht<-Heatmap(mat,name='ROC-AUC',col = col_fun,
            #top_annotation = column_ha,
            #right_annotation = row_ha,
            bottom_annotation = botom_ha,
            height = unit(15,'cm'),
            width = unit(12,'cm'))
draw(ht,padding=unit(c(5,5,5,5),'cm'),heatmap_legend_side='left')

# saving the plot 
tiff('plots/figure3.tiff',units = 'in',width = 15,height = 12,res = 300)
draw(ht,padding=unit(c(10,4,10,10),'cm'),heatmap_legend_side='left')
dev.off()
save(significant_endpoints,file='outputData/TOXcast_endpoints_with_significant_permutation_ROC_moa.RData') # endpoints with p-permutations=0



# 4. Integration of harmonic and average  of ToxCast endpoints vs VAM scores ROC analysis ------------------------------------------
load('outputData/permutation_ROC_hs_VAM_score_moa.RData')
colnames(ROC_res)<-c('harmonic_edc_score_ROC','harmonic_score_permutation');harmonic_roc<-ROC_res;harmonic_roc$endpoints<-rownames(harmonic_roc)
load('outputData/permutation_ROC_average_VAM_score_moa.RData')
colnames(ROC_res)<-c('average_edc_score_ROC','average_score_permutation');average_roc<-ROC_res;average_roc$endpoints<-rownames(average_roc)
p_test<-merge(average_roc,harmonic_roc)
load('outputData/most_informative_layers_final_moa.RData')
roc_result<-merge(p_test,data)
write.csv(roc_result,'outputData/excel_files/ROC_final_results_ToxCast_moa.csv')





