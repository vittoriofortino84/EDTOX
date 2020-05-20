
# 4.1   EDC and decoy selection -----------------------------------------------
load("outputData/toxcast_direct_based_NR_coreg_from_TOXCAST_binary.RData")   #list of toxcast compounds for NR receptors
rm(hitc_mat_MIE);rm(toxcast_3_1)
edcs<-read.csv('inputData/DEDuCT_ChemicalBasicInformation.csv',header = TRUE) #list of edcs from DEDuCT
edcs<-unique(edcs$CAS.Number)
hitc_mat_NR_Coreg<-as.matrix(hitc_mat_NR_Coreg)
edc_ind<-which(rownames(hitc_mat_NR_Coreg) %in% edcs)
edc_hitc_NR<-hitc_mat_NR_Coreg[edc_ind,]
selectedendponints<-colnames(edc_hitc_NR)
save(selectedendponints,file='outputData/265_endpoints_significant_LR.RData') # The list of 265 NR receptors 

# to get the cutoff for the optimum number of EDCs  tested by TOXCAST assay endpoints with no NA values
number_compounds<-number_assays<-rep(NA,nrow(edc_hitc_NR))
for (i in 1:length(number_compounds)){
  selected_NR_endpoints<-which(colSums(!is.na(edc_hitc_NR))>i) # select endpoints with more than i compounds tested
  edc_hitc_NR2<-edc_hitc_NR[,selected_NR_endpoints]            # subseting these endpoints
  edc_hitc_NR2<-na.omit(edc_hitc_NR2)                          # removing compounds with NA values
  number_assays[i]<-ncol(edc_hitc_NR2)                         # Number of assay endpoints
  number_compounds[i]<-nrow(edc_hitc_NR2)                      # number of compounds after removing NAs
}
rm(edc_hitc_NR2)
datam<-as.data.frame(list(number_assays,number_compounds))
colnames(datam)<-c('number_assays','number_compounds')
datam$label<-paste(datam$number_compounds,'Compounds',datam$number_assays,'assays',sep = '_')
datam<-unique(datam)
datam$label[!datam$number_compounds==304]<-''
library(ggplot2)
library(ggrepel)
ggplot(datam,aes(x=number_assays,y=number_compounds))+
  geom_point(size=ifelse(datam$number_compounds==304,4,2),color=ifelse(datam$number_compounds==304,'black','red'))+ geom_line()+ 
  ggtitle('Subsets with none NA values')+ xlab('Number of assay endpoints')+ylab('Number of compounds')+
  scale_y_continuous(breaks = seq(0, nrow(edc_hitc_NR), 20)) +
  scale_x_continuous(breaks = seq(0, ncol(edc_hitc_NR), 20))+
  theme_minimal() +
  geom_text_repel(label=datam$label,size=4,segment.size = 0.8, segment.color = "grey") +
  theme(text = element_text(size=20))

# data_pareto<-as.data.frame(list(number_assays,number_compounds))
# colnames(data_pareto)<-c('assay','edcs')
# ps<-rPref::psel(data_pareto,high(edcs)*high(assay))

# After getting the cutoff we select 116 to cover 304 compounds
selected_NR_endpoints<-which(colSums(!is.na(edc_hitc_NR))>300) # we select 116 NRs with none NA values for at least 300 EDCs (between Deduct and ToxCast)
edc_hitc_NR<-edc_hitc_NR[,selected_NR_endpoints]
edc_hitc_NR<-na.omit(edc_hitc_NR)
less_than_5_percent<-which(colSums(edc_hitc_NR)<round((dim(edc_hitc_NR)[1]*5)/100))                    
edc_hitc_NR<-edc_hitc_NR[,-less_than_5_percent]  # if an endpoint is active for less then 5% of edcs is elimiated

library(magrittr)
selected_NR_endpoints<-colnames(edc_hitc_NR)
types<-selected_NR_endpoints %>% strsplit(split='_',fixed = T) %>% lapply("[", 2) %>% unique() %>% unlist() # Different types of endpoints
# performing fisher exact test to see if the proportion of edcs are different for endpoints directions
z_tst<-function(x,var){
  var<-paste('_',var,'_',sep = '')
  sel<-grep(var,colnames(x))
  if (var =='_CAR_'){sel<-sel[-c(1,2)]}
  if (var =='_ERR_'){sel<-sel[-c(3,4)]}
  if (length(sel)==2) { # the test is done only if two endpoints are present with different directions
    x<-x[,sel]
    m<-colSums(x)
    m<-which(m==min(m))
    m<-names(m)
    model.1 <- prop.test(x=c(colSums(x)[1],colSums(x)[2]),c(nrow(x),nrow(x)),alternative = 'two.sided')
    q<-model.1$p.value
    print(paste(colnames(x)[1],'(',colSums(x)[1],') ',colnames(x)[2],'(',colSums(x)[2],')','p_value=',q))
    names(q)<-m
  }else{
    q<-'not'
    names(q)<-var
  }
  return(q)
} # function for z fisher test
ad_l<-c('AR','FXR','PPARd','VDR','ERa','GR','PPARg')
types<-types[-which(types %in% ad_l)]
ad_l<-data.frame(expand.grid(c('AR','FXR','PPARd','VDR','ERa','GR','PPARg'),c('_TRANS','_BLA')))
colnames(ad_l)<-c('v1','v2')
ad_l<-paste(ad_l$v1,ad_l$v2,sep = '')
types<-c(ad_l,types)
less_significant_mies<-sapply(types, function(x)z_tst(edc_hitc_NR,x))  #function z_tst call
less_significant_mies<-less_significant_mies[-which(less_significant_mies=='not')]
less_significant_mies<-less_significant_mies[which(as.numeric(less_significant_mies)<0.05)] # Those with p.values >0.05 in Fisher test
to_be_removed<-unlist(sapply(colnames(edc_hitc_NR), function(x)grep(x,names(less_significant_mies))))
p_to_be_removed<-unlist(sapply(colnames(edc_hitc_NR), function(x)less_significant_mies[grep(x,names(less_significant_mies))]))

ind<-which(colnames(edc_hitc_NR) %in% names(to_be_removed))
edc_hitc_NR<-edc_hitc_NR[,-ind]                                    #Removing un_significant endpoints based on Logistic regression 
selected_NR_endpoints<-colnames(edc_hitc_NR)                     
print(paste(length(which(rowSums(edc_hitc_NR) %in% 0)),'will be removed from EDCs list'))
edc_hitc_NR<-edc_hitc_NR[-which(rowSums(edc_hitc_NR) %in% 0),]     #EDCs which are inactive for all significant endpoints are eliminated
data<- as.data.frame(colSums(edc_hitc_NR))
colnames(data)<-'Active_Compounds'
data$Assay_Endpoints<-rownames(data)

library(ggplot2)

#data$Active_Compounds<-factor(data$Active_Compounds,levels = sort(unique(data$Active_Compounds)))
data$Active_Compounds<-as.numeric(as.character(data$Active_Compounds))
ggplot(data,aes(x=reorder(data$Assay_Endpoints,-data$Active_Compounds),y=Active_Compounds,fill=Active_Compounds))+geom_bar(stat = 'identity')+
  theme_classic()+
  scale_y_continuous()+
  ggtitle('Significant TOXCAST endpoints based on Fisher test')+ylab('Number of Active EDCs')+
  #labs(x='',fill=NULL)+
  theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1,size = 8),    #labels in x axis
        plot.margin = unit(c(2,2,2,2),'cm'))


# integrating with the compounds in CTD with calculated MIEs
ixns<-read.csv('inputData/CTD_chem_gene_ixns.csv.gz',comment.char = c("#"),stringsAsFactors = F) #all compounds from CTD
library(data.table)
ixns<-as.data.table(ixns)
library(dplyr)
ixns<-ixns %>% distinct(ixns$CasRN,ixns$ChemicalID)
chemnames<-rownames(edc_hitc_NR)
edc_hitc_NR<-as.data.table(edc_hitc_NR)
edc_hitc_NR$cas<-chemnames
colnames(ixns)<-c('cas','chem')
edc_in_ctd<-merge(edc_hitc_NR,ixns)               #228 edcs
load('outputData/chem2gene_no_out.RData')
ed_ind<-which(rownames(chem2gene) %in% edc_in_ctd$chem)

edc199<-chem2gene[ed_ind]


#### NEGATIVE CONTROL SELECTION
# selection based on the most disimilar compound compared to edcs in terms of MIEs
# calculation of jaccard dissimilarity between edcs and OTHER COMPOUNDS IN CTD
jan<-matrix(NA, nrow = length(edc199), ncol = length(chem2gene))
for (i in 1:dim(jan)[1]){
  for (j in 1:dim(jan)[2]){
    l1<-length(intersect(edc199[[i]],chem2gene[[j]]))
    l2<-length(union(edc199[[i]],chem2gene[[j]]))
    jan[i,j]<-(1-(l1/l2))
  }
}

ind_disimilar<-which(apply(jan, 2,mean)==1)  #selecting those compounds as decoys with average jaccard dissimilarity equal to 1 
decoys<-chem2gene[ind_disimilar]
new_edc199_decoy_1462<-c(edc199,decoys)
save(new_edc199_decoy_1462,file='outputData/new199edc_1462dec.RData')

# 4.2   Visualization of EDCs and decoys dissimililarities as heatmap PLOT---------

load('outputData/new199edc_1462dec.RData')
input<-new_edc199_decoy_1462              #decoys based on CTD
nedc<-199
ndecoy<-1462
# 
jan<-matrix(NA, nrow = length(input), ncol = length(input))
for (i in 1:(length(input)-1)){
  er<-i+1
  for (j in er:length(input)){
    
    l1<-length(intersect(input[[i]],input[[j]]))
    l2<-length(union(input[[i]],input[[j]]))
    
    jan[i,j]<-jan[j,i]<-(1-(l1/l2))
  }
}
diag(jan)<-0
#Heatmap Plot
melted_cormat <- reshape::melt(jan)
source("functions/plot_functions.R") #heatmap
plot(heat_map(melted_cormat,melted_cormat$value,melted_cormat$X1,melted_cormat$X2))


# 4.3   Calculation of silhouette score between edcs based on jaccard distance TABLE S3:supplementary--------

source('functions/general_functions.R')                # jaccrd_dissimilarity function for binary matrices
load('outputData/chem2gene_no_out.RData')              # All mies in CTD
load('outputData/new199edc_1462dec.RData')             # MIEs for edcs and decoys
all_genenes<-unique(unlist(sapply(chem2gene,function(x)x)))
binary_MIE<-matrix(0, nrow = length(chem2gene), ncol = length(all_genenes))
for  (i in 1:length(chem2gene)){
  binary_MIE[i,which(all_genenes %in% chem2gene[[i]])]<-1
}
colnames(binary_MIE)<-all_genenes
rownames(binary_MIE)<-names(chem2gene)
training_set<-binary_MIE[names(new_edc199_decoy_1462),]
jac<-jaccrd_dissimilarity(training_set)
km <- kmeans(jac, centers = 2, nstart=25)
km$cluster<-c(rep(1,199),rep(2,1462))
library(cluster)
ss<-silhouette(km$cluster, jac)
mean(ss[1:199,3 ]) 
sil_MIE_edcs<-mean(ss[1:199,3 ])
save(sil_MIE_edcs,file='outputData/silhouettescores_MIE_edcs.RData')



#OLD functions 
# LR_tst<-function(x,var){ #logistic regression
#require(car)
#   var<-paste('_',var,'_',sep = '')
#   sel<-grep(var,colnames(x))
#   if (var =='_CAR_'){sel<-sel[-c(1,2)]}
#   if (var =='_ERR_'){sel<-sel[-c(3,4)]}
#   
#   if (length(sel)>1) {
#     x<-x[,sel]
#     print(colnames(x))
#     m<-apply(x, 2, mean)
#     m<-which(m==min(m))
#     m<-names(m)
#     x<-stack(as.data.frame(x))
#     model.1 <- glm(ind~values, data=x, family=binomial(link="logit"))
#     q<-car::Anova(model.1, test="LR", type="III")
#     q<-q$`Pr(>Chisq)`
#     print(paste(length(sel),'   ',q))
#     names(q)<-m
#   }else{
#     
#     q<-'not'
#     names(q)<-var
#   }
#   return(q)
# } # function for logisitc regression
# t_tst<-function(x,var){ #t-test
#   var<-paste('_',var,'_',sep = '')
#   sel<-grep(var,colnames(x))
#   if (var =='_CAR_'){sel<-sel[-c(1,2)]}
#   if (var =='_ERR_'){sel<-sel[-c(3,4)]}
#   if (length(sel)==2) {
#     x<-x[,sel]
#     print(colnames(x))
#     m<-apply(x, 2, mean)
#     m<-which(m==min(m))
#     m<-names(m)
#     
#     model.1 <- t.test(x[,1],x[,2])
#     
#     q<-model.1$p.value
#     #print(paste(length(sel),'   ',q))
#     names(q)<-m
#   }else{
#     
#     q<-'not'
#     names(q)<-var
#   }
#   return(q)
# } # function for t-test
# 

