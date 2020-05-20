# 16 april 2020 we searched the term endocrine disruption pubmed
library(pubmed.mineR)
abstracts <- readabs('inputData/text_mining/pubmed_result_endocrine_disruption.txt')
pmids<-abstracts@PMID
res<-list()
for (i in 1:length(pmids)){res[[i]]<-pubtator_function(pmids[i]);print(i)}
res[which(sapply(res, length)==1)]<-c()
res_genes<-lapply(res, function(x)x$Genes) # all genes reported in the papers
res_genes<-unlist(res_genes)
library(magrittr)
entrez<-res_genes %>% strsplit('>') %>% lapply('[[', 2) %>% unlist() %>% strsplit(';') %>% unlist()
gene_freq<-as.data.frame(table(entrez))
# saving
save(gene_freq,res,file = 'inputData/text_mining/endocrine_disruption_result.RData')
