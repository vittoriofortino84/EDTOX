library(STRINGdb)
library(igraph)
library(org.Hs.eg.db)

combscore<-0.85           #based on pareto
string_db<-STRINGdb$new(version="10",species=9606)
entz_id<-data.frame(entz_id=keys(org.Hs.egGENENAME),stringsAsFactors = F)
string_ids<-string_db$map(my_data_frame=entz_id,my_data_frame_id_col_names='entz_id',removeUnmappedRows=T)
string_inter<-string_db$get_interactions(string_ids$STRING_id)
string_inter<-string_inter[,-grep('coexpression',names(string_inter))]
sel_inter<-1-(string_inter[,3:(ncol(string_inter)-1)]/1000)
string_inter$combined_score<-1-Reduce('*',sel_inter)
edge_list<-string_inter[string_inter$combined_score>=combscore,1:2]
edge_list<-merge(edge_list,string_ids,by.x='from',by.y='STRING_id')
edge_list<-merge(edge_list,string_ids,by.x='to',by.y='STRING_id')[,3:4]
colnames(edge_list)<-c('Node.1','Node.2')
gr<-graph_from_edgelist(as.matrix(edge_list),directed=F)
save(gr,file='outputData/network/ppi_network.RData')