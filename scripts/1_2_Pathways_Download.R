# 2.   Retrieving pathways from KEGG,REACTOME,WIKI,msigdb and GO for FGSEA ----------------------------------------------------------------------
# Pathways related to wikiaop (GO and GO with offsprings), KEGG, REACTOME,msigdb and wiki pathways are retrieved
# output :outputData/toxdb2gene_final.RData

library(XML)
#####USING aopwiki
AOPwiki<-xmlParse('inputData/aop-wiki-xml-2019-01-01')
AOPwiki<-xmlToList(AOPwiki)                                                  # Wiki adverse outcome pathways
bp<-AOPwiki[names(AOPwiki)=="biological-process"]                            # Seperate  "biological process" list from other processes like stressor,chemical so on  
bp<-as.data.frame(t(sapply(bp,unlist)),stringsAsFactors=F)                   # Change it into transposed data Frame
goid<-bp$`source-id`[bp$source=='GO']                                        # Seperate those biological processes with GO term from others like MESH 
library(GO.db)
library(org.Hs.eg.db)
ks<-keys(GO.db)                                                             # keys : all GO pathways in GO database
GO_ann<-select(GO.db,keys=ks,columns=c('TERM',"ONTOLOGY"),keytype="GOID")   # Mapping terms and ontology from all keys in a dataframe (GOID,Pathway_TERM,ONTOLOGY)
GO_sel<-GO_ann[GO_ann$GOID%in%goid,]                                        # selecting only those GO terms among the list of all GOs that we have GOID for them in wiki_AOP
orgdb_goid<-unique(as.data.frame(org.Hs.egGO)$go_id)                        # GO terms in Genome wide annotation for Human (org.Hs.egGO)
go2gene_noDesc<-sapply(intersect(GO_sel$GOID,orgdb_goid),function(x)        # Search intersection of selected GO terms and those in org.Hs.egGO and reurn>>
  unique(unlist(mget(x,revmap(org.Hs.egGO)))))                              # GO terms with their related genes in "go2gene_noDesc"

##### Finding the descendants of the pathways
get_go_descendants<-function(x,onto,goid){                                   # Adds offsppring for the specific ontology
  if(onto=='BP')off<-GOBPOFFSPRING[[x]]                                      # of the given GO term and 
  else if (onto=='MF')off<-GOMFOFFSPRING[[x]]
  else off<-GOCCOFFSPRING[[x]]
  desc<-intersect(c(x,off),goid)                                             # returns the intersection of GO+offsprings with annotaions of org.Hs.egGO
  names(desc)<-NULL
  return(desc)
}
go_desc<-apply(GO_sel,1,function(x)get_go_descendants(x[1],x[3],orgdb_goid)) #For each selected GO term returns a list of its offspring GO terms
names(go_desc)<-GO_sel$GOID
go_desc<-go_desc[sapply(go_desc,length)>0]
go2gene<-lapply(go_desc,function(x)unique(unlist(mget(x,revmap(org.Hs.egGO))))) # GO terms with their corresponding genes in "go2gene"

#####USING KEGG  pathways from CTD genes pathways
ptw<-read.csv('inputData/CTD_genes_pathways.csv.gz',stringsAsFactors = F,comment.char ="#",header = F)  
kegg<-ptw[grepl('KEGG',ptw$V4),]                                                    #Getting Column V4 with KEGG terms as (pathways)
kegg2gene<-sapply(unique(kegg$V4),function(x)as.character(kegg$V2[kegg$V4==x]))     #For each KEGG term the related genes are assigned and saved in "kegg2gene"
#####REACTOME
#https://reactome.org/download-data: 
library(GSA)
human_rec_names<-read.delim('inputData/NCBI2Reactome_PE_All_Levels.txt',header = FALSE)
human_rec_names<-human_rec_names[human_rec_names$V8=='Homo sapiens',]
reactom<-GSA.read.gmt('inputData/ReactomePathways.gmt')
react_anntoations<-cbind(reactom$geneset.names,reactom$geneset.descriptions)
write.csv(react_anntoations,'outputData/reactome.csv')
human_ind<-which(as.character(reactom$geneset.descriptions) %in% as.character(human_rec_names$V4))
react<-reactom$genesets[human_ind]
names(react)<-reactom$geneset.descriptions[human_ind]
library(org.Hs.eg.db)
all_genes <- org.Hs.egSYMBOL # Get the gene symbol that are mapped to an entrez gene identifiers 
mapped_genes <- mappedkeys(all_genes) # Convert to a list 
all_genes <- as.list(all_genes[mapped_genes])
reactome2gene<-lapply(react,function(x)names(all_genes[all_genes %in% as.character(x)]))
to_remove<-which((unlist(lapply(react,function(x)length(x)))-unlist(lapply(reactome2gene,function(x)length(x))))>1) #Omit if there was even a single gene with no ENTREZ for the symbol
reactome2gene<-reactome2gene[-to_remove]
rm(list=c('all_genes','human_rec_names','react','reactom','to_remove','mapped_genes'))  
# msigdb pathways
library(msigdbr)
m_df = msigdbr(species = "Homo sapiens")
m_df_ids = m_df$gs_subcat=='BP'| m_df$gs_subcat=='MF'
m_df<-m_df[m_df_ids,]                         
pat_nam=unique(m_df$gs_name)
msi<-lapply(pat_nam,function(x)m_df$entrez_gene[m_df$gs_name==as.character(x)])
names(msi)=pat_nam
msigdb2gene<-lapply(msi, function(x)as.character(x))
load("inputData/c5.go.mapping.rda")
q<-unlist(lapply(names(msigdb2gene),function(x)gsub(pattern = 'GO_',replacement = '',x)))
goids<-c5.go.mapping$goid[which(c5.go.mapping$description %in% q)]
go2gene<-go2gene[-which(names(go2gene)%in% goids)]
go2gene_noDesc<-go2gene_noDesc[-which(names(go2gene_noDesc)%in% goids)]
# wiki pathways
library(rWikiPathways)
library(GSA)
downloadPathwayArchive(destpath ='inputData/', date=20190510,organism="Homo sapiens", format='gmt')
wikiptw<-GSA.read.gmt('inputData/wikipathways-20190510-gmt-Homo_sapiens.gmt')
wikiptw2gene<-wikiptw$genesets                                                      #saving related genes to wiki pathways in wikiptw2gene
names(wikiptw2gene)<-wikiptw$geneset.names
#finalization
g_max_size<-200
go2gene<-go2gene[-which(sapply(go2gene, function(x)length(x)>g_max_size))]
go2gene_noDesc<-go2gene_noDesc[-which(sapply(go2gene_noDesc, function(x)length(x)>g_max_size))]
kegg2gene<-kegg2gene[-which(sapply(kegg2gene, function(x)length(x)>g_max_size))]
msigdb2gene<-msigdb2gene[-which(sapply(msigdb2gene, function(x)length(x)>g_max_size))]
reactome2gene<-reactome2gene[-which(sapply(reactome2gene, function(x)length(x)>g_max_size))]
wikiptw2gene<-wikiptw2gene[-which(sapply(wikiptw2gene, function(x)length(x)>g_max_size))]
toxdb_names<-names(c(go2gene,wikiptw2gene,msigdb2gene,reactome2gene,kegg2gene))
save(toxdb_names,msigdb2gene,go2gene,go2gene_noDesc,kegg2gene,reactome2gene,wikiptw2gene,file='outputData/toxdb2gene_final.RData')


# 2.1. Go to AOP annotations  ------------------------------------------------
# converting aop_ke_ec.tsv from wiki_aop and their  into ta binary matrix each AO-pathway and GO terms are correlated
# the output is aop_mat.csv and is needed for the function pathway_annotate
#output:outputData/aop_mat.csv
#  https://aopwiki.org/downloads/aop_ke_ec.tsv

aop_data<-read.delim(file='inputData/aop_ke_ec.tsv',sep='\t',header = F)
source_go<-aop_data$V8[which(aop_data$V7=='GO')]
object_go<-aop_data$V5[which(aop_data$V4=='GO')]
all_go<-unique(c(as.character(object_go),as.character(source_go)))
aops<-unique(aop_data$V1)
fm<-matrix(0, nrow = length(aops), ncol = length(all_go))
colnames(fm)<-all_go
rownames(fm)<-aops
for (i in 1:nrow(aop_data)){
  if (any(aop_data[i,5] %in% all_go)==T){
    ind<-which(all_go %in% aop_data[i,5])
    fm[which(rownames(fm) %in% aop_data[i,1]),ind]=1
  }
  if (any(aop_data[i,8] %in% all_go)==T){
    ind<-which(all_go %in% aop_data[i,8])
    fm[which(rownames(fm) %in% aop_data[i,1]),ind]=1
  }
}

annot_correction<-function(inp){
  library(GO.db)
  Go_annotations<-sapply(GOTERM,function(x)x@Term)
  nn<-inp
  for (i in 1:length(inp)){
    # GO annotaions
    if (any(grep(pattern = 'GO:',x = inp[i] ))==TRUE){
      nam<-Go_annotations[inp[i]]
      nam<-gsub(' ','_',nam)
      nn[i]<-paste('GO_',nam,sep = '')
    } #end if
    nn[i]<-toupper(nn[i])
  }#end for
  return(nn)
} #end function

colnames(fm)<-annot_correction(colnames(fm))
write.csv(fm,file='outputData/aop_mat.csv')







# 2.2  kegg_categorization --------------------------------------------------
system('cp ./scripts/bash/keg2tag.sh   ./inputData/annotaion/kegg')
setwd('inputData/annotaion/kegg')
system('./keg2tag.sh')
main_cat_files<-list.files('.',pattern = "\\.tag$")
main_tags<-lapply(main_cat_files, function(x)unlist(read.table(x,header=F,stringsAsFactors = F,
                                                               colClasses = 'character')))
library(magrittr)
names(main_tags)<-main_cat_files %>% strsplit('.tag') %>% lapply('[[', 1) %>% unlist()
setwd('Metabolism/')
met_cat_files<-list.files('.',pattern = "\\.tag$")
met_tags<-lapply(met_cat_files, function(x)unlist(read.table(x,header=F,stringsAsFactors = F,
                                                             colClasses = 'character')))
names(met_tags)<-met_cat_files %>% strsplit('.tag') %>% lapply('[[', 1) %>% unlist()
setwd('../../../..')
save(main_tags,met_tags,file = 'inputData/annotaion/keg_from_web.RData')
####
# ## DICTIONARY FOR keg2category function
rm(list=ls())
# library(keggorthology) ; library(igraph)
# data(KOgraph) # because the object Kograph was not updated we decided to retrieve the tags from internet as aobve
# g<-igraph::igraph.from.graphNEL(KOgraph)
# df<-as_data_frame(g,'both')
# v<-df$vertices
# edg<-df$edges
# colnames(edg)<-c('class','name')
# merged<-merge(edg,v,all=F)
# # extracting pathway  tags from the mappings
# main_class<-merged$name[merged$depth==1]
# kegg_category_lib<-lapply(main_class, function(x)
#   c(merged$tag[merged$class %in% merged$name[merged$class==x]],
#     merged[merged$class==x,'tag'],
#     merged[merged$name==x,'tag']))
# kegg_category_lib$Metabolism<-c()
load('inputData/annotaion/keg_from_web.RData')
kegg_category_lib<-main_tags
names(kegg_category_lib)<-gsub(x=names(kegg_category_lib),pattern="Cellular_Processes",
                       replacement = 'Cellular Processes')
names(kegg_category_lib)<-gsub(x=names(kegg_category_lib),pattern="Environmental_Information_Processing",
                       replacement = 'Environmental Information Processing')
names(kegg_category_lib)<-gsub(x=names(kegg_category_lib),pattern="Human_Diseases",
                       replacement = 'Human Diseases')
names(kegg_category_lib)<-gsub(x=names(kegg_category_lib),pattern="Genetic_Information_Processing",
                       replacement = 'Genetic Information Processing')
names(kegg_category_lib)<-gsub(x=names(kegg_category_lib),pattern="Organismal_Systems",
                       replacement = 'Organismal Systems')
names(kegg_category_lib)<-gsub(x=names(kegg_category_lib),pattern="Drug_Development",
                       replacement = 'Drug Development')
kegg_category_lib$Metabolism<-c()
# 
# 
# 
# metabolism_categories<-merged$name[merged$class=='Metabolism']  
# metabolism_sub_lib<- lapply(metabolism_categories, function(x) 
#   c(merged[merged$class==x,'tag'],
#     merged[merged$name==x,'tag']))             
# names(metabolism_sub_lib)<-tolower(metabolism_categories)

metabolism_sub_lib<-met_tags
names(metabolism_sub_lib)<-gsub(x=names(metabolism_sub_lib),pattern="Aminoacidmetabolism",
                       replacement = 'amino acid metabolism')
names(metabolism_sub_lib)<-gsub(x=names(metabolism_sub_lib),pattern="Carbohydratemetabolism",
                                replacement = 'carbohydrate metabolism')
names(metabolism_sub_lib)<-gsub(x=names(metabolism_sub_lib),pattern="Energymetabolism",
                                replacement = 'energy metabolism')
names(metabolism_sub_lib)<-gsub(x=names(metabolism_sub_lib),pattern="Glycanbiosynthesisandmetabolism",
                                replacement = 'glycan biosynthesis and metabolism')
names(metabolism_sub_lib)<-gsub(x=names(metabolism_sub_lib),pattern="Metabolismofcofactorsand",
                                replacement = 'metabolism of cofactors and vitamins')
names(metabolism_sub_lib)<-gsub(x=names(metabolism_sub_lib),pattern="Metabolismofterpenoidsand",
                                replacement = 'biosynthesis of polyketides and nonribosomal peptides')
names(metabolism_sub_lib)<-gsub(x=names(metabolism_sub_lib),pattern="Xenobioticsbiodegradationandmetabolism",
                                replacement = 'xenobiotics biodegradation and metabolism')
names(metabolism_sub_lib)<-gsub(x=names(metabolism_sub_lib),pattern="Biosynthesisofothersecondary",
                                replacement = 'biosynthesis of secondary metabolites')
names(metabolism_sub_lib)<-gsub(x=names(metabolism_sub_lib),pattern="Lipidmetabolism",
                                replacement = 'lipid metabolism')
names(metabolism_sub_lib)<-gsub(x=names(metabolism_sub_lib),pattern="Metabolismofotheramino",
                                replacement = 'metabolism of other amino acids')
names(metabolism_sub_lib)<-gsub(x=names(metabolism_sub_lib),pattern="Nucleotidemetabolism",
                                replacement = 'nucleotide metabolism')
metabolism_sub_lib$Chemicalstructuretransformationmaps<-c()
metabolism_sub_lib$Globalandoverviewmaps<-c()

# Combining with kegg modules
library(rjson)
result <- as.matrix(as.data.frame(fromJSON(file = "inputData/annotaion/kegg.json")))

class_ind<-grep(x=colnames(result),pattern = '^children.children.name') # indexes for the classes
#grep(x=colnames(result),pattern = '^children.children.children.name') # indexes for the subclasses
ind_mat<-matrix(NA, nrow = 2, ncol = length(class_ind))
ind_mat[1,]<-class_ind+2
ind_mat[2,]<-c(class_ind[2:length(class_ind)]-1,length(result))
colnames(ind_mat)<-result[1,class_ind]

modules_tags<-list()
for (i in 1:ncol(ind_mat)){
  modules_tags[[i]]<-sapply(sapply(result[1,ind_mat[1,i]:ind_mat[2,i]], 
                                   function(x)strsplit(x,' ',fixed=T)),
                            function(y)y[[1]])
  
}
names(modules_tags)<-tolower(colnames(ind_mat))
modules_tags$`module set`<-c()
modules_tags$`gene set`<-c()
save(modules_tags,file='inputData/annotaion/kegg_modules.RData')
inters<-intersect(names(modules_tags),names(metabolism_sub_lib))
intersected_met<-mapply(function(x,y)c(x,y),modules_tags[inters],metabolism_sub_lib[inters])

all_met<-c(intersected_met,modules_tags[setdiff(names(modules_tags),inters)],
           metabolism_sub_lib[setdiff(names(metabolism_sub_lib),inters)])

all_met$`amino acid metabolism`<-c(all_met$`amino acid metabolism`,all_met$`metabolism of other amino acids`)
all_met$`metabolism of other amino acids`<-c()
all_met$`xenobiotics biodegradation`<-c(all_met$`xenobiotics biodegradation`,all_met$`xenobiotics biodegradation and metabolism`)
all_met$`xenobiotics biodegradation and metabolism`<-c()
all_met$`biosynthesis polyketides,terpenoids`<-c(all_met$`biosynthesis of polyketides and nonribosomal peptides`,
                                                                           all_met$`biosynthesis of terpenoids and polyketides`)


all_met$`biosynthesis of polyketides and nonribosomal peptides`<-c()
all_met$`biosynthesis of terpenoids and polyketides`<-c()

all_met$`biosynthesis of secondary metabolites`<-c(all_met$`biosynthesis of secondary metabolites`,
                                                   all_met$`biosynthesis of other secondary metabolites`)
all_met$`biosynthesis of other secondary metabolites`<-c()

all_met$`glycan metabolism`<-c(all_met$`glycan metabolism`,
                               all_met$`glycan biosynthesis and metabolism`)

all_met$`glycan biosynthesis and metabolism`<-c()

names(all_met)<-gsub(x=names(all_met),pattern = 'metabolism',replacement = '')
names(all_met)<-gsub(x=names(all_met),pattern = ' of ',replacement = ' ')

names(all_met)<-paste('Metabolism',names(all_met),sep = '-')
names(all_met)<-gsub(x=names(all_met),pattern = 'Metabolism- ',replacement ='Metabolism-'  )

kegg_category_lib<-c(kegg_category_lib,all_met)

save(kegg_category_lib,file = 'inputData/annotaion/kegg_category_lib.RData')

# OLD # DICTIONARY FOR keg2category function
# library(keggorthology) ; library(igraph)
# data(KOgraph)
# g<-igraph::igraph.from.graphNEL(KOgraph)
# df<-as_data_frame(g,'both')
# v<-df$vertices
# edg<-df$edges
# colnames(edg)<-c('class','name')
# merged<-merge(edg,v,all=F)
# # extracting pathway  tags from the mappings
# kegg_category_lib<-list(c(merged[merged$class=='Cell Growth and Death','tag'],
#                           merged[merged$name== 'Cell Growth and Death','tag']),              #1."Cellular Processes, Cell Growth and Death"
#                         
#                         c(merged[merged$class=='Transport and Catabolism','tag'],
#                           merged[merged$name== 'Transport and Catabolism','tag']),           #2. "Cellular Processes, Transport and catabolism"
#                         
#                         c(merged[merged$class=='Translation','tag'],
#                           merged[merged$name== 'Translation','tag']),                        #3. "Genetic Information Processing, Translation"
#                         
#                         c(merged[merged$class=='Folding, Sorting and Degradation','tag'],
#                           merged[merged$name== 'Folding, Sorting and Degradation','tag']),   #4."Genetic Information Processing, Folding, sorting and degradation"
#                         
#                         c(merged[merged$class=='Replication and Repair','tag'],
#                           merged[merged$name== 'Replication and Repair','tag']),             #5. "Genetic Information Processing, Replication and repair"
#                         
#                         c(merged[merged$class=='Transcription','tag'],
#                           merged[merged$name== 'Transcription','tag']),                      #6. "Genetic Information Processing, Transcription"
#                         
#                         c(merged$tag[merged$class %in% merged$name[merged$class=='Metabolism']],
#                           merged[merged$class=='Metabolism','tag'],
#                           merged[merged$name=='Metabolism','tag']),                          #7. metabolism
#                         
#                         c(merged$tag[merged$class %in% merged$name[merged$class=="Environmental Information Processing"]],
#                           merged[merged$class=='Environmental Information Processing','tag'],
#                           merged[merged$name=='Environmental Information Processing','tag']),#8. Environmental Information Processing
#                         
#                         c(merged$tag[merged$class %in% merged$name[merged$class=="Human Diseases"]],           
#                           merged[merged$class=='Human Diseases','tag'],
#                           merged[merged$name=='Human Diseases','tag']),                      #9.  Human Diseases
#                         
#                         c(merged[merged$class=='Cell Communication','tag'],
#                           merged[merged$name== 'Cell Communication','tag']),                #10. Cell communication
#                         
#                         c(merged$tag[merged$class %in% merged$name[merged$class=="Organismal Systems"]],
#                           merged[merged$class=='Organismal Systems','tag'],
#                           merged[merged$name=='Organismal Systems','tag'])                  #11. Organismal Systems
#                         
# )
# names(kegg_category_lib)<-c("Cellular Processes, Cell Growth and Death",                         #1
#                             "Cellular Processes, Transport and catabolism",                      #2
#                             "Genetic Information Processing, Translation",                       #3
#                             "Genetic Information Processing, Folding, sorting and degradation",  #4
#                             "Genetic Information Processing, Replication and repair",            #5
#                             "Genetic Information Processing, Transcription",                     #6
#                             "Metabolism",                                                        #7
#                             "Environmental Information Processing",                              #8
#                             "Human Diseases",                                                    #9
#                             "Cell Communication",                                                #10
#                             "Organismal Systems"                                                 #11
# )
# save(kegg_category_lib,file = 'inputData/annotaion/kegg_category_lib.RData')
# 

# 2.3. Reactome categorization -------------------------------------------------------------
# Dictionary for reactome2category function
library(data.table)
reactome = data.frame(fread('https://reactome.org/download/current/ReactomePathways.txt'))
rownames(reactome) = reactome[,1]; reactome = reactome[grep("HSA", rownames(reactome)),]
colnames(reactome)<-c('ID','term','Organism')
react_hie = data.frame(fread('https://reactome.org/download/current/ReactomePathwaysRelation.txt'))
react_hie = react_hie[grep("HSA", react_hie[,1]),]
react_top = data.frame(fread('https://reactome.org/download/current/Complex_2_Pathway_human.txt'))
react_top = react_top[,c(2,3)]
react_top = react_top[which(!duplicated(react_top$pathway)),]
save(reactome,react_hie,react_top,file = 'inputData/annotaion/reactome.RData')


