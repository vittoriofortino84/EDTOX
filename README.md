# EDTOX

Please click the follwoing link to download the [Supplementary Data Tables](https://github.com/vittoriofortino84/EDTOX/blob/master/Supplementary_Data_Tables.xlsx) for

*A toxicogenomic platform for system-level understanding and prediction of EDC-induced toxicity*
 
Article by *Amirhossein Sakhteman, Mario Failli, Jenni Kublbeck, Anna-Liisa Levonen and Vittorio Fortino\* (\*Corresponding author)*

## Description of the R scripts used in the pipeline

## Part I: Development of toxicogenomics-driven predictive model to identify and prioritize compounds with ED activity
### 1. Preparation of MIES, pathways and training datasets

|**R Script**|[MIEs_from_CTD.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_1_MIEs_from_CTD.R) |
| ------------- |--------------|
|**Input**|  http://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz tested  for the release of june 2020|
| |  http://ctdbase.org/reports/CTD_chemicals.csv.gz tested  for the release of june 2020|
|**Output**| A list of chemical-gene associations|  
|**Dependencies**| data.table, FactoMineR, factoextra|
|**Summary**|Preparation of a binary data matrix indicating molecular initiating events (MIEs) of compounds annotated in CTD. Performing multiple correspondence analysis on the resulting matrix uisng FactoMineR and factoextra to determine the most distant gene-compound association types. Selection of five types of interactions: reaction, binding, activity, expression and metabolic processing. Removing compounds with too many MIEs as outliars.|

|**R Script**| [Pathways_Download.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_2_Pathways_Download.R)
| ------------- |--------------|
|**Input**| https://aopwiki.org/downloads/aop-wiki-xml-2019-01-01.gz|
||https://aopwiki.org/downloads/aop_ke_ec.tsv|
|| https://reactome.org/download/current/ReactomePathways.txt|
|| http://ctdbase.org/reports/CTD_genes_pathways.csv.gz|
|| https://reactome.org/download/current/miRBase2Reactome_PE_All_Levels.txt|
|**Output**| A list of pathways extracted from GO, KEGG, wiki and Reactome|
|**Dependencies**| XML, GO.db, org.Hs.eg.db, GSA, msigdbr, rWikiPathways, magrittr, rjson, data.table|
|**Summary**| Pathways related to KEGG, REACTOME, MSIGDB, GO and WIKI with less than 200 genes are retrieved. GO terms are systeamtically linked to known Wiki-AOPs. 

|**R Script**| [TOXCAST_nuclear_receptors_coregulators.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_3_TOXCAST_nuclear_receptors_coregulators.R)|
| ------------- |--------------|
|**Input**| hitc_Matrix_190226.csv from ToxCast 3.1 https://www.epa.gov/chemical-research/exploring-toxcast-data-downloadable-data |
||Assay_Summary_190226.csv' from ToxCast 3.1 https://www.epa.gov/chemical-research/exploring-toxcast-data-downloadable-data |
||DSSTox_Identifiers_and_CASRN.xlsx from ToxCast 3.1 https://www.epa.gov/chemical-research/exploring-toxcast-data-downloadable-data |
|| Result of the script [MIEs_from_CTD.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_1_MIEs_from_CTD.R)|
|| http://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz |
|**Output**| A binary matrix of chemicals and their corresponding hitcalls for the assay endpoitnts related to nuclear receptors and their co-regulators|
|**Dependencies**|tidyr, dplyr,org.Hs.eg.db, readxl, data.table|
|**Summary**|The genes related to nuclear receptors and their co-regulators from domain experts and the NURSA are merged. Then, ToxCast assays targeting these genes are selected. 

|**R Script**|[EDC_Decoy_selection.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_4_EDC_Decoy_selection.R)|
| ------------- |--------------|
|**Input**| https://cb.imsc.res.in/deduct/images/Batch_Download/DEDuCT_ChemicalBasicInformation.csv|
|| http://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz |
|| Result of the script [MIEs_from_CTD.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_1_MIEs_from_CTD.R)|
|| Result of the script [TOXCAST_nuclear_receptors_coregulators.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_3_TOXCAST_nuclear_receptors_coregulators.R)|
|**Output**| A list of EDCs and decoys, and their corresponding MIEs|
|**Dependencies**| ggplot2, ggrepel, magrittr, data.table, dplyr, reshape, cluster|
|**Summary**| An initial list of EDCs was retrieved from DEDuCT. Then, ToxCast assays targeting ED-reated nuclear receptor and co-regulators are used to refine the initial selection of EDCs. The most significat in vitro assay endpoints for identifying compounds with ED activity are selected based on a statistical test aiming at identifying the most relevant directions in ED-gene interactions (e.g. up/down regulation). Moreover, EDCs that are inactive for all the significant assay endpoints are systemaatically removed from the inital list of EDCs. After selecting a reliable set of EDCs, pairwise jaccard distance between the MIEs of the identified EDCs and other compounds annotated in CTD was calculated.The compounds with the maximum Jaccard distance towards EDCs will be seleceted as negative controls (decoys).|

|**R Script**| [ToxCast_dictionaries.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_5_ToxCast_dictionaries.R) |
| ------------- |--------------|
|**Input**|hitc_Matrix_190226.csv from ToxCast 3.1 https://www.epa.gov/chemical-research/exploring-toxcast-data-downloadable-data|
|| DSSTox_Identifiers_and_CASRN.xlsx from ToxCast 3.1 https://www.epa.gov/chemical-research/exploring-toxcast-data-downloadable-data|
||Assay_Summary_190226.csv' from ToxCast 3.1 https://www.epa.gov/chemical-research/exploring-toxcast-data-downloadable-data |
|**Output**| Hitcall matrix retried from ToxCast |
|**Dependencies**| tidyr, dplyr, readxl|
|**Summary**|Preparation of a dictionary for ToxCast target genes and their corresponding endpoints. Conversion of ToxCast DSSTox_Identifiers to CAS registry identifiers and preparation of the final Hitcall matrix with all  ToxCast endpoints.|
<br/>
<br/>
<br/>



### 2. Compiling gene co-expression networks 


|**R Script**|[Drug_matrix_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_1_Drug_matrix_wTO.R)  |
| ------------- |--------------|
|**Input**| LFCs and annottations from https://www.ebi.ac.uk/biostudies/studies/S-DIXA-AN-009?query=S-DIXA-AN-009|
|| journal.pcbi.1004847.s026.XLS from https://pubmed.ncbi.nlm.nih.gov/27028627/|
|**Output**| A set of four gene networks derived from DrugMatrix|
|**Dependencies**| xlsx, doParallel, wTO|
|**Summary**|The control samples from the preprocessed and normalized LFC values of Drug Matrix for rat in vitro hepatocytes and rat in vivo were removed. Three exposure time points 1,3 and 5 days for in vivo and 1 day for in vitro were selected as data frames. The genes expressed in liver were selected and orthology mapping of the probe IDs to entrez gene values was done. Four gene co-expression networks from the data frames using wTO package with bootstrap resampling method were selected.|

|**R Script**| [TG_Gates_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_2_TG_Gates_wTO.R)|
| ------------- |--------------|
| **Input**| LFCs and annottations from https://www.ebi.ac.uk/biostudies/studies/S-DIXA-AN-005?query=S-DIXA-AN-005|
|| LFCs and annottations from https://www.ebi.ac.uk/biostudies/studies/S-DIXA-AN-004?query=S-DIXA-AN-004|
|| LFCs and annottations from https://www.ebi.ac.uk/biostudies/studies/S-DIXA-AN-007?query=S-DIXA-AN-007|
|| journal.pcbi.1004847.s026.XLS from https://pubmed.ncbi.nlm.nih.gov/27028627/|
|**Output**|A set of eight gene networks derived from TG-GATEs|
|**Dependencies**| XLSX, doParallel, wTO|
| **Summary**|Removing control samples from the preprocesses and normalized LFC values related to TG-Gates data source for rat in vitro, human invitro and rat in vivo. Selection of three dose levels (high, middle and low) and three time points (8, 15 and 29 days) from the LFC values related to TG-Gates rat in vivo.(6 data frames) Selection of 1 day time exposure related to human and rat in vitro LFC values. (two data frames)Selection of the genes expressed in the liver from each data frame and orthology mapping of probe IDs to gene entrez IDs. Compiling 8 gene co-expression networks from the resulting data frames using wTO package with bootstrap resampling method. |

|**R Script** |[LINCS_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_3_LINCS_wTO.R)
| ------------- |--------------|
|**Input**| LFCs and annottations for level 5 studies  from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742 (phase I)|
||LFCs and annottations for level 5 studies from from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742 (phase II)|
|| journal.pcbi.1004847.s026.XLS from https://pubmed.ncbi.nlm.nih.gov/27028627/|
|**Output**| 1 consensus gene network for LINCS from phase 1 and phase 2 studies|
|**Dependencies** |cmapR, doParallel, wTO|
| **Summary**|Normalized and preproessed LFC values from the level 5 of phase 1 and phase 2 LINCS data source will be used. Selection of cell line HEPG2 with expousre time of 24 hours from phase1 and phase 2 gene expression data in LINCS.Selection of the gene IDS which are expressed in the liver.Compiling 2 gene networks from phase 1 and phase 2 using wTO package with bootstrapping resampling method.Compiling one consensus gene co-expression netowrk from the overlapping genes of the two networks using wTO package.|

|**R Script**|[Consensus_Rat_in_vitro_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_4_Consensus_Rat_in_vitro_wTO.R)|
| ------------- |--------------|
|**Input**| Results of the scripts [Drug_matrix_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_1_Drug_matrix_wTO.R) and [TG_Gates_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_2_TG_Gates_wTO.R)|
|**Output**| 1 consensus gene network for Drug matrix and TG-Gates hepatocytes after 1 day treatment|
|**Dependencies**|wTO|
|**Summary**| Compiling a consensus gene co-expression network from the overlapping genes of the two networks related to in vitro rat from drug matrix and TG-GATEs using wTO package|

|**R Script**|[PPI_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_5_PPI_wTO.R)|  
| ------------- |--------------|
|**Input**| No input is needed|
|**Output**| PPI network with new combined score|
|**Dependencies**|STRINGdb, igraph, org.Hs.eg.db |
|**Summary**|Retrieving a protein-protein interaction network from stringDB. Mapping nodes to entreg gene IDs.|
<br/>
<br/>
<br/>

### 3. Optimization of toxicogenomics-driven gene networks network 

|**R Script**|[Drug_matrix_tuner.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/3_1_Drug_matrix_tuner.R)|
| ------------- |--------------|
|**Input**| result of the scripts [Drug_matrix_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_1_Drug_matrix_wTO.R), [EDC_Decoy_selection.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_4_EDC_Decoy_selection.R) |
|**Output**|A matrix of silhouette scores for different combination of Edges percentiles and sorted genes for Drug Matrix networks |
|**Dependencies**| dnet,igraph,cluster |
|**Summary**|The top %2, %3, %5 and  %10 edge portions are extracted from each network. Each network will be subjected to random walk with restart uisng the seeds related to the MIEs of EDCs and decoys (benchmark set). The 200,500,700 and 1000 top most visited genes will be extracted after the random walk. A binary vector will be genrated for each compound (EDCs and decoys) with 1 representing the gene in the list of top visited gene.The pairwise jaccard distance between the binary vector of each EDCs and decoy will be calculated. Using the jaccard dismilarity matrix and a vector representing the class of each componud as EDC or Decoy average silhouette score was calculated for EDCs.|

|**R Script**|[TG_GATEs_tuner.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/3_2_TG_GATEs_tuner.R)|
| ------------- |--------------|
|**Input**|result of the scripts [TG_Gates_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_2_TG_Gates_wTO.R), [EDC_Decoy_selection.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_4_EDC_Decoy_selection.R)  |
|**Output**|A matrix of silhouette scores for different combination of Edges percentiles and sorted genes for  TG-GATEs networks |
|**Dependencies**| dnet,igraph,cluster |
|**Summary**|The top %2, %3, %5 and  %10 edge portions  are extracted from each network.Each network will be subjected to random walk with restart uisng the seeds related to the MIEs of EDCs and decoys (benchmark set). The 200,500,700 and 1000 top most visited genes will be extracted after the random walk. A binary vector will be genrated for each compound (EDCs and decoys) with 1 representing the gene in the list of top visited gene.The pairwise jaccard distance beween the binary vector of each EDCs and decoy will be calculated. Using the jaccard dismilarity matrix and a vector representing the class of each componud as EDC or Decoy average silhouette score was calculated for EDCs.|

|**R Script**|[Consensus_tuner.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/3_3_Consensus_tuner.R)|
| ------------- |--------------|
|**Input**|result of the scripts [Consensus_Rat_in_vitro_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_4_Consensus_Rat_in_vitro_wTO.R), [LINCS_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_3_LINCS_wTO.R), [EDC_Decoy_selection.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_4_EDC_Decoy_selection.R) |
|**Output**|A matrix of silhouette scores for different combination of Edges percentiles and sorted genes for  consensus networks|
|**Dependencies**| dnet,igraph,cluster |
|**Summary**|The top %2, %3, %5 and  %10 edge portions  are extracted from each network.Each network will be subjected to random walk with restart uisng the seeds related to the MIEs of EDCs and decoys (benchmark set). The 200,500,700 and 1000 top most visited genes will be extracted after the random walk. A binary vector will be genrated for each compound (EDCs and decoys) with 1 representing the gene in the list of top visited gene.The pairwise jaccard distance beween the binary vector of each EDCs and decoy will be calculated. Using the jaccard dismilarity matrix and a vector representing the class of each componud as EDC or Decoy average silhouette score was calculated for EDCs.|

|**R Script**|[PPI_tuner.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/3_4_PPI_tuner.R)|
| ------------- |--------------|
|**Input**|result of the script [EDC_Decoy_selection.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_4_EDC_Decoy_selection.R) |
|**Output**|A matrix of silhouette scores for different combination of combined scores and sorted genes for PPI network |
|**Dependencies**|dnet,igraph,cluster,STRINGdb,org.Hs.eg.db |
|**Summary**| New ppi networks were compiled using the 0.6,0.65,0.7,0.75,0.8,0.85 values as the cutoffs for the combined score. All networks will be subjected to random walk with restart uisng the seeds related to the MIEs of EDCs and decoys (benchmark set). The 200,500,700 and 1000 top most visited genes will be extracted after the random walk. A binary vector will be genrated for each compound (EDCs and decoys) with 1 representing the gene in the list of top visited gene. The pairwise jaccard distance beween the binary vector of each EDCs and decoy will be calculated. Using the jaccard dismilarity matrix and a vector representing the class of each componud as EDC or Decoy average silhouette score was calculated for EDCs.|

|**R Script**|[pareto_solution_on_tuning_results.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/3_5_pareto_solution_on_tuning_results.R)| 
| ------------- |--------------|
|**Input**|Results of the scripts [PPI_tuner.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/3_4_PPI_tuner.R), [Consensus_tuner.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/3_3_Consensus_tuner.R), [TG_GATEs_tuner.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/3_2_TG_GATEs_tuner.R) ,  [TG_GATEs_tuner.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/3_2_TG_GATEs_tuner.R), [Drug_matrix_tuner.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/3_1_Drug_matrix_tuner.R)  |
|**Output**|The optimized genes and edges solutions for all networks based on pareto solution |
|**Dependencies**|rPref,knitr |
|**Summary**| Using pareto to select the subnetworks that (1) maiximize the distance between EDCs and decoys,(2) minimize the selected edges and maximize the combined score. |
<br/>
<br/>
<br/>

### 4. Random walk with restart on gene networks and fast gene set enrichment analysis (RWR-FGSEA) 

|**R Script**|[RWR_FGSEA_for_edc_decoys.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/4_1_RWR_FGSEA_for_edc_decoys.R)|
| ------------- |--------------|
|**Input**|The results of the scripts [Drug_matrix_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_1_Drug_matrix_wTO.R), [TG_Gates_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_2_TG_Gates_wTO.R), [LINCS_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_3_LINCS_wTO.R), [Consensus_Rat_in_vitro_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_4_Consensus_Rat_in_vitro_wTO.R), [PPI_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_5_PPI_wTO.R), [Pathways_Download.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_2_Pathways_Download.R), [EDC_Decoy_selection.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_4_EDC_Decoy_selection.R),[pareto_solution_on_tuning_results.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/3_5_pareto_solution_on_tuning_results.R) |
|**Output**|A matrix of pathway scores for EDCs and decoys for each network|
|**Dependencies**|fgsea,dnet,igraph,BiocParallel|
|**Summary**|Random walk with restart was performed in order to extend the initial set of EDC-MIEs retrieved from CTD. This operation is repeaatd for eaach compound and network. Then Fast Gene Set Enrichment Analysis was performed to identify the most responsive pathways of the genes selected with the RWR algorithm.|

|**R Script**|[RWR_FGSEA_for_all_compounds_in_CTD.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/4_2_RWR_FGSEA_for_all_compounds_in_CTD.R)|
| ------------- |--------------|
|**Input**|The results of the scripts [Drug_matrix_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_1_Drug_matrix_wTO.R), [TG_Gates_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_2_TG_Gates_wTO.R), [LINCS_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_3_LINCS_wTO.R), [Consensus_Rat_in_vitro_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_4_Consensus_Rat_in_vitro_wTO.R), [PPI_wTO.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/2_5_PPI_wTO.R), [Pathways_Download.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_2_Pathways_Download.R),[MIEs_from_CTD.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_1_MIEs_from_CTD.R)  ,[pareto_solution_on_tuning_results.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/3_5_pareto_solution_on_tuning_results.R)|
|**Output**|A matrix of pathway scores for 12k chemical in CTD for each network|
|**Dependencies**|fgsea,dnet,igraph,BiocParallel|
|**Summary**|Random walk with restart was performed in order to extend the initial set of EDC-MIEs retrieved from CTD. This operation is repeaatd for eaach compound and network. Then Fast Gene Set Enrichment Analysis was performed to identify the most responsive pathways of the genes selected with the RWR algorithm.|
<br/>
<br/>
<br/>

### 5. Training of elastic-net-based classifiers for the claassification of EDCs and decoys

|**R Script**|[manual_curation_of_pathways_as_features.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_1_manual_curation_of_pathways_as_features.R)| 
| ------------- |--------------|
|**Input**| The result of the script [Pathways_Download.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_2_Pathways_Download.R)|
|**Output**|A list of pathways to be used in machine learning|
|**Dependencies**| No dependencies|
|**Summary**| Elastic-net based classifiers were trained on a more restricted set of pathways, which excludes pathways related to viral, bacterial and radiation. Duplicated pathways or pathways with no genes expressed in liver were also removed.|


|**R Script**|[Preparation_of_training_datasets.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_2_Preparation_of_training_datasets.R)| 
| ------------- |--------------|
|**Input**|The results of the scripts [manual_curation_of_pathways_as_features.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_1_manual_curation_of_pathways_as_features.R),[RWR_FGSEA_for_all_compounds_in_CTD.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/4_2_RWR_FGSEA_for_all_compounds_in_CTD.R), [RWR_FGSEA_for_edc_decoys.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/4_1_RWR_FGSEA_for_edc_decoys.R),[EDC_Decoy_selection.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_4_EDC_Decoy_selection.R)   |
|**Output**| A list containing training set for each network. A list containing test set for each network|
|**Dependencies**|No dependencies|
|**Summary**|Preparation of a training dataset for each toxicogenomics network. Each training dataset consists of NES scores or simply patwahy activation scores, which were obtained by using the RWR-FGSEA pipeline; and a vector labels indicateing EDCs and decoys.|

|**R Script**|[glm_modeling.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_3_glm_modeling.R)| 
| ------------- |--------------|
|**Input**|The output of the scripts [Preparation_of_training_datasets.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_2_Preparation_of_training_datasets.R), [MIEs_from_CTD.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_1_MIEs_from_CTD.R), [EDC_Decoy_selection.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_4_EDC_Decoy_selection.R)   |
|**Output**|List objects with GLM coefficients and models for each network and MIEs level|
|**Dependencies**|caret, doParallel|
|**Summary**|Training elastic-net-based classifiers on the computed training dataset by using cross-fold validation .|

|**R Script**|[k_fold_cross_validation.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_4_k_fold_cross_validation.R)| 
| ------------- |--------------|
|**Input**|The output of the scripts [Preparation_of_training_datasets.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_2_Preparation_of_training_datasets.R), [MIEs_from_CTD.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_1_MIEs_from_CTD.R), [EDC_Decoy_selection.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_4_EDC_Decoy_selection.R)   |
|**Output**|List of objects including Accuracy, F1-scores, confusion matrix, specificity and sensitivity for each trained classifier|
|**Dependencies**|caret, doParallel|
|**Summary**|Repeated 5-fold cross-validation was performed for  be performed on all 15 models. (Pathway level). Repeated 5-fold cross validation will be performed on a binary data matrix of the genes as columns. Compounds of benchmark (EDC,decoy) as rows. (The genes related to MIEs are characterized as 1 in the binary matrix. (Gene level))|

|**R Script**|[comparing_cross_validation_across_all_layers_ANOVA.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_5_comparing_cross_validation_across_all_layers_ANOVA.R)| 
| ------------- |--------------|
|**Input**|The output of the script [k_fold_cross_validation.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_4_k_fold_cross_validation.R)|
|**Output**|list of ANOVA results between F1-scores and the boxplot of F1-scores|
|**Dependencies**|ggplot2|
|**Summary**|The F1 scores of the k-fold-cross validation will be compared using ANOVA. Boxplot will be used to represent the obtained F1 scors across all GLM models.|
<br/>
<br/>
<br/>

### 6. Selection of the most informative pathways by using the coeffient estimates and ROC analysis based on individual pathways.

|**R Script**|[Integration_of_coefficients_stabilties_NES_scores.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_6_Integration_of_coefficients_stabilties_NES_scores.R)| 
| ------------- |--------------|
|**Input**|The outputs of the scripts [k_fold_cross_validation.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_4_k_fold_cross_validation.R), [glm_modeling.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_3_glm_modeling.R), [Preparation_of_training_datasets.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_2_Preparation_of_training_datasets.R), [manual_curation_of_pathways_as_features.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_1_manual_curation_of_pathways_as_features.R) |
|**Output**|CSV and excel files with GLM coefs, ROC-AUCS, mean NES scores for each pathway and each network |
|**Dependencies**|PRROC|
|**Summary**|ROC analysis was used to verify whether single pathways, which are selected by the elastic-net process, can distinguish alone EDCs from decoys. The NES scores, elastic-net coefficients, ROC-AUCs, average of NES scores will be integrated across all data layers (suppl. data).|

|**R Script**|[NES_bubble_plot_MOA.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_7_NES_bubble_plot_MOA.R)| 
| ------------- |--------------|
|**Input**|The output of the script [Integration_of_coefficients_stabilties_NES_scores.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_6_Integration_of_coefficients_stabilties_NES_scores.R) |
|**Output**| graphical represntation of the results of [Integration_of_coefficients_stabilties_NES_scores.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_6_Integration_of_coefficients_stabilties_NES_scores.R) as bubble plot|
|**Dependencies**|ggplot2, RColorBrewer|
|**Summary**|The map of pathway activation scores and GLM coefficient for each pathway and network will be represented as bubble plot.The generated plot can be used to indicate the putative pathways as mode of action for the EDCs.|
<br/>
<br/>
<br/>

### 7. Compiling EDC class probilities for 12k compounds annotated in CTD 

|**R Script**|[prediction_all_compounds_class_probilities.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/6_1_prediction_all_compounds_class_probilities.R)|
| ------------- |--------------|
|**Input**| Output of the scripts [glm_modeling.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_3_glm_modeling.R), [RWR_FGSEA_for_all_compounds_in_CTD.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/4_2_RWR_FGSEA_for_all_compounds_in_CTD.R)|
|**Output**|Class probability of a compound to be EDC for each data layer|
|**Dependencies**| caret|
|**Summary**| GLM coefs for each network will be used to predict class probability of all compounds in CTD using their NES scores across different pathways.|

|**R Script**|[ROC_TOXCAST_vs_class_probabilities.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/6_2_ROC_TOXCAST_vs_class_probabilities.R)|  
| ------------- |--------------|
|**Input**|Output of the scripts [TOXCAST_nuclear_receptors_coregulators.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_3_TOXCAST_nuclear_receptors_coregulators.R), [prediction_all_compounds_class_probilities.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/6_1_prediction_all_compounds_class_probilities.R)|
|**Output**|A matrix of area under the curve for ROC between class probability of each data layer and Hitcall result of ToxCast endpoints |
|**Dependencies**|PRROC|
|**Summary**| ROC analysis will be used between the class probabilty of each network and the binary in vitro experimental hitcall ToxCast assays related to nuclear receptors for each network.- Selection of most informative networks based on the results of ROC curve analysis.|


|**R Script**|[Developing_Harmonic_and_average_EDC_scores.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/6_3_Developing_Harmonic_and_average_EDC_scores.R)|
| ------------- |--------------|
|**Input**| Output of the scripts [EDC_Decoy_selection.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_4_EDC_Decoy_selection.R), [ROC_TOXCAST_vs_class_probabilities.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/6_2_ROC_TOXCAST_vs_class_probabilities.R), [prediction_all_compounds_class_probilities.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/6_1_prediction_all_compounds_class_probilities.R)   |
|**Output**|vectors of harmonic scores and average scores for 12K compounds in CTD|
|**Dependencies**|ggplot2, ggrepel|
|**Summary**|Average EDC score is defined as the average of class probabilties across selected networks for each compound.Harmonic EDC score is defined as the harmonic sum class probabilties across selected networks for each compound.|
<br/>
<br/>
<br/>


### 8. Prediction of EDC class probilities and EDC scores for a new compound using its MIEs

|**R Script**|[edc_score_new_comps.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/edc_score_new_comps.R)|
| ------------- |--------------|
|**Input**|A list of compounds and their corresponding MIEs |
|**Output**|Class probability of a compound to be EDC for each data layer|
|**Dependencies**| caret, fgsea, dnet, igraph|
|**Summary**| The MIEs of the compounds will be used as seeds to start random walk with restart and enrichment analysis on  different networks. The resulting pathway scores will be used to calculate the EDC scores and the class probabilities for the compounds. |
<br/>
<br/>
<br/>

# Part II: Evaluation and validation of EDC scores
|**R Script**|[ROC_analysis_EDC_scores_vs_TOXCAST_endpoints.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/6_4_ROC_analysis_EDC_scores_vs_TOXCAST_endpoints.R)|  
| ------------- |--------------|
|**Input**| Output of the scripts Output of the scripts [TOXCAST_nuclear_receptors_coregulators.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_3_TOXCAST_nuclear_receptors_coregulators.R), [Developing_Harmonic_and_average_EDC_scores.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/6_3_Developing_Harmonic_and_average_EDC_scores.R), [ROC_TOXCAST_vs_class_probabilities.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/6_2_ROC_TOXCAST_vs_class_probabilities.R)|
|**Output**|Matrix of Area under the curve of ROC between harmonic and average EDC scores with ToxCast Hitcall, P_values of ROC proportion test |
|**Dependencies**|PRROC, gplots |
|**Summary**|Validating the compiled EDC scores by ROC analysis between  EDC scores and ToxCast hitcall data for each assay endpoint.|

|**R Script**|[Comparison_of_VAM_on_DeDuCt.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/7_1_Comparison_of_VAM_on_DeDuCt.R)| 
| ------------- |--------------|
|**Input**|https://cb.imsc.res.in/deduct/images/Batch_Download/DEDuCT_ChemicalBasicInformation.csv, Output of the scripts [EDC_Decoy_selection.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_4_EDC_Decoy_selection.R), [prediction_all_compounds_class_probilities.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/6_1_prediction_all_compounds_class_probilities.R) |
|**Output**| Excel table with harmonic scores and the compounds in DeDUCT |
|**Dependencies**| xlsx|
|**Summary**|Evaluation of the scores for the compounds in DeDUCT |



|**R Script**|[Comparison_VAM_scores_vs_TOXPI_scores.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/7_2_Comparison_VAM_scores_vs_TOXPI_scores.R)| 
| ------------- |--------------|
|**Input**| Suppl. file of http://dx.doi.org/10.1016/j.coph.2014.09.021, https://cb.imsc.res.in/deduct/images/Batch_Download/DEDuCT_ChemicalBasicInformation.csv, Output of the script [prediction_all_compounds_class_probilities.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/6_1_prediction_all_compounds_class_probilities.R)|
|**Output**| plot of the Toxpi scores VS EDC scores, excel files of toxpi and edc scores |
|**Dependencies**| ggplot2, ggrepel|
|**Summary**|Evaluation of the EDC scores with ToxPi scores (The scores developed from ToxCast assay endpoints) |


|**R Script**|[Comparison_pathway_scores_vs_TOXDB_pathway_scores.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/7_3_Comparison_pathway_scores_vs_TOXDB_pathway_scores.R)|  
| ------------- |--------------|
|**Input**| ToxDB pathway scores for pathways aryl hydrocarbon receptor,  Breaset cancer and estrogen receptor from http://toxdb.molgen.mpg.de/ , Output of the script [RWR_FGSEA_for_all_compounds_in_CTD.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/4_2_RWR_FGSEA_for_all_compounds_in_CTD.R), |
|**Output**| CSV file with pathway scores of ToxDB and the pathway scores resutling from RWR-FGSEA|
|**Dependencies**|webchem, ggplot2, ggrepel |
|**Summary**| CHEMBL ids will be mapped to CAS ids using webchem. Comparison between EDC scores with pathway scores from TOXDB for the pathways aryl hydrocarbon receptor, Breaset cancer and estrogen receptor|




|**R Script**|[Validation_VAM_scores_Eurion_External_set_compounds.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/8_1_Validation_VAM_scores_Eurion_External_set_compounds.R)|
| ------------- |--------------|
|**Input**|https://cb.imsc.res.in/deduct/images/Batch_Download/DEDuCT_ChemicalBasicInformation.csv, output of the script [prediction_all_compounds_class_probilities.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/6_1_prediction_all_compounds_class_probilities.R) |
|**Output**| Excel file and heatmap plot of EDC scores for external validation set|
|**Dependencies**|reshape2,RColorBrewer,ggplot,magrittr |
|**Summary**|Evaluation of EDC scores using known EDCs from expert domain and DeDuCT. Calculation of accuaracy for the EDC scores.Representing the validation set compounds as heatmap plot| 





|**R Script**|[Validation_with_ToxCast_mies.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/8_3_Validation_with_ToxCast_mies.R)| 
| ------------- |--------------|
|**Input**|Output of the script  [ToxCast_dictionaries.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_5_ToxCast_dictionaries.R), [MIEs_from_CTD.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_1_MIEs_from_CTD.R),  |
|**Output**|Plot of predicted EDC scores from ToxCast MIES VS EDC scores predicted from CTD MIES|
|**Dependencies**|ggplot2, ggrepel, knitr |
|**Summary**|The target genes for the compounds with positive assay endpoint in ToxCast wil be considered as the MIEs of the compounds.Random walk with resetart and fast gene set enrichment analysis pipeline will be repeated starting with the MIEs from ToxCast. The obtained NES scores will be subjected to GLM models and the resulting class probabilities will be compared with the class probabilities of the MIEs from CTD. |

<br/>
<br/>
<br/>

# Part III: Sensitivity and enrichment analysis
|**R Script**|[dose_response_TG_GATEs_single.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/8_4_dose_response_TG_GATEs_single.R)| 
| ------------- |--------------|
|**Input**|Output of the scripts [Developing_Harmonic_and_average_EDC_scores.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/6_3_Developing_Harmonic_and_average_EDC_scores.R), [RWR_FGSEA_for_all_compounds_in_CTD.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/4_2_RWR_FGSEA_for_all_compounds_in_CTD.R),[Integration_of_coefficients_stabilties_NES_scores.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_6_Integration_of_coefficients_stabilties_NES_scores.R)|
|**Output**|List of the compounds based on the patterns of increase and decrease, List of ANOVA results for pathways in different dose response scenarios, plot depictions|
|**Dependencies**|dplyr,ggplot2|
|**Summary**|TG-Gates class probability profile of the compounds will be categorized as 4 different groups (based on increase or decrease pattern at different doses Low, Middle and High) The significant pathways based on increases or decrease of pathway activation scores will be determined using ANOVA followed by post tukey test with bootstrap sampling method.|

|**R Script**|[time_exposure_response_TG_GATEs_repeated.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/8_5_time_exposure_response_TG_GATEs_repeated.R)|
| ------------- |--------------|
|**Input**|Output of the scripts [Developing_Harmonic_and_average_EDC_scores.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/6_3_Developing_Harmonic_and_average_EDC_scores.R), [RWR_FGSEA_for_all_compounds_in_CTD.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/4_2_RWR_FGSEA_for_all_compounds_in_CTD.R),[Integration_of_coefficients_stabilties_NES_scores.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_6_Integration_of_coefficients_stabilties_NES_scores.R)|
|**Output**|List of the compounds based on the patterns of increase and decrease, List of ANOVA results for pathways in different time of exposure scenarios, Plot depictions|
|**Dependencies**|dplyr,ggplot2|
|**Summary**|TG-Gates class probability profile of the compounds will be categorized as 4 different patterns (based on increase or decrease pattern at different Time of exposure 8 days, 15 days and 29 days). The significant pathways based on increases or decrease of pathway activation scores will be determined using ANOVA followed by post tukey test with bootstrap sampling method|

|**R Script**|[textmining_16_april_2020_endocrine_disruption.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/8_6_textmining_16_april_2020_endocrine_disruption.R)|
| ------------- |--------------|
|**Input**|The abstract results of pubmed search for endocrine disruption |
|**Output**|A data frame with entrez gene names and their frequencies across different papers |
|**Dependencies**|pubmed.mineR, magrittr|
|**Summary**|After the search in pubmed for the term endocrine disruption (ED) the genes related to ED will be determined using text mining with pubmed.mineR.|

|**R Script**|[time_dose_plots_heatplot.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/8_7_time_dose_plots_heatplot.R)|
| ------------- |--------------|
|**Input**|Output of the scripts [time_exposure_response_TG_GATEs_repeated.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/8_5_time_exposure_response_TG_GATEs_repeated.R), [dose_response_TG_GATEs_single.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/8_4_dose_response_TG_GATEs_single.R),[Pathways_Download.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_2_Pathways_Download.R) |
|**Output**|Plots of time of exposure and dose for sensitivity analysis, heatplot eseential genes and pathways|
|**Dependencies**|ggpubr, RColorBrewer, ggplot2|
|**Summary**|The pathways versus thier genes will be represented as two heatplots (dose and time point).The highly cited gened will be determined in the plot.The essential genes will be highlighted in the plot.|

|**R Script**|[categorization_enrichment_result.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/8_8_categorization_enrichment_result.R)|
| ------------- |--------------|
|**Input**|Output of the script [Integration_of_coefficients_stabilties_NES_scores.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/5_6_Integration_of_coefficients_stabilties_NES_scores.R), |
|**Output**|plots of categorization for specific pathways related to KEGG, REACTOME and wiki AOP|
|**Dependencies**|UpSetR, ggplot2, RColorBrewer |
|**Summary**|For each data layer (network) the number selected KEGG pathways by elastic GLM will be categorized and represented in plots.For each data layer (network) the number selected REACTOME pathways by elastic GLM will be categorized and represented in plots.|
<br/>
<br/>
<br/>


# Part IV: Linking MIEs to Adverse outcomes
### 1. Preparation of the training benchmark set 
|**R Script**|[preparation_data_set_disease_biomarker.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_1_preparation_data_set_disease_biomarker.R)|  
| ------------- |--------------|
|**Input**|Output of the scripts [Pathways_Download.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_2_Pathways_Download.R), [EDC_Decoy_selection.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_4_EDC_Decoy_selection.R),[RWR_FGSEA_for_all_compounds_in_CTD.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/4_2_RWR_FGSEA_for_all_compounds_in_CTD.R)  |
|**Output**|A list of pathway scores for each data layer with specefic pathways|
|**Dependencies**|No dependencies|
|**Summary**|Retaining the pathways with max length of 50 genes for biomarker discovery study. Preparation of the NES scores for each data layer based on the selected pathways. |

|**R Script**|[chem2disease_CTD.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_2_chem2disease_CTD.R)| 
| ------------- |--------------|
|**Input**|http://ctdbase.org/reports/CTD_chemicals_diseases.csv.gz, output of the script  [EDC_Decoy_selection.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/1_4_EDC_Decoy_selection.R),  |
|**Output**|A binary matrix of chemicals and diseases with 1 representing the association of a disease and the chemical|
|**Dependencies**|xlsx|
|**Summary**|Compiling a binary matrix from chem-disease associations in CTD. The rows of the matrix are the compound names and the columns are the disease names. The compound disease associations are represented by 1 in the matrix.|


|**R Script** |[preparation_two_class_training_set_by_disease_name.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_3_preparation_two_class_training_set_by_disease_name.R) |
| ------------- |--------------|
|**Input**|Output of the scripts [preparation_data_set_disease_biomarker.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_1_preparation_data_set_disease_biomarker.R) |
|**Output**|A list object with training x data and y vector for each data layer|
|**Dependencies**|No dependencies|
|**Summary**|Preparation of training set (NES scores and class labels) for atherosclerosis, metabolic syndrome and diabetes type 2|
<br/>
<br/>
<br/>


### 2. Training elastic-net GLM classifier for the metabolic syndrome, diabetes type2 and atherosclerosis 
|**R Script** |[metabolic_syndrome_two_class_traning_models_glm.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_4_1_metabolic_syndrome_two_class_traning_models_glm.R) |
| ------------- |--------------|
|**Input**|Output of the script [preparation_two_class_training_set_by_disease_name.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_3_preparation_two_class_training_set_by_disease_name.R) |
|**Output**|A list classification models and coeffecients for each data layer|
|**Dependencies**|caret, doParallel|
|**Summary**|Training GLM  models of metabolic syndrome for all 15 data layers|

|**R Script** | [artherosclerosis_two_class_training_models_glm.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_4_2_artherosclerosis_two_class_training_models_glm.R)|
| ------------- |--------------|
|**Input**|Output of the script [preparation_two_class_training_set_by_disease_name.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_3_preparation_two_class_training_set_by_disease_name.R) |
|**Output**|A list classification models and coeffecients for each data layer|
|**Dependencies**|caret, doParallel|
|**Summary**|Training GLM  models of atherosclerosis for all 15 data layers|

|**R Script** | [diabetes_2_two_class_traning_models_glm.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_4_3_diabetes_2_two_class_traning_models_glm.R)|
| ------------- |--------------|
|**Input**|Output of the script [preparation_two_class_training_set_by_disease_name.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_3_preparation_two_class_training_set_by_disease_name.R) |
|**Output**|A list classification models and coeffecients for each data layer|
|**Dependencies**|caret, doParallel|
|**Summary**|Training GLM  models of diabetes T2 for all 15 data layers|


|**R Script**| [kfold_CV_two_class_metabolic_syndrome_glm.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_5_1_kfold_CV_two_class_metabolic_syndrome_glm.R)|
| ------------- |--------------|
|**Input**|Output of the script [preparation_two_class_training_set_by_disease_name.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_3_preparation_two_class_training_set_by_disease_name.R) |
|**Output**|List of accuaracy metrics for the classification task for each data layer|
|**Dependencies**|caret, doParallel|
|**Summary**|Performing K-fold-Cross validation on models of metabolic syndrome for all 15 data layers|

|**R Script** |[kfold_CV_two_class_artherosclerosis_glm.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_5_2_kfold_CV_two_class_artherosclerosis_glm.R) |
| ------------- |--------------|
|**Input**|Output of the script [preparation_two_class_training_set_by_disease_name.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_3_preparation_two_class_training_set_by_disease_name.R) |
|**Output**|List of accuaracy metrics for the classification task for each data layer|
|**Dependencies**|caret, doParallel|
|**Summary**|Performing K-fold-Cross validation on models of atherosclerosis for all 15 data layers|

|**R Script** |[kfold_CV_two_class_diabetes_2_glm.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_5_3_kfold_CV_two_class_diabetes_2_glm.R) |
| ------------- |--------------|
|**Input**|Output of the script [preparation_two_class_training_set_by_disease_name.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_3_preparation_two_class_training_set_by_disease_name.R) |
|**Output**|List of accuaracy metrics for the classification task for each data layer|
|**Dependencies**|caret, doParallel|
|**Summary**|Performing K-fold-Cross validation on models of diabetes type 2 for all 15 data layers|

|**R Script** |[comparison_CV_F1_scores_all_layers_ANOVA_metabolic_syndrome.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_6_1_comparison_CV_F1_scores_all_layers_ANOVA_metabolic_syndrome.R) |
| ------------- |--------------|
|**Input**|Output of the script [kfold_CV_two_class_metabolic_syndrome_glm.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_5_1_kfold_CV_two_class_metabolic_syndrome_glm.R) |
|**Output**|A list of ANOVA results between different data layers and boxplot representation|
|**Dependencies**|plotly, dplyr, ggplot2|
|**Summary**| The F1 scores of the k-fold-cross validation will be compared using ANOVA across all 15 data layer for metabolic syndrome.Boxplot will be used to represent the obtained F1 scors across all GLM models for metabolic syndrome.|

|**R Script**|[comparison_CV_F1_scores_all_layers_ANOVA_atherosclerosis.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_6_2_comparison_CV_F1_scores_all_layers_ANOVA_atherosclerosis.R) |
| ------------- |--------------|
|**Input**|Output of the script [kfold_CV_two_class_artherosclerosis_glm.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_5_2_kfold_CV_two_class_artherosclerosis_glm.R)|
|**Output**|A list of ANOVA results between different data layers and boxplot representation|
|**Dependencies**|plotly, dplyr, ggplot2|
|**Summary**|The F1 scores of the k-fold-cross validation will be compared using ANOVA across all 15 data layer for atherosclerosis. Boxplot will be used to represent the obtained F1 scors across all GLM models for atherosclerosis.|

|**R Script** |[comparison_CV_F1_scores_all_layers_ANOVA_diabetes.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_6_3_comparison_CV_F1_scores_all_layers_ANOVA_diabetes.R)  |
| ------------- |--------------|
|**Input**|Output of the script [kfold_CV_two_class_diabetes_2_glm.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_5_3_kfold_CV_two_class_diabetes_2_glm.R)|
|**Output**|A list of ANOVA results between different data layers and boxplot representation|
|**Dependencies**|plotly, dplyr, ggplot2|
|**Summary**|The F1 scores of the k-fold-cross validation will be compared using ANOVA across all 15 data layer for diabetes type 2. Boxplot will be used to represent the obtained F1 scors across all GLM models for diyabetes type 2.|
<br/>
<br/>
<br/>

### 3. Univariate analysis and integration with GLM resutls
|**R Script** |[Integration_of_glm_coefs_stabilities_NES_ROC_metabolic_syndrome.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_7_1_Integration_of_glm_coefs_stabilities_NES_ROC_metabolic_syndrome.R) |
| ------------- |--------------|
|**Input**|Output of the scripts [preparation_two_class_training_set_by_disease_name.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_3_preparation_two_class_training_set_by_disease_name.R), [kfold_CV_two_class_metabolic_syndrome_glm.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_5_1_kfold_CV_two_class_metabolic_syndrome_glm.R), [metabolic_syndrome_two_class_traning_models_glm.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_4_1_metabolic_syndrome_two_class_traning_models_glm.R)   |
|**Output**|A data frame with stability index, AUC-ROC and GLM coef for each pathway across all data layers|
|**Dependencies**|magrittr, PRROC|
|**Summary**|For each network, ROC analysis will be used as a univariate method to evaluate the NES scores for each pathway for EDCs leading to metabolic synndrome vs edcs not leading to metabolic syndrome. The NES scores, glm coefficients ROC-AUCs, average of NES scores will be integrated across all data layers (suppl. data).|

|**R Script** |[Integration_of_glm_coefs_stabilities_NES_ROC_atherosclerosis.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_7_2_Integration_of_glm_coefs_stabilities_NES_ROC_atherosclerosis.R) |
| ------------- |--------------|
|**Input**|Output of the scripts [preparation_two_class_training_set_by_disease_name.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_3_preparation_two_class_training_set_by_disease_name.R), [kfold_CV_two_class_artherosclerosis_glm.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_5_2_kfold_CV_two_class_artherosclerosis_glm.R), [artherosclerosis_two_class_training_models_glm.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_4_2_artherosclerosis_two_class_training_models_glm.R)|
|**Output**|A data frame with stability index, AUC-ROC and GLM coef for each pathway across all data layers|
|**Dependencies**|magrittr, PRROC|
|**Summary**|For each network, ROC analysis will be used as a univariate method to evaluate the NES scores for each pathway for EDCs leading to atherosclerosis vs edcs not leading to atherosclerosis. The NES scores, glm coefficients ROC-AUCs, average of NES scores will be integrated across all data layers (suppl. data).|

|**R Script** | [Integration_of_glm_coefs_stabilities_NES_ROC_diabetes_2.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_7_3_Integration_of_glm_coefs_stabilities_NES_ROC_diabetes_2.R)|
| ------------- |--------------|
|**Input**|Output of the scripts [preparation_two_class_training_set_by_disease_name.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_3_preparation_two_class_training_set_by_disease_name.R), [kfold_CV_two_class_diabetes_2_glm.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_5_3_kfold_CV_two_class_diabetes_2_glm.R), [diabetes_2_two_class_traning_models_glm.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_4_3_diabetes_2_two_class_traning_models_glm.R)|
|**Output**|A data frame with stability index, AUC-ROC and GLM coef for each pathway across all data layers|
|**Dependencies**|magrittr, PRROC|
|**Summary**|For each network, ROC analysis will be used as a univariate method to evaluate the NES scores for each pathway for EDCs leading to diabetes type 2 vs edcs not leading to diabetes type 2.The NES scores, glm coefficients ROC-AUCs, average of NES scores will be integrated across all data layers (suppl. data).|


|**R Script** |[bubble_plots_pathways.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_8_bubble_plots_pathways.R)   |
| ------------- |--------------|
|**Input**|The output of the scripts [Integration_of_glm_coefs_stabilities_NES_ROC_diabetes_2.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_7_3_Integration_of_glm_coefs_stabilities_NES_ROC_diabetes_2.R), [Integration_of_glm_coefs_stabilities_NES_ROC_atherosclerosis.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_7_2_Integration_of_glm_coefs_stabilities_NES_ROC_atherosclerosis.R), [Integration_of_glm_coefs_stabilities_NES_ROC_metabolic_syndrome.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_7_1_Integration_of_glm_coefs_stabilities_NES_ROC_metabolic_syndrome.R)   |
|**Output**|Bubble plot for each disease across all data layers representing the pathways and their coeffeicients|
|**Dependencies**|ggplot2, ggarrange|
|**Summary**|The map of pathway activation scores and GLM coefficient for each pathway and network will be represented as bubble plots for the three diseases metabolic syndrome, atherosclerosis and diabetes type 2. The generated plots can be used to indicate the putative pathways linking EDCs to adverse outcome.|
<br/>
<br/>
<br/>


### 4. Compiling and validation of  disease scores for atherosclerosis, metabolic syndrome and diabetes type 2 
|**R Script** |[disease_scores.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_9_disease_scores.R) |
| ------------- |--------------|
|**Input**|Output of the scripts [preparation_data_set_disease_biomarker.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_1_preparation_data_set_disease_biomarker.R), [chem2disease_CTD.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_2_chem2disease_CTD.R), [metabolic_syndrome_two_class_traning_models_glm.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_4_1_metabolic_syndrome_two_class_traning_models_glm.R),  [artherosclerosis_two_class_training_models_glm.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_4_2_artherosclerosis_two_class_training_models_glm.R), [diabetes_2_two_class_traning_models_glm.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_4_3_diabetes_2_two_class_traning_models_glm.R) |
|**Output**|CSV files with the disease scores related to chemical in CTD for each data layer, ROC aucs of the predicted scores with disease vectors in CTD|
|**Dependencies**|PRROC|
|**Summary**|For each disease (metabolic syndrome, atherosclerosis and diabetes type 2), harmonic sum and average score will be compiled across all 15 networks for 12k compounds in CTD. A ROC analysis will be performed between the predicted scores and the results in CTD|

|**R Script**|[disease_scores_pie_chart.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_10_disease_scores_pie_chart.R) |
| ------------- |--------------|
|**Input**|Output of the script [disease_scores.R](https://github.com/vittoriofortino84/EDTOX/blob/master/scripts/9_9_disease_scores.R) |
|**Output**|Piechart of disease scores|
|**Dependencies**|ggplot2|
|**Summary**|For each disease (metabolic syndrome, atherosclerosis and diabetes type 2) a pie chart will be depicted to reveal distribution of disease scores across 12k compounds in CTD.|

