# EDTOX

Please click the follwoing link to download the [Supplementary Data Tables](https://github.com/vittoriofortino84/EDTOX/blob/master/Supplementary_Data_Tables.xlsx) for

*A toxicogenomic platform for system-level understanding and prediction of EDC-induced toxicity*
 
Article by *Amirhossein Sakhteman, Mario Failli, Jenni Kublbeck, Anna-Liisa Levonen and Vittorio Fortino\* (\*Corresponding author)*


## Description of the Supplementary data Tables 

* Suppl. Tab. 1 – List of ToxCast assay endpoints related to nuclear receptors.
* Suppl. Tab. 2 – Testing of proportions to identify ToxCast assay endpoints for the selection of EDCs.
* Suppl. Tab. 3 – The result of ANOVA and Tukey post hoc tests for the comparison of test accuracy scores of EDC-based classifiers.
* Suppl. Tab. 4 – Predicted EDC probability scores for more than 10K compounds annotated in CTD.
* Suppl. Tab. 5 – Predicted EDC-Atherosclerosis scores for compounds annotated in CTD.
* Suppl. Tab. 6 – Predicted EDC-Metabolic-Syndrome scores for compounds annotated in CTD.
* Suppl. Tab. 7 – Predicted EDC-Type2-Diabetes scores for compounds annotated in CTD.
* Suppl. Tab. 8 – List of relevant EDC-MIEs and related pathways.
* Suppl. Tab. 9 – List of informative genes that are currently not annotated as EDC-MIEs in CTD. 
* Suppl. Tab. 10 – Pathway activation scores generated for the classification between EDCs and decoys (or negative controls).
* Suppl. Tab. 11 – Pathway activation scores generated for the classification of EDCs leading to Atherosclerosis.
* Suppl. Tab. 12 – Pathway activation scores generated for the classification of EDCs leading to Metabolic Syndrome.
* Suppl. Tab. 13 – Pathway activation scores generated for the classification of EDCs leading to Type2 Diabetes.
* Suppl. Tab. 14– Pathway selected based on AUC values greater than 0.7. 
* Suppl. Tab. 15 – ROC-curve-based analysis between EDC-class probabilities across different data layers and ToxCast assay endpoints.
* Suppl. Tab. 16 – Predicted EDC scores for the compounds in DEDuCT.
* Suppl. Tab. 17 – Predicted EDC scores for the compounds selected by domain experts.
* Suppl. Tab. 18 – Comparison between predicted EDC scores and ToxPi scores.
* Suppl. Tab. 19 – EDC Class probabilities of in vivo models with MIEs from in vitro ToxCast endpoints.
* Suppl. Tab. 20 – Comparison between predicted EDC scores and ToxDB scores.


## Description of the R scripts used in the pipeline
# Part I: Developing EDC scores
### 1. Preparation of the MIES, pathways and training benchmark set 

#### Script '1_1_MIEs_from_CTD.R'
-  input:  http://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz 
-  output: A list object for the chemical and their related MIEs (genes)  
-  needed libraries data.table, FactoMineR, factoextra
- Preparation of a binary data matrix for molecular initiating events (MIEs) from compound-gene interactions in CTD. 
- The interactions subtypes related to metabolism were grouped as metabolism and the interaction types related to transport are grouped as transport.
- Performing multiple correspondence analysis on the resulting matrix uisng FactoMineR and factoextra. 
- Selection of reaction,binding,activity,expression,metabolic processing as the more distant types of the interaction based on the plot of MCA.
- For the compounds with more than 50 gene interactions the less informative gene interactions will be removed.

#### Script 1_2_Pathways_Download.R 
- input: https://aopwiki.org/downloads/aop-wiki-xml-2019-01-01.gz
-        https://aopwiki.org/downloads/aop_ke_ec.tsv
-        https://reactome.org/download/current/ReactomePathways.txt
-        http://ctdbase.org/reports/CTD_genes_pathways.csv.gz
-        https://reactome.org/download/current/miRBase2Reactome_PE_All_Levels.txt
- Pathways related to KEGG, REACTOME,MSIGDB, GO and WIKI with the size of less than 200 will be retrieved.
- A binary dictionary to link the GO terms with Wiki-AOPs will be generated.
- The classifications tags for the pathways related to KEGG and REACTOME pathways will be downloaded and preprocessed for enrichment analysis.

#### Script 1_3_TOXCAST_nuclear_receptors_coregulators.R 
- The genes related to nuclear receptors and their co-regulators from experts domain and NURSA will be merged.
- The target gene ids from ToxCast will be extracted.
- The ToxCast assay endpoints which their target genes are in the list of nuclear receptor genes will be saved as endpoints related to nuclear receptor.

#### Script 1_4_EDC_Decoy_selection.R 
- The list of EDCs will be retrieved from DEDuCT as CAS ids.
- The ToxCast assay endpoints related to nuclear receptor and co-regulators of EDCs will be extracted. 
- The most significat in vitro assay endpoints for the mechanism of EDCs will be characterized using statistical proportion test.
- EDCs (DEDuCT list) which are incative for all the significant assay endpoints will be removed from the final list of EDcs.
- Pairwise jaccard distance between the MIEs related to remaining EDCs and other compounds in CTD will be calculated.
- The compounds with the maximum Jaccard distance with EDCs will be seleceted as negative controls (decoys).

#### 1_5_ToxCast_dictionaries.R
- Preparation of a dictionary for ToxCast target genes and their corresponding endpoints.
- Conversion of ToxCast DSSTox_Identifiers to CAS registry identifiers and preparation of the final Hitcall list for ToxCast.

### 2. Generating gene co-expression networks 


#### Script 2_1_Drug_matrix_wTO.R                                    
- Removing the control samples from the preprocessed and normalized LFC values related to Drug Matrix data source for rat in vitro hepatocytes and rat in vivo.
- Selection of the three exposure time points 1,3 and 5 days for in vivo and 1 day for in vitro and splitting the data as four data frames.
- Selection of the genes expressed in liver and orthology mapping of the probe IDs to entrez gene values
- Compiling 4 gene co-expression networks from the data frames using wTO package with bootstrap resampling method.

#### Script 2_2_TG_Gates_wTO.R  
- Removing control samples from the preprocesses and normalized LFC values related to TG-Gates data source for rat in vitro, human invitro and rat in vivo.
- Selection of three dose levels (high, middle and low) and three time points (8, 15 and 29 days) from the LFC values related to TG-Gates rat in vivo. (6 data frames)
- Selection of 1 day time exposure related to human and rat in vitro LFC values. (two data frames)
- Selection of the genes expressed in the liver from each data frame and orthology mapping of probe IDs to gene entrez IDs.
- Compiling 8 gene co-expression networks from the resulting data frames using wTO package with bootstrap resampling method. 

#### Script 2_3_LINCS_wTO.R
- Normalized and preproessed LFC values from the level 5 of phase 1 and phase 2 LINCS data source will be used.
- Selection of cell line HEPG2 with expousre time of 24 hours from phase1 and phase 2 gene expression data in LINCS.
- Selection of the gene IDS which are expressed in the liver.
- Compiling 2 gene networks from phase 1 and phase 2 using wTO package with bootstrapping resampling method.
- Compiling one consensus gene co-expression netowrk from the overlapping genes of the two networks using wTO package.

#### Script 2_4_Consensus_Rat_in_vitro_wTO.R 
- Compiling one consensus gene co-expression netowrk from the overlapping genes of the two networks related to in vitro rat from drug matrix and TG-GATEs using wTO package.

#### Script 2_5_PPI_wTO.R   
- Retrieving protein protein interaction network from stringDB.
- Mapping nodes to entreg gene IDs.
- Recompiling a new combined score after elimination of coexpression from the network.
- Recompiling the final ppi network.

### 3. Intra tuning and optimization of the pipeline based on  different genesets from Random walk and network edges  
- The top %2, %3, %5 and  %10 edge portions  are extracted from each network. 
- Each network will be subjected to random walk with restart uisng the seeds related to the MIEs of EDCs and decoys (benchmark set). 
- The 200,500,700 and 1000 top most visited genes will be extracted after the random walk.
- A binary vector will be genrated for each compound (EDCs and decoys) with 1 representing the gene in the list of top visited gene.
- The pairwise jaccard distance beween the binary vector of each EDCs and decoy will be calculated. 
- Using the jaccard dismilarity matrix and a vector representing the class of each componud as EDC or Decoy average silhouette score was calculated for EDCs. 

#### Script 3_1_Drug_matrix_tuner.R                                 
- Intra tuning and optimization of 4 Drug Matrix networks

#### Script 3_2_TG_GATEs_tuner.R                                
- Intra tuning and optimization of 8 TG-GATEs networks

#### Script 3_3_Consensus_tuner.R                                   
- Intra tuning and optimization of 2 consensus networks related to LINCS and DrugMatrix-TG-Gates

#### Script 3_4_PPI_tuner.R                                          
- New ppi networks were compiled using the 0.6,0.65,0.7,0.75,0.8,0.85 values as the cutoffs for the combined score. 
- The network will be subjected to random walk with restart uisng the seeds related to the MIEs of EDCs and decoys (benchmark set). 
- The 200,500,700 and 1000 top most visited genes will be extracted after the random walk.
- A binary vector will be genrated for each compound (EDCs and decoys) with 1 representing the gene in the list of top visited gene.
- The pairwise jaccard distance beween the binary vector of each EDCs and decoy will be calculated. 
- Using the jaccard dismilarity matrix and a vector representing the class of each componud as EDC or Decoy average silhouette score was calculated for EDCs. 

#### Script 3_5_pareto_solution_on_tuning_results.R                
- Using pareto solution to obtain final genesets size and edge percents among none dominant solutions. 
- The pareto solution is used to maximize the silhouette score, minimize the gene and edge percent for the networks. 
- In PPI  pareto is being used to maximize the silhouette score, minimize the edge percent and maximize the combined score. 


### 4. Random walk with restart and fgsea (RWR-FGSEA) 

#### Script 4_1_RWR_FGSEA_for_edc_decoys.R                        
- Random walk with restart will be performed for all 15 networks using the optimized edge percent/combined scores from the the MIEs of each EDCs and decoys as seeds 
- Fast gene set enrichment analysis is performed using the retrieved pathways as gene sets.

#### Script 4_2_RWR_FGSEA_for_all_compounds_in_CTD.R               
- Random walk with restart will be performed for all 15 networks using the optimized edge percent/combined scores from the the MIEs of each compounds in CTD as seeds 
- Fast gene set enrichment analysis is performed using the retrieved pathways as gene sets.
- Normalized enrichment scores for all netwworks will be saved.

### 5. GLM modeling of training set (Pathway scores of EDCs and decoys  and labels), accuracy tests and Visualization 

#### Script 5_1_manual_curation_of_pathways_as_features.R          
- The pathways related to viral,bacterial, radiation will be removed. 
- The duplicated pathways based on jaccard similarity will be removed.
- The pathway with no genes expressed in liver are being removed.

#### Script 5_2_Preparation_of_training_datasets.R                   
- Preparation of a list for each network x= the matrix of NES scores from FGSEA and Y = labels as EDC and decoy, n-edc=number of EDCs for each layer and n-decoy= number of decoys for each layer
- The pathways with non significant values will be removed.

#### Script 5_3_glm_modeling.R                                  
- Performing elastic net GLM on training set related to each network using 5 fold cross-validation as tuning method for the parameters.
- Saving the GLM model coefficients for eahc network.

#### Script 5_4_k_fold_cross_validation.R                          
- Repeated 5_fold_cross_validation will be performed on all 15 models. (Pathway level)
- Repeated 5_fold cross_validation will be performed on a binary data matrix of the genes as columns. Compounds of benchmark (EDC,decoy) as rows.
(The genes related to MIEs are characterized as 1 in the binary matrix. (Gene level))

#### Script 5_5_comparing_cross_validation_across_all_layers_ANOVA.R
- The F1 scores of the k-fold-cross validation will be compared using ANOVA.
- Boxplot will be used to represent the obtained F1 scors across all GLM models.

#### Script 5_6_Integration_of_coefficients_stabilties_NES_scores.R
- For each network, ROC analysis will be used as a univariate method to evaluate the NES scores for each pathway for EDCs and decoys.
- The NES scores, glm coefficients ROC-AUCs, average of NES scores will be integrated across all data layers (suppl. data).

#### Script 5_7_NES_bubble_plot_MOA.R                               
- The map of pathway activation scores and GLM coefficient for each pathway and network will be represented as bubble plot.
- The generated plot can be used to indicate the putative pathways as mode of action for the EDCs.


### 6. Prediction of all compounds toxicity class probility and developing EDC score

#### Script 6_1_prediction_all_compounds_class_probilities.R       
- GLM coefs for each network will be used to predict class probability of all compounds in CTD using their NES scores across different pathways.

#### Script 6_2_ROC_TOXCAST_vs_class_probabilities.R  
- ROC analysis will be used between the class probabilty of each network and the binary in vitro experimental hitcall ToxCast assays related to nuclear receptors for each network.
- Selection of most informative networks based on the results of ROC curve analysis.

#### Script 6_3_Developing_Harmonic_and_average_EDC_scores.R         
- Average EDC score is defined as the average of class probabilties across selected networks for each compound.
- Harmonic EDC score is Defining defined as the harmonic sum class probabilties across selected networks for each compound.

#### Script 6_4_ROC_analysis_EDC_scores_vs_TOXCAST_endpoints.R      
- Validating the compiled EDC scores by ROC analysis between  EDC scores and ToxCast hitcall data for each assay endpoint.


### 7. Comparison and evaluation of EDC scores with other Toxicity Scores and tools 

#### Script 7_1_Comparison_of_VAM_on_DeDuCt.R                      
- Representation of the compounds with evidence as EDC in DEDuCT using their EDC scores.

#### Script 7_2_Comparison_VAM_scores_vs_TOXPI_scores.R           
- Evaluation of the EDC scores by comparing them with ToxPi scores (The scores developed from ToxCast assay endpoints) 

#### Script 7_3_Comparison_pathway_scores_vs_TOXDB_pathway_scores.R  
- Comparison between EDC scores with pathway scores from TOXDB for the pathways aryl hydrocarbon receptor, Breaset cancer and estrogen receptor


### 8. Validation of EDC scores with external data 

#### Script 8_1_Validation_VAM_scores_Eurion_External_set_compounds.R
- Evaluation of EDC scores using known EDCs from expert domain and DeDuCT 
- Calculation of accuaracy for the EDC scores
- Representing the validation set compounds as heatmap plot 

#### Script 8_2_validation_with_disease_score.R  
- For each disease (atherosclerosis, diabetes type2 and metabolic syndrome) a disease score was developed as described in section 9.
- The disease scores for the known EDCs are represented for each disease as a heatmap plot.

#### Script 8_3_Validation_with_ToxCast_mies.R  
- The target genes for the compounds with positive assay endpoint in ToxCast wil be considered as the MIEs of the compounds. 
- Random walk with resetart and fast gene set enrichment analysis pipeline will be repeated starting with the MIEs from ToxCast.
- The obtained NES scores will be subjected to GLM models and the resulting class probabilities will be compared with the class probabilities of the MIEs from CTD.

# Sensitivity and enrichment analysis
#### Script 8_4_dose_response_TG_GATEs_single.R 
- TG-Gates class probability profile of the compounds will be categorized as 4 different patterns (based on increase or decrease pattern at different doses Low, Middle and High)  
-  The significant pathways based on increases or decrease of pathway activation scores will be determined using ANOVA followed by post tukey test with bootstrap sampling method.

#### Script 8_5_time_exposure_response_TG_GATEs_repeated.R          
- TG-Gates class probability profile of the compounds will be categorized as 4 different patterns (based on increase or decrease pattern at different Time of exposure 8 days, 15 days and 29 days)  
-  The significant pathways based on increases or decrease of pathway activation scores will be determined using ANOVA followed by post tukey test with bootstrap sampling method.

#### Script 8_6_textmining_16_april_2020_endocrine_disruption.R
- After the search in pubmed for the term endocrine disruption (ED) the genes related to ED will be determined using text mining with pubmed.mineR.

#### Script 8_7_time_dose_plots_heatplot.R                      
- The pathways versus thier genes will be represented as two heatplots (dose and time point).
- The highly cited gened will be determined in the plot.
- The essential genes will be highlighted in the plot.

## 8_8_categorization_enrichment_result.R  
- For each data layer (network) the number selected KEGG pathways by elastic GLM will be categorized and represented in plots 
- For each data layer (network) the number selected REACTOME pathways by elastic GLM will be categorized and represented in plots 

### 9. Linking ED MIEs to AOPs  

#### Script 9_1_preparation_data_set_disease_biomarker.R           
- Retaining the pathways with length of 50 genes for biomarker discovery.
- Preparation of the NES scores for each data layer based on the selected pathways. 

#### Script 9_2_chem2disease_CTD.R                               
- Compiling a binary matrix from chem-disease associations in CTD.
- The rows of the matrix are the compound names and the columns are the disease names.
- The compound disease associations are represented by 1 in the matrix.

#### Script 9_3_preparation_two_class_training_set_by_disease_name.R 
- Preparation of training set (NES scores and class labels) for atherosclerosis, metabolic syndrome, diabetes type 2 and coronary artery disease

#### Script 9_4_1_metabolic_syndrome_two_class_traning_models_glm.R   
- Training GLM  models of metabolic syndrome for all 15 data layers

#### Script 9_4_2_artherosclerosis_two_class_training_models_glm.R  
- Training GLM  models of atherosclerosis for all 15 data layers

#### Script 9_4_3_diabetes_2_two_class_traning_models_glm.R         
- Training GLM  models of diabetes T2 for all 15 data layers

#### Script 9_5_1_kfold_CV_two_class_metabolic_syndrome_glm.R        
- Performing K-fold-Cross validation on models of metabolic syndrome for all 15 data layers

#### Script 9_5_2_kfold_CV_two_class_artherosclerosis_glm.R         
- Performing K-fold-Cross validation on models of atherosclerosis for all 15 data layers

#### Script 9_5_3_kfold_CV_two_class_diabetes_2_glm.R                
- Performing K-fold-Cross validation on models of metabolic syndrome for all 15 data layers

#### Script 9_6_1_comparison_CV_F1_scores_all_layers_ANOVA_metabolic_syndrome.R
- The F1 scores of the k-fold-cross validation will be compared using ANOVA across all 15 data layer for metabolic syndrome.
- Boxplot will be used to represent the obtained F1 scors across all GLM models for metabolic syndrome.

#### Script 9_6_2_comparison_CV_F1_scores_all_layers_ANOVA_atherosclerosis.R
- The F1 scores of the k-fold-cross validation will be compared using ANOVA across all 15 data layer for atherosclerosis.
- Boxplot will be used to represent the obtained F1 scors across all GLM models for atherosclerosis.

#### Script 9_6_3_comparison_CV_F1_scores_all_layers_ANOVA_diabetes.R 
- The F1 scores of the k-fold-cross validation will be compared using ANOVA across all 15 data layer for diabetes type 2.
- Boxplot will be used to represent the obtained F1 scors across all GLM models for dyabetes type 2.

#### Script 9_7_1_Integration_of_glm_coefs_stabilities_NES_ROC_metabolic_syndrome.R
- For each network, ROC analysis will be used as a univariate method to evaluate the NES scores for each pathway for EDCs leading to metabolic synndrome vs edcs not leading to metabolic syndrome.
- The NES scores, glm coefficients ROC-AUCs, average of NES scores will be integrated across all data layers (suppl. data).

#### Script 9_7_2_Integration_of_glm_coefs_stabilities_NES_ROC_atherosclerosis.R
- For each network, ROC analysis will be used as a univariate method to evaluate the NES scores for each pathway for EDCs leading to atherosclerosis vs edcs not leading to atherosclerosis.
- The NES scores, glm coefficients ROC-AUCs, average of NES scores will be integrated across all data layers (suppl. data).

#### Script 9_7_3_Integration_of_glm_coefs_stabilities_NES_ROC_diabetes_2.R
- For each network, ROC analysis will be used as a univariate method to evaluate the NES scores for each pathway for EDCs leading to diabetes type 2 vs edcs not leading to diabetes type 2.
- The NES scores, glm coefficients ROC-AUCs, average of NES scores will be integrated across all data layers (suppl. data).

#### Script 9_8_bubble_plots_pathways.R                           
- The map of pathway activation scores and GLM coefficient for each pathway and network will be represented as bubble plots for the three diseases metabolic syndrome,
atherosclerosis and diabetes type 2.
- The generated plots can be used to indicate the putative pathways linking EDCs to adverse outcome.

#### Script 9_9_disease_scores.R
- For each disease (metabolic syndrome, atherosclerosis and diabetes type 2), harmonic sum and average score will be compiled across all 15 networks for 12k compounds in CTD.

#### Script 9_10_disease_scores_pie_chart.R                         
- For each disease (metabolic syndrome, atherosclerosis and diabetes type 2) a pie chart will be depicted to reveal distribution of disease scores across 12k compounds in CTD.
