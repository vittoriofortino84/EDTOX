# EDTOX
A toxicogenomics data space for system-level understanding and prediction of EDC-induced toxicity

# Description of the R scripts used in the pipeline
# 1. Preparation of the MIES, pathways and training benchmark set 

## 1_1_MIEs_from_CTD.R
- Preparation of a binary data matrix for molecular initiating events (MIEs) from compound-gene interactions in CTD. The interactions subtypes related to metabolism
were grouped as metabolism and the interaction types related to transport are grouped as transport.
- Performing multiple correspondence analysis on the resulting matrix uisng FactoMineR and factoextra. 
- Selection of reaction,binding,activity,expression,metabolic processing as the more distant types of the interaction based on the plot of MCA.
- For the compounds with more than 50 gene interactions the less informative gene interactions will be removed.

## 1_2_Pathways_Download.R 
- Pathways related to KEGG, REACTOME,MSIGDB, GO and WIKI with the size of less than 200 will be retrieved.
- A binary dictionary to link the GO terms with Wiki-AOPs will be generated.
- The classifications tags for the pathways related to KEGG and REACTOME pathways will be downloaded and preprocessed for enrichment analysis.

## 1_3_TOXCAST_nuclear_receptors_coregulators.R 
- The genes related to nuclear receptors and their co-regulators from experts domain and NURSA will be merged.
- The target gene ids from ToxCast will be extracted.
- The ToxCast assay endpoints which their target genes are in the list of nuclear receptor genes will be saved as endpoints related to nuclear receptor.

## 1_4_EDC_Decoy_selection.R 
- The list of EDCs will be retrieved from DEDuCT as CAS ids.
- The ToxCast assay endpoints related to nuclear receptor and co-regulators of EDCs will be extracted. 
- The most significat in vitro assay endpoints for the mechanism of EDCs will be characterized using statistical proportion test.
- EDCs (DEDuCT list) which are incative for all the significant assay endpoints will be removed from the final list of EDcs.
- Pairwise jaccard distance between the MIEs related to remaining EDCs and other compounds in CTD will be calculated.
- The compounds with the maximum Jaccard distance with EDCs will be seleceted as negative controls (decoys).

## 1_5_ToxCast_dictionaries.R
- Preparation of a dictionary for ToxCast target genes and their corresponding endpoints.
- Conversion of ToxCast DSSTox_Identifiers to CAS registry identifiers and preparation of the final Hitcall list for ToxCast.

# 2. Generating gene co-expression networks 


## 2_1_Drug_matrix_wTO.R                                    
- Removing the control samples from the preprocessed and normalized LFC values related to Drug Matrix data source for rat in vitro hepatocytes and rat in vivo.
- Selection of the three exposure time points 1,3 and 5 days for in vivo and 1 day for in vitro and splitting the data as four data frames.
- Selection of the genes expressed in liver and orthology mapping of the probe IDs to entrez gene values
- Compiling 4 gene co-expression networks from the data frames using wTO package with bootstrap resampling method.

## 2_2_TG_Gates_wTO.R  
- Removing control samples from the preprocesses and normalized LFC values related to TG-Gates data source for rat in vitro, human invitro and rat in vivo.
- Selection of three dose levels (high, middle and low) and three time points (8, 15 and 29 days) from the LFC values related to TG-Gates rat in vivo. (6 data frames)
- Selection of 1 day time exposure related to human and rat in vitro LFC values. (two data frames)
- Selection of the genes expressed in the liver from each data frame and orthology mapping of probe IDs to gene entrez IDs.
- Compiling 8 gene co-expression networks from the resulting data frames using wTO package with bootstrap resampling method. 


## 2_3_LINCS_wTO.R
- Normalized and preproessed LFC values from the level 5 of phase 1 and phase 2 LINCS data source will be used.
- Selection of cell line HEPG2 with expousre time of 24 hours from phase1 and phase 2 gene expression data in LINCS.
- Selection of the gene IDS which are expressed in the liver.
- Compiling 2 gene networks from phase 1 and phase 2 using wTO package with bootstrapping resampling method.
- Compiling one consensus gene co-expression netowrk from the overlapping genes of the two networks using wTO package.

## 2_4_Consensus_Rat_in_vitro_wTO.R 
- Compiling one consensus gene co-expression netowrk from the overlapping genes of the two networks related to in vitro rat from drug matrix and TG-GATEs using wTO package.

## 2_5_PPI_wTO.R   
- Retrieving protein protein interaction network from stringDB.
- Mapping nodes to entreg gene IDs.
- Recompiling a new combined score after elimination of coexpression from the network.
- Recompiling the final ppi network.

# 3. Intra tuning and optimization of the pipeline based on  different genesets from Random walk and network edges  
- The top %2, %3, %5 and  %10 edge portions  are extracted from each network. 
- Each network will be subjected to random walk with restart uisng the seeds related to the MIEs of EDCs and decoys (benchmark set). 
- The 200,500,700 and 1000 top most visited genes will be extracted after the random walk.
- A binary vector will be genrated for each compound (EDCs and decoys) with 1 representing the gene in the list of top visited gene.
- The pairwise jaccard distance beween the binary vector of each EDCs and decoy will be calculated. 
- Using the jaccard dismilarity matrix and a vector representing the class of each componud as EDC or Decoy average silhouette score was calculated for EDCs. 

## 3_1_Drug_matrix_tuner.R                                 
- Intra tuning and optimization of 4 Drug Matrix networks

## 3_2_TG_GATEs_tuner.R                                
- Intra tuning and optimization of 8 TG-GATEs networks

## 3_3_Consensus_tuner.R                                   
- Intra tuning and optimization of 2 consensus networks related to LINCS and DrugMatrix-TG-Gates

## 3_4_PPI_tuner.R                                          
- New ppi networks were compiled using the 0.6,0.65,0.7,0.75,0.8,0.85 values as the cutoffs for the combined score. 
- The network will be subjected to random walk with restart uisng the seeds related to the MIEs of EDCs and decoys (benchmark set). 
- The 200,500,700 and 1000 top most visited genes will be extracted after the random walk.
- A binary vector will be genrated for each compound (EDCs and decoys) with 1 representing the gene in the list of top visited gene.
- The pairwise jaccard distance beween the binary vector of each EDCs and decoy will be calculated. 
- Using the jaccard dismilarity matrix and a vector representing the class of each componud as EDC or Decoy average silhouette score was calculated for EDCs. 

## 3_5_pareto_solution_on_tuning_results.R                
- Using pareto solution to obtain final genesets size and edge percents among none dominant solutions. 
- The pareto solution is used to maximize the silhouette score, minimize the gene and edge percent for the networks. 
- In PPI  pareto is being used to maximize the silhouette score, minimize the edge percent and maximize the combined score. 


# 4. Random walk with restart and fgsea (RWR-FGSEA) 
## 4_1_RWR_FGSEA_for_edc_decoys.R                        
- Random walk with restart will be performed for all 15 networks using the optimized edge percent/combined scores from the the MIEs of each EDCs and decoys as seeds 
- Fast gene set enrichment analysis is performed using the retrieved pathways as gene sets.

## 4_2_RWR_FGSEA_for_all_compounds_in_CTD.R               
- Random walk with restart will be performed for all 15 networks using the optimized edge percent/combined scores from the the MIEs of each compounds in CTD as seeds 
- Fast gene set enrichment analysis is performed using the retrieved pathways as gene sets.
- Normalized enrichment scores for all netwworks will be saved.

# 5. GLM modeling of training set (Pathway scores of EDCs and decoys  and labels), accuracy tests and Visualization 
## 5_1_manual_curation_of_pathways_as_features.R          
- The pathways related to viral,bacterial, radiation will be removed. 
- The duplicated pathways based on jaccard similarity will be removed.
- The pathway with no genes expressed in liver are being removed.

## 5_2_Preparation_of_training_datasets.R                   
- Preparation of a list for each network x= the matrix of NES scores from FGSEA and Y = labels as EDC and decoy, n-edc=number of EDCs for each layer and n-decoy= number of decoys for each layer
- The pathways with non significant values will be removed.

## 5_3_glm_modeling.R                                  
- Performing elastic net GLM on training set related to each network using 5 fold cross-validation as tuning method for the parameters.
- Saving the GLM model coefficients for eahc network.

## 5_4_k_fold_cross_validation.R                          
- Repeated 5_fold_cross_validation will be performed on all 15 models. (Pathway level)
- Repeated 5_fold cross_validation will be performed on a binary data matrix of the genes as columns. Compounds of benchmark (EDC,decoy) as rows.
The genes related to MIEs are characterized as 1 in the binary matrix. (Gene level)

## 5_5_comparing_cross_validation_across_all_layers_ANOVA.R
- The F1 scores of the k-fold-cross validation will be compared using ANOVA.
- Boxplot will be used to represent the obtained F1 scors across all GLM models.

## 5_6_Integration_of_coefficients_stabilties_NES_scores.R
- For each network, ROC analysis will be used as a univariate method to evaluate the NES scores for each pathway for EDCs and decoys.
- The NES scores, glm coefficients ROC-AUCs, average of NES scores will be integrated across all data layers (suppl. data).

## 5_7_NES_bubble_plot_MOA.R                               
- The map of pathway activation scores for each pathway and network will be represented as bubble plot.


# 6. Prediction of all compounds toxicity class probility and developing EDC score

## 6_1_prediction_all_compounds_class_probilities.R       
- GLM coefs for each network will be used to predict class probability of all compounds in CTD using their NES scores across different pathways.

## 6_2_ROC_TOXCAST_vs_class_probabilities.R  
- ROC analysis will be used between the class probabilty of each network and the binary in vitro experimental hitcall ToxCast assays related to nuclear receptors for each network.
- Selection of most informative networks based on the results of ROC curve analysis.

## 6_3_Developing_Harmonic_and_average_EDC_scores.R         
- Average EDC score is defined as the average of class probabilties across selected networks for each compound.
- Harmonic EDC score is Defining defined as the harmonic sum class probabilties across selected networks for each compound.

## 6_4_ROC_analysis_EDC_scores_vs_TOXCAST_endpoints.R      
- Validating the compiled EDC scores by ROC analysis between  EDC scores and ToxCast hitcall data for each assay endpoint.


# 7. Comparison and evaluation of EDC scores with other Toxicity Scores and tools 

## 7_1_Comparison_of_VAM_on_DeDuCt.R                      
- Representation of the compounds with evidence as EDC in DEDuCT using their EDC scores.

## 7_2_Comparison_VAM_scores_vs_TOXPI_scores.R           
- Evaluation of the EDC scores by comparing them with ToxPi scores (The scores developed from ToxCast assay endpoints) 

## 7_3_Comparison_pathway_scores_vs_TOXDB_pathway_scores.R  
- Comparison between EDC scores with pathway scores from TOXDB for the pathways aryl hydrocarbon receptor, Breaset cancer and estrogen receptor


# 8. Validation of EDC scores with external data 
## 8_1_Validation_VAM_scores_Eurion_External_set_compounds.R
- Evaluation of EDC scores using known EDCs from expert domain and DeDuCT 
- Calculation of accuaracy for the EDC scores
- Representing the validation set compounds as heatmap plot 

## 8_2_validation_with_disease_score.R  
- For each disease (atherosclerosis, diabetes type2 and metabolic syndrome) a disease score was developed as described in section 9.
- The disease scores for the known EDCs are represented for each disease as a heatmap plot.

## 8_3_Validation_with_ToxCast_mies.R  
- The target genes for the compounds with positive assay endpoint in ToxCast wil be considered as the MIEs of the compounds. 
- Random walk with resetart and fast gene set enrichment analysis pipeline will be repeated starting with the MIEs from ToxCast.
- The obtained NES scores will be subjected to GLM models and the resulting class probabilities will be compared with the class probabilities of the MIEs from CTD.

# Sensitivity and enrichment categorization analysis
## 8_4_dose_response_TG_GATEs_single.R 
- TG-Gates class probability profile of the compounds will be categorized as 4 different patterns (based on increase or decrease pattern at different doses Low, Middle and High)  
-  The significant pathways based on increases or decrease of pathway activation scores will be determined using ANOVA followed by post tukey test with bootstrap sampling method.

## 8_5_time_exposure_response_TG_GATEs_repeated.R          
- TG-Gates class probability profile of the compounds will be categorized as 4 different patterns (based on increase or decrease pattern at different Time of exposure 8 days, 15 days and 29 days)  
-  The significant pathways based on increases or decrease of pathway activation scores will be determined using ANOVA followed by post tukey test with bootstrap sampling method.

## 8_6_textmining_16_april_2020_endocrine_disruption.R
- After the search in pubmed for the term endocrine disruption (ED) the genes related to ED will be determined using text mining with pubmineR.


## 8_7_time_dose_plots_heatplot.R                      
- The pathways versus thier genes will be represented as two heatplots (dose and time point).
- The highly cited gened will be determined in the plot.
- The essential genes will be highlighted in the plot.

## 8_8_categorization_enrichment_result.R  
- For each data layer (network) the number selected KEGG pathways by elastic GLM will be categorized and represented in plots 
- For each data layer (network) the number selected REACTOME pathways by elastic GLM will be categorized and represented in plots 

# 9. Linking ED MIEs to AOPs  
## 9_1_preparation_data_set_disease_biomarker.R           
- Retaining the pathways with length of 50 genes for biomarker discovery.
- Preparation of the NES scores for each data layer based on the selected pathways. 

## 9_2_chem2disease_CTD.R                               
Generation of binary matrix for componds and disease

## 9_3_preparation_two_class_training_set_by_disease_name.R 
Preparation of training set for each disease

## 9_4_1_metabolic_syndrome_two_class_traning_models_glm.R   
GLM  model metabolic syndrome all data layers

## 9_4_2_artherosclerosis_two_class_training_models_glm.R  
GLM  model atherosclerosis all data layers

## 9_4_3_diabetes_2_two_class_traning_models_glm.R         
GLM  model diabetes T2 all data layers

## 9_5_1_kfold_CV_two_class_metabolic_syndrome_glm.R        
K-fold-CV metabolic syndrome all data layers

## 9_5_2_kfold_CV_two_class_artherosclerosis_glm.R         
K-fold-CV atherosclerosis all data layers

## 9_5_3_kfold_CV_two_class_diabetes_2_glm.R                
K-fold-CV diabetes T2 all data layers

## 9_6_1_comparison_CV_F1_scores_all_layers_ANOVA_metabolic_syndrome.R
Boxplot F1-scores metabolic syndrome

## 9_6_2_comparison_CV_F1_scores_all_layers_ANOVA_atherosclerosis.R
Boxplot F1-scores atherosclerosis

## 9_6_3_comparison_CV_F1_scores_all_layers_ANOVA_diabetes.R 
Boxplot F1-scores diabets t2

## 9_7_1_Integration_of_glm_coefs_stabilities_NES_ROC_metabolic_syndrome.R
ROC analysis metabolic syndrome

## 9_7_2_Integration_of_glm_coefs_stabilities_NES_ROC_atherosclerosis.R
ROC analysis atherosclerosis

## 9_7_3_Integration_of_glm_coefs_stabilities_NES_ROC_diabetes_2.R
ROC analysis diaetes T2

## 9_8_bubble_plots_pathways.R                           
Bubble plot for the pathways related to AOPs

## 9_9_disease_scores.R                                   
Class probability and harmonic and average score for disease

## 9_10_disease_scores_pie_chart.R                         
Pie chart for the disease and EDC scores
