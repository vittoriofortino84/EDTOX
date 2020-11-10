# EDTOX
A toxicogenomics data space for system-level understanding and prediction of EDC-induced toxicity

# Description of the R scripts used in the pipeline
# 1. Preparation of the MIES, pathways and training benchmark set 

## 1_1_MIEs_from_CTD.R
- Preparation of a binary data matrix of molecular initiating events (MIEs) from compound-gene interactions in CTD. The interactions subtypes related to metabolism
were grouped as metabolism and the interaction types related to transport are grouped as transport.
- Performing multiple correspondence analysis on the resulting matrix uisng FactoMineR and factoextra. 
- Selection of reaction,binding,activity,expression,metabolic processing as the more distant types of the interaction based on the plot of MCA.
- For the compounds with more than 50 gene interactions the less informative gene interactions will be removed.

## 1_2_Pathways_Download.R 
- Pathways related to KEGG, REACTOME,MSIGDB, GO and WIKI with the size of less than 200 will be retrieved.
- A binary dictionary to link the GO terms with Wiki-AOPs will be generated.
- The classifications tags for the pathways related to KEGG and REACTOME pathways will be downloaded and preprocessed for enrichment analysis.

## 1_3_TOXCAST_nuclear_receptors_coregulators.R 
- The genes related to nuclear receptors and their co-regulators from our experts and NURSA data source will be merged.
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

# 3. Intra tuning and optimization of the pipeline based on  different genesets from Random walk and network edges  
The top edge portion 0.02,0.03,0.05,0.1 are extracted from different networks after sorting the weighted topological overlap of the edges. 
Each network was subjected to random walk with restart starting from the MIEs of EDCs and decoys. The top most visited genes
200,500,700 and 1000 were extracted and the jaccard distance beween the EDCs and decoys is calculated. The average of silhouette score 
for the edcs are calculated. 

## 3_1_Drug_matrix_tuner.R                                 
Optimization of Drug Matrix

## 3_2_TG_GATEs_tuner.R                                  
Optimization of TG-GATEs

## 3_3_Consensus_tuner.R                                   
Optimization of Consensus LINCS hepatocytes

## 3_4_PPI_tuner.R                                          
Optimization of PPI

## 3_5_pareto_solution_on_tuning_results.R                
Using pareto solution to obtain final genesets size and edge percents. The pareto solution is used to maximize the silhouette score, minimize the gene and edge percent for the networks except PPI where pareto is being used to maximize the siilhouette score, minimize the edge percent and maximize the combined score. 


# 4. Random walk with restart and fgsea (RWR-FGSEA) 
## 4_1_RWR_FGSEA_for_edc_decoys.R                        
Random walk and restart with restart and gene set enrichment analysis is performed by the MIEs of EDCs and decoys.

## 4_2_RWR_FGSEA_for_all_compounds_in_CTD.R               
Random walk and restart with restart and gene set enrichment analysis is performed by the MIEs of compounds in CTD.


# 5. GLM modeling of training set (Pathway scores of EDCs and decoys  and labels), accuracy tests and Visualization 
## 5_1_manual_curation_of_pathways_as_features.R          
The pathways related to viral,bacterial, radiation were removed. The duplicated pathways based on jaccard similarit were removed.
The pathway with no genes expressed in liver are being removed.

## 5_2_Preparation_of_training_datasets.R                   
Preparation of training set x= the matrix of NES scores from FGSEA and Y the vector of labels as EDC and decoy for the compounds in each data layer for elastic net generalized linear models for different data layers. 

## 5_3_glm_modeling.R                                  
Performing elastic net GLM on training set

## 5_4_k_fold_cross_validation.R                          
k_fold_cross_validation of the models and the gene level models (using MIEs)

## 5_5_comparing_cross_validation_across_all_layers_ANOVA.R
ANOVA on cross-validation results of pathway and MIE based models 

## 5_6_Integration_of_coefficients_stabilties_NES_scores.R
Integration of NES score, glm coefficients and performing ROC analysis on pathways scores for EDCs and decoys.

## 5_7_NES_bubble_plot_MOA.R                               
Bubble plot of pathways to charecterize the possible mechanism of action (MOA) for EDCs


# 6. Prediction of all compounds toxicity class probility and developing EDC score

## 6_1_prediction_all_compounds_class_probilities.R       
Prediction of class probilities for compounds in CTD 

## 6_2_ROC_TOXCAST_vs_class_probabilities.R                
Selection of networks based on ROC curve analysis of the class probabilities from different networks with TOXCAST hitc matrix.

## 6_3_Developing_Harmonic_and_average_EDC_scores.R         
Defining EDC scores as the harmonic sum of class probilites for most informative layers

## 6_4_ROC_analysis_EDC_scores_vs_TOXCAST_endpoints.R      
ROC analysis between final EDC scores and ToxCast Data and heatmap analysis


# 7. Comparison and evaluation of EDC scores with other Toxicity Scores and tools 

## 7_1_Comparison_of_VAM_on_DeDuCt.R                      
Calculation of EDC scores for compounds in DEDuCT.

## 7_2_Comparison_VAM_scores_vs_TOXPI_scores.R           
Comparison between EDC scores with scores from TOXPI  

## 7_3_Comparison_pathway_scores_vs_TOXDB_pathway_scores.R  
Comparison between EDC scores with pathway scores from TOXDB


# 8. Validation of EDC scores with external data 
## 8_1_Validation_VAM_scores_Eurion_External_set_compounds.R
Validation of EDC scores with known chemicals from expert domain which are considered as EDC together with the compounds from DEDuCT list and the heatmap analysis and accuracy.

## 8_2_validation_with_disease_score.R                       
Validation of EDC scores Using disease labels in CTD with the class probabilty calculated from the compounds and the heatmap plot

## 8_3_Validation_with_ToxCast_mies.R                        
Validation of EDC scores by repeating the pipeline with MIES from ToxCast and comparing the result with the class probabilities calculated from the  MIEs from CTD to infer the overall efficiency of the pipeline in toxicity prediction

# Sensitivity and enrichment categorization analysis
## 8_4_dose_response_TG_GATEs_single.R                
Class probability vs low middle and high TG-GATEs are categorized in 4 patterns based on increase at different doses. ANOVA followed by Tukey post test are performed on the pathway activation scores at different doses for each pattern group and the 

## 8_5_time_exposure_response_TG_GATEs_repeated.R          
Class probability vs 8,15,29 days TG-GATEs

## 8_6_textmining_16_april_2020_endocrine_disruption.R      
Text mining for genes related to endocrine disruption

## 8_7_time_dose_plots_heatplot.R                      
Heatplots for gene biomarkers

## 8_8_categorization_enrichment_result.R                 
Categorization of pathways in KEGG and REACTOME

# 9. Linking ED MIEs to AOPs  
## 9_1_preparation_data_set_disease_biomarker.R           
Preapration of the pathways for AOP modeling

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
