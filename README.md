# EDTOX
A toxicogenomics data space for system-level understanding and prediction of EDC-induced toxicity

# Description of the R files used in the pipeline
# 1. Scripts used for the prepration of the MIES, pathways and training benchmark set 

# 1_1_MIEs_from_CTD.R
Preparation of MIEs from compound-gene interactions in CTD

# 1_2_Pathways_Download.R 
Pathways related to KEGG, REACTOME,MSIGDB, GO and WIKI and annotation dictionary  for the pathways

# 1_3_TOXCAST_nuclear_receptors_coregulators.R 
To get TOXCAST nuclear receptors assay endpoints using NURSA and nuclear receptor genes

# 1_4_EDC_Decoy_selection.R 
Selection of EDCs and decoys based on proportion test of ToxCast assays and jaccard dissimilarity of MIEs

# 1_5_ToxCast_dictionaries.R
Preparation of ToxCast target genes and endpoints dictionaries

# 2. Scripts used for generating gene co-expression networks 

# 2_1_Drug_matrix_wTO.R                                    
wTO network from Drug Matrix LFCs

# 2_2_TG_Gates_wTO.R                                       
wTO network from TG-GATEs LFCs

# 2_3_LINCS_wTO.R
wTO network from LINCS LFCs

# 2_4_Consensus_Rat_in_vitro_wTO.R 
wTO network from CONSENSUS rat hepatocytes from DrugMatrix and TG-GATEs

# 2_5_PPI_wTO.R                                             
ppi Network from StringDB for Homo Sapiens

# 3. Intra Tuning and optimization of the pipeline based on  different genesets from Random walk and network edges  

# 3_1_Drug_matrix_tuner.R                                 
optmization of Drug Matrix

# 3_2_TG_GATEs_tuner.R)                                  
optimization of TG-GATEs

# 3_3_Consensus_tuner.R                                   
optimization of Consensus LINCS hepatocytes

# 3_4_PPI_tuner.R                                          
optimization of PPI

# 3_5_pareto_solution_on_tuning_results.R                
Using pareto solution to obtain final genesets size and edge percents


# 4. Random walk with restart and fgsea (RWR-FGSEA) 
# 4_1_RWR_FGSEA_for_edc_decoys.R                        
RWR-FGSEA on EDC and decoys set

# 4_2_RWR_FGSEA_for_all_compounds_in_CTD.R               
RWR-FGSEA on all compounds in CTD


# 5. GLM modeling of training set (Pathway scores of EDCs and decoys  and labels), accuracy tests and Visualization 
# 5_1_manual_curation_of_pathways_as_features.R          
preprocessing of training set :Curation of specific patheway

# 5_2_Preparation_of_training_datasets.R                   
Preparation of training set for elastic net generalize linear models 

# 5_3_glm_modeling.R                                  
Performing elastic net GLM on training set

# 5_4_k_fold_cross_validation.R                          
k_fold_cross_validation of the models and the gene level models (using MIEs)

# 5_5_comparing_cross_validation_across_all_layers_ANOVA.R
ANOVA on cross-validation results of pathway and MIE based models 

# 5_6_Integration_of_coefficients_stabilties_NES_scores.R
Integration of NES score, glm coefficients and performing ROC analysis on pathways scores 

# 5_7_NES_bubble_plot_MOA.R                               
Bubble plot of pathways to predict MOA for EDCs


# 6. Prediction of all compounds toxicity class probility and developing EDC score

# 6_1_prediction_all_compounds_class_probilities.R       
Prediction of class probilities for compounds in CTD 

# 6_2_ROC_TOXCAST_vs_class_probabilities.R                
Selection of networks based on ROC curve analysis with TOXCAST hitc matrix

# 6_3_Developing_Harmonic_and_average_EDC_scores.R         
Defining EDC scores as the harmonic sum of class probilites for most informative layers

# 6_4_ROC_analysis_EDC_scores_vs_TOXCAST_endpoints.R      
ROC analysis between final EDC scores and ToxCast Data and heatmap analysis


# 7. Comparison and evaluation of EDC scores with other Toxicity Scores and tools -------------------------

source('scripts/7_1_Comparison_of_VAM_on_DeDuCt.R')                        # 7_1. Calculation of EDC scores for compounds in DeDUCT list
source('scripts/7_2_Comparison_VAM_scores_vs_TOXPI_scores.R')              # 7_2. Comparison between EDC scores with scores from TOXPI      
source('scripts/7_3_Comparison_pathway_scores_vs_TOXDB_pathway_scores.R')  # 7_3. Comparison between EDC scores with pathway scores from TOXDB


# 8. Validation of EDC scores with external data ------------------------
source('scripts/8_1_Validation_VAM_scores_Eurion_External_set_compounds.R') # 8_1. Validation of EDC scores with known chemical considered as EDC by: Nick Plant 
source('scripts/8_2_validation_with_disease_score.R')                       # 8_2. Validation of EDC scores Using disease labels in CTD
source('scripts/8_3_Validation_with_ToxCast_mies.R')                        # 8_3. Validation of EDC scores by repeating the pipeline with MIES from ToxCast

# sensitivity and enrichment categorization analysis
source('scripts/8_4_dose_response_TG_GATEs_single.R')                       # 8_4. Class probability vs low middle and high TG-GATEs
source('scripts/8_5_time_exposure_response_TG_GATEs_repeated.R')            # 8_5. Class probability vs 8,15,29 days TG-GATEs
source('scripts/8_6_textmining_16_april_2020_endocrine_disruption.R')       # 8_6. Text mining for genes related to endocrine disruption
source('scripts/8_7_time_dose_plots_heatplot.R')                            # 8_7. Heatplots for gene biomarkers
source('scripts/8_8_categorization_enrichment_result.R')                    # 8_8. Categorization of pathways in KEGG and REACTOME


# 9. Linkin MIEs to AOPs which are related to EDCs --------------------------------------------------  
source('scripts/9_1_preparation_data_set_disease_biomarker.R')              # 9_1. Preapration of the pathways for AOP modeling          
source('scripts/9_2_chem2disease_CTD.R')                                    # 9_2. Generation of binary matrix for componds and disease
source('scripts/9_3_preparation_two_class_training_set_by_disease_name.R')  # 9_3. Preparation of training set for each disease
source('scripts/9_4_1_metabolic_syndrome_two_class_traning_models_glm.R')   # 9_4_1. GLM  model metabolic syndrome all data layers
source('scripts/9_4_2_artherosclerosis_two_class_training_models_glm.R')    # 9_4_2. GLM  model atherosclerosis all data layers
source('scripts/9_4_3_diabetes_2_two_class_traning_models_glm.R')           # 9_4_3. GLM  model diabetes T2 all data layers
source('scripts/9_5_1_kfold_CV_two_class_metabolic_syndrome_glm.R')         # 9_5_1. K-fold-CV metabolic syndrome all data layers
source('scripts/9_5_2_kfold_CV_two_class_artherosclerosis_glm.R')           # 9_5_2. K-fold-CV atherosclerosis all data layers
source('scripts/9_5_3_kfold_CV_two_class_diabetes_2_glm.R')                 # 9_5_3. K-fold-CV diabetes T2 all data layers
source('scripts/9_6_1_comparison_CV_F1_scores_all_layers_ANOVA_metabolic_syndrome.R')#9_6_1. Boxplot F1-scores metabolic syndrome
source('scripts/9_6_2_comparison_CV_F1_scores_all_layers_ANOVA_atherosclerosis.R')#9_6_2.Boxplot F1-scores atherosclerosis
source('scripts/9_6_3_comparison_CV_F1_scores_all_layers_ANOVA_diabetes.R') #9_6_3. Boxplot F1-scores diabets t2
source('scripts/9_7_1_Integration_of_glm_coefs_stabilities_NES_ROC_metabolic_syndrome.R')# 9_7_1. ROC analysis metabolic syndrome
source('scripts/9_7_2_Integration_of_glm_coefs_stabilities_NES_ROC_atherosclerosis.R')# 9_7_2. ROC analysis atherosclerosis
source('scripts/9_7_3_Integration_of_glm_coefs_stabilities_NES_ROC_diabetes_2.R')# 9_7_3. ROC analysis diaetes T2
source('scripts/9_8_bubble_plots_pathways.R')                               # 9_8. Bubble plot for the pathways related to AOPs
source('scripts/9_9_disease_scores.R')                                      # 9_9. Class probability and harmonic and average score for disease
source('scripts/9_10_disease_scores_pie_chart.R')                           # 9_10.Pie chart for the disease and EDC scores



