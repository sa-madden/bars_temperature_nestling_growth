# The effects of temperature on nestling growth in a songbird depend on developmental and social context

### Author names, affiliations, and contact information 
XXX 

### Paper summary

Climate change can adversely impact animals, especially those that cannot independently thermoregulate or avoid exposure. Cold, hot, and variable temperatures may impede nestling songbird growth due to increased thermoregulatory costs and reduced food delivery by parents. At a broad scale, temperature effects on nestling growth vary across climatic zones, but how temperature effects vary with early-life developmental constraints imposed by the timing of thermoregulatory development, competition with siblings, and the amount of parental care has received less attention. We investigated whether the effect of temperature on the mass of wild barn swallow (Hirundo rustica erythrogaster) nestlings (n = 113) in Boulder County, CO depends on timing of exposure during development, relative size within the brood, or level of parental feeding. We found that the temperature effects differed before versus after putative development of thermoregulatory independence and may be more pronounced among the smallest nestlings within each nest and for nestlings that received less parental feeding. These findings indicate the existence of fine-scale heterogeneity in which the effects of temperature on nestling development are sensitive to metabolic constraints and early-life social environment, providing key insight into the factors that may ameliorate or exacerbate climate impacts on individual birds.   

### File organization and description

This repository contains the scripts, data, and output necessary to reproduce analyses for this paper. An overview of the repository organization and file descriptions are provided below. 

### There are 3 subdirectories in this repository.
* _Data_ - this subdirectory contains the raw data and the processed data (.RData, contained within the _Tidy_ subfolder) for all analyses. The .RData follow the modular scripts such that the processed data saved at the end of one script is loaded in the sequential script, thus allowing the user to enter the analysis pipeline at an stage without having to re-run all previous steps. 
* _Scripts_ - this subdirectory contains the sequentially numbered R scripts for load and cleaning the data and for all downstream analyses. Given that the scripts are modular, they can be run alone. The exception to this is that "00_calculate_parental_care.R" contains functions necessary to run "01_load_clean_parental_care_data.R."  For complete reproducibility run the scripts in sequential order starting with script 00.
* _Output_ - this subdirectory contains summary tables and figures generated as part of the analysis.

#### Scripts description
* _"00_calculate_parental_care.R"_ - Create functions used to summarize parental care data (e.g., total number of visits) from raw Animal Behavior Pro .csv logs.
* _"01_load_clean_parental_care_data.R"_ - Run functions to calculate parental care variables from raw Animal Behavior Pro logs and join in with parental care observation metadata. 
* _"01_load_clean_govee_data.R*_ - Combine data from all thermometer .csv data logs into a single dataframe and tidy it. 
* _"02_tidy_summarize_data.R*_ - Join, tidy, and summarize data in preparation for downstream analyses.
* _"03_parental_care_indices.R"_ - Extract Best Linear Unbiased Predictions for parental feeding rate using a generalized linear mixed model.
* _"04_descriptive_stats_vis.R"_ - Run descriptive statistics and create visualizations for the nestling, parental care, and temperature data.
* _*05_covariate_models.R"_ - Models bivariate associations of potential covariates with nestling mass, parental care, and temperature variables.
* _*06_hypothesis_testing_models.R"_ - Run analyses corresponding the the three questions in the paper and create tables and visualizations of model results.
