# The effects of temperature on nestling growth in a songbird depend on developmental and social context

## Author names, affiliations, and contact information 
XXX 

----
## Paper summary

Climate change can adversely impact animals, especially those that cannot independently thermoregulate or avoid exposure. Cold, hot, and variable temperatures may impede nestling songbird growth due to increased thermoregulatory costs and reduced food delivery by parents. At a broad scale, temperature effects on nestling growth vary across climatic zones, but how temperature effects vary with early-life developmental constraints imposed by the timing of thermoregulatory development, competition with siblings, and the amount of parental care has received less attention. We investigated whether the effect of temperature on the mass of wild barn swallow (_Hirundo rustica erythrogaster_) nestlings (n = 113) in Boulder County, CO depends on timing of exposure during development, relative size within the brood, or level of parental feeding. We found that the temperature effects differed before versus after putative development of thermoregulatory independence and may be more pronounced among the smallest nestlings within each nest and for nestlings that received less parental feeding. These findings indicate the existence of fine-scale heterogeneity in which the effects of temperature on nestling development are sensitive to metabolic constraints and early-life social environment, providing key insight into the factors that may ameliorate or exacerbate climate impacts on individual birds.   

----
## File organization and description

This repository contains the scripts, data, and output necessary to reproduce analyses for this paper. An overview of the repository organization and file descriptions are provided below. 

### There are 3 subdirectories in this repository.
* _Data_ - this subdirectory contains the raw data and the processed data (.RData, contained within the _Tidy_ subfolder) for all analyses. The .RData follow the modular scripts such that the processed data saved at the end of one script is loaded in the sequential script, thus allowing the user to enter the analysis pipeline at an stage without having to re-run all previous steps. 
* _Scripts_ - this subdirectory contains the sequentially numbered R scripts for load and cleaning the data and for all downstream analyses. Given that the scripts are modular, they can be run alone. The exception to this is that "00_calculate_parental_care.R" contains functions necessary to run "01_load_clean_parental_care_data.R."  For complete reproducibility run the scripts in sequential order starting with script 00.
* _Output_ - this subdirectory contains summary tables and figures generated as part of the analysis.

#### Raw data description
* _nestling_data.csv -_ This file has 354 rows and 29 columns that contain the raw data measured for each nestling barn swallow. Repeated measures for each individual nestling are stored in rows in long format (i.e., there are multiple rows for each individual). The columns include four data types, including 11 characater variables, 13 double numeric variables, 1 logical variable and 4 time variables.
	-	female_band = col_character(), the mother's ID
	-	male_band = col_character(), the father's ID
	-	nestling_band = col_character(), the nestling ID
	-	hatch_order = col_double(), numbering of hatching order
	-	nest = col_double(), the nest ID
	-	site = col_character(), the name of the breeding site location
	-	brood = col_double(), the brood attempt number (either first or second brood)
	-	hatch_date = col_character(), date when the first nestling hatched
	-	sample_date = col_character(), date when a samples were collected from nestlings
	-	nestling_age = col_double(), estimated nestling age (days) on date when samples were collected
	-	nestling_number = col_double(), total number of nestlings in a nest on the date when samples were collected
	-	extract_time = col_time(format = ""), 24 hr time when first nestling was extracted from a nest for sample collection before the behavioral observations
	-	nos_mites = col_double(), estimated number of mites counted on nestling
	-	mites_tp = col_logical(), an unused variable, number of mites counted on paper put in a nest
	-	rt_wing_length = col_double(), the length in mm of the nestlings right wing
	-	mass_pre_obs = col_double(), the nestlings mass before behavioral observations
	-	mass_post_obs = col_double(), the nestlings mass before behavioral observations
	-	post_obs_extract_time = col_time(format = ""), 24 hr time when first nestling was extracted from a nest for sample collection after the behavioral observations
	-	`3min_glucose` = col_double(), the first/baseline blood glucose measurement (mg/dl)
	-	`3min_glucose_time` = col_time(format = ""), the time (minutes:seconds) between when the first nestling was extracted and the first/baseline blood glucose reading
	-	`15min_glucose` = col_double(), the second/baseline blood glucose measurement (mg/dl)
	-	`15min_glucose_time` = col_time(format = ""), the time (minutes:seconds) between when the first nestling was extracted and the first/baseline blood glucose reading
	-	blood_amount_lysis = col_double(), an estimated volume of blood (uL) collected in lysis buffer
	-	blood_amount_rna = col_double(), an estimated volume of blood (uL) collected in RNAlater buffer
	-	lysis_sample = col_character(), a binary indicating if blood collection in lysis buffer was successful
	-	rna_sample = col_character(), a binary indicating if blood collection in lysis buffer was successful
	-	feathers = col_character(), a binary indicating if feather samples were collected
	-	survive_at_sampling = col_character(), a binary indicating if the individual nestling was alive on the sample collection date
	-	notes = col_character(), unstructured notes about the data

* _parent_care_data.csv -_ This file has 92 rows and 18 columns that contain the raw data measured for each nest. Repeated measures for each nest are stored in rows in long format (i.e., there are multiple rows for each nest; up to three observation trials). The columns include three data types, including 8 character variables, 10 double numeric variables, and 2 time variables.
	-	female_band = col_character(), the mother's ID
	-	male_band = col_character(), the father's ID
	-	nest = col_double(), the nest ID
	-	site = col_character(), the name of the breeding site location
	-	obs_date = col_character(), date when each behavior observation and sample collection occurred
  -	observer = col_character(), the initials of the observer
  -	obs_method = col_character(), whether the behavior observation happened in person or was scored from a video recording
  -	nestling_number = col_double(), total number of nestlings in a nest on the date when samples were collected
  -	nestling_age = col_double(), estimated nestling age (days) on date when samples were collected (based on the date when the first nestling hatched)
  -	blind_camera_distance = col_double(), distance (m) between the observer's blind or the camera from the nest
  -	obs_start_time = col_time(format = ""), 24 hr time when behavioral observation trial began
  -	obs_end_time = col_time(format = ""), 24 hr time when behavioral observation trial ended
  -	trial_temp = col_double(), the temperature (deg.C) recorded from the Govee by the observer at the approximate time of the behavior observation
  -	min_temp = col_double(), the 24 hour min temperature (deg.C) recorded from the Govee by the observer at the approximate time of the behavior observation
  -	max_temp = col_double(), the 24 hour max temperature (deg.C) recorded from the Govee by the observer at the approximate time of the behavior observation
  -	trial_humidity = col_double(), the percent humidity recorded from the Govee by the observer at the approximate time of the behavior observation
  -	wind_speed = col_double(), the wind speed recorded from local weather stations by the observer at the approximate time of the behavior observation
  -	wind_gust = col_double(), the wind speed gusts recorded from local weather stations by the observer at the approximate time of the behavior observation
  -	cloud_cover = col_character(), a categorical description of the amount of cloud cover recorded by the observer at the approximate time of the behavior observation
  -	notes = col_character(), unstructured notes about the data

 * _nestling_parental_care_summary.csv -_ This file has 33 rows and 26 columns that contain summary raw data measured for each nest. Data for each nest are stored in rows. The columns include two data types, including22 character variables and 4 double numeric variables.
	-	site = col_character(), the name of the breeding site location
	-	`Nest #` = col_double(), a numeric ID for each nest
   -	`Govee date` = col_character(), date when Govee thermometers were deployed near each nest
    -	`Govee initials` = col_character(), initials of the person who deployed the Govee
    -	CI = col_character(), date when the clutch was initiated at each nest
    -	`No. eggs` = col_double(), number of eggs laid
    -	`No. hatchling` = col_double(), number of nestlings that hatched
    -	`No. survive day 12` = col_double(), number of nestlings that survived to ~ day 12
    -	EHD = col_character(), estimated hatch date based on clutch initiation
    -	HD = col_character(), actual hatch date when first nestling hatches
    -	`Day 3` = col_character(), date when early-development parental care observation occurred
    -	`Day 3 Observer` = col_character(), initials of the parental care observer for early-development
    -	`Day 8` = col_character(), date when mid-development parental care observation occurred
    -	`Day 8 Observer` = col_character(), initials of the parental care observer for mid-development
    -	`Day 12` = col_character(), date when late-development parental care observation occurred
    -	`Day 12 Observer` = col_character(), initials of the parental care observer for late-development
    -	`All data sheets` = col_character(), a binary for if all nestling and parental care data sheets were accounted for
    -	`Data entered` = col_character(), initials of who entered data
    -	Complete_obs = col_character(), a binary for if parental care data was collected at early, mid and late development, i.e., complete data set for each nest
    -	Notes = col_character(), unstructured notes about the data
    -	nestling_measures = col_character(), a binary for if nestling measures were collected for each nest
    -	ab_pro_3 = col_character(), a binary for if the early-development parental care data were summarized in animal behavior pro app.
    -	ab_pro_8 = col_character(), a binary for if the mid-development parental care data were summarized in animal behavior pro app.
    -	ab_pro_12 = col_character(), a binary for if the late-development parental care data were summarized in animal behavior pro app.
    -	parent_care_data_sheet = col_character(), a binary for if all parental care data sheets were accounted for
    -	govee_data = col_character() a binary for if Govee temperature data were collected for each nest

* _revised_abpro_logs-_ This subfolder contains the raw logs of parental care behavior exported for each parental care observation from the Animal Behavior Pro mobile app.

* _used_nests_15_min-_ This subfolder contains the raw logs of temperature data exported from the thermometer placed near each nest.

#### Scripts description
* _"00_calculate_parental_care.R"_ - Create functions used to summarize parental care data (e.g., total number of visits) from raw Animal Behavior Pro .csv logs.
* _"01_load_clean_parental_care_data.R"_ - Run functions to calculate parental care variables from raw Animal Behavior Pro logs and join in with parental care observation metadata. 
* _"01_load_clean_govee_data.R*_ - Combine data from all thermometer .csv data logs into a single dataframe and tidy it. 
* _"02_tidy_summarize_data.R*_ - Join, tidy, and summarize data in preparation for downstream analyses.
* _"03_parental_care_indices.R"_ - Extract Best Linear Unbiased Predictions for parental feeding rate using a generalized linear mixed model.
* _"04_descriptive_stats_vis.R"_ - Run descriptive statistics and create visualizations for the nestling, parental care, and temperature data.
* _*05_covariate_models.R"_ - Models bivariate associations of potential covariates with nestling mass, parental care, and temperature variables.
* _*06_hypothesis_testing_models.R"_ - Run analyses corresponding the the three questions in the paper and create tables and visualizations of model results.


----
### Software and version information

#### This analysis was performed in R
*	R version 4.4.1 (2024-06-14)
*	Platform: x86_64-w64-mingw32 
*	Running under: Windows 11

#### R package versions
attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] DHARMa_0.4.7      broom_1.0.6       lme4_1.1-35.5     Matrix_1.7-0      ggeffects_1.7.2  
 [6] ggpubr_0.6.0      gridExtra_2.3     viridis_0.6.5     viridisLite_0.4.2 Hmisc_5.2-0      
[11] lubridate_1.9.3   forcats_1.0.0     stringr_1.5.1     dplyr_1.1.4       purrr_1.0.2      
[16] readr_2.1.5       tidyr_1.3.1       tibble_3.2.1      ggplot2_3.5.1     tidyverse_2.0.0  
