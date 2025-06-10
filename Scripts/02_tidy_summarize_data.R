####### Purpose: data tidying and summary for barn swallow weather/
####### nestling growth project
####### By: Sage Madden
####### Created: 5/11/2022
####### Last modified: 6/10/2025

# Code Blocks
# 1: Configure work space
# 2: Load data
# 3: Tidy nestling data
# 4: Tidy + summarize parental care data
# 5: Join nestling + parental care data
# 6: Tidy + summarize govee data
# 7: Add colony size
# 8: Merge datasets
# 9: Additional data tidying
# 10: Save data

###############################################################################
##############                Configure work space              ##############
###############################################################################

### Global options
# clear global environment
rm(list = ls())

# prevent R from automatically reading character strings as factors
options(stringsAsFactors = FALSE)


### Load relevant packages 
library('ggplot2')
library('tidyverse')
library('lubridate')


### Get Version and Session Info
R.Version()
sessionInfo()

###############################################################################
##############                      Load data                   ##############
###############################################################################

## import data files
parent <- read.csv("./Data/parental_care_data.csv")
nestling <- read.csv("./Data/nestling_data.csv")
ab_pro <- read.csv("./Data/Tidy/ab_pro_summary_data.csv")
govee <- read.csv("./Data/Tidy/govee_used_nests_filt.csv")



###############################################################################
##############                Tidy nestling data                ##############
###############################################################################

### Tidy nestling data
## Format site names to lower case and remove white space
nestling$site <- tolower(nestling$site)
nestling$site <- gsub(" ", "", nestling$site)

## Rename variables
nestling <- nestling %>%
  rename(base_gluc = X3min_glucose,
         base_gluc_time = X3min_glucose_time,
         stress_gluc = X15min_glucose,
         stress_gluc_time = X15min_glucose_time)

## Format dates
nestling$hatch_date <- as.Date(nestling$hatch_date, '%m/%d/%y')

nestling$sample_date <- as.Date(nestling$sample_date, '%m/%d/%y')

## Concatenate sample_date and extract_time and format for lubridate
nestling$sample_date_time <- paste(nestling$sample_date, 
                                   nestling$extract_time)

nestling$sample_date_time <- gsub(' ','-', nestling$sample_date_time)

nestling$sample_date_time <- gsub(':','-', nestling$sample_date_time)

## Convert sample_date_time to date time class
nestling$sample_date_time <- ymd_hm(nestling$sample_date_time)

## Format nest ID
nestling$nest <- as.integer(nestling$nest)


nestling$nest_id <- as.factor(paste(nestling$site, 
                                    nestling$nest))


## Format site
nestling$site <- as.factor(nestling$site)

## Create a factor to identify developmental stage at sampling
nestling <- nestling %>%
  group_by(nest_id) %>%
  mutate(sample_state =
           case_when(nestling_age <= 6
                     ~ c('early'),
                     nestling_age > 6 & nestling_age <= 10
                     ~ c('mid'),
                     nestling_age > 10 
                     ~ c('late'))) %>%
  ungroup()

# Re-code *nominal* factor (with ordered levels)
# Set levels (odering) of state variable 
nestling <- transform(nestling, 
                      sample_state = factor(sample_state,
                                            levels = c("early", "mid", 
                                                       "late")))
levels(nestling$sample_state)

## Classify nestling size smallest vs other by nest by develop. state
nestling <- nestling  %>%
  group_by(nest_id, sample_state) %>%
  mutate(size_order = as.factor(case_when(rt_wing_length == 
                                            min(rt_wing_length, 
                                                na.rm = T)
                                          ~ 'min',
                                          rt_wing_length != 
                                            min(rt_wing_length, 
                                                na.rm = T)
                                          ~ 'other'))) %>%
  ungroup()

## Classify nestling size above and below average size by nest by
# develop. state
avg_size_nest <- nestling  %>%
  group_by(nest_id, sample_state) %>%
  summarise(avg_size_by_nest = round(mean(rt_wing_length, 
                                          na.rm = T),2)) %>%
  ungroup()

nestling <- nestling  %>%
  left_join(avg_size_nest, by = c('nest_id' = 'nest_id', 
                                  'sample_state' = 'sample_state'), 
            copy = F)

nestling <- nestling  %>%
  group_by(nest_id, sample_state) %>%
  mutate(size_by_avg = as.factor(case_when(rt_wing_length < 
                                             avg_size_by_nest
                                           ~ 'sm',
                                           rt_wing_length >= 
                                             avg_size_by_nest
                                           ~ 'lrg'))) %>%
  ungroup()


## Make a binary variable for some or no mites
nestling <- nestling %>% 
  mutate(mite_bin = ifelse(nestling$nos_mites >= 1, 'yes', 'no'))

# Remove hayes 7 and schaaps 131 from both data sets
#* Hayes 7 no mid and late nestling measurements (depredation)
#* Schaaps 131.2 second brood

nestling <- nestling  %>%
  filter(nest_id != 'hayes 7') %>%
  filter(nest_id != 'schaaps 131') 

## Reorder variables
nestling_col <- colnames(nestling)

nestling <- nestling[, c("female_band", "male_band","nestling_band",
                         "hatch_order", "site", "brood", "nest_id",
                         "hatch_date", "sample_date", "sample_date_time",
                         "nestling_age", "nestling_number", 
                         "sample_state", "extract_time", "rt_wing_length", 
                         "mass_pre_obs", "size_order", "avg_size_by_nest",      
                         "size_by_avg",
                         "base_gluc", "base_gluc_time", 
                         "stress_gluc","stress_gluc_time",  
                         "blood_amount_lysis", "lysis_sample", 
                         "blood_amount_rna", "rna_sample", "feathers", 
                         "nos_mites", "mite_bin", "mites_tp", 
                         "survive_at_sampling", "post_obs_extract_time", 
                         "mass_post_obs", "notes" )] 


###############################################################################
##############      Tidy + summarize parental care data         ##############
###############################################################################

### Tidy parental_care_data
## Format data to all lower case
parent$site <- tolower(parent$site)
parent$site <- gsub(" ", "", parent$site)

## Format dates
parent$obs_date <- as.Date(parent$obs_date, '%m/%d/%y')

## Format nest ID
parent$nest <- as.integer(parent$nest)

parent$nest_id <- as.factor(paste(parent$site, parent$nest))

## Rename notes column
parent <- parent %>% rename(notes_obs = notes)

## Create a factor to identify developmental stage at obs
parent <- parent %>%
  group_by(nest_id) %>%
  mutate(obs_state =
           case_when(nestling_age <= 6
                     ~ c('early'),
                     nestling_age > 6 & nestling_age <= 10
                     ~ c('mid'),
                     nestling_age > 10 
                     ~ c('late'))) %>%
  ungroup()

# Re-code *nominal* factor (with ordered levels)
# Set levels (odering) of state variable 
parent <- transform(parent, obs_state = factor(obs_state,
                                               levels = c("early", "mid", 
                                                             "late")))
levels(parent$obs_state)


### Tidy animal behavior pro data
## Format site names to lower case and remove white space
ab_pro$site <- tolower(ab_pro$site)
ab_pro$site <- gsub(" ", "", ab_pro$site)

## Format dates
ab_pro$obs_date <- as.Date(ab_pro$obs_date, 
                                    '%m/%d/%y')

## d) Format nest ID
ab_pro$nest <- as.integer(ab_pro$nest)

ab_pro$nest_id <- as.factor(paste(ab_pro$site, 
                                  ab_pro$nest))

## Left join parent_care_obs to parent_care_trial 
parent_care <- parent  %>%
  left_join(select(ab_pro, -c(X, female_band, nest, site, 
                                       nestling_age)),
            by = c('nest_id' = 'nest_id', 
                   'obs_date' = 'obs_date'), 
            copy = F)

# Fix 2012s in datasheet
parent_care$obs_date[parent_care$obs_date == ymd("2012-06-22")] <- ymd(	
  "2021-06-22")
parent_care$obs_date[parent_care$obs_date == ymd("2012-07-02")] <- ymd(	
  "2021-07-02")

# Remove hayes 7 and schaaps 131 from both data sets
#* Hayes 7 no mid and late nestling measurements (depredation)
#* Schaaps 131.2 second brood
parent_care <- parent_care  %>%
  filter(nest_id != 'hayes 7') %>%
  filter(nest_id != 'schaaps 131') 

###############################################################################
##############      Join nestling + parental care data         ##############
###############################################################################

### Join and tidy nestling and parental care data
## Left join parent_care to nestling 
nestling_parent_care <- nestling  %>%
  left_join(select(parent_care, -c(female_band, male_band, site, 
                                   nestling_age, nestling_number, 
                                   obs_state)),
            by = c('nest_id' = 'nest_id', 
                   'sample_date' = 'obs_date'), 
            copy = F)


###############################################################################
##############           Tidy + summarize govee data             ##############
###############################################################################

### Tidy govee data
## Format data to all lower case
govee$site <- tolower(govee$site)
govee$site <- gsub(" ", "", govee$site)

## Format dates
govee$date <- as.Date(govee$date)
govee$ymd_hms <- ymd_hms(govee$ymd_hms)

## Format nest ID
govee$nest <- as.integer(govee$nest)

govee$nest_id <- as.factor(paste(govee$site, govee$nest))

### Summarize govee data across obs periods, nestling stages, etc.
## Summarize govee data by day
govee_daily <- govee %>% group_by(nest_id, date) %>% 
  summarize(daily_max_temp = max(temp_c), 
            daily_avg_temp = mean(temp_c),
            daily_min_temp = min(temp_c),
            daily_med_temp = median(temp_c),
            daily_iqr_temp = IQR(temp_c),
            daily_max_humid = max(humid_perc),
            daily_avg_humid = mean(humid_perc), 
            daily_min_humid = min(humid_perc)) %>%
  ungroup()



### Summarize govee data over the behavioral observation period 
parent_care_govee <- parent_care %>% 
  select(-c("nest", "site")) %>% 
  left_join(govee, by = c("nest_id" = "nest_id", 
                              "obs_date" = "date"))
parent_care_govee$obs_start_time <- hm(parent_care_govee$obs_start_time)
parent_care_govee$obs_end_time <- hm(parent_care_govee$obs_end_time)

# Create variable to indicate whether a govee time is within 
# the obs time
parent_care_govee$time_within_obs <- NA
for(i in 1:length(parent_care_govee$female_band)){
  if(is.na(parent_care_govee$obs_start_time[i]) == TRUE | 
     is.na(parent_care_govee$obs_end_time[i]) == TRUE |
     is.na(parent_care_govee$time[i]) == TRUE){
    parent_care_govee$time_within_obs[i] <- NA
  }
  else if(parent_care_govee$time[i] >= parent_care_govee$obs_start_time[i] & 
          parent_care_govee$time[i] <= parent_care_govee$obs_end_time[i]){
    parent_care_govee$time_within_obs[i] <- "Y"
  } 
  else{
    parent_care_govee$time_within_obs[i] <- "N"
  }
}

test <- filter(parent_care_govee, is.na(time_within_obs) == TRUE)

# Summarize govee data by obs period
govee_obs <- parent_care_govee %>% 
  filter(time_within_obs == "Y") %>%
  group_by(nest_id, obs_state) %>% 
  summarize(obs_max_temp = max(temp_c), 
            obs_avg_temp = mean(temp_c),
            obs_min_temp = min(temp_c),
            obs_med_temp = median(temp_c),
            obs_sd_temp = sd(temp_c),
            obs_iqr_temp = IQR(temp_c),
            obs_max_humid = max(humid_perc),
            obs_avg_humid = mean(humid_perc), 
            obs_min_humid = min(humid_perc),
            obs_sd_humid = sd(humid_perc)) %>% 
  ungroup()


## Summarize govee data before and after thermoregulatory indepdendence
## Age of thermoregulatory indpendence estimated at 6 days based
## on Dunn 1979 and avg brood size of 3.68

govee <- rename(govee, measure_date_time = ymd_hms)

# Get variables indicating the day at which nestlings were 
# 6 days old on average 
early_samples <- nestling_parent_care %>% filter(nestling_age <= 5) %>%
  select(nest_id, sample_date, nestling_age) %>% 
  group_by(nest_id) %>%
  mutate(nest_age = as.integer(mean(nestling_age))) %>% ungroup()

# Create column to populate with dates
early_samples$thermo_date <- ymd("1990-01-01")

# Determine date of expected thermoregulatory ind. (day 6) based on age and 
# date of first sample
for(i in 1:length(early_samples$nest_id)){
  if(early_samples$nest_age[i] == 2){
    early_samples$thermo_date[i] <- ymd(early_samples$sample_date[i] 
                                        + 4)
  } else if(early_samples$nest_age[i] == 3){
    early_samples$thermo_date[i] <- ymd(early_samples$sample_date[i] 
                                        + 3)
  } else{
    early_samples$thermo_date[i] <- ymd(early_samples$sample_date[i] 
                                        + 2)
  }
}

# Select on the columns I need -- nest_id and date of thermoreg ind.
thermo_meas <- early_samples %>% 
  select(nest_id, thermo_date)
colnames(thermo_meas) <- c("nest_id", "thermo_meas")

# Remove duplicates due to presence of multiple nestlings per nest
thermo_meas$duplicate <- duplicated(thermo_meas)
thermo_meas <- thermo_meas %>% filter(duplicate == FALSE & 
                                      is.na(thermo_meas) == FALSE) %>%
  select(-duplicate)

# Join date of thermoreg ind. with govee data
govee2 <- left_join(govee, thermo_meas, by = c("nest_id" = "nest_id"))

# For loop creating variable to indicate whether time and date are 
# before or after thermoreg ind.
govee2$thermo_state <- NA
for(i in 1:length(govee2$nest_id)){
  if(is.na(govee2$thermo_meas[i]) == FALSE &
     govee2$measure_date[i] <= govee2$thermo_meas[i]){
    govee2$thermo_state[i] <- "before"
  }
  else if(is.na(govee2$thermo_meas[i]) == FALSE & 
          govee2$measure_date[i] > govee2$thermo_meas[i]){
    govee2$thermo_state[i] <- "after"
  }
  else{
    govee2$thermo_state[i] <- NA
  }
}
unknown <- filter(govee2, is.na(thermo_state) == TRUE)
unique(unknown$nest_id)

# Summarize 
govee_thermo_before <- govee2 %>% 
  filter(thermo_state == "before") %>%
  group_by(nest_id) %>% 
  summarize(thermo_bef_max_temp = max(temp_c), 
            thermo_bef_avg_temp = mean(temp_c),
            thermo_bef_min_temp = min(temp_c),
            thermo_bef_med_temp = median(temp_c),
            thermo_bef_sd_temp = sd(temp_c),
            thermo_bef_iqr_temp = IQR(temp_c),
            thermo_bef_max_humid = max(humid_perc),
            thermo_bef_avg_humid = mean(humid_perc),
            thermo_bef_min_humid = min(humid_perc),
            thermo_bef_sd_humid = sd(humid_perc)) %>%
  ungroup()

govee_thermo_after <- govee2 %>% 
  filter(thermo_state == "after") %>%
  group_by(nest_id) %>% 
  summarize(thermo_aft_max_temp = max(temp_c), 
            thermo_aft_avg_temp = mean(temp_c),
            thermo_aft_min_temp = min(temp_c),
            thermo_aft_med_temp = median(temp_c),
            thermo_aft_sd_temp = sd(temp_c),
            thermo_aft_iqr_temp = IQR(temp_c),
            thermo_aft_max_humid = max(humid_perc),
            thermo_aft_avg_humid = mean(humid_perc),
            thermo_aft_min_humid = min(humid_perc),
            thermo_aft_sd_humid = sd(humid_perc)) %>%
  ungroup()

# Repeat the process with 5 days as the cutoff (for sensitivity analyses)
# Get variables indicating the day at which nestlings were 
# 5 days old on average 
early_samples_5 <- nestling_parent_care %>% filter(nestling_age <= 5) %>%
  select(nest_id, sample_date, nestling_age) %>% 
  group_by(nest_id) %>%
  mutate(nest_age = as.integer(mean(nestling_age))) %>% ungroup()

# Create column to populate with dates
early_samples_5$thermo_date <- ymd("1990-01-01")

# Determine date of expected thermoregulatory ind. (day 6) based on age and 
# date of first sample
for(i in 1:length(early_samples_5$nest_id)){
  if(early_samples_5$nest_age[i] == 2){
    early_samples_5$thermo_date[i] <- ymd(early_samples_5$sample_date[i] 
                                        + 3)
  } else if(early_samples_5$nest_age[i] == 3){
    early_samples_5$thermo_date[i] <- ymd(early_samples_5$sample_date[i] 
                                        + 2)
  } else{
    early_samples_5$thermo_date[i] <- ymd(early_samples_5$sample_date[i] 
                                        + 1)
  }
}

# Select on the columns I need -- nest_id and date of thermoreg ind.
thermo_meas_5 <- early_samples_5 %>% 
  select(nest_id, thermo_date)
colnames(thermo_meas_5) <- c("nest_id", "thermo_meas_5")

# Remove duplicates due to presence of multiple nestlings per nest
thermo_meas_5$duplicate <- duplicated(thermo_meas_5)
thermo_meas_5 <- thermo_meas_5 %>% filter(duplicate == FALSE & 
                                        is.na(thermo_meas_5) == FALSE) %>%
  select(-duplicate)

# Join date of thermoreg ind. with govee data
govee2 <- left_join(govee2, thermo_meas_5, by = c("nest_id" = "nest_id"))

# For loop creating variable to indicate whether time and date are 
# before or after thermoreg ind.
govee2$thermo_state_5 <- NA
for(i in 1:length(govee2$nest_id)){
  if(is.na(govee2$thermo_meas_5[i]) == FALSE &
     govee2$measure_date[i] <= govee2$thermo_meas_5[i]){
    govee2$thermo_state_5[i] <- "before"
  }
  else if(is.na(govee2$thermo_meas_5[i]) == FALSE & 
          govee2$measure_date[i] > govee2$thermo_meas_5[i]){
    govee2$thermo_state_5[i] <- "after"
  }
  else{
    govee2$thermo_state_5[i] <- NA
  }
}
unknown <- filter(govee2, is.na(thermo_state_5) == TRUE)
unique(unknown$nest_id)

# Summarize 
govee_thermo_before_5 <- govee2 %>% 
  filter(thermo_state_5 == "before") %>%
  group_by(nest_id) %>% 
  summarize(thermo_bef_max_temp_5 = max(temp_c), 
            thermo_bef_avg_temp_5 = mean(temp_c),
            thermo_bef_min_temp_5 = min(temp_c),
            thermo_bef_med_temp_5 = median(temp_c),
            thermo_bef_sd_temp_5 = sd(temp_c),
            thermo_bef_iqr_temp_5 = IQR(temp_c),
            thermo_bef_max_humid_5 = max(humid_perc),
            thermo_bef_avg_humid_5 = mean(humid_perc),
            thermo_bef_min_humid_5 = min(humid_perc),
            thermo_bef_sd_humid_5 = sd(humid_perc)) %>%
  ungroup()

govee_thermo_after_5 <- govee2 %>% 
  filter(thermo_state_5 == "after") %>%
  group_by(nest_id) %>% 
  summarize(thermo_aft_max_temp_5 = max(temp_c), 
            thermo_aft_avg_temp_5 = mean(temp_c),
            thermo_aft_min_temp_5 = min(temp_c),
            thermo_aft_med_temp_5 = median(temp_c),
            thermo_aft_sd_temp_5 = sd(temp_c),
            thermo_aft_iqr_temp_5 = IQR(temp_c),
            thermo_aft_max_humid_5 = max(humid_perc),
            thermo_aft_avg_humid_5 = mean(humid_perc),
            thermo_aft_min_humid_5 = min(humid_perc),
            thermo_aft_sd_humid_5 = sd(humid_perc)) %>%
  ungroup()

# Repeat the process with 5 days as the cutoff
# Get variables indicating the day at which nestlings were 
# 7 days old on average 
early_samples_7 <- nestling_parent_care %>% filter(nestling_age <= 5) %>%
  select(nest_id, sample_date, nestling_age) %>% 
  group_by(nest_id) %>%
  mutate(nest_age = as.integer(mean(nestling_age))) %>% ungroup()

# Create column to populate with dates
early_samples_7$thermo_date <- ymd("1990-01-01")

# Determine date of expected thermoregulatory ind. (day 6) based on age and 
# date of first sample
for(i in 1:length(early_samples_7$nest_id)){
  if(early_samples_7$nest_age[i] == 2){
    early_samples_7$thermo_date[i] <- ymd(early_samples_7$sample_date[i] 
                                          + 5)
  } else if(early_samples_7$nest_age[i] == 3){
    early_samples_7$thermo_date[i] <- ymd(early_samples_7$sample_date[i] 
                                          + 4)
  } else{
    early_samples_7$thermo_date[i] <- ymd(early_samples_7$sample_date[i] 
                                          + 3)
  }
}

# Select on the columns I need -- nest_id and date of thermoreg ind.
thermo_meas_7 <- early_samples_7 %>% 
  select(nest_id, thermo_date)
colnames(thermo_meas_7) <- c("nest_id", "thermo_meas_7")

# Remove duplicates due to presence of multiple nestlings per nest
thermo_meas_7$duplicate <- duplicated(thermo_meas_7)
thermo_meas_7 <- thermo_meas_7 %>% filter(duplicate == FALSE & 
                                            is.na(thermo_meas_7) == FALSE) %>%
  select(-duplicate)

# Join date of thermoreg ind. with govee data
govee2 <- left_join(govee2, thermo_meas_7, by = c("nest_id" = "nest_id"))

# For loop creating variable to indicate whether time and date are 
# wihtin each devel stage
govee2$thermo_state_7 <- NA
for(i in 1:length(govee2$nest_id)){
  if(is.na(govee2$thermo_meas_7[i]) == FALSE &
     govee2$measure_date[i] <= govee2$thermo_meas_7[i]){
    govee2$thermo_state_7[i] <- "before"
  }
  else if(is.na(govee2$thermo_meas_7[i]) == FALSE & 
          govee2$measure_date[i] > govee2$thermo_meas_7[i]){
    govee2$thermo_state_7[i] <- "after"
  }
  else{
    govee2$thermo_state_7[i] <- NA
  }
}
unknown <- filter(govee2, is.na(thermo_state_7) == TRUE)
unique(unknown$nest_id)

# Summarize 
govee_thermo_before_7 <- govee2 %>% 
  filter(thermo_state_7 == "before") %>%
  group_by(nest_id) %>% 
  summarize(thermo_bef_max_temp_7 = max(temp_c), 
            thermo_bef_avg_temp_7 = mean(temp_c),
            thermo_bef_min_temp_7 = min(temp_c),
            thermo_bef_med_temp_7 = median(temp_c),
            thermo_bef_sd_temp_7 = sd(temp_c),
            thermo_bef_iqr_temp_7= IQR(temp_c),
            thermo_bef_max_humid_7 = max(humid_perc),
            thermo_bef_avg_humid_7 = mean(humid_perc),
            thermo_bef_min_humid_7 = min(humid_perc),
            thermo_bef_sd_humid_7 = sd(humid_perc)) %>%
  ungroup()

govee_thermo_after_7 <- govee2 %>% 
  filter(thermo_state_7 == "after") %>%
  group_by(nest_id) %>% 
  summarize(thermo_aft_max_temp_7 = max(temp_c), 
            thermo_aft_avg_temp_7 = mean(temp_c),
            thermo_aft_min_temp_7 = min(temp_c),
            thermo_aft_med_temp_7 = median(temp_c),
            thermo_aft_sd_temp_7 = sd(temp_c),
            thermo_aft_iqr_temp_7 = IQR(temp_c),
            thermo_aft_max_humid_7 = max(humid_perc),
            thermo_aft_avg_humid_7 = mean(humid_perc),
            thermo_aft_min_humid_7 = min(humid_perc),
            thermo_aft_sd_humid_7 = sd(humid_perc)) %>%
  ungroup()


###  Summarize govee data over the entire nestling period
govee_nest <- govee %>% 
  group_by(nest_id) %>% 
  summarize(nest_max_temp = max(temp_c), 
            nest_avg_temp = mean(temp_c),
            nest_min_temp = min(temp_c),
            nest_med_temp = median(temp_c),
            nest_sd_temp = sd(temp_c),
            nest_iqr_temp = IQR(temp_c),
            nest_max_humid = max(humid_perc),
            nest_avg_humid = mean(humid_perc),
            nest_min_humid = min(humid_perc), 
            nest_sd_humid = sd(humid_perc)) %>%
  ungroup()

###############################################################################
##############                Add colony size                   ##############
###############################################################################

### Add colony size variable
site <- c("grizz", "vanloon", "schaaps", "hayes", "bluecloud", 
          "urbanfarm", "cooks", "hoops")
num_pairs <- c(17, 3, 11, 1, 5, 3, 5, 4)
colony <- data.frame(cbind(site, num_pairs))


###############################################################################
##############                Merge datasets                     ##############
###############################################################################

### Merge all of the variables into one data frame at nest level

prim_merged <- left_join(parent_care, colony, 
                         by = c("site" = "site"))

prim_merged <- left_join(prim_merged, govee_obs, 
                         by = c("nest_id" = "nest_id", 
                                "obs_state" = "obs_state"))

prim_merged <- left_join(prim_merged, govee_thermo_before, 
                         by = c("nest_id" = "nest_id"))

prim_merged <- left_join(prim_merged, govee_thermo_after, 
                         by = c("nest_id" = "nest_id"))

prim_merged <- left_join(prim_merged, govee_thermo_before_5, 
                         by = c("nest_id" = "nest_id"))

prim_merged <- left_join(prim_merged, govee_thermo_after_5, 
                         by = c("nest_id" = "nest_id"))

prim_merged <- left_join(prim_merged, govee_thermo_before_7, 
                         by = c("nest_id" = "nest_id"))

prim_merged <- left_join(prim_merged, govee_thermo_after_7, 
                         by = c("nest_id" = "nest_id"))

prim_merged <- left_join(prim_merged, govee_nest, 
                         by = c("nest_id" = "nest_id"))


str(prim_merged)

# Also one with individual rows for each nestling (for glucose analyses)
nestl_merged <- left_join(nestling_parent_care, colony, 
                          by = c("site" = "site"))

nestl_merged <- left_join(nestl_merged, govee_thermo_before, 
                          by = c("nest_id" = "nest_id"))

nestl_merged <- left_join(nestl_merged, govee_thermo_after, 
                          by = c("nest_id" = "nest_id"))

nestl_merged <- left_join(nestl_merged, govee_thermo_before_5, 
                          by = c("nest_id" = "nest_id"))

nestl_merged <- left_join(nestl_merged, govee_thermo_after_5, 
                          by = c("nest_id" = "nest_id"))

nestl_merged <- left_join(nestl_merged, govee_thermo_before_7, 
                          by = c("nest_id" = "nest_id"))

nestl_merged <- left_join(nestl_merged, govee_thermo_after_7, 
                          by = c("nest_id" = "nest_id"))

nestl_merged <- left_join(nestl_merged, govee_nest, 
                         by = c("nest_id" = "nest_id"))


str(nestl_merged)
colnames(nestl_merged)

###############################################################################
##############             Additional data tidying              ##############
###############################################################################
## Subset the data to include only the late (~day 12) measurements
## mid and early development nestling measures are not used in this project
late_nestling_parent_care <- nestl_merged %>%
  filter(sample_state == 'late') 

## Extract early development behavioral obs duration for each nest
early_nest_dur <- nestl_merged %>%
  select(nest_id, sample_state, obs_duration) %>%
  filter(sample_state == 'early') %>%
  distinct(nest_id, .keep_all = T) %>%
  pivot_wider(names_from = sample_state,
              values_from = c(obs_duration))

## Extract mid development obs duration for each nest
mid_nest_dur <- nestl_merged %>%
  select(nest_id, sample_state, obs_duration) %>%
  filter(sample_state == 'mid') %>%
  distinct(nest_id, .keep_all = T) %>%
  pivot_wider(names_from = sample_state,
              values_from = c(obs_duration))

## Extract late development obs duration for each nest
late_nest_dur <- nestl_merged %>%
  select(nest_id, sample_state, obs_duration) %>%
  filter(sample_state == 'late') %>%
  distinct(nest_id, .keep_all = T) %>%
  pivot_wider(names_from = sample_state,
              values_from = c(obs_duration))

## Left join early, mid, and late duration
duration <- late_nest_dur %>%
  left_join(mid_nest_dur, by = c('nest_id' = 'nest_id'),
            copy = F)

duration <- duration %>%
  left_join(early_nest_dur, by = c('nest_id' = 'nest_id'),
            copy = F)

## Left join duration to late_nestling_parent_care data frame 
late_nestling_parent_care <- late_nestling_parent_care %>%
  left_join(duration, by = c('nest_id' = 'nest_id'),
            copy = F)

### Add mid development size
mid_nest_relative_size <- nestl_merged %>%
  select(nestling_band, sample_state, size_by_avg) %>%
  filter(sample_state == 'mid') %>%
  filter(!is.na(nestling_band)) %>%
  distinct(nestling_band, .keep_all = T) %>%
  pivot_wider(names_from = sample_state,
              values_from = c(size_by_avg))

## Rename
mid_nest_relative_size <- mid_nest_relative_size %>%
  rename(mid_relative_size = mid)

## Left join mid_nest_relative_size to late_nestling_parent_care
late_nestling_parent_care <- late_nestling_parent_care %>%
  left_join(mid_nest_relative_size, by = c('nestling_band' = 'nestling_band'), 
            copy = F)

## Add mid size order
mid_size_order <- nestl_merged %>%
  select(nestling_band, sample_state, size_order) %>%
  filter(sample_state == 'mid') %>%
  filter(!is.na(nestling_band)) %>%
  distinct(nestling_band, .keep_all = T) %>%
  pivot_wider(names_from = sample_state,
              values_from = c(size_order))

## Rename
mid_size_order <- mid_size_order %>%
  rename(mid_size_order = mid)

## Left join mid_nest_relative_size to late_nestling_parent_care
late_nestling_parent_care <- late_nestling_parent_care %>%
  left_join(mid_size_order, by = c('nestling_band' = 'nestling_band'), 
            copy = F)

### Add mid development mite bin
mid_nest_mite_bin <- nestl_merged %>%
  select(nestling_band, sample_state, mite_bin) %>%
  filter(sample_state == 'mid') %>%
  distinct(nestling_band, .keep_all = T) %>%
  pivot_wider(names_from = sample_state,
              values_from = c(mite_bin))

## Rename
mid_nest_mite_bin <- mid_nest_mite_bin %>%
  rename(mid_mite_bin = mid)

## Left join mite bin to late_nestling_parent_care
late_nestling_parent_care <- late_nestling_parent_care %>%
  left_join(mid_nest_mite_bin, by = c('nestling_band' = 'nestling_band'), 
            copy = F)

### Add mid development brood size
mid_nest_brood_size <- nestl_merged %>%
  select(nest_id, sample_state, nestling_number) %>%
  filter(sample_state == 'mid') %>%
  distinct(nest_id, .keep_all = T) %>%
  pivot_wider(names_from = sample_state,
              values_from = c(nestling_number))

## Rename
mid_nest_brood_size <- mid_nest_brood_size %>%
  rename(mid_brood_size = mid)

## Left join brood size to late_nestling_parent_care
late_nestling_parent_care <- late_nestling_parent_care %>%
  left_join(mid_nest_brood_size, by = c('nest_id' = 'nest_id'), 
            copy = F)


# Create a variable that is the time from initial disturbance (i.e.,
# time of nestling extraction) until observation start time for all devel stages

# Get observation start time and extract time in usable format in nestling data
nestl_merged$obs_start_time_t <- hm(nestl_merged$obs_start_time)

nestl_merged$extract_time_t <- hm(nestl_merged$extract_time)

test <- filter(nestl_merged, is.na(obs_start_time_t) == TRUE)

# Calculate difference between extract time and obs start time
nestl_merged <- nestl_merged %>% 
  mutate(disturb_time = obs_start_time_t - extract_time_t) %>%
  mutate(disturb_min = hour(disturb_time)*60 + minute(disturb_time))

## e) Subset data to assess if care is influenced by nestling removal
disturb_data <- nestl_merged %>%
  group_by(nest_id, sample_date) %>%
  arrange(disturb_min) %>%
  dplyr::slice_head() %>%
  select(nest_id, sample_date, disturb_min) %>%
  ungroup()

## f) Left join the disturb.min to parent_care data 
prim_merged <- prim_merged %>%
  left_join(disturb_data, by = c('nest_id' = 'nest_id', 
                                 'obs_date' = 'sample_date'), 
            copy = F)


# Get observation start time in usable format parental care data
prim_merged$obs_start_time_split <- prim_merged$obs_start_time
prim_merged

prim_merged <- separate(prim_merged, col = obs_start_time_split,
                        into = c("hour", "minute"), sep = ":")
prim_merged <- prim_merged %>% mutate(obs_start_24hr = 
                                        as.numeric(hour) + 
                                        as.numeric(minute)/60)
prim_merged$obs_start_24hr <- as.numeric(prim_merged$obs_start_24hr)

# Make site a factor
prim_merged$fsite <- as.factor(prim_merged$site)

# Make nestID a factor
prim_merged$fnest_id <- as.factor(prim_merged$nest_id)

# Convert hatch date into days since season start 
# Function to convert dates into days since June 1
calculate_days_of_summer <- function(hatch) {
  splitDate <- str_split(hatch, "-")
  splitDate <- unlist(splitDate)
  day <- as.numeric(splitDate[3])
  month <- as.numeric(splitDate[2])
  if(month == 6){
    dos <- day
  }
  else{
    dos <- day + 30
  }
  return(as.character(dos))
}

# Run function to get days since June 1 for each dataset
late_nestling_parent_care$days_summer <- NA
for(i in 1:length(late_nestling_parent_care$hatch_date)){
  late_nestling_parent_care$days_summer[i] <- 
    calculate_days_of_summer(late_nestling_parent_care$hatch_date[i])
}

late_nestling_parent_care$days_summer <- as.numeric(late_nestling_parent_care$days_summer)


###############################################################################
##############                      Save data                    ##############
###############################################################################

# Save
save(file = 'Data/Tidy/tidy_parent_nestl_weather_data_10-4.RData', 
     list = c('prim_merged', 'nestl_merged', 'govee_daily', 'late_nestling_parent_care'))
write.csv(prim_merged, file = "Data/Tidy/tidy_parent_nestl_weather_data_10-4.csv")


