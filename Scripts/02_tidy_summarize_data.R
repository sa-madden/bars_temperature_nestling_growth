####### Purpose: data tidying and summary for barn swallow weather/
####### nestling growth project
####### By: Sage Madden
####### Created: 5/11/2022
####### Last modified: 12/15/2022

### 1.1 Global options
## a) clear global environment
rm(list = ls())

## b) prevent R from automatically reading character strings as factors
options(stringsAsFactors = FALSE)


### 1.2 Load relevant packages 
library('ggplot2')
library('tidyverse')
library('lubridate')


### 1.3 Get Version and Session Info
R.Version()
sessionInfo()


## 2.1 import data files
parent <- read.csv("./Data/parental_care_data.csv")
nestling <- read.csv("./Data/nestling_data.csv")
ab_pro <- read.csv("./Output/ab_pro_summary_data.csv")
govee <- read.csv("./Output/govee_used_nests_filt.csv")
noaa <- read.csv("./Output/noaa_data_cleaned.csv")


### 3.1 Tidy nestling data
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

## Format glucose times 
nestling$base_gluc_time <-substr(nestling$base_gluc_time, 1, 5)
nestling$base_gluc_time <-ifelse(nestling$base_gluc_time > 0, 
                                 paste0('00:',  nestling$base_gluc_time),
                                 '')

nestling$stress_gluc_time <-substr(nestling$stress_gluc_time, 1, 5)
nestling$stress_gluc_time <-ifelse(nestling$stress_gluc_time>0, 
                                   paste0('00:', nestling$stress_gluc_time),
                                   '')

## Convert blood times to seconds
nestling$base_gluc_s <- 
  period_to_seconds(hms(nestling$base_gluc_time))

nestling$stress_gluc_s <- 
  period_to_seconds(hms(nestling$stress_gluc_time))

## Format nest ID
nestling$nest <- as.integer(nestling$nest)


nestling$nest_id <- as.factor(paste(nestling$site, 
                                    nestling$nest))

#nestling <- nestling %>%
#  select(-c (nest))

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

#Re-code *nominal* factor (with ordered levels)
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

## Calculate mass to wing length ratio (like human BMI)
nestling <- nestling %>%  
  mutate(mass_wing_index = (mass_pre_obs/rt_wing_length)) 

## Make a binary variable for some or no mites
nestling <- nestling %>% 
  mutate(mite_bin = ifelse(nestling$nos_mites >= 1, 'yes', 'no'))

## Calculate difference second and first glucose (mg/dl)
nestling <- nestling %>% 
  mutate(gluc_diff = (stress_gluc 
                      - base_gluc)) 

## Create a variable to indicate positive or negative differences
# between stress and baseline glucose
nestling$diff_dir = ifelse(nestling$gluc_diff <= 0, 
                           'negative', 'positive')

## Reorder variables
nestling_col <- colnames(nestling)

nestling <- nestling[, c("female_band", "male_band","nestling_band",
                         "hatch_order", "site", "brood", "nest_id",
                         "hatch_date", "sample_date", "sample_date_time",
                         "nestling_age", "nestling_number", 
                         "sample_state", "extract_time", "rt_wing_length", 
                         "mass_pre_obs", "size_order", "avg_size_by_nest",      
                         "size_by_avg", "mass_wing_index", 
                         "base_gluc", "base_gluc_time", "base_gluc_s",
                         "stress_gluc","stress_gluc_time", 
                         "stress_gluc_s", "gluc_diff", "diff_dir", 
                         "blood_amount_lysis", "lysis_sample", 
                         "blood_amount_rna", "rna_sample", "feathers", 
                         "nos_mites", "mite_bin", "mites_tp", 
                         "survive_at_sampling", "post_obs_extract_time", 
                         "mass_post_obs", "notes" )] 


### 3.2 Tidy parental_care_data
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


## f) Create a factor to identify developmental stage at obs
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

#Re-code *nominal* factor (with ordered levels)
# Set levels (odering) of state variable 
parent <- transform(parent, obs_state = factor(obs_state,
                                               levels = c("early", "mid", 
                                                             "late")))
levels(parent$obs_state)


### 3.3 Tidy animal behavior pro data
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


### 3.4 Join and tidy nestling and parental care data
## Left join parent_care to nestling 
nestling_parent_care <- nestling  %>%
  left_join(select(parent_care, -c(female_band, male_band, site, 
                                   nestling_age, nestling_number, 
                                   obs_state)),
            by = c('nest_id' = 'nest_id', 
                   'sample_date' = 'obs_date'), 
            copy = F)



## b) Create tertiles of parental care behaviors within each developmental
# state
nestling_parent_care <- nestling_parent_care %>%
  group_by(sample_state)  %>%
  mutate(tert_tot_visit = ntile(total_visits, 3)) %>%
  mutate(tert_tot_feed = ntile(total_feeding_visits, 3)) %>%
  mutate(tert_tot_time = ntile(total_an_duration, 3)) %>%
  mutate(tert_tot_brood = ntile(total_brooding_duration, 3)) %>%
  ungroup()

parent_care <- parent_care %>%
  group_by(obs_state)  %>%
  mutate(tert_tot_visit = ntile(total_visits, 3)) %>%
  mutate(tert_tot_feed = ntile(total_feeding_visits, 3)) %>%
  mutate(tert_tot_time = ntile(total_an_duration, 3)) %>%
  mutate(tert_tot_brood = ntile(total_brooding_duration, 3)) %>%
  ungroup()


### 3.5 Tidy govee data
## Format data to all lower case
govee$site <- tolower(govee$site)
govee$site <- gsub(" ", "", govee$site)

## Format dates
govee$date <- as.Date(govee$date)
govee$ymd_hms <- ymd_hms(govee$ymd_hms)

## Format nest ID
govee$nest <- as.integer(govee$nest)

govee$nest_id <- as.factor(paste(govee$site, govee$nest))

### 3.5 Summarize govee data across obs periods, nestling stages, etc.
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



### 3.6 Summarize govee data over the behavioral observation period 
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


### 3.7 Create variable encoding number of hours over various 
## threshold temps on the day before obs
# Create thresholds for low, intermed, high, and very high temps 
govee_thresholds <- quantile(govee$temp_c, c(0.25, 0.50, 0.75)) 

parent_care_thresh <- parent_care
parent_care_thresh <- mutate(parent_care, day_bef_obs = obs_date - 1)


parent_care_govee_thresh <- parent_care_thresh %>% 
  select(-c("nest", "site")) %>% 
  left_join(govee, by = c("nest_id" = "nest_id", 
                          "day_bef_obs" = "date"))

parent_care_govee_thresh$low <- NA
parent_care_govee_thresh$med <- NA
parent_care_govee_thresh$high <- NA
parent_care_govee_thresh$very_high <- NA

for(i in 1:length(parent_care_govee_thresh$female_band)){
  if(is.na(parent_care_govee_thresh$temp_c[i]) == TRUE){
    parent_care_govee_thresh$low[i] <- NA
    parent_care_govee_thresh$med[i] <- NA
    parent_care_govee_thresh$high[i] <- NA
    parent_care_govee_thresh$very_high[i] <- NA
  }
  else if(parent_care_govee_thresh$temp_c[i] < 
     as.numeric(govee_thresholds[1])){
    parent_care_govee_thresh$low[i] <- 0.25
  }
  else if(parent_care_govee_thresh$temp_c[i] < 
          as.numeric(govee_thresholds[2])){
    parent_care_govee_thresh$med[i] <- 0.25
  } 
  else if(parent_care_govee_thresh$temp_c[i] < 
          as.numeric(govee_thresholds[3])){
    parent_care_govee_thresh$high[i] <- 0.25
  }
  else if(parent_care_govee_thresh$temp_c[i] >= 
          as.numeric(govee_thresholds[3])){
    parent_care_govee_thresh$very_high[i] <- 0.25
  }
}

govee_thresh_summ <- parent_care_govee_thresh %>% 
  group_by(nest_id, obs_state) %>% 
  summarize(hours_low = sum(low, na.rm = T), 
            hours_med = sum(med, na.rm = T), 
            hours_high = sum(high, na.rm = T),
            hours_very_high = sum(very_high, na.rm = T)) %>% 
  ungroup() %>%
  mutate(hours_in_day = hours_low + hours_med + hours_high + hours_very_high)

### 3.8 Summarize govee data by developmental period
# Create variable to indicate whether a govee time is within
# dates for devel state
govee <- rename(govee, measure_date_time = ymd_hms)


# Variables for first and second measure date to be used in loop
first_meas <- nestling_parent_care %>% filter(nestling_age <= 5) %>%
  select(nest_id, sample_date_time)
colnames(first_meas) <- c("nest_id", "first_meas")
second_meas <- nestling_parent_care %>% filter(nestling_age > 5 & 
                                   nestling_age <= 10) %>%
  select(nest_id, sample_date_time)
colnames(second_meas) <- c("nest_id", "second_meas")
first_meas$duplicate <- duplicated(first_meas)
second_meas$duplicate <- duplicated(second_meas)
first_meas <- first_meas %>% filter(duplicate == FALSE & 
                                     is.na(first_meas) == FALSE) %>%
  select(-duplicate)
second_meas <- second_meas %>% filter(duplicate == FALSE & 
  is.na(second_meas) == FALSE) %>%
  select(-duplicate)

govee2 <- left_join(govee, first_meas, by = c("nest_id" = "nest_id"))
govee2 <- left_join(govee2, second_meas, by = c("nest_id" = "nest_id"))

# For loop creating variable to indicate whether time and date are 
# wihtin each devel stage
for(i in 1:length(govee2$nest_id)){
  if(is.na(govee2$first_meas[i]) == FALSE &
    govee2$measure_date_time[i] <= govee2$first_meas[i]){
    govee2$measure_state[i] <- "early"
  }
  else if(is.na(govee2$first_meas[i]) == FALSE & 
          is.na(govee2$second_meas[i]) == FALSE & 
          govee2$measure_date_time[i] > govee2$first_meas[i] &
          govee2$measure_date_time[i] <= govee2$second_meas[i]){
    govee2$measure_state[i] <- "mid"
  }
  else if(is.na(govee2$second_meas[i]) == FALSE &
          govee2$measure_date_time[i] > govee2$second_meas[i]){
    govee2$measure_state[i] <- "late"
  }
  else{
    govee2$measure_state[i] <- NA
  }
}
unknown <- filter(govee2, is.na(measure_state) == TRUE)
unique(unknown$nest_id)

# Summarize 
govee_devel <- govee2 %>% 
  filter(is.na(measure_state) == FALSE) %>%
  group_by(nest_id, measure_state) %>% 
  summarize(devel_max_temp = max(temp_c), 
            devel_avg_temp = mean(temp_c),
            devel_min_temp = min(temp_c),
            devel_med_temp = median(temp_c),
            devel_sd_temp = sd(temp_c),
            devel_iqr_temp = IQR(temp_c),
            devel_max_humid = max(humid_perc),
            devel_avg_humid = mean(humid_perc),
            devel_min_humid = min(humid_perc),
            devel_sd_humid = sd(humid_perc)) %>%
  ungroup()

## Add hours above threshold for devel stages
devel_thresh <- govee2

devel_thresh$low <- NA
devel_thresh$med <- NA
devel_thresh$high <- NA
devel_thresh$very_high <- NA

for(i in 1:length(devel_thresh$nest_id)){
  if(is.na(devel_thresh$temp_c[i]) == TRUE){
    devel_thresh$low[i] <- NA
    devel_thresh$med[i] <- NA
    devel_thresh$high[i] <- NA
    devel_thresh$very_high[i] <- NA
  }
  else if(devel_thresh$temp_c[i] < 
          as.numeric(govee_thresholds[1])){
    devel_thresh$low[i] <- 0.25
  }
  else if(devel_thresh$temp_c[i] < 
          as.numeric(govee_thresholds[2])){
    devel_thresh$med[i] <- 0.25
  } 
  else if(devel_thresh$temp_c[i] < 
          as.numeric(govee_thresholds[3])){
    devel_thresh$high[i] <- 0.25
  }
  else if(devel_thresh$temp_c[i] >= 
          as.numeric(govee_thresholds[3])){
    devel_thresh$very_high[i] <- 0.25
  }
}

devel_thresh_summ <- devel_thresh %>% 
  group_by(nest_id, measure_state) %>% 
  summarize(devel_hours_low = sum(low, na.rm = T), 
            devel_hours_med = sum(med, na.rm = T), 
            devel_hours_high = sum(high, na.rm = T),
            devel_hours_very_high = sum(very_high, na.rm = T)) %>% 
  ungroup() %>%
  mutate(devel_hours_in_stage = devel_hours_low + devel_hours_med + 
           devel_hours_high + devel_hours_very_high) %>%
  filter(is.na(measure_state) == FALSE)

max(devel_thresh_summ$devel_hours_in_stage)
min(devel_thresh_summ$devel_hours_in_stage)


## Summarizer govee data before and after thermoregulatory indepdendence
## Age of thermoregulatory indpendence estimated at 6 days based
## on Dunn 1979 and avg brood size of 3.68
# Get variables indicating the day at which nestlings were 
# 6 days old on average 
early_samples <- nestling_parent_care %>% filter(nestling_age <= 5) %>%
  select(nest_id, sample_date, nestling_age) %>% 
  group_by(nest_id) %>%
  mutate(nest_age = as.integer(mean(nestling_age))) %>% ungroup()

early_samples$thermo_date <- ymd("1990-01-01")

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


thermo_meas <- early_samples %>% 
  select(nest_id, thermo_date)
colnames(thermo_meas) <- c("nest_id", "thermo_meas")

thermo_meas$duplicate <- duplicated(thermo_meas)
thermo_meas <- thermo_meas %>% filter(duplicate == FALSE & 
                                      is.na(thermo_meas) == FALSE) %>%
  select(-duplicate)

govee2 <- left_join(govee, thermo_meas, by = c("nest_id" = "nest_id"))

# For loop creating variable to indicate whether time and date are 
# wihtin each devel stage
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

# Repeat the process with 5 days as the cutoff
# Get variables indicating the day at which nestlings were 
# 6 days old on average 
early_samples_5 <- nestling_parent_care %>% filter(nestling_age <= 5) %>%
  select(nest_id, sample_date, nestling_age) %>% 
  group_by(nest_id) %>%
  mutate(nest_age = as.integer(mean(nestling_age))) %>% ungroup()

early_samples_5$thermo_date <- ymd("1990-01-01")

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


thermo_meas_5 <- early_samples_5 %>% 
  select(nest_id, thermo_date)
colnames(thermo_meas_5) <- c("nest_id", "thermo_meas_5")

thermo_meas_5$duplicate <- duplicated(thermo_meas_5)
thermo_meas_5 <- thermo_meas_5 %>% filter(duplicate == FALSE & 
                                        is.na(thermo_meas_5) == FALSE) %>%
  select(-duplicate)

govee2 <- left_join(govee2, thermo_meas_5, by = c("nest_id" = "nest_id"))

# For loop creating variable to indicate whether time and date are 
# wihtin each devel stage
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
# 6 days old on average 
early_samples_7 <- nestling_parent_care %>% filter(nestling_age <= 5) %>%
  select(nest_id, sample_date, nestling_age) %>% 
  group_by(nest_id) %>%
  mutate(nest_age = as.integer(mean(nestling_age))) %>% ungroup()

early_samples_7$thermo_date <- ymd("1990-01-01")

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


thermo_meas_7 <- early_samples_7 %>% 
  select(nest_id, thermo_date)
colnames(thermo_meas_7) <- c("nest_id", "thermo_meas_7")

thermo_meas_7$duplicate <- duplicated(thermo_meas_7)
thermo_meas_7 <- thermo_meas_7 %>% filter(duplicate == FALSE & 
                                            is.na(thermo_meas_7) == FALSE) %>%
  select(-duplicate)

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


### 3.9 Summarize govee data over the entire nestling period
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

# Calculate hours over threshold for entire nestling period 
nest_thresh_summ <- devel_thresh %>% 
  filter(is.na(measure_state) == FALSE) %>%
  group_by(nest_id) %>% 
  summarize(nest_hours_low = sum(low, na.rm = T), 
            nest_hours_med = sum(med, na.rm = T), 
            nest_hours_high = sum(high, na.rm = T),
            nest_hours_very_high = sum(very_high, na.rm = T)) %>% 
  ungroup() %>%
  mutate(nest_hours_in_stage = nest_hours_low + nest_hours_med + 
           nest_hours_high + nest_hours_very_high) 
  

### 3.10 Calculate max outdoor temperature for each obs day
noaa <- noaa %>% mutate(date = paste(month, day, year, sep = "/"))
noaa$date <- mdy(noaa$date)

daily_temps <- left_join(govee_daily, noaa, 
                         by = "date")

daily_temps <- daily_temps %>% 
  mutate(outdoor_max_c = ((max_temp_f-32)*5/9), 
         outdoor_min_c = ((min_temp_f-32)*5/9))

### 3.11 Calculate avg mass of brood, wing chord for each developmental stage
## Also sum number of mites found

brood_summ <- nestling_parent_care %>% group_by(nest_id, sample_state) %>% 
  summarize(avg_nestl_mass = mean(mass_pre_obs, na.rm = T),
            avg_wing_chord = mean(rt_wing_length, na.rm = T),
            brood_nos_mites = sum(nos_mites, na.rm = T)) %>%
  ungroup()

brood_summ <- brood_summ %>% mutate(brood_mite_bin = ifelse(brood_summ$brood_nos_mites >= 1, 'yes', 'no'))

## And calculate number surviving to day 12
nestling_parent_care$survived <- NA

for(i in 1:length(nestling_parent_care$nest_id)){
  if(nestling_parent_care$survive_at_sampling[i] == "Y"){
    nestling_parent_care$survived[i] <- 1
  } else{
    nestling_parent_care$survived[i] <- 0
  }
}

survive_summ <- nestling_parent_care %>% 
  filter(sample_state == "late") %>%
  group_by(nest_id) %>%
  summarize(num_survived = sum(survived),
            prop_survived = sum(survived)/n()) %>%
  ungroup %>%
  mutate(
    all_survived = case_when(
      prop_survived == 1 ~ "Y",
      prop_survived != 1 ~ "N",
    )
  )



### 3.12 Calculate change in brood mass and wing chord over stage
brood_summ$avg_nestl_growth <- NA
brood_summ$avg_wing_growth <- NA

for(i in 1:length(brood_summ$nest_id)){
  if(brood_summ$sample_state[i] == "early"){
    brood_summ$avg_nestl_growth[i] <- brood_summ$avg_nestl_mass[i]
    brood_summ$avg_wing_growth[i] <- brood_summ$avg_wing_chord[i]
  } 
  else if(brood_summ$sample_state[i] == "mid" &
          brood_summ$sample_state[i-1] == "early"){
    brood_summ$avg_nestl_growth[i] <- brood_summ$avg_nestl_mass[i] -
      brood_summ$avg_nestl_mass[i-1]
    brood_summ$avg_wing_growth[i] <- brood_summ$avg_wing_chord[i] -
      brood_summ$avg_wing_chord[i-1]
  } 
  else if(brood_summ$sample_state[i] == "late" &
          brood_summ$sample_state[i-1] == "mid"){
    brood_summ$avg_nestl_growth[i] <- brood_summ$avg_nestl_mass[i] -
      brood_summ$avg_nestl_mass[i-1]
    brood_summ$avg_wing_growth[i] <- brood_summ$avg_wing_chord[i] -
      brood_summ$avg_wing_chord[i-1]
  }
  else{
    brood_summ$avg_nestl_growth[i] <- NA
    brood_summ$avg_wing_growth[i] <- NA
  }
}

### 3.13 Add colony size variable
site <- c("grizz", "vanloon", "schaaps", "hayes", "bluecloud", 
          "urbanfarm", "cooks", "hoops")
num_pairs <- c(17, 3, 11, 1, 5, 3, 5, 4)
colony <- data.frame(cbind(site, num_pairs))


### 4.0 Merge all of the variables into one data frame at nest level
prim_merged <- left_join(parent_care, govee_obs, 
                         by = c("nest_id" = "nest_id", 
                                "obs_state" = "obs_state"))
prim_merged <- left_join(prim_merged, govee_thresh_summ, 
                         by = c("nest_id" = "nest_id", 
                                "obs_state" = "obs_state"))
prim_merged <- left_join(prim_merged, govee_devel, 
                         by = c("nest_id" = "nest_id", 
                                "obs_state" = "measure_state"))

prim_merged <- left_join(prim_merged, devel_thresh_summ, 
                         by = c("nest_id" = "nest_id", 
                                "obs_state" = "measure_state"))

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

prim_merged <- left_join(prim_merged, nest_thresh_summ, 
                         by = c("nest_id" = "nest_id"))

prim_merged <- left_join(prim_merged, brood_summ, 
                         by = c("nest_id" = "nest_id", 
                                "obs_state" = "sample_state"))

prim_merged <- left_join(prim_merged, colony, 
                         by = c("site" = "site"))

str(prim_merged)

# Also one with individual rows for each nestling (for glucose analyses)
nestl_merged <- left_join(nestling_parent_care, govee_devel, 
                         by = c("nest_id" = "nest_id", 
                                "sample_state" = "measure_state"))

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

nestl_merged <- left_join(nestl_merged, nest_thresh_summ, 
                         by = c("nest_id" = "nest_id"))

nestl_merged <- left_join(nestl_merged, survive_summ, 
                          by = c("nest_id" = "nest_id"))

nestl_merged <- left_join(nestl_merged, brood_summ, 
                         by = c("nest_id" = "nest_id", 
                                "sample_state" = "sample_state"))

nestl_merged <- left_join(nestl_merged, colony, 
                          by = c("site" = "site"))
str(nestl_merged)
colnames(nestl_merged)

# Remove columns I don't need
# colnames(nestl_merged)
# nestl_merged_sel <- select(nestl_merged, 1:16, 21:28, 34:35, 37:41, 47:50, 75:95)
nestl_merged_sel <- select(nestl_merged, -41)

### 4.1 Reorder columns 
colnames(prim_merged)
prim_merged <- prim_merged[, c("female_band", "male_band","nest", "site", 
                          "num_pairs", "nest_id",
                         "obs_date", "observer", "obs_method",
                         "nestling_age", "nestling_number", 
                         "obs_state", "blind_camera_distance", "obs_start_time", 
                         "obs_end_time", "obs_duration", "total_visits",
                         "female_visits", "male_visits",
                         "total_feeding_visits", "female_feeding_visits", 
                         "male_feeding_visits", "total_sanitizing_visits", 
                         "female_sanitizing_visits", "male_sanitizing_visits", 
                         "total_an_duration",
                         "female_an_duration", "male_an_duration", 
                         "total_brooding_duration", "female_brooding_duration", 
                         "male_brooding_duration", "tert_tot_visit", 
                         "tert_tot_feed", "tert_tot_time", 
                         "tert_tot_brood", "trial_temp", "min_temp",      
                         "max_temp", "trial_humidity", 
                         "wind_speed", "wind_gust", "cloud_cover",
                         "obs_max_temp","obs_avg_temp", 
                         "obs_min_temp", "obs_med_temp", "obs_sd_temp", 
                         "obs_iqr_temp", "obs_max_humid", "obs_avg_humid", 
                         "obs_min_humid", "obs_sd_humid", "hours_low", 
                         "hours_med", "hours_high", 
                         "hours_very_high", "hours_in_day", 
                         "devel_max_temp", 
                         "devel_avg_temp", "devel_min_temp",
                         "devel_med_temp", "devel_sd_temp",
                         "devel_iqr_temp",
                         "devel_max_humid", "devel_avg_humid",
                         "devel_min_humid", "devel_sd_humid", 
                         "devel_hours_low", "devel_hours_med", 
                         "devel_hours_high", 
                         "devel_hours_very_high", "devel_hours_in_stage",
                         "thermo_bef_max_temp", "thermo_bef_avg_temp",
                         "thermo_bef_min_temp", "thermo_bef_med_temp",
                         "thermo_bef_sd_temp", "thermo_bef_iqr_temp",
                         "thermo_bef_max_humid", "thermo_bef_avg_humid",
                         "thermo_bef_min_humid", "thermo_bef_sd_humid",
                         "thermo_aft_max_temp", "thermo_aft_avg_temp",
                         "thermo_aft_min_temp", "thermo_aft_med_temp",
                         "thermo_aft_sd_temp", "thermo_aft_iqr_temp",
                         "thermo_aft_max_humid", "thermo_aft_avg_humid", 
                         "thermo_aft_min_humid", "thermo_aft_sd_humid",
                         "thermo_bef_max_temp_5", "thermo_bef_avg_temp_5",
                         "thermo_bef_min_temp_5", "thermo_bef_med_temp_5",
                         "thermo_bef_sd_temp_5", "thermo_bef_iqr_temp_5",
                         "thermo_bef_max_humid_5", "thermo_bef_avg_humid_5",
                         "thermo_bef_min_humid_5", "thermo_bef_sd_humid_5",
                         "thermo_aft_max_temp_5", "thermo_aft_avg_temp_5",
                         "thermo_aft_min_temp_5", "thermo_aft_med_temp_5",
                         "thermo_aft_sd_temp_5", "thermo_aft_iqr_temp_5",
                         "thermo_aft_max_humid_5", "thermo_aft_avg_humid_5", 
                         "thermo_aft_min_humid_5", "thermo_aft_sd_humid_5",
                         "thermo_bef_max_temp_7", "thermo_bef_avg_temp_7",
                         "thermo_bef_min_temp_7", "thermo_bef_med_temp_7",
                         "thermo_bef_sd_temp_7", "thermo_bef_iqr_temp_7",
                         "thermo_bef_max_humid_7", "thermo_bef_avg_humid_7",
                         "thermo_bef_min_humid_7", "thermo_bef_sd_humid_7",
                         "thermo_aft_max_temp_7", "thermo_aft_avg_temp_7",
                         "thermo_aft_min_temp_7", "thermo_aft_med_temp_7",
                         "thermo_aft_sd_temp_7", "thermo_aft_iqr_temp_7",
                         "thermo_aft_max_humid_7", "thermo_aft_avg_humid_7", 
                         "thermo_aft_min_humid_7", "thermo_aft_sd_humid_7",
                         "nest_max_temp", 
                         "nest_avg_temp", "nest_min_temp",
                         "nest_med_temp","nest_sd_temp","nest_iqr_temp",
                         "nest_max_humid", "nest_avg_humid", 
                         "nest_min_humid", "nest_sd_humid",
                         "nest_hours_low", 
                         "nest_hours_med", "nest_hours_high", 
                         "nest_hours_very_high", "nest_hours_in_stage",
                         "avg_nestl_mass", "avg_wing_chord",
                         "avg_nestl_growth", "avg_wing_growth",
                         "brood_nos_mites", 
                         "brood_mite_bin", "notes_obs")] 

colnames(nestl_merged_sel)
nestl_merged <- nestl_merged[, c("female_band", "male_band", "nestling_band",
                               "nest", "site", "num_pairs", 
                               "nest_id", "brood",
                               "hatch_date", "sample_date", "sample_date_time",
                               "nestling_age", "nestling_number", 
                               "sample_state", "extract_time", 
                               "rt_wing_length", "mass_pre_obs", 
                               "size_order", "avg_size_by_nest", 
                               "size_by_avg", "mass_wing_index", "base_gluc",
                               "base_gluc_time", "base_gluc_s",
                               "stress_gluc", "stress_gluc_time", 
                               "stress_gluc_s", "gluc_diff", 
                               "diff_dir", "blood_amount_lysis", 
                               "lysis_sample", "blood_amount_rna", 
                               "rna_sample", "feathers", "nos_mites", "mite_bin",
                               "mites_tp", 
                               "survive_at_sampling", 
                               "survived", "post_obs_extract_time",
                               "mass_post_obs",
                               "avg_nestl_mass", "avg_wing_chord", 
                               "avg_nestl_growth", "avg_wing_growth",
                               "brood_nos_mites", 
                               "brood_mite_bin", "num_survived",
                               "prop_survived", "all_survived", 
                               "observer", "obs_method",
                               "blind_camera_distance", "obs_start_time",
                               "obs_end_time",
                               "obs_duration", "total_visits",
                               "female_visits", "male_visits",
                               "total_feeding_visits", "female_feeding_visits", 
                               "male_feeding_visits", "total_sanitizing_visits", 
                               "female_sanitizing_visits", "male_sanitizing_visits", 
                               "total_an_duration",
                               "female_an_duration", "male_an_duration", 
                               "total_brooding_duration", "female_brooding_duration", 
                               "male_brooding_duration", "tert_tot_visit", 
                               "tert_tot_feed", "tert_tot_time", 
                               "tert_tot_brood", "trial_temp", "min_temp",      
                               "max_temp", "trial_humidity", 
                               "wind_speed", "wind_gust", "cloud_cover",
                               "devel_max_temp", 
                               "devel_avg_temp", "devel_min_temp",
                               "devel_med_temp", "devel_sd_temp",
                               "devel_iqr_temp",
                               "devel_max_humid", "devel_avg_humid",
                               "devel_min_humid", "devel_sd_humid",
                               "thermo_bef_max_temp", "thermo_bef_avg_temp",
                               "thermo_bef_min_temp", "thermo_bef_med_temp",
                               "thermo_bef_sd_temp", "thermo_bef_iqr_temp",
                               "thermo_bef_max_humid", "thermo_bef_avg_humid",
                               "thermo_bef_min_humid", "thermo_bef_sd_humid",
                               "thermo_aft_max_temp", "thermo_aft_avg_temp",
                               "thermo_aft_min_temp", "thermo_aft_med_temp",
                               "thermo_aft_sd_temp", "thermo_aft_iqr_temp",
                               "thermo_aft_max_humid", "thermo_aft_avg_humid", 
                               "thermo_aft_min_humid", "thermo_aft_sd_humid",
                               "thermo_bef_max_temp_5", "thermo_bef_avg_temp_5",
                               "thermo_bef_min_temp_5", "thermo_bef_med_temp_5",
                               "thermo_bef_sd_temp_5", "thermo_bef_iqr_temp_5",
                               "thermo_bef_max_humid_5", "thermo_bef_avg_humid_5",
                               "thermo_bef_min_humid_5", "thermo_bef_sd_humid_5",
                               "thermo_aft_max_temp_5", "thermo_aft_avg_temp_5",
                               "thermo_aft_min_temp_5", "thermo_aft_med_temp_5",
                               "thermo_aft_sd_temp_5", "thermo_aft_iqr_temp_5",
                               "thermo_aft_max_humid_5", "thermo_aft_avg_humid_5", 
                               "thermo_aft_min_humid_5", "thermo_aft_sd_humid_5",
                               "thermo_bef_max_temp_7", "thermo_bef_avg_temp_7",
                               "thermo_bef_min_temp_7", "thermo_bef_med_temp_7",
                               "thermo_bef_sd_temp_7", "thermo_bef_iqr_temp_7",
                               "thermo_bef_max_humid_7", "thermo_bef_avg_humid_7",
                               "thermo_bef_min_humid_7", "thermo_bef_sd_humid_7",
                               "thermo_aft_max_temp_7", "thermo_aft_avg_temp_7",
                               "thermo_aft_min_temp_7", "thermo_aft_med_temp_7",
                               "thermo_aft_sd_temp_7", "thermo_aft_iqr_temp_7",
                               "thermo_aft_max_humid_7", "thermo_aft_avg_humid_7", 
                               "thermo_aft_min_humid_7", "thermo_aft_sd_humid_7",
                               "nest_max_temp", 
                               "nest_avg_temp", "nest_min_temp", 
                               "nest_med_temp","nest_sd_temp",
                               "nest_iqr_temp",
                               "nest_max_humid", "nest_avg_humid", 
                               "nest_min_humid", "nest_sd_humid", 
                               "nest_hours_low", 
                               "nest_hours_med", "nest_hours_high", 
                               "nest_hours_very_high", "nest_hours_in_stage",
                               "notes", "notes_obs")] 


# Save
#save(file = 'Data/tidy_parent_nestl_weather_data_8-23.RData', 
#     list = c('prim_merged', 'nestl_merged', 'govee_daily',
#              'noaa'))
#write.csv(prim_merged, file = "Data/tidy_parent_nestl_weather_data_8-23.csv")


