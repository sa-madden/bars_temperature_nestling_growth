####### Purpose: calculate parental care indeces and format final
####### dataset for barn swallow nest mircoclimate and nestling 
####### growth dataset at nestling level
####### By: Sage Madden
####### Created: 12/16/2022
####### Last modified: 12/16/2022


###############################################################################
##############             1.  Configure work space              ##############
###############################################################################

### 1.1 Global options
## clear global environment
rm(list = ls())

## b) prevent R from automatically reading charater strins as factors
options(stringsAsFactors = FALSE)


### 1.2 Install and load CRAN packages   
## Data Manipulation and Descriptive Stats Packages
library('tidyverse')

## Graph Plotting and Visualization Packages
library ('ggplot2')
library('hrbrthemes')
library('viridis')
library ('gridExtra')

## Modelling Packages
library('lme4')


### 1.3 Get Version and Session Info
R.Version()
sessionInfo()


###############################################################################
##############                    2. Load RData                  ##############
###############################################################################  

### 2.1 Load RData
## Load RData tidy barn swallow data
load('Data/tidy_parent_nestl_weather_data_8-23.RData')

### 2.2 Create new dataframe with the specific variables for analyses:
### growth from mid to late development, parental care index, 
### mid development nestling size, mid development mites, mid 
### development brood size, colony size
### temp mid to late development

### Add parental care index
## Extract early development total brooding duration for each nest
early_nest_brood <- nestl_merged %>%
  select(nest_id, sample_state, total_brooding_duration, 
         tert_tot_brood) %>%
  filter(sample_state == 'early') %>%
  distinct(nest_id, .keep_all = T) %>%
  #mutate(row = row_number()) %>% # used to create unique identifier
  pivot_wider(names_from = sample_state,
              values_from = c(total_brooding_duration,
                              tert_tot_brood))
#select(-row)

## Extract mid development total brooding duration for each nest
mid_nest_brood <- nestl_merged %>%
  select(nest_id, sample_state, total_brooding_duration, 
         tert_tot_brood) %>%
  filter(sample_state == 'mid') %>%
  distinct(nest_id, .keep_all = T) %>%
  pivot_wider(names_from = sample_state,
              values_from = c(total_brooding_duration,
                              tert_tot_brood))

## Left join early and mid brooding
brood <- early_nest_brood %>%
  left_join(mid_nest_brood, by = c('nest_id' = 'nest_id'),
            copy = F)

## Extract mid development total feeding visits for each nest
mid_nest_feed <- nestl_merged %>%
  select(nest_id, sample_state, total_feeding_visits, 
         tert_tot_feed) %>%
  filter(sample_state == 'mid') %>%
  distinct(nest_id, .keep_all = T) %>%
  pivot_wider(names_from = sample_state,
              values_from = c(total_feeding_visits,
                              tert_tot_feed))

# f) Extract late development feeding visits for each nest
late_nest_feed <- nestl_merged %>%
  select(nest_id, sample_state, total_feeding_visits, 
         tert_tot_feed) %>%
  filter(sample_state == 'late') %>%
  distinct(nest_id, .keep_all = T) %>%
  pivot_wider(names_from = sample_state,
              values_from = c(total_feeding_visits,
                              tert_tot_feed))

## Left join mid and late feeding
feed <- late_nest_feed %>%
  left_join(mid_nest_feed, by = c('nest_id' = 'nest_id'),
            copy = F)

## Subset the data to include only the late (~day 12) measurements
late_nestling_parent_care <- nestl_merged %>%
  filter(sample_state == 'late') %>%
  select(-c(total_brooding_duration, total_feeding_visits))

## Left join brood to late_nestling_parent_care data frame
late_nestling_parent_care <- late_nestling_parent_care %>%
  left_join(brood, by = c('nest_id' = 'nest_id'),
            copy = F)

## Left join feed to late_nestling_parent_care data frame 
late_nestling_parent_care <- late_nestling_parent_care %>%
  left_join(feed, by = c('nest_id' = 'nest_id'),
            copy = F)

## Extract early development obs duration for each nest
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

## Left join durationto late_nestling_parent_care data frame 
late_nestling_parent_care <- late_nestling_parent_care %>%
  left_join(duration, by = c('nest_id' = 'nest_id'),
            copy = F)

## Create parental care behavior sum adding tertile values for total 
## early brooding duration and total mid and late feeding visits
late_nestling_parent_care <- late_nestling_parent_care %>%
  rowwise() %>%
  mutate(care_sum = sum(tert_tot_brood_early, tert_tot_feed_mid, 
                        tert_tot_feed_late, na.rm = T))  %>%
  mutate(care_sum = na_if(care_sum, 0)) %>%
  ungroup()

## Create parental care behavior possible total to account for 
# missing data
late_nestling_parent_care <- late_nestling_parent_care %>%
  mutate(care_tot = case_when(
    !is.na(tert_tot_brood_early) &
      !is.na(tert_tot_feed_mid) &
      !is.na(tert_tot_feed_late)
    ~ 9,
    is.na(tert_tot_brood_early) &
      !is.na(tert_tot_feed_mid) &
      !is.na(tert_tot_feed_late)
    ~ 6,
    !is.na(tert_tot_brood_early) &
      is.na(tert_tot_feed_mid) &
      !is.na(tert_tot_feed_late)
    ~ 6,
    !is.na(tert_tot_brood_early) &
      !is.na(tert_tot_feed_mid) &
      is.na(tert_tot_feed_late)
    ~ 6,
    is.na(tert_tot_brood_early) &
      is.na(tert_tot_feed_mid) &
      !is.na(tert_tot_feed_late)
    ~ 3,
    is.na(tert_tot_brood_early) &
      !is.na(tert_tot_feed_mid) &
      is.na(tert_tot_feed_late)
    ~ 3,
    !is.na(tert_tot_brood_early) &
      is.na(tert_tot_feed_mid) &
      is.na(tert_tot_feed_late)
    ~ 3)) 

## Create continuous score of parental care behaviors 
# sum of behavior tertile/total possible tertile sum
late_nestling_parent_care <- late_nestling_parent_care %>%
  mutate(care_indx_cont = care_sum/care_tot) %>%
  mutate(care_indx = as.factor(ntile(care_indx_cont, 3))) 

## Re-label the care_indx factor levels
late_nestling_parent_care <- late_nestling_parent_care %>%
  mutate(care_indx =
           case_when(care_indx == 1
                     ~ c('low'),
                     care_indx == 2
                     ~ c('avg'),
                     care_indx == 3
                     ~ c('hi')))

## e) Re-code *nominal* factor (with ordered levels)
# Set levels (odering) of care_indx variable 
late_nestling_parent_care <- transform(late_nestling_parent_care, 
                                       care_indx = factor(care_indx,
                                                          levels = c("low", "avg", 
                                                                     "hi")))

## Re-label the tert_brood_early factor levels
late_nestling_parent_care <- late_nestling_parent_care %>%
  mutate(tert_tot_brood_early =
           case_when(tert_tot_brood_early == 1
                     ~ c('low'),
                     tert_tot_brood_early == 2
                     ~ c('avg'),
                     tert_tot_brood_early == 3
                     ~ c('hi')))

## Re-code *nominal* factor (with ordered levels)
# Set levels (odering) of tert_brood_early variable 
late_nestling_parent_care <- transform(late_nestling_parent_care, 
                                       tert_tot_brood_early = 
                                         factor(tert_tot_brood_early,
                                                levels = c("low", "avg",
                                                           "hi")))

## h) Re-label the tert_tot_feed_mid factor levels
late_nestling_parent_care <- late_nestling_parent_care %>%
  mutate(tert_tot_feed_mid =
           case_when(tert_tot_feed_mid == 1
                     ~ c('low'),
                     tert_tot_feed_mid == 2
                     ~ c('avg'),
                     tert_tot_feed_mid == 3
                     ~ c('hi')))

## Re-code *nominal* factor (with ordered levels)
# Set levels (odering) of tert_tot_feed_mid variable 
late_nestling_parent_care <- transform(late_nestling_parent_care, 
                                       tert_tot_feed_mid = 
                                         factor(tert_tot_feed_mid,
                                                levels = c("low", "avg",
                                                           "hi")))

## j) Re-label the tert_tot_feed_late factor levels
late_nestling_parent_care <- late_nestling_parent_care %>%
  mutate(tert_tot_feed_late =
           case_when(tert_tot_feed_late == 1
                     ~ c('low'),
                     tert_tot_feed_late == 2
                     ~ c('avg'),
                     tert_tot_feed_late == 3
                     ~ c('hi')))

## k) Re-code *nominal* factor (with ordered levels)
# Set levels (odering) of tert_tot_feed_late variable 
late_nestling_parent_care <- transform(late_nestling_parent_care, 
                                       tert_tot_feed_late = 
                                         factor(tert_tot_feed_late,
                                                levels = c("low", "avg",
                                                           "hi")))

## Calculate nestling growth from mid to late development 
diff_wing <- nestl_merged %>%
  select(nestling_band, sample_state, rt_wing_length) %>%
  filter(sample_state != 'early') %>%
  filter(!is.na(rt_wing_length)) %>%
  pivot_wider(id_cols = nestling_band, names_from = sample_state, 
              values_from = rt_wing_length) %>% 
  mutate(rt_wing_diff = (late - mid))

## Rename variables
diff_wing <- diff_wing %>%
  rename(mid_rt_wing_length = mid,
         late_rt_wing_length = late)

## Left join diff_size to late_nestling_parent_care
late_nestling_parent_care <- late_nestling_parent_care %>%
  left_join(diff_wing, by = c('nestling_band' = 'nestling_band'), 
            copy = F)

## Calculate nestling growth from mid to late development 
diff_mass <- nestl_merged %>%
  select(nestling_band, sample_state, mass_pre_obs) %>%
  filter(sample_state != 'early') %>%
  filter(!is.na(mass_pre_obs)) %>%
  pivot_wider(id_cols = nestling_band, names_from = sample_state, 
              values_from = mass_pre_obs) %>% 
  mutate(mass_diff = (late - mid))

## Rename variables
diff_mass <- diff_mass %>%
  rename(mid_mass = mid,
         late_mass = late)

## Left join diff_size to late_nestling_parent_care
late_nestling_parent_care <- late_nestling_parent_care %>%
  left_join(diff_mass, by = c('nestling_band' = 'nestling_band'), 
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

### Add parental care BLUPs
## First, calculate them from parental care dataset (prim_merged)
### 5.3 Close proximity BLUP extractions

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

# Get observation start time and extract time in usable format in nestling data
late_nestling_parent_care$obs_start_time_t <- hm(late_nestling_parent_care$obs_start_time)

late_nestling_parent_care$extract_time_t <- hm(late_nestling_parent_care$extract_time)


test <- filter(late_nestling_parent_care, is.na(obs_start_time_t) == TRUE)


## Create a variable that is the time from initial disturbance (i.e.,
# time of nestling extraction) until observation start time
late_nestling_parent_care <- late_nestling_parent_care %>% 
  mutate(disturb_time = obs_start_time_t - extract_time_t) %>%
  mutate(disturb_min = hour(disturb_time)*60 + minute(disturb_time))

## e) Subset data to assess if care is influenced by nestling removal
disturb_data <- late_nestling_parent_care %>%
  group_by(nest_id, obs_date) %>%
  arrange(disturb_min) %>%
  dplyr::slice_head() %>%
  select(nest_id, sample_date, disturb_min) %>%
  ungroup()

## f) Left join the disturb.min to parent_care data 
prim_merged <- prim_merged %>%
  left_join(disturb_data, by = c('nest_id' = 'nest_id', 
                                 'obs_date' = 'sample_date'), 
            copy = F)

#***NOTE: Remove Hayes 7 and Schaaps 131.2 because not used in main analysis
#* Hayes 7 no mid and late nestling measurements (depredation)
#* Schaaps 131.2 second brood
prim_merged <- prim_merged  %>%
  filter(nest_id != 'hayes 7') %>%
  filter(nest_id != 'schaaps 131') 


### Feeding BLUPs


# NOTE: Use when there are repeated measuresments for a variable
# that is to be used as an explanatory variable in another analysis.
# Can control for other variables that bias estimates of explanatory
# variable

# NOTE: BLUPs are conditional modes from a generalized linear model
# (according to Doug Bates). 
# BLUP = fixef(intrcpt) + ranef
# Model  using glmmTMB (since these were zero inflated behavior data)
feeding_blups_lmm <- lmer(total_feeding_visits ~ scale(nestling_age) + 
                             scale(nestling_number) + scale(obs_med_temp) +
                            scale(disturb_min) + 
                             offset(obs_duration/3600) +
                             (1|fnest_id),
                           data = subset(prim_merged,
                                         is.na(total_feeding_visits) == F &
                                           is.na(nestling_age) == F &
                                           is.na(nestling_number) == F &
                                           is.na(obs_med_temp) == F  &
                                           is.na(disturb_min) == F)
                           )

plot(feeding_blups_lmm)

test <- subset(prim_merged,
               is.na(total_feeding_visits) == F &
                 is.na(nestling_age) == F &
                 is.na(nestling_number) == F &
                 is.na(obs_med_temp) == F  &
                 is.na(disturb_min) == F)
# Residuals look a bit cone shaped

# Try poisson
feeding_blups_glmm_poss <- glmer(total_feeding_visits ~ scale(nestling_age) + 
                            scale(nestling_number) + scale(obs_med_temp) + 
                            offset(log(obs_duration/3600)) +
                            (1|fnest_id),
                          data = subset(prim_merged,
                                        is.na(total_feeding_visits) == F &
                                          is.na(nestling_age) == F &
                                          is.na(nestling_number) == F &
                                          is.na(obs_med_temp) == F),
                          family = poisson()
                          )

simulation_output <- simulateResiduals(fittedModel = feeding_blups_glmm_poss, 
                                       plot = T)
# Not terrible but some issues

# Negative binomial
feeding_blups_glmm_nb <- glmer.nb(total_feeding_visits ~ scale(nestling_age) + 
                                   scale(nestling_number) + scale(obs_med_temp) + 
                                   offset(log(obs_duration/3600)) +
                                    (1|fnest_id),
                                 data = subset(prim_merged,
                                               is.na(total_feeding_visits) == F &
                                                 is.na(nestling_age) == F &
                                                 is.na(nestling_number) == F &
                                                 is.na(obs_med_temp) == F))

simulation_output <- simulateResiduals(fittedModel = feeding_blups_glmm_nb, 
                                       plot = T)
# Looks pretty good
summary(feeding_blups_glmm_nb)


# Calucate the BLUPs, individual variation in feeding visits
# from the best fitting model (above)

## a) Generate best fit model summary
ranef(feeding_blups_glmm_nb) # random effect
fixef(feeding_blups_glmm_nb) # fixed effect
coef(feeding_blups_glmm_nb) # fixed effect

## b) extract BLUPs from mixed model object
feeding_blups <- as.data.frame(ranef(feeding_blups_glmm_nb)) # extract ranef as
# a dataframe, BLUPs = rand effects + intercept (from poiss/neg. binom)

## c) Rename variables in blups table
feeding_blups <- feeding_blups %>%
  rename('fnest_id' = 'grp') %>%
  rename('feeding_ranef' = 'condval') %>%
  rename('feeding_ranef_sd' = 'condsd') %>%
  select(c('fnest_id', 'feeding_ranef', 'feeding_ranef_sd'))


## d) extract fixed effect (intercept) from poisson/neg. binomial model
feeding_intrcpt <- (fixef(feeding_blups_glmm_nb)[[1]])[[1]] # fixed effect
#* INTERPETATION ****#  
# Interpet as expected log count of behavior when controlling for /
# holding constant the effects of covariates (or if exponentiated
# the intercept is the incident rate of the expected count
# of behavior) as a proporiton of time overlap (the offset)

## f) Create a new variable that is ranef plus both poisson/neg. binomial
# model intercept. This provides estimates of individual level
# variation in proportion of time spent of obs vs feeding visits
feeding_blups <-  feeding_blups  %>%
  mutate(feeding_blups = feeding_intrcpt + feeding_ranef) %>%
  mutate(feeding_expontd_blups = exp(feeding_blups))
#* INTERPETATION ****#  
# Interpret as the conditional mode or the log counts / incident rate
# of expected counts as a proportion of the overlap time (offset) for
# each individual...while holding constant effect of other covariates


### Brooding BLUPs
# Model  using glmmTMB (since these were zero inflated behavior data)
brooding_blups_lmm <- lmer((total_brooding_duration/60) ~ scale(nestling_age) + 
                            scale(nestling_number) + scale(obs_med_temp) + 
                            offset(obs_duration/3600) +
                            (1|fnest_id),
                          data = subset(prim_merged,
                                        is.na(total_brooding_duration) == F &
                                          is.na(nestling_age) == F &
                                          is.na(nestling_number) == F &
                                          is.na(obs_med_temp) == F),
)

plot(brooding_blups_lmm)
# Doesn't look great, but not sure what to do with it...

summary(brooding_blups_lmm)


# Calucate the BLUPs, individual variation in brooding duration
# from the best fitting model (above)

## a) Generate best fit model summary
ranef(brooding_blups_lmm) # random effect
fixef(brooding_blups_lmm) # fixed effect
coef(brooding_blups_lmm) # fixed effect

## b) extract BLUPs from mixed model object
brooding_blups <- as.data.frame(ranef(brooding_blups_lmm)) # extract ranef as
# a dataframe, BLUPs = rand effects + intercept (from poiss/neg. binom)

## c) Rename variables in blups table
brooding_blups <- brooding_blups %>%
  rename('fnest_id' = 'grp') %>%
  rename('brooding_ranef' = 'condval') %>%
  rename('brooding_ranef_sd' = 'condsd') %>%
  select(c('fnest_id', 'brooding_ranef', 'brooding_ranef_sd'))


## d) extract fixed effect (intercept) from poisson/neg. binomial model
brooding_intrcpt <- (fixef(brooding_blups_lmm)[[1]])[[1]] # fixed effect
#* INTERPETATION ****#  
# Interpet as expected log count of behavior when controlling for /
# holding constant the effects of covariates (or if exponentiated
# the intercept is the incident rate of the expected count
# of behavior) as a proporiton of time overlap (the offset)

## f) Create a new variable that is ranef plus both poisson/neg. binomial
# model intercept. This provides estimates of individual level
# variation in proportion of time spent of obs vs feeding visits
brooding_blups <-  brooding_blups  %>%
  mutate(brooding_blups = brooding_intrcpt + brooding_ranef)
#* INTERPETATION ****#  
# Interpret as the conditional mode or the log counts / incident rate
# of expected counts as a proportion of the overlap time (offset) for
# each individual...while holding constant effect of other covariates

## Left join feeding BLUPs to late_nestling_parent_care
late_nestling_parent_care <- late_nestling_parent_care %>%
  left_join(feeding_blups, by = c('nest_id' = 'fnest_id'), 
            copy = F)

## Left join brooding BLUPs to late_nestling_parent_care
late_nestling_parent_care <- late_nestling_parent_care %>%
  left_join(brooding_blups, by = c('nest_id' = 'fnest_id'), 
            copy = F)


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
for(i in 1:length(late_nestling_parent_care$hatch_date)){
  late_nestling_parent_care$days_summer[i] <- 
    calculate_days_of_summer(late_nestling_parent_care$hatch_date[i])
}

late_nestling_parent_care$days_summer <- as.numeric(late_nestling_parent_care$days_summer)


save(file = 'Data/tidy_parent_nestl_weather_data_8-23_with_pci.RData', 
     list = c('prim_merged', 'nestl_merged', 'govee_daily',
              'noaa', 'late_nestling_parent_care'))



