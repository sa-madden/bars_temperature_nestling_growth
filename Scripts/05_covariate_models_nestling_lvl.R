####### Purpose: run covariate models for barn swallow
####### nest microclimate and nestling growth dataset
####### By: Sage Madden
####### Created: 1/16/2023
####### Last modified: 1/16/2023

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
library('lubridate')
library("Hmisc")

## Graph Plotting and Visualization Packages
library ('ggplot2')
library('hrbrthemes')
library('viridis')
library ('gridExtra')
library ('dotwhisker')

## Modelling Packages
library('lme4')
library ('emmeans')
library('MuMIn')
library('car')

# load pbkrtest and lmertest (emmeans dependency)
library('pbkrtest')
library('lmerTest')
# prevent lmerTest from masking the lme4 function for lmer
lmer <- lme4::lmer

## Broom packages
library('broom')
library('broom.mixed')

## Model checking packages
library('DHARMa')


### 1.3 Get Version and Session Info
R.Version()
sessionInfo()

###############################################################################
##############                    2. Load RData                  ##############
###############################################################################  

### 2.1 Load RData
## Load RData tidy barn swallow data
load('Data/tidy_parent_nestl_weather_data_with_pci.RData')

# Make site a factor
late_nestling_parent_care$fsite <- as.factor(late_nestling_parent_care$site)
late_nestling_parent_care$fnest_id <- as.factor(late_nestling_parent_care$nest_id)

late_nestling_parent_care$num_pairs <- as.numeric(late_nestling_parent_care$num_pairs)


###############################################################################
##############      3. Growth                                   ##############
###############################################################################

## Mass and hatch date -- not significant, but marginal  
mass_hatch_lmer <- lmer(mass_pre_obs ~ scale(days_summer) + (1|fsite) +
                                   (1|fnest_id), 
                                 data = subset(late_nestling_parent_care,
                                               !is.na(x = mass_pre_obs) & 
                                                 !is.na(x = days_summer)))

## Check diagnostics for the full model
plot(mass_hatch_lmer)
# Normal QQplot
{qqnorm(resid(mass_hatch_lmer))
  qqline(resid(mass_hatch_lmer))}
# Histogram of residuals
hist(resid(mass_hatch_lmer))
# Checking for influential outliers
infIndexPlot(mass_hatch_lmer, vars=c("Cook"))
infIndexPlot(mass_hatch_lmer, vars=c("Studentized"))

summary(mass_hatch_lmer)
confint(mass_hatch_lmer)

## Mass and colony size -- not significant 
mass_colony_lmer <- lmer(mass_pre_obs ~ scale(num_pairs) + (1|fsite) +
                          (1|fnest_id), 
                        data = subset(late_nestling_parent_care,
                                      !is.na(x = mass_pre_obs) & 
                                        !is.na(x = num_pairs)))

## Check diagnostics for the full model
plot(mass_colony_lmer)
# Normal QQplot
{qqnorm(resid(mass_colony_lmer))
  qqline(resid(mass_colony_lmer))}
# Histogram of residuals
hist(resid(mass_colony_lmer))
# Checking for influential outliers
infIndexPlot(mass_colony_lmer, vars=c("Cook"))
infIndexPlot(mass_colony_lmer, vars=c("Studentized"))

summary(mass_colony_lmer)
confint(mass_colony_lmer)

## Mass and brood size -- not significant, but marginal
mass_brood_lmer <- lmer(mass_pre_obs ~ scale(nestling_number) + (1|fsite) +
                           (1|fnest_id), 
                         data = subset(late_nestling_parent_care,
                                       !is.na(x = mass_pre_obs) & 
                                         !is.na(x = nestling_number)))

## Check diagnostics for the full model
plot(mass_brood_lmer)
# Normal QQplot
{qqnorm(resid(mass_brood_lmer))
  qqline(resid(mass_brood_lmer))}
# Histogram of residuals
hist(resid(mass_brood_lmer))
# Checking for influential outliers
infIndexPlot(mass_brood_lmer, vars=c("Cook"))
infIndexPlot(mass_brood_lmer, vars=c("Studentized"))

summary(mass_brood_lmer)
confint(mass_brood_lmer)

## Mass and nestling age -- not significant 
mass_age_lmer <- lmer(mass_pre_obs ~ scale(nestling_age) + (1|fsite) +
                          (1|fnest_id), 
                        data = subset(late_nestling_parent_care,
                                      !is.na(x = mass_pre_obs) & 
                                        !is.na(x = nestling_age)))

## Check diagnostics for the full model
plot(mass_age_lmer)
# Normal QQplot
{qqnorm(resid(mass_age_lmer))
  qqline(resid(mass_age_lmer))}
# Histogram of residuals
hist(resid(mass_age_lmer))
# Checking for influential outliers
infIndexPlot(mass_age_lmer, vars=c("Cook"))
infIndexPlot(mass_age_lmer, vars=c("Studentized"))

summary(mass_age_lmer)
confint(mass_age_lmer)


### Wing chord

## wing and hatch date -- not significant, but marginal  
wing_hatch_lmer <- lmer(rt_wing_length ~ scale(days_summer) + (1|fsite) +
                          (1|fnest_id), 
                        data = subset(late_nestling_parent_care,
                                      !is.na(x = rt_wing_length) & 
                                        !is.na(x = days_summer)))

## Check diagnostics for the full model
plot(wing_hatch_lmer)
# Normal QQplot
{qqnorm(resid(wing_hatch_lmer))
  qqline(resid(wing_hatch_lmer))}
# Histogram of residuals
hist(resid(wing_hatch_lmer))
# Checking for influential outliers
infIndexPlot(wing_hatch_lmer, vars=c("Cook"))
infIndexPlot(wing_hatch_lmer, vars=c("Studentized"))

summary(wing_hatch_lmer)
confint(wing_hatch_lmer)

## wing and colony size -- not significant 
wing_colony_lmer <- lmer(rt_wing_length ~ scale(num_pairs) + (1|fsite) +
                           (1|fnest_id), 
                         data = subset(late_nestling_parent_care,
                                       !is.na(x = rt_wing_length) & 
                                         !is.na(x = num_pairs)))

## Check diagnostics for the full model
plot(wing_colony_lmer)
# Normal QQplot
{qqnorm(resid(wing_colony_lmer))
  qqline(resid(wing_colony_lmer))}
# Histogram of residuals
hist(resid(wing_colony_lmer))
# Checking for influential outliers
infIndexPlot(wing_colony_lmer, vars=c("Cook"))
infIndexPlot(wing_colony_lmer, vars=c("Studentized"))

summary(wing_colony_lmer)
confint(wing_colony_lmer)

## wing and brood size -- not significant, but marginal
wing_brood_lmer <- lmer(rt_wing_length ~ scale(nestling_number) + (1|fsite) +
                          (1|fnest_id), 
                        data = subset(late_nestling_parent_care,
                                      !is.na(x = rt_wing_length) & 
                                        !is.na(x = nestling_number)))

## Check diagnostics for the full model
plot(wing_brood_lmer)
# Normal QQplot
{qqnorm(resid(wing_brood_lmer))
  qqline(resid(wing_brood_lmer))}
# Histogram of residuals
hist(resid(wing_brood_lmer))
# Checking for influential outliers
infIndexPlot(wing_brood_lmer, vars=c("Cook"))
infIndexPlot(wing_brood_lmer, vars=c("Studentized"))

summary(wing_brood_lmer)
confint(wing_brood_lmer)

## wing and nestling age -- not significant 
wing_age_lmer <- lmer(rt_wing_length ~ scale(nestling_age) + (1|fsite) +
                        (1|fnest_id), 
                      data = subset(late_nestling_parent_care,
                                    !is.na(x = rt_wing_length) & 
                                      !is.na(x = nestling_age)))

## Check diagnostics for the full model
plot(wing_age_lmer)
# Normal QQplot
{qqnorm(resid(wing_age_lmer))
  qqline(resid(wing_age_lmer))}
# Histogram of residuals
hist(resid(wing_age_lmer))
# Checking for influential outliers
infIndexPlot(wing_age_lmer, vars=c("Cook"))
infIndexPlot(wing_age_lmer, vars=c("Studentized"))

summary(wing_age_lmer)
confint(wing_age_lmer)



###############################################################################
##############        Temperature                               ##############
###############################################################################
# Remove duplicates for nest level measures
one_measure <- distinct(late_nestling_parent_care, nest_id, .keep_all = TRUE)


## max_temp and hatch date -- not significant, but marginal  
max_temp_hatch_lmer <- lmer(nest_max_temp ~ scale(days_summer) + (1|fsite), 
                            data = subset(one_measure,
                                          !is.na(x = nest_max_temp) & 
                                            !is.na(x = days_summer)))

## Check diagnostics for the full model
plot(max_temp_hatch_lmer)
# Normal QQplot
{qqnorm(resid(max_temp_hatch_lmer))
  qqline(resid(max_temp_hatch_lmer))}
# Histogram of residuals
hist(resid(max_temp_hatch_lmer))
# Checking for influential outliers
infIndexPlot(max_temp_hatch_lmer, vars=c("Cook"))
infIndexPlot(max_temp_hatch_lmer, vars=c("Studentized"))

summary(max_temp_hatch_lmer)
confint(max_temp_hatch_lmer)


## min_temp and hatch date -- not significant, but marginal  
min_temp_hatch_lmer <- lmer(nest_min_temp ~ scale(days_summer) + (1|fsite), 
                            data = subset(one_measure,
                                          !is.na(x = nest_min_temp) & 
                                            !is.na(x = days_summer)))

## Check diagnostics for the full model
plot(min_temp_hatch_lmer)
# Normal QQplot
{qqnorm(resid(min_temp_hatch_lmer))
  qqline(resid(min_temp_hatch_lmer))}
# Histogram of residuals
hist(resid(min_temp_hatch_lmer))
# Checking for influential outliers
infIndexPlot(min_temp_hatch_lmer, vars=c("Cook"))
infIndexPlot(min_temp_hatch_lmer, vars=c("Studentized"))

summary(min_temp_hatch_lmer)
confint(min_temp_hatch_lmer)


## iqr_temp and hatch date -- not significant, but marginal  
iqr_temp_hatch_lmer <- lmer(nest_iqr_temp ~ scale(days_summer) + (1|fsite), 
                            data = subset(one_measure,
                                          !is.na(x = nest_iqr_temp) & 
                                            !is.na(x = days_summer)))

## Check diagnostics for the full model
plot(iqr_temp_hatch_lmer)
# Normal QQplot
{qqnorm(resid(iqr_temp_hatch_lmer))
  qqline(resid(iqr_temp_hatch_lmer))}
# Histogram of residuals
hist(resid(iqr_temp_hatch_lmer))
# Checking for influential outliers
infIndexPlot(iqr_temp_hatch_lmer, vars=c("Cook"))
infIndexPlot(iqr_temp_hatch_lmer, vars=c("Studentized"))

summary(iqr_temp_hatch_lmer)
confint(iqr_temp_hatch_lmer)

###############################################################################
##############      Parental care                                ##############
###############################################################################
## feeding_blups and hatch date -- not significant, but marginal  
feeding_blups_hatch_lmer <- lmer(feeding_expontd_blups ~ scale(days_summer) + (1|fsite), 
                        data = subset(one_measure,
                                      !is.na(x = feeding_expontd_blups) & 
                                        !is.na(x = days_summer)))

## Check diagnostics for the full model
plot(feeding_blups_hatch_lmer)
# Normal QQplot
{qqnorm(resid(feeding_blups_hatch_lmer))
  qqline(resid(feeding_blups_hatch_lmer))}
# Histogram of residuals
hist(resid(feeding_blups_hatch_lmer))
# Checking for influential outliers
infIndexPlot(feeding_blups_hatch_lmer, vars=c("Cook"))
infIndexPlot(feeding_blups_hatch_lmer, vars=c("Studentized"))

summary(feeding_blups_hatch_lmer)
confint(feeding_blups_hatch_lmer)

## feeding_blups and colony size -- not significant 
feeding_blups_colony_lmer <- lmer(feeding_expontd_blups ~ scale(num_pairs) + (1|fsite), 
                         data = subset(one_measure,
                                       !is.na(x = feeding_expontd_blups) & 
                                         !is.na(x = num_pairs)))

## Check diagnostics for the full model
plot(feeding_blups_colony_lmer)
# Normal QQplot
{qqnorm(resid(feeding_blups_colony_lmer))
  qqline(resid(feeding_blups_colony_lmer))}
# Histogram of residuals
hist(resid(feeding_blups_colony_lmer))
# Checking for influential outliers
infIndexPlot(feeding_blups_colony_lmer, vars=c("Cook"))
infIndexPlot(feeding_blups_colony_lmer, vars=c("Studentized"))

summary(feeding_blups_colony_lmer)
confint(feeding_blups_colony_lmer)

## feeding_blups and brood size -- not significant, but marginal
feeding_blups_brood_lmer <- lmer(feeding_expontd_blups ~ scale(nestling_number) + (1|fsite), 
                        data = subset(one_measure,
                                      !is.na(x = feeding_expontd_blups) & 
                                        !is.na(x = nestling_number)))

## Check diagnostics for the full model
plot(feeding_blups_brood_lmer)
# Normal QQplot
{qqnorm(resid(feeding_blups_brood_lmer))
  qqline(resid(feeding_blups_brood_lmer))}
# Histogram of residuals
hist(resid(feeding_blups_brood_lmer))
# Checking for influential outliers
infIndexPlot(feeding_blups_brood_lmer, vars=c("Cook"))
infIndexPlot(feeding_blups_brood_lmer, vars=c("Studentized"))

summary(feeding_blups_brood_lmer)
confint(feeding_blups_brood_lmer)

## feeding_blups and nestling age -- not significant 
feeding_blups_age_lmer <- lmer(feeding_expontd_blups ~ scale(nestling_age) + (1|fsite), 
                      data = subset(one_measure,
                                    !is.na(x = feeding_expontd_blups) & 
                                      !is.na(x = nestling_age)))

## Check diagnostics for the full model
plot(feeding_blups_age_lmer)
# Normal QQplot
{qqnorm(resid(feeding_blups_age_lmer))
  qqline(resid(feeding_blups_age_lmer))}
# Histogram of residuals
hist(resid(feeding_blups_age_lmer))
# Checking for influential outliers
infIndexPlot(feeding_blups_age_lmer, vars=c("Cook"))
infIndexPlot(feeding_blups_age_lmer, vars=c("Studentized"))

summary(feeding_blups_age_lmer)
confint(feeding_blups_age_lmer)


## brooding_blups and hatch date -- not significant, but marginal  
brooding_blups_hatch_lmer <- lmer(brooding_blups ~ scale(days_summer) + (1|fsite), 
                                 data = subset(one_measure,
                                               !is.na(x = brooding_blups) & 
                                                 !is.na(x = days_summer)))

## Check diagnostics for the full model
plot(brooding_blups_hatch_lmer)
# Normal QQplot
{qqnorm(resid(brooding_blups_hatch_lmer))
  qqline(resid(brooding_blups_hatch_lmer))}
# Histogram of residuals
hist(resid(brooding_blups_hatch_lmer))
# Checking for influential outliers
infIndexPlot(brooding_blups_hatch_lmer, vars=c("Cook"))
infIndexPlot(brooding_blups_hatch_lmer, vars=c("Studentized"))

summary(brooding_blups_hatch_lmer)
confint(brooding_blups_hatch_lmer)

## brooding_blupsand colony size -- not significant 
brooding_blups_colony_lmer <- lmer(brooding_blups ~ scale(num_pairs) + (1|fsite), 
                                  data = subset(one_measure,
                                                !is.na(x = brooding_blups) & 
                                                  !is.na(x = num_pairs)))

## Check diagnostics for the full model
plot(brooding_blups_colony_lmer)
# Normal QQplot
{qqnorm(resid(brooding_blups_colony_lmer))
  qqline(resid(brooding_blups_colony_lmer))}
# Histogram of residuals
hist(resid(brooding_blups_colony_lmer))
# Checking for influential outliers
infIndexPlot(brooding_blups_colony_lmer, vars=c("Cook"))
infIndexPlot(brooding_blups_colony_lmer, vars=c("Studentized"))

summary(brooding_blups_colony_lmer)
confint(brooding_blups_colony_lmer)

## brooding_blupsand brood size -- not significant, but marginal
brooding_blups_brood_lmer <- lmer(brooding_blups ~ scale(nestling_number) + (1|fsite), 
                                 data = subset(one_measure,
                                               !is.na(x = brooding_blups) & 
                                                 !is.na(x = nestling_number)))

## Check diagnostics for the full model
plot(brooding_blups_brood_lmer)
# Normal QQplot
{qqnorm(resid(brooding_blups_brood_lmer))
  qqline(resid(brooding_blups_brood_lmer))}
# Histogram of residuals
hist(resid(brooding_blups_brood_lmer))
# Checking for influential outliers
infIndexPlot(brooding_blups_brood_lmer, vars=c("Cook"))
infIndexPlot(brooding_blups_brood_lmer, vars=c("Studentized"))

summary(brooding_blups_brood_lmer)
confint(brooding_blups_brood_lmer)

## brooding_blupsand nestling age -- not significant 
brooding_blups_age_lmer <- lmer(brooding_blups ~ scale(nestling_age) + (1|fsite), 
                               data = subset(one_measure,
                                             !is.na(x = brooding_blups) & 
                                               !is.na(x = nestling_age)))

## Check diagnostics for the full model
plot(brooding_blups_age_lmer)
# Normal QQplot
{qqnorm(resid(brooding_blups_age_lmer))
  qqline(resid(brooding_blups_age_lmer))}
# Histogram of residuals
hist(resid(brooding_blups_age_lmer))
# Checking for influential outliers
infIndexPlot(brooding_blups_age_lmer, vars=c("Cook"))
infIndexPlot(brooding_blups_age_lmer, vars=c("Studentized"))

summary(brooding_blups_age_lmer)
confint(brooding_blups_age_lmer)





###############################################################################
##############        Growth + parental care                     ##############
###############################################################################

## Mass and feeding BLUPs
mass_feeding_lmer <- lmer(mass_pre_obs ~ feeding_expontd_blups + (1|fsite) +
                          (1|fnest_id), 
                        data = subset(late_nestling_parent_care,
                                      !is.na(x = mass_pre_obs) & 
                                        !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_feeding_lmer)
# Normal QQplot
{qqnorm(resid(mass_feeding_lmer))
  qqline(resid(mass_feeding_lmer))}
# Histogram of residuals
hist(resid(mass_feeding_lmer))
# Checking for influential outliers
infIndexPlot(mass_feeding_lmer, vars=c("Cook"))
infIndexPlot(mass_feeding_lmer, vars=c("Studentized"))

summary(mass_feeding_lmer)
confint(mass_feeding_lmer)

## Mass and brooding BLUPs
mass_brooding_lmer <- lmer(mass_pre_obs ~ brooding_blups + (1|fsite) +
                            (1|fnest_id), 
                          data = subset(late_nestling_parent_care,
                                        !is.na(x = mass_pre_obs) & 
                                          !is.na(x = brooding_blups)))

## Check diagnostics for the full model
plot(mass_brooding_lmer)
# Normal QQplot
{qqnorm(resid(mass_brooding_lmer))
  qqline(resid(mass_brooding_lmer))}
# Histogram of residuals
hist(resid(mass_brooding_lmer))
# Checking for influential outliers
infIndexPlot(mass_brooding_lmer, vars=c("Cook"))
infIndexPlot(mass_brooding_lmer, vars=c("Studentized"))

summary(mass_brooding_lmer)
confint(mass_brooding_lmer)

## Wing and feeding BLUPs
wing_feeding_lmer <- lmer(rt_wing_length ~ feeding_expontd_blups + (1|fsite) +
                            (1|fnest_id), 
                          data = subset(late_nestling_parent_care,
                                        !is.na(x = rt_wing_length) & 
                                          !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(wing_feeding_lmer)
# Normal QQplot
{qqnorm(resid(wing_feeding_lmer))
  qqline(resid(wing_feeding_lmer))}
# Histogram of residuals
hist(resid(wing_feeding_lmer))
# Checking for influential outliers
infIndexPlot(wing_feeding_lmer, vars=c("Cook"))
infIndexPlot(wing_feeding_lmer, vars=c("Studentized"))

summary(wing_feeding_lmer)
confint(wing_feeding_lmer)

## wing and brooding BLUPs
wing_brooding_lmer <- lmer(rt_wing_length ~ brooding_blups + (1|fsite) +
                             (1|fnest_id), 
                           data = subset(late_nestling_parent_care,
                                         !is.na(x = rt_wing_length) & 
                                           !is.na(x = brooding_blups)))

## Check diagnostics for the full model
plot(wing_brooding_lmer)
# Normal QQplot
{qqnorm(resid(wing_brooding_lmer))
  qqline(resid(wing_brooding_lmer))}
# Histogram of residuals
hist(resid(wing_brooding_lmer))
# Checking for influential outliers
infIndexPlot(wing_brooding_lmer, vars=c("Cook"))
infIndexPlot(wing_brooding_lmer, vars=c("Studentized"))

summary(wing_brooding_lmer)
confint(wing_brooding_lmer)


avg_parent_care <- prim_merged %>% group_by(nest_id) %>%
  summarise(avg_feeding = mean(total_feeding_visits),
            avg_brooding = mean(total_brooding_duration))

joined_care <- left_join(late_nestling_parent_care, avg_parent_care,
                         by = "nest_id")

## Mass and average feeding rate
mass_feeding_avg_lmer <- lmer(mass_pre_obs ~ avg_feeding + (1|fsite) +
                            (1|fnest_id), 
                          data = subset(joined_care,
                                        !is.na(x = mass_pre_obs) & 
                                          !is.na(x = avg_feeding)))

## Check diagnostics for the full model
plot(mass_feeding_avg_lmer)
# Normal QQplot
{qqnorm(resid(mass_feeding_avg_lmer))
  qqline(resid(mass_feeding_avg_lmer))}
# Histogram of residuals
hist(resid(mass_feeding_avg_lmer))
# Checking for influential outliers
infIndexPlot(mass_feeding_avg_lmer, vars=c("Cook"))
infIndexPlot(mass_feeding_avg_lmer, vars=c("Studentized"))

summary(mass_feeding_avg_lmer)
confint(mass_feeding_avg_lmer)

## Mass and avg brooding duration
joined_care$avg_brooding_min <- joined_care$avg_brooding/3600
mass_brooding_avg_lmer <- lmer(mass_pre_obs ~ avg_brooding_min + (1|fsite) +
                                (1|fnest_id), 
                              data = subset(joined_care,
                                            !is.na(x = mass_pre_obs) & 
                                              !is.na(x = avg_brooding_min)))

## Check diagnostics for the full model
plot(mass_brooding_avg_lmer)
# Normal QQplot
{qqnorm(resid(mass_brooding_avg_lmer))
  qqline(resid(mass_brooding_avg_lmer))}
# Histogram of residuals
hist(resid(mass_brooding_avg_lmer))
# Checking for influential outliers
infIndexPlot(mass_brooding_avg_lmer, vars=c("Cook"))
infIndexPlot(mass_brooding_avg_lmer, vars=c("Studentized"))

summary(mass_brooding_avg_lmer)
confint(mass_brooding_avg_lmer)


## Wing and average feeding rate
wing_feeding_avg_lmer <- lmer(rt_wing_length ~ avg_feeding + (1|fsite) +
                                (1|fnest_id), 
                              data = subset(joined_care,
                                            !is.na(x = rt_wing_length) & 
                                              !is.na(x = avg_feeding)))

## Check diagnostics for the full model
plot(wing_feeding_avg_lmer)
# Normal QQplot
{qqnorm(resid(wing_feeding_avg_lmer))
  qqline(resid(wing_feeding_avg_lmer))}
# Histogram of residuals
hist(resid(wing_feeding_avg_lmer))
# Checking for influential outliers
infIndexPlot(wing_feeding_avg_lmer, vars=c("Cook"))
infIndexPlot(wing_feeding_avg_lmer, vars=c("Studentized"))

summary(wing_feeding_avg_lmer)
confint(wing_feeding_avg_lmer)

## wing and avg brooding duration
wing_brooding_avg_lmer <- lmer(rt_wing_length ~ avg_brooding_min + (1|fsite) +
                                 (1|fnest_id), 
                               data = subset(joined_care,
                                             !is.na(x = rt_wing_length) & 
                                               !is.na(x = avg_brooding_min)))

## Check diagnostics for the full model
plot(wing_brooding_avg_lmer)
# Normal QQplot
{qqnorm(resid(wing_brooding_avg_lmer))
  qqline(resid(wing_brooding_avg_lmer))}
# Histogram of residuals
hist(resid(wing_brooding_avg_lmer))
# Checking for influential outliers
infIndexPlot(wing_brooding_avg_lmer, vars=c("Cook"))
infIndexPlot(wing_brooding_avg_lmer, vars=c("Studentized"))

summary(wing_brooding_avg_lmer)
confint(wing_brooding_avg_lmer)

