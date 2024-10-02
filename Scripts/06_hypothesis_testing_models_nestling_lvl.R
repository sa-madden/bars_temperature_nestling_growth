####### Purpose: run models to test hypotheses for barn swallow
####### nest mircoclimate and nestling growth dataset
####### By: Sage Madden
####### Created: 12/19/2022
####### Last modified: 1/16/2023

###############################################################################
##############             1.  Configure work space              ##############
###############################################################################

### 1Global options
## clear global environment
rm(list = ls())

## bprevent R from automatically reading charater strins as factors
options(stringsAsFactors = FALSE)


### Install and load CRAN packages   
## Data Manipulation and Descriptive Stats Packages
library('tidyverse')
library('lubridate')
library("Hmisc")
library("dplyr")

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


### Get Version and Session Info
R.Version()
sessionInfo()

###############################################################################
##############                    2. Load RData                  ##############
###############################################################################  

### Load RData
## Load RData tidy barn swallow data
load('Data/tidy_parent_nestl_weather_data_with_pci.RData')

# Make site a factor
late_nestling_parent_care$fsite <- as.factor(late_nestling_parent_care$site)
late_nestling_parent_care$fnest_id <- as.factor(late_nestling_parent_care$nest_id)

late_nestling_parent_care$num_pairs <- as.numeric(late_nestling_parent_care$num_pairs)


###############################################################################
##############      3. Growth, temp, and parental care      ##############
###############################################################################
##### FEEDING 
#### Minimum temp
### Mass unadjusted
mass_min_temp_blups_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) *
                                   scale(feeding_expontd_blups) + (1|fsite) +
                                   (1|fnest_id), 
                                 data = subset(late_nestling_parent_care,
                                               !is.na(x = mass_pre_obs) & 
                                                 !is.na(x = nest_min_temp) &
                                                 !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_min_temp_blups_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_blups_lmer))
  qqline(resid(mass_min_temp_blups_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_blups_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_temp_blups_lmer, vars=c("Cook"))
infIndexPlot(mass_min_temp_blups_lmer, vars=c("Studentized"))

summary(mass_min_temp_blups_lmer)
confint(mass_min_temp_blups_lmer)  

## Create three categories for BLUPs to allow stratification 
late_nestling_parent_care <- late_nestling_parent_care %>%
  mutate(feeding_expontd_blups_strat =  as.integer(Hmisc::cut2(feeding_expontd_blups, g=3)))

# Hmisc alternative using cut and quantile   
#    (cut(care.indx.cont.3, 
#         quantile(care.indx.cont.3, probs=0:3/3,
#                  na.rm = T), 
#         include.lowest = T)))
# don't use ntile, which allows the same value to occur in mult. quantile         
#as.factor(ntile(care.indx.cont.3, 3))) 

## LOW parental care model
mass_min_temp_blups_low_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) + 
                                       (1|fsite) +
                                   (1|fnest_id), 
                                 data = subset(late_nestling_parent_care,
                                               feeding_expontd_blups_strat == 1,
                                               !is.na(x = mass_pre_obs) & 
                                                 !is.na(x = nest_min_temp) &
                                                 !is.na(x = feeding_expontd_blups)))

summary(mass_min_temp_blups_low_lmer)
confint(mass_min_temp_blups_low_lmer) 

## MED parental care
mass_min_temp_blups_med_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) + 
                                       (1|fsite) +
                                       (1|fnest_id), 
                                     data = subset(late_nestling_parent_care,
                                                   feeding_expontd_blups_strat == 2,
                                                   !is.na(x = mass_pre_obs) & 
                                                     !is.na(x = nest_min_temp) &
                                                     !is.na(x = feeding_expontd_blups)))

summary(mass_min_temp_blups_med_lmer)
confint(mass_min_temp_blups_med_lmer) 

## HIGH parental care
mass_min_temp_blups_high_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) + 
                                       (1|fsite) +
                                       (1|fnest_id), 
                                     data = subset(late_nestling_parent_care,
                                                   feeding_expontd_blups_strat == 3,
                                                   !is.na(x = mass_pre_obs) & 
                                                     !is.na(x = nest_min_temp) &
                                                     !is.na(x = feeding_expontd_blups)))

summary(mass_min_temp_blups_high_lmer)
confint(mass_min_temp_blups_high_lmer) 


### Mass adjusted
mass_min_temp_blups_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) *
                                   scale(feeding_expontd_blups) + scale(nestling_number) + 
                                     scale(days_summer) + (1|fsite) +
                                   (1|fnest_id), 
                                 data = subset(late_nestling_parent_care,
                                               !is.na(x = mass_pre_obs) & 
                                                 !is.na(x = nest_min_temp) &
                                                 !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_min_temp_blups_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_blups_adj_lmer))
  qqline(resid(mass_min_temp_blups_adj_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_blups_adj_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_temp_blups_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_min_temp_blups_adj_lmer, vars=c("Studentized"))

summary(mass_min_temp_blups_adj_lmer)
confint(mass_min_temp_blups_adj_lmer) 

## LOW parental care model - adjusted
mass_min_temp_blups_low_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                           (1|fsite) +
                                           (1|fnest_id), 
                                         data = subset(late_nestling_parent_care,
                                                       feeding_expontd_blups_strat == 1,
                                                       !is.na(x = mass_pre_obs) & 
                                                         !is.na(x = nest_min_temp) &
                                                         !is.na(x = feeding_expontd_blups)))

summary(mass_min_temp_blups_low_adj_lmer)
confint(mass_min_temp_blups_low_adj_lmer) 

## MED parental care - adjusted
mass_min_temp_blups_med_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                         (1|fsite) +
                                       (1|fnest_id), 
                                     data = subset(late_nestling_parent_care,
                                                   feeding_expontd_blups_strat == 2,
                                                   !is.na(x = mass_pre_obs) & 
                                                     !is.na(x = nest_min_temp) &
                                                     !is.na(x = feeding_expontd_blups)))

summary(mass_min_temp_blups_med_adj_lmer)
confint(mass_min_temp_blups_med_adj_lmer) 

## HIGH parental care - adjusted 
mass_min_temp_blups_high_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) +
                                            scale(nestling_number) + scale(days_summer) + 
                                        (1|fsite) +
                                        (1|fnest_id), 
                                      data = subset(late_nestling_parent_care,
                                                    feeding_expontd_blups_strat == 3,
                                                    !is.na(x = mass_pre_obs) & 
                                                      !is.na(x = nest_min_temp) &
                                                      !is.na(x = feeding_expontd_blups)))

summary(mass_min_temp_blups_high_adj_lmer)
confint(mass_min_temp_blups_high_adj_lmer) 

### Mass adjusted with outliers removed
mass_outliers_removed <- late_nestling_parent_care[-c(4, 7, 44), ]

mass_min_temp_blups_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) *
                                       scale(feeding_expontd_blups) + scale(nestling_number) + 
                                       scale(days_summer) + (1|fsite) +
                                       (1|fnest_id), 
                                     data = subset(mass_outliers_removed,
                                                   !is.na(x = mass_pre_obs) & 
                                                     !is.na(x = nest_min_temp) &
                                                     !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_min_temp_blups_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_blups_adj_noout_lmer))
  qqline(resid(mass_min_temp_blups_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_blups_adj_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_temp_blups_adj_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_min_temp_blups_adj_noout_lmer, vars=c("Studentized"))

summary(mass_min_temp_blups_adj_noout_lmer)
confint(mass_min_temp_blups_adj_noout_lmer) 

## LOW parental care model - adjusted
mass_min_temp_blups_low_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                           (1|fsite) +
                                           (1|fnest_id), 
                                         data = subset(mass_outliers_removed,
                                                       feeding_expontd_blups_strat == 1,
                                                       !is.na(x = mass_pre_obs) & 
                                                         !is.na(x = nest_min_temp) &
                                                         !is.na(x = feeding_expontd_blups)))

summary(mass_min_temp_blups_low_adj_noout_lmer)
confint(mass_min_temp_blups_low_adj_noout_lmer) 

## MED parental care - adjusted
mass_min_temp_blups_med_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                           (1|fsite) +
                                           (1|fnest_id), 
                                         data = subset(mass_outliers_removed,
                                                       feeding_expontd_blups_strat == 2,
                                                       !is.na(x = mass_pre_obs) & 
                                                         !is.na(x = nest_min_temp) &
                                                         !is.na(x = feeding_expontd_blups)))

summary(mass_min_temp_blups_med_adj_noout_lmer)
confint(mass_min_temp_blups_med_adj_noout_lmer) 

## HIGH parental care - adjusted 
mass_min_temp_blups_high_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) +
                                            scale(nestling_number) + scale(days_summer) + 
                                            (1|fsite) +
                                            (1|fnest_id), 
                                          data = subset(mass_outliers_removed,
                                                        feeding_expontd_blups_strat == 3,
                                                        !is.na(x = mass_pre_obs) & 
                                                          !is.na(x = nest_min_temp) &
                                                          !is.na(x = feeding_expontd_blups)))

summary(mass_min_temp_blups_high_adj_noout_lmer)
confint(mass_min_temp_blups_high_adj_noout_lmer) 

### Wing unadjusted
wing_min_temp_blups_lmer <- lmer(rt_wing_length ~ scale(nest_min_temp) *
                                   scale(feeding_expontd_blups) + (1|fsite) +
                                   (1|fnest_id), 
                                 data = subset(mass_outliers_removed,
                                               !is.na(x = rt_wing_length) & 
                                                 !is.na(x = nest_min_temp) &
                                                 !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(wing_min_temp_blups_lmer)
# Normal QQplot
{qqnorm(resid(wing_min_temp_blups_lmer))
  qqline(resid(wing_min_temp_blups_lmer))}
# Histogram of residuals
hist(resid(wing_min_temp_blups_lmer))
# Checking for influential outliers
infIndexPlot(wing_min_temp_blups_lmer, vars=c("Cook"))
infIndexPlot(wing_min_temp_blups_lmer, vars=c("Studentized"))

summary(wing_min_temp_blups_lmer)
confint(wing_min_temp_blups_lmer)

# Model without interaction 
wing_min_temp_blups_noint_lmer <- lmer(rt_wing_length ~ scale(nest_min_temp) +
                                   (1|fsite) +
                                   (1|fnest_id), 
                                 data = subset(late_nestling_parent_care,
                                               !is.na(x = rt_wing_length) & 
                                                 !is.na(x = nest_min_temp) &
                                                 !is.na(x = feeding_expontd_blups)))

summary(wing_min_temp_blups_noint_lmer)
confint(wing_min_temp_blups_noint_lmer)


### Wing adjusted
wing_min_temp_blups_adj_lmer <- lmer(rt_wing_length ~ scale(nest_min_temp) *
                                       scale(feeding_expontd_blups) + scale(nestling_number) + 
                                       scale(days_summer) + (1|fsite) +
                                       (1|fnest_id), 
                                     data = subset(late_nestling_parent_care,
                                                   !is.na(x = rt_wing_length) & 
                                                     !is.na(x = nest_min_temp) &
                                                     !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(wing_min_temp_blups_adj_lmer)
# Normal QQplot
{qqnorm(resid(wing_min_temp_blups_adj_lmer))
  qqline(resid(wing_min_temp_blups_adj_lmer))}
# Histogram of residuals
hist(resid(wing_min_temp_blups_adj_lmer))
# Checking for influential outliers
infIndexPlot(wing_min_temp_blups_adj_lmer, vars=c("Cook"))
infIndexPlot(wing_min_temp_blups_adj_lmer, vars=c("Studentized"))

summary(wing_min_temp_blups_adj_lmer)
confint(wing_min_temp_blups_adj_lmer) 


#### Maximum temp
## Mass unadjusted
mass_max_temp_blups_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) *
                                   scale(feeding_expontd_blups) + (1|fsite) +
                                   (1|fnest_id), 
                                 data = subset(late_nestling_parent_care,
                                               !is.na(x = mass_pre_obs) & 
                                                 !is.na(x = nest_max_temp) &
                                                 !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_max_temp_blups_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_blups_lmer))
  qqline(resid(mass_max_temp_blups_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_blups_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_temp_blups_lmer, vars=c("Cook"))
infIndexPlot(mass_max_temp_blups_lmer, vars=c("Studentized"))

summary(mass_max_temp_blups_lmer)
confint(mass_max_temp_blups_lmer)  

# Model without interaction
mass_max_temp_blups_noint_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                   (1|fsite) +
                                   (1|fnest_id), 
                                 data = subset(late_nestling_parent_care,
                                               !is.na(x = mass_pre_obs) & 
                                                 !is.na(x = nest_max_temp) &
                                                 !is.na(x = feeding_expontd_blups)))

summary(mass_max_temp_blups_noint_lmer)
confint(mass_max_temp_blups_noint_lmer)

### Mass adjusted
mass_max_temp_blups_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) *
                                   scale(feeding_expontd_blups) + 
                                     scale(nestling_number) + scale(days_summer) +
                                     (1|fsite) +
                                   (1|fnest_id), 
                                 data = subset(late_nestling_parent_care,
                                               !is.na(x = mass_pre_obs) & 
                                                 !is.na(x = nest_max_temp) &
                                                 !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_max_temp_blups_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_blups_adj_lmer))
  qqline(resid(mass_max_temp_blups_adj_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_blups_adj_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_temp_blups_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_max_temp_blups_adj_lmer, vars=c("Studentized"))

summary(mass_max_temp_blups_adj_lmer)
confint(mass_max_temp_blups_adj_lmer)  

# Model without interaction
mass_max_temp_blups_noint_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                          scale(nestling_number) + scale(days_summer) +
                                         (1|fsite) +
                                         (1|fnest_id), 
                                       data = subset(late_nestling_parent_care,
                                                     !is.na(x = mass_pre_obs) & 
                                                       !is.na(x = nest_max_temp) &
                                                       !is.na(x = feeding_expontd_blups)))

summary(mass_max_temp_blups_noint_adj_lmer)
confint(mass_max_temp_blups_noint_adj_lmer)

### Mass adjusted -- outliers removed
mass_max_temp_blups_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) *
                                       scale(feeding_expontd_blups) + 
                                       scale(nestling_number) + scale(days_summer) +
                                       (1|fsite) +
                                       (1|fnest_id), 
                                     data = subset(mass_outliers_removed,
                                                   !is.na(x = mass_pre_obs) & 
                                                     !is.na(x = nest_max_temp) &
                                                     !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_max_temp_blups_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_blups_adj_noout_lmer))
  qqline(resid(mass_max_temp_blups_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_blups_adj_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_temp_blups_adj_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_max_temp_blups_adj_noout_lmer, vars=c("Studentized"))

summary(mass_max_temp_blups_adj_noout_lmer)
confint(mass_max_temp_blups_adj_noout_lmer)

# Model without interaction
mass_max_temp_blups_noint_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                             scale(nestling_number) + scale(days_summer) +
                                             (1|fsite) +
                                             (1|fnest_id), 
                                           data = subset(mass_outliers_removed,
                                                         !is.na(x = mass_pre_obs) & 
                                                           !is.na(x = nest_max_temp) &
                                                           !is.na(x = feeding_expontd_blups)))

summary(mass_max_temp_blups_noint_adj_noout_lmer)
confint(mass_max_temp_blups_noint_adj_noout_lmer)

### Wing unadjusted
wing_max_temp_blups_lmer <- lmer(rt_wing_length ~ scale(nest_max_temp) *
                                   scale(feeding_expontd_blups) + (1|fsite) +
                                   (1|fnest_id), 
                                 data = subset(late_nestling_parent_care,
                                               !is.na(x = rt_wing_length) & 
                                                 !is.na(x = nest_max_temp) &
                                                 !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(wing_max_temp_blups_lmer)
# Normal QQplot
{qqnorm(resid(wing_max_temp_blups_lmer))
  qqline(resid(wing_max_temp_blups_lmer))}
# Histogram of residuals
hist(resid(wing_max_temp_blups_lmer))
# Checking for influential outliers
infIndexPlot(wing_max_temp_blups_lmer, vars=c("Cook"))
infIndexPlot(wing_max_temp_blups_lmer, vars=c("Studentized"))

summary(wing_max_temp_blups_lmer)
confint(wing_max_temp_blups_lmer)

# Model with no interaction 
wing_max_temp_blups_noint_lmer <- lmer(rt_wing_length ~ scale(nest_max_temp) +
                                   (1|fsite) +
                                   (1|fnest_id), 
                                 data = subset(late_nestling_parent_care,
                                               !is.na(x = rt_wing_length) & 
                                                 !is.na(x = nest_max_temp) &
                                                 !is.na(x = feeding_expontd_blups)))

summary(wing_max_temp_blups_noint_lmer)
confint(wing_max_temp_blups_noint_lmer)

### Wing adjusted
wing_max_temp_blups_adj_lmer <- lmer(rt_wing_length ~ scale(nest_max_temp) *
                                       scale(feeding_expontd_blups) + scale(nestling_number) + 
                                       scale(days_summer) + (1|fsite) +
                                       (1|fnest_id), 
                                     data = subset(late_nestling_parent_care,
                                                   !is.na(x = rt_wing_length) & 
                                                     !is.na(x = nest_max_temp) &
                                                     !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(wing_max_temp_blups_adj_lmer)
# Normal QQplot
{qqnorm(resid(wing_max_temp_blups_adj_lmer))
  qqline(resid(wing_max_temp_blups_adj_lmer))}
# Histogram of residuals
hist(resid(wing_max_temp_blups_adj_lmer))
# Checking for influential outliers
infIndexPlot(wing_max_temp_blups_adj_lmer, vars=c("Cook"))
infIndexPlot(wing_max_temp_blups_adj_lmer, vars=c("Studentized"))

summary(wing_max_temp_blups_adj_lmer)
confint(wing_max_temp_blups_adj_lmer) 


#### IQR of temp
### Mass unadjusted
mass_iqr_temp_blups_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) *
                                   scale(feeding_expontd_blups) + (1|fsite) +
                                   (1|fnest_id), 
                                 data = subset(late_nestling_parent_care,
                                               !is.na(x = mass_pre_obs) & 
                                                 !is.na(x = nest_iqr_temp) &
                                                 !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_iqr_temp_blups_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_blups_lmer))
  qqline(resid(mass_iqr_temp_blups_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_blups_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_temp_blups_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_temp_blups_lmer, vars=c("Studentized"))

summary(mass_iqr_temp_blups_lmer)
confint(mass_iqr_temp_blups_lmer)  

# Model without interaction
mass_iqr_temp_blups_noint_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) +
                                   (1|fsite) +
                                   (1|fnest_id), 
                                 data = subset(late_nestling_parent_care,
                                               !is.na(x = mass_pre_obs) & 
                                                 !is.na(x = nest_iqr_temp) &
                                                 !is.na(x = feeding_expontd_blups)))

summary(mass_iqr_temp_blups_noint_lmer)
confint(mass_iqr_temp_blups_noint_lmer)


### Mass adjusted
mass_iqr_temp_blups_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) *
                                   scale(feeding_expontd_blups) + 
                                     scale(nestling_number) + scale(days_summer) +
                                     (1|fsite) +
                                   (1|fnest_id), 
                                 data = subset(late_nestling_parent_care,
                                               !is.na(x = mass_pre_obs) & 
                                                 !is.na(x = nest_iqr_temp) &
                                                 !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_iqr_temp_blups_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_blups_adj_lmer))
  qqline(resid(mass_iqr_temp_blups_adj_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_blups_adj_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_temp_blups_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_temp_blups_adj_lmer, vars=c("Studentized"))

summary(mass_iqr_temp_blups_adj_lmer)
confint(mass_iqr_temp_blups_adj_lmer) 

# Model without interaction
mass_iqr_temp_blups_noint_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) +
                                         scale(nestling_number) + scale(days_summer) + 
                                         (1|fsite) +
                                         (1|fnest_id), 
                                       data = subset(late_nestling_parent_care,
                                                     !is.na(x = mass_pre_obs) & 
                                                       !is.na(x = nest_iqr_temp) &
                                                       !is.na(x = feeding_expontd_blups)))

summary(mass_iqr_temp_blups_noint_adj_lmer)
confint(mass_iqr_temp_blups_noint_adj_lmer)

### Mass adjusted -- outliers removed 
mass_iqr_temp_blups_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) *
                                       scale(feeding_expontd_blups) + 
                                       scale(nestling_number) + scale(days_summer) +
                                       (1|fsite) +
                                       (1|fnest_id), 
                                     data = subset(mass_outliers_removed,
                                                   !is.na(x = mass_pre_obs) & 
                                                     !is.na(x = nest_iqr_temp) &
                                                     !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_iqr_temp_blups_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_blups_adj_noout_lmer))
  qqline(resid(mass_iqr_temp_blups_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_blups_adj_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_temp_blups_adj_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_temp_blups_adj_noout_lmer, vars=c("Studentized"))

summary(mass_iqr_temp_blups_adj_noout_lmer)
confint(mass_iqr_temp_blups_adj_noout_lmer)

# Model without interaction
mass_iqr_temp_blups_noint_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) +
                                             scale(nestling_number) + scale(days_summer) + 
                                             (1|fsite) +
                                             (1|fnest_id), 
                                           data = subset( mass_outliers_removed,
                                                         !is.na(x = mass_pre_obs) & 
                                                           !is.na(x = nest_iqr_temp) &
                                                           !is.na(x = feeding_expontd_blups)))

summary(mass_iqr_temp_blups_noint_adj_noout_lmer)
confint(mass_iqr_temp_blups_noint_adj_noout_lmer)

### Wing unadjusted
wing_iqr_temp_blups_lmer <- lmer(rt_wing_length ~ scale(nest_iqr_temp) *
                                   scale(feeding_expontd_blups) + (1|fsite) +
                                   (1|fnest_id), 
                                 data = subset(late_nestling_parent_care,
                                               !is.na(x = rt_wing_length) & 
                                                 !is.na(x = nest_iqr_temp) &
                                                 !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(wing_iqr_temp_blups_lmer)
# Normal QQplot
{qqnorm(resid(wing_iqr_temp_blups_lmer))
  qqline(resid(wing_iqr_temp_blups_lmer))}
# Histogram of residuals
hist(resid(wing_iqr_temp_blups_lmer))
# Checking for influential outliers
infIndexPlot(wing_iqr_temp_blups_lmer, vars=c("Cook"))
infIndexPlot(wing_iqr_temp_blups_lmer, vars=c("Studentized"))

summary(wing_iqr_temp_blups_lmer)
confint(wing_iqr_temp_blups_lmer)

# Model without interaction 
wing_iqr_temp_blups_noint_lmer <- lmer(rt_wing_length ~ scale(nest_iqr_temp) +
                                   (1|fsite) +
                                   (1|fnest_id), 
                                 data = subset(late_nestling_parent_care,
                                               !is.na(x = rt_wing_length) & 
                                                 !is.na(x = nest_iqr_temp) &
                                                 !is.na(x = feeding_expontd_blups)))

summary(wing_iqr_temp_blups_noint_lmer)
confint(wing_iqr_temp_blups_noint_lmer)

### Wing adjusted
wing_iqr_temp_blups_adj_lmer <- lmer(rt_wing_length ~ scale(nest_iqr_temp) *
                                       scale(feeding_expontd_blups) + scale(nestling_number) + 
                                       scale(days_summer) + (1|fsite) +
                                       (1|fnest_id), 
                                     data = subset(late_nestling_parent_care,
                                                   !is.na(x = rt_wing_length) & 
                                                     !is.na(x = nest_iqr_temp) &
                                                     !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(wing_iqr_temp_blups_adj_lmer)
# Normal QQplot
{qqnorm(resid(wing_iqr_temp_blups_adj_lmer))
  qqline(resid(wing_iqr_temp_blups_adj_lmer))}
# Histogram of residuals
hist(resid(wing_iqr_temp_blups_adj_lmer))
# Checking for influential outliers
infIndexPlot(wing_iqr_temp_blups_adj_lmer, vars=c("Cook"))
infIndexPlot(wing_iqr_temp_blups_adj_lmer, vars=c("Studentized"))

summary(wing_iqr_temp_blups_adj_lmer)
confint(wing_iqr_temp_blups_adj_lmer)


##### BROODING
#### Minimum temp
### Mass unadjusted
mass_min_temp_brooding_blups_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) *
                                   scale(brooding_blups) + (1|fsite) +
                                   (1|fnest_id), 
                                 data = subset(late_nestling_parent_care,
                                               !is.na(x = mass_pre_obs) & 
                                                 !is.na(x = nest_min_temp) &
                                                 !is.na(x = brooding_blups)))

## Check diagnostics for the full model
plot(mass_min_temp_brooding_blups_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_brooding_blups_lmer))
  qqline(resid(mass_min_temp_brooding_blups_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_brooding_blups_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_temp_brooding_blups_lmer, vars=c("Cook"))
infIndexPlot(mass_min_temp_brooding_blups_lmer, vars=c("Studentized"))

summary(mass_min_temp_brooding_blups_lmer)
confint(mass_min_temp_brooding_blups_lmer) 

# Model without interaction 
mass_min_temp_brooding_blups_noint_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) +
                                            (1|fsite) +
                                            (1|fnest_id), 
                                          data = subset(late_nestling_parent_care,
                                                        !is.na(x = mass_pre_obs) & 
                                                          !is.na(x = nest_min_temp) &
                                                          !is.na(x = brooding_blups)))

summary(mass_min_temp_brooding_blups_noint_lmer)
confint(mass_min_temp_brooding_blups_noint_lmer)

### Mass adjusted
mass_min_temp_brooding_blups_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) *
                                       scale(brooding_blups) + scale(nestling_number) + 
                                       scale(days_summer) + (1|fsite) +
                                       (1|fnest_id), 
                                     data = subset(late_nestling_parent_care,
                                                   !is.na(x = mass_pre_obs) & 
                                                     !is.na(x = nest_min_temp) &
                                                     !is.na(x = brooding_blups)))

## Check diagnostics for the full model
plot(mass_min_temp_brooding_blups_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_brooding_blups_adj_lmer))
  qqline(resid(mass_min_temp_brooding_blups_adj_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_brooding_blups_adj_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_temp_brooding_blups_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_min_temp_brooding_blups_adj_lmer, vars=c("Studentized"))

summary(mass_min_temp_brooding_blups_adj_lmer)
confint(mass_min_temp_brooding_blups_adj_lmer)


mass_brooding_outliers_removed <- late_nestling_parent_care[-c(4, 7), ]

### Mass adjusted -- outliers removed 
mass_min_temp_brooding_blups_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) *
                                             scale(brooding_blups) + 
                                             scale(nestling_number) + scale(days_summer) +
                                             (1|fsite) +
                                             (1|fnest_id), 
                                           data = subset(mass_brooding_outliers_removed,
                                                         !is.na(x = mass_pre_obs) & 
                                                           !is.na(x = nest_min_temp) &
                                                           !is.na(x = brooding_blups)))

## Check diagnostics for the full model
plot(mass_min_temp_brooding_blups_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_brooding_blups_adj_noout_lmer))
  qqline(resid(mass_min_temp_brooding_blups_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_brooding_blups_adj_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_temp_brooding_blups_adj_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_min_temp_brooding_blups_adj_noout_lmer, vars=c("Studentized"))

summary(mass_min_temp_brooding_blups_adj_noout_lmer)
confint(mass_min_temp_brooding_blups_adj_noout_lmer)


### Wing unadjusted
wing_min_temp_brooding_blups_lmer <- lmer(rt_wing_length ~ scale(nest_min_temp) *
                                   scale(brooding_blups) + (1|fsite) +
                                   (1|fnest_id), 
                                 data = subset(late_nestling_parent_care,
                                               !is.na(x = rt_wing_length) & 
                                                 !is.na(x = nest_min_temp) &
                                                 !is.na(x = brooding_blups)))

## Check diagnostics for the full model
plot(wing_min_temp_brooding_blups_lmer)
# Normal QQplot
{qqnorm(resid(wing_min_temp_brooding_blups_lmer))
  qqline(resid(wing_min_temp_brooding_blups_lmer))}
# Histogram of residuals
hist(resid(wing_min_temp_brooding_blups_lmer))
# Checking for influential outliers
infIndexPlot(wing_min_temp_brooding_blups_lmer, vars=c("Cook"))
infIndexPlot(wing_min_temp_brooding_blups_lmer, vars=c("Studentized"))

summary(wing_min_temp_brooding_blups_lmer)
confint(wing_min_temp_brooding_blups_lmer)

## Create three categories for BLUPs to allow stratification 
late_nestling_parent_care <- late_nestling_parent_care %>%
  mutate(brooding_blups_strat =  as.integer(Hmisc::cut2(brooding_blups, g=3)))

# LOW parental care 
wing_min_temp_brooding_blups_low_lmer <- lmer(rt_wing_length ~ scale(nest_min_temp) +
                                                (1|fsite) +
                                                (1|fnest_id), 
                                              data = subset(late_nestling_parent_care,
                                                            brooding_blups_strat == 1,
                                                            !is.na(x = rt_wing_length) & 
                                                              !is.na(x = nest_min_temp) &
                                                              !is.na(x = brooding_blups)))

summary(wing_min_temp_brooding_blups_low_lmer)
confint(wing_min_temp_brooding_blups_low_lmer)

# MED parental care
wing_min_temp_brooding_blups_med_lmer <- lmer(rt_wing_length ~ scale(nest_min_temp) +
                                                (1|fsite) +
                                            (1|fnest_id), 
                                          data = subset(late_nestling_parent_care,
                                                        brooding_blups_strat == 2,
                                                        !is.na(x = rt_wing_length) & 
                                                          !is.na(x = nest_min_temp) &
                                                          !is.na(x = brooding_blups)))

summary(wing_min_temp_brooding_blups_med_lmer)
confint(wing_min_temp_brooding_blups_med_lmer)

# HIGH parental care
wing_min_temp_brooding_blups_high_lmer <- lmer(rt_wing_length ~ scale(nest_min_temp) +
                                                (1|fsite) +
                                                (1|fnest_id), 
                                              data = subset(late_nestling_parent_care,
                                                            brooding_blups_strat == 3,
                                                            !is.na(x = rt_wing_length) & 
                                                              !is.na(x = nest_min_temp) &
                                                              !is.na(x = brooding_blups)))

summary(wing_min_temp_brooding_blups_high_lmer)
confint(wing_min_temp_brooding_blups_high_lmer)


### Wing adjusted
wing_min_temp_brooding_blups_adj_lmer <- lmer(rt_wing_length ~ scale(nest_min_temp) *
                                                scale(brooding_blups) + scale(nestling_number) + 
                                                scale(days_summer) + (1|fsite) +
                                                (1|fnest_id), 
                                              data = subset(late_nestling_parent_care,
                                                            !is.na(x = rt_wing_length) & 
                                                              !is.na(x = nest_min_temp) &
                                                              !is.na(x = brooding_blups)))

## Check diagnostics for the full model
plot(wing_min_temp_brooding_blups_adj_lmer)
# Normal QQplot
{qqnorm(resid(wing_min_temp_brooding_blups_adj_lmer))
  qqline(resid(wing_min_temp_brooding_blups_adj_lmer))}
# Histogram of residuals
hist(resid(wing_min_temp_brooding_blups_adj_lmer))
# Checking for influential outliers
infIndexPlot(wing_min_temp_brooding_blups_adj_lmer, vars=c("Cook"))
infIndexPlot(wing_min_temp_brooding_blups_adj_lmer, vars=c("Studentized"))

summary(wing_min_temp_brooding_blups_adj_lmer)
confint(wing_min_temp_brooding_blups_adj_lmer)


#### Maximum temp
### Mass unadjusted
mass_max_temp_brooding_blups_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) *
                                            scale(brooding_blups) + (1|fsite) +
                                            (1|fnest_id), 
                                          data = subset(late_nestling_parent_care,
                                                        !is.na(x = mass_pre_obs) & 
                                                          !is.na(x = nest_max_temp) &
                                                          !is.na(x = brooding_blups)))

## Check diagnostics for the full model
plot(mass_max_temp_brooding_blups_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_brooding_blups_lmer))
  qqline(resid(mass_max_temp_brooding_blups_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_brooding_blups_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_temp_brooding_blups_lmer, vars=c("Cook"))
infIndexPlot(mass_max_temp_brooding_blups_lmer, vars=c("Studentized"))

summary(mass_max_temp_brooding_blups_lmer)
confint(mass_max_temp_brooding_blups_lmer)

# Model without interaction 
mass_max_temp_brooding_blups_noint_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                            (1|fsite) +
                                            (1|fnest_id), 
                                          data = subset(late_nestling_parent_care,
                                                        !is.na(x = mass_pre_obs) & 
                                                          !is.na(x = nest_max_temp) &
                                                          !is.na(x = brooding_blups)))

summary(mass_max_temp_brooding_blups_noint_lmer)
confint(mass_max_temp_brooding_blups_noint_lmer)

### Mass adjusted
mass_max_temp_brooding_blups_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) *
                                                scale(brooding_blups) + scale(nestling_number) + 
                                                scale(days_summer) + (1|fsite) +
                                                (1|fnest_id), 
                                              data = subset(late_nestling_parent_care,
                                                            !is.na(x = mass_pre_obs) & 
                                                              !is.na(x = nest_max_temp) &
                                                              !is.na(x = brooding_blups)))

## Check diagnostics for the full model
plot(mass_max_temp_brooding_blups_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_brooding_blups_adj_lmer))
  qqline(resid(mass_max_temp_brooding_blups_adj_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_brooding_blups_adj_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_temp_brooding_blups_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_max_temp_brooding_blups_adj_lmer, vars=c("Studentized"))

summary(mass_max_temp_brooding_blups_adj_lmer)
confint(mass_max_temp_brooding_blups_adj_lmer)

### Mass adjusted -- outliers removed 
mass_max_temp_brooding_blups_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) *
                                                      scale(brooding_blups) + 
                                                      scale(nestling_number) + scale(days_summer) +
                                                      (1|fsite) +
                                                      (1|fnest_id), 
                                                    data = subset(mass_brooding_outliers_removed,
                                                                  !is.na(x = mass_pre_obs) & 
                                                                    !is.na(x = nest_max_temp) &
                                                                    !is.na(x = brooding_blups)))

## Check diagnostics for the full model
plot(mass_max_temp_brooding_blups_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_brooding_blups_adj_noout_lmer))
  qqline(resid(mass_max_temp_brooding_blups_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_brooding_blups_adj_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_temp_brooding_blups_adj_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_max_temp_brooding_blups_adj_noout_lmer, vars=c("Studentized"))

summary(mass_max_temp_brooding_blups_adj_noout_lmer)
confint(mass_max_temp_brooding_blups_adj_noout_lmer)


### Wing unadjusted
wing_max_temp_brooding_blups_lmer <- lmer(rt_wing_length ~ scale(nest_max_temp) *
                                            scale(brooding_blups) + (1|fsite) +
                                            (1|fnest_id), 
                                          data = subset(late_nestling_parent_care,
                                                        !is.na(x = rt_wing_length) & 
                                                          !is.na(x = nest_max_temp) &
                                                          !is.na(x = brooding_blups)))

## Check diagnostics for the full model
plot(wing_max_temp_brooding_blups_lmer)
# Normal QQplot
{qqnorm(resid(wing_max_temp_brooding_blups_lmer))
  qqline(resid(wing_max_temp_brooding_blups_lmer))}
# Histogram of residuals
hist(resid(wing_max_temp_brooding_blups_lmer))
# Checking for influential outliers
infIndexPlot(wing_max_temp_brooding_blups_lmer, vars=c("Cook"))
infIndexPlot(wing_max_temp_brooding_blups_lmer, vars=c("Studentized"))

summary(wing_max_temp_brooding_blups_lmer)
confint(wing_max_temp_brooding_blups_lmer)

# Model without interaction 
wing_max_temp_brooding_blups_noint_lmer <- lmer(rt_wing_length ~ scale(nest_max_temp) +
                                            (1|fsite) +
                                            (1|fnest_id), 
                                          data = subset(late_nestling_parent_care,
                                                        !is.na(x = rt_wing_length) & 
                                                          !is.na(x = nest_max_temp) &
                                                          !is.na(x = brooding_blups)))

summary(wing_max_temp_brooding_blups_noint_lmer)
confint(wing_max_temp_brooding_blups_noint_lmer)

### Wing adjusted
wing_max_temp_brooding_blups_adj_lmer <- lmer(rt_wing_length ~ scale(nest_max_temp) *
                                                scale(brooding_blups) + scale(nestling_number) + 
                                                scale(days_summer) + (1|fsite) +
                                                (1|fnest_id), 
                                              data = subset(late_nestling_parent_care,
                                                            !is.na(x = rt_wing_length) & 
                                                              !is.na(x = nest_max_temp) &
                                                              !is.na(x = brooding_blups)))

## Check diagnostics for the full model
plot(wing_max_temp_brooding_blups_adj_lmer)
# Normal QQplot
{qqnorm(resid(wing_max_temp_brooding_blups_adj_lmer))
  qqline(resid(wing_max_temp_brooding_blups_adj_lmer))}
# Histogram of residuals
hist(resid(wing_max_temp_brooding_blups_adj_lmer))
# Checking for influential outliers
infIndexPlot(wing_max_temp_brooding_blups_adj_lmer, vars=c("Cook"))
infIndexPlot(wing_max_temp_brooding_blups_adj_lmer, vars=c("Studentized"))

summary(wing_max_temp_brooding_blups_adj_lmer)
confint(wing_max_temp_brooding_blups_adj_lmer)

#### IQR of temp
### Mass unadjusted
mass_iqr_temp_brooding_blups_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) *
                                            scale(brooding_blups) + (1|fsite) +
                                            (1|fnest_id), 
                                          data = subset(late_nestling_parent_care,
                                                        !is.na(x = mass_pre_obs) & 
                                                          !is.na(x = nest_iqr_temp) &
                                                          !is.na(x = brooding_blups)))

## Check diagnostics for the full model
plot(mass_iqr_temp_brooding_blups_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_brooding_blups_lmer))
  qqline(resid(mass_iqr_temp_brooding_blups_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_brooding_blups_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_temp_brooding_blups_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_temp_brooding_blups_lmer, vars=c("Studentized"))

summary(mass_iqr_temp_brooding_blups_lmer)
confint(mass_iqr_temp_brooding_blups_lmer)  

# Model without interaction  
mass_iqr_temp_brooding_blups_noint_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) +
                                            (1|fsite) +
                                            (1|fnest_id), 
                                          data = subset(late_nestling_parent_care,
                                                        !is.na(x = mass_pre_obs) & 
                                                          !is.na(x = nest_iqr_temp) &
                                                          !is.na(x = brooding_blups)))
summary(mass_iqr_temp_brooding_blups_noint_lmer)
confint(mass_iqr_temp_brooding_blups_noint_lmer)

### Mass adjusted
mass_iqr_temp_brooding_blups_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) *
                                                scale(brooding_blups) + scale(nestling_number) + 
                                                scale(days_summer) + (1|fsite) +
                                                (1|fnest_id), 
                                              data = subset(late_nestling_parent_care,
                                                            !is.na(x = mass_pre_obs) & 
                                                              !is.na(x = nest_iqr_temp) &
                                                              !is.na(x = brooding_blups)))

## Check diagnostics for the full model
plot(mass_iqr_temp_brooding_blups_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_brooding_blups_adj_lmer))
  qqline(resid(mass_iqr_temp_brooding_blups_adj_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_brooding_blups_adj_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_temp_brooding_blups_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_temp_brooding_blups_adj_lmer, vars=c("Studentized"))

summary(mass_iqr_temp_brooding_blups_adj_lmer)
confint(mass_iqr_temp_brooding_blups_adj_lmer)

### Mass adjusted -- outliers removed 
mass_iqr_temp_brooding_blups_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) *
                                                      scale(brooding_blups) + 
                                                      scale(nestling_number) + scale(days_summer) +
                                                      (1|fsite) +
                                                      (1|fnest_id), 
                                                    data = subset(mass_brooding_outliers_removed,
                                                                  !is.na(x = mass_pre_obs) & 
                                                                    !is.na(x = nest_iqr_temp) &
                                                                    !is.na(x = brooding_blups)))

## Check diagnostics for the full model
plot(mass_iqr_temp_brooding_blups_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_brooding_blups_adj_noout_lmer))
  qqline(resid(mass_iqr_temp_brooding_blups_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_brooding_blups_adj_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_temp_brooding_blups_adj_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_temp_brooding_blups_adj_noout_lmer, vars=c("Studentized"))

summary(mass_iqr_temp_brooding_blups_adj_noout_lmer)
confint(mass_iqr_temp_brooding_blups_adj_noout_lmer)


### Wing unadjusted
wing_iqr_temp_brooding_blups_lmer <- lmer(rt_wing_length ~ scale(nest_iqr_temp) *
                                            scale(brooding_blups) + (1|fsite) +
                                            (1|fnest_id), 
                                          data = subset(late_nestling_parent_care,
                                                        !is.na(x = rt_wing_length) & 
                                                          !is.na(x = nest_iqr_temp) &
                                                          !is.na(x = brooding_blups)))

## Check diagnostics for the full model
plot(wing_iqr_temp_brooding_blups_lmer)
# Normal QQplot
{qqnorm(resid(wing_iqr_temp_brooding_blups_lmer))
  qqline(resid(wing_iqr_temp_brooding_blups_lmer))}
# Histogram of residuals
hist(resid(wing_iqr_temp_brooding_blups_lmer))
# Checking for influential outliers
infIndexPlot(wing_iqr_temp_brooding_blups_lmer, vars=c("Cook"))
infIndexPlot(wing_iqr_temp_brooding_blups_lmer, vars=c("Studentized"))

summary(wing_iqr_temp_brooding_blups_lmer)
confint(wing_iqr_temp_brooding_blups_lmer)

# Model without interaction 
wing_iqr_temp_brooding_blups_noint_lmer <- lmer(rt_wing_length ~ scale(nest_iqr_temp) +
                                            (1|fsite) +
                                            (1|fnest_id), 
                                          data = subset(late_nestling_parent_care,
                                                        !is.na(x = rt_wing_length) & 
                                                          !is.na(x = nest_iqr_temp) &
                                                          !is.na(x = brooding_blups)))

summary(wing_iqr_temp_brooding_blups_noint_lmer)
confint(wing_iqr_temp_brooding_blups_noint_lmer)

### Wing adjusted
wing_iqr_temp_brooding_blups_adj_lmer <- lmer(rt_wing_length ~ scale(nest_iqr_temp) *
                                                scale(brooding_blups) + scale(nestling_number) + 
                                                scale(days_summer) + (1|fsite) +
                                                (1|fnest_id), 
                                              data = subset(late_nestling_parent_care,
                                                            !is.na(x = rt_wing_length) & 
                                                              !is.na(x = nest_iqr_temp) &
                                                              !is.na(x = brooding_blups)))

## Check diagnostics for the full model
plot(wing_iqr_temp_brooding_blups_adj_lmer)
# Normal QQplot
{qqnorm(resid(wing_iqr_temp_brooding_blups_adj_lmer))
  qqline(resid(wing_iqr_temp_brooding_blups_adj_lmer))}
# Histogram of residuals
hist(resid(wing_iqr_temp_brooding_blups_adj_lmer))
# Checking for influential outliers
infIndexPlot(wing_iqr_temp_brooding_blups_adj_lmer, vars=c("Cook"))
infIndexPlot(wing_iqr_temp_brooding_blups_adj_lmer, vars=c("Studentized"))

summary(wing_iqr_temp_brooding_blups_adj_lmer)
confint(wing_iqr_temp_brooding_blups_adj_lmer)

###############################################################################
##############      3. Growth, temp, and relative size     ##############
###############################################################################

########################### MID SIZE #######################################

#### Minimum temp
### Mass unadjusted
mass_min_temp_mid_size_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) *
                                  mid_size_order + (1|fsite) +
                                  (1|fnest_id), 
                                data = subset(late_nestling_parent_care,
                                              !is.na(x = mid_size_order) & 
                                                !is.na(x = nest_min_temp)&
                                                !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_temp_mid_size_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_mid_size_lmer))
  qqline(resid(mass_min_temp_mid_size_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_mid_size_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_temp_mid_size_lmer, vars=c("Cook"))
infIndexPlot(mass_min_temp_mid_size_lmer, vars=c("Studentized"))

summary(mass_min_temp_mid_size_lmer)
confint(mass_min_temp_mid_size_lmer)


### Wing unadjusted
wing_min_temp_mid_size_lmer <- lmer(rt_wing_length ~ scale(nest_min_temp) *
                                  mid_size_order + (1|fsite) +
                                  (1|fnest_id), 
                                data = subset(late_nestling_parent_care,
                                              !is.na(x = rt_wing_length) & 
                                                !is.na(x = mid_size_order)&
                                                !is.na(x = nest_min_temp)))

## Check diagnostics for the full model
plot(wing_min_temp_mid_size_lmer)
# Normal QQplot
{qqnorm(resid(wing_min_temp_mid_size_lmer))
  qqline(resid(wing_min_temp_mid_size_lmer))}
# Histogram of residuals
hist(resid(wing_min_temp_mid_size_lmer))
# Checking for influential outliers
infIndexPlot(wing_min_temp_mid_size_lmer, vars=c("Cook"))
infIndexPlot(wing_min_temp_mid_size_lmer, vars=c("Studentized"))

summary(wing_min_temp_mid_size_lmer)
confint(wing_min_temp_mid_size_lmer)


### Mass adjusted
mass_min_temp_mid_size_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) *
                                      mid_size_order + scale(nestling_number) +
                                      scale(days_summer) + (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  !is.na(x = mid_size_order) & 
                                                    !is.na(x = nest_min_temp)&
                                                    !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_temp_mid_size_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_mid_size_adj_lmer))
  qqline(resid(mass_min_temp_mid_size_adj_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_mid_size_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_temp_mid_size_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_min_temp_mid_size_adj_lmer, vars=c("Studentized"))

summary(mass_min_temp_mid_size_adj_lmer)
confint(mass_min_temp_mid_size_adj_lmer)


#### Maximum temp
### Mass unadjusted
mass_max_temp_mid_size_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) *
                                  mid_size_order + (1|fsite) +
                                  (1|fnest_id), 
                                data = subset(late_nestling_parent_care,
                                              !is.na(x = mass_pre_obs) & 
                                                !is.na(x = nest_max_temp) &
                                                !is.na(x = mid_size_order)))

## Check diagnostics for the full model
plot(mass_max_temp_mid_size_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_mid_size_lmer))
  qqline(resid(mass_max_temp_mid_size_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_mid_size_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_temp_mid_size_lmer, vars=c("Cook"))
infIndexPlot(mass_max_temp_mid_size_lmer, vars=c("Studentized"))

summary(mass_max_temp_mid_size_lmer)
confint(mass_max_temp_mid_size_lmer) 

# SMALL nestlings
mass_max_temp_mid_size_small_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                        (1|fsite), 
                                      data = subset(late_nestling_parent_care,
                                                    mid_size_order == "min",
                                                    !is.na(x = mass_pre_obs) & 
                                                      !is.na(x = nest_max_temp) &
                                                      !is.na(x = mid_size_order)))

summary(mass_max_temp_mid_size_small_lmer)
confint(mass_max_temp_mid_size_small_lmer)

# BIG nestlings
mass_max_temp_mid_size_big_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                      (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  mid_size_order == "other",
                                                  !is.na(x = mass_pre_obs) & 
                                                    !is.na(x = nest_max_temp) &
                                                    !is.na(x = mid_size_order)))

summary(mass_max_temp_mid_size_big_lmer)
confint(mass_max_temp_mid_size_big_lmer)


### Wing unadjusted
wing_max_temp_mid_size_lmer <- lmer(rt_wing_length ~ scale(nest_max_temp) *
                                  mid_size_order + (1|fsite) +
                                  (1|fnest_id), 
                                data = subset(late_nestling_parent_care,
                                              !is.na(x = rt_wing_length) & 
                                                !is.na(x = mid_size_order)&
                                                !is.na(x = nest_max_temp)))

## Check diagnostics for the full model
plot(wing_max_temp_mid_size_lmer)
# Normal QQplot
{qqnorm(resid(wing_max_temp_mid_size_lmer))
  qqline(resid(wing_max_temp_mid_size_lmer))}
# Histogram of residuals
hist(resid(wing_max_temp_mid_size_lmer))
# Checking for influential outliers
infIndexPlot(wing_max_temp_mid_size_lmer, vars=c("Cook"))
infIndexPlot(wing_max_temp_mid_size_lmer, vars=c("Studentized"))

summary(wing_max_temp_mid_size_lmer)
confint(wing_max_temp_mid_size_lmer)

# SMALL nestlings 
wing_max_temp_mid_size_small_lmer <- lmer(rt_wing_length ~ scale(nest_max_temp) +
                                        (1|fsite), 
                                      data = subset(late_nestling_parent_care,
                                                    mid_size_order == "min",
                                                    !is.na(x = rt_wing_length) & 
                                                      !is.na(x = mid_size_order)&
                                                      !is.na(x = nest_max_temp)))

summary(wing_max_temp_mid_size_small_lmer)
confint(wing_max_temp_mid_size_small_lmer)

# LARGE nestlings 
wing_max_temp_mid_size_big_lmer <- lmer(rt_wing_length ~ scale(nest_max_temp) +
                                      (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  mid_size_order == "other",
                                                  !is.na(x = rt_wing_length) & 
                                                    !is.na(x = mid_size_order)&
                                                    !is.na(x = nest_max_temp)))

summary(wing_max_temp_mid_size_big_lmer)
confint(wing_max_temp_mid_size_big_lmer)


### Mass adjusted
mass_max_temp_mid_size_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) *
                                      mid_size_order + scale(nestling_number) +
                                      scale(days_summer) + (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  !is.na(x = mass_pre_obs) & 
                                                    !is.na(x = nest_max_temp) &
                                                    !is.na(x = mid_size_order)))

## Check diagnostics for the full model
plot(mass_max_temp_mid_size_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_mid_size_adj_lmer))
  qqline(resid(mass_max_temp_mid_size_adj_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_mid_size_adj_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_temp_mid_size_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_max_temp_mid_size_adj_lmer, vars=c("Studentized"))

summary(mass_max_temp_mid_size_adj_lmer)
confint(mass_max_temp_mid_size_adj_lmer) 

# SMALL nestlings
mass_max_temp_mid_size_small_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                            scale(nestling_number) + scale(days_summer) +
                                            (1|fsite), 
                                          data = subset(late_nestling_parent_care,
                                                        mid_size_order == "min",
                                                        !is.na(x = mass_pre_obs) & 
                                                          !is.na(x = nest_max_temp) &
                                                          !is.na(x = mid_size_order)))

summary(mass_max_temp_mid_size_small_adj_lmer)
confint(mass_max_temp_mid_size_small_adj_lmer)

# BIG nestlings
mass_max_temp_mid_size_big_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                          scale(nestling_number) + scale(days_summer) +
                                          (1|fsite) +
                                          (1|fnest_id), 
                                        data = subset(late_nestling_parent_care,
                                                      mid_size_order == "other",
                                                      !is.na(x = mass_pre_obs) & 
                                                        !is.na(x = nest_max_temp) &
                                                        !is.na(x = mid_size_order)))

summary(mass_max_temp_mid_size_big_adj_lmer)
confint(mass_max_temp_mid_size_big_adj_lmer)


mid_size_outliers_removed <- late_nestling_parent_care[-c(44), ]

### Mass adjusted outliers removed
mass_max_temp_mid_size_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) *
                                            mid_size_order + scale(nestling_number) +
                                            scale(days_summer) + (1|fsite) +
                                            (1|fnest_id), 
                                          data = subset(mid_size_outliers_removed,
                                                        !is.na(x = mass_pre_obs) & 
                                                          !is.na(x = nest_max_temp) &
                                                          !is.na(x = mid_size_order)))

## Check diagnostics for the full model
plot(mass_max_temp_mid_size_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_mid_size_adj_noout_lmer))
  qqline(resid(mass_max_temp_mid_size_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_mid_size_adj_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_temp_mid_size_adj_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_max_temp_mid_size_adj_noout_lmer, vars=c("Studentized"))

summary(mass_max_temp_mid_size_adj_noout_lmer)
confint(mass_max_temp_mid_size_adj_noout_lmer) 

# SMALL nestlings
mass_max_temp_mid_size_small_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                                  scale(nestling_number) + scale(days_summer) +
                                                  (1|fsite), 
                                                data = subset(mid_size_outliers_removed,
                                                              mid_size_order == "min",
                                                              !is.na(x = mass_pre_obs) & 
                                                                !is.na(x = nest_max_temp) &
                                                                !is.na(x = mid_size_order)))

summary(mass_max_temp_mid_size_small_adj_noout_lmer)
confint(mass_max_temp_mid_size_small_adj_noout_lmer)

# BIG nestlings
mass_max_temp_mid_size_big_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                                scale(nestling_number) + scale(days_summer) +
                                                (1|fsite) +
                                                (1|fnest_id), 
                                              data = subset(mid_size_outliers_removed,
                                                            mid_size_order == "other",
                                                            !is.na(x = mass_pre_obs) & 
                                                              !is.na(x = nest_max_temp) &
                                                              !is.na(x = mid_size_order)))

summary(mass_max_temp_mid_size_big_adj_noout_lmer)
confint(mass_max_temp_mid_size_big_adj_noout_lmer)


#### IQR temp
### Mass unadjusted
mass_iqr_temp_mid_size_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) *
                                  mid_size_order + (1|fsite) +
                                  (1|fnest_id), 
                                data = subset(late_nestling_parent_care,
                                              !is.na(x = mass_pre_obs) & 
                                                !is.na(x = nest_iqr_temp) &
                                                !is.na(x = mid_size_order)))

## Check diagnostics for the full model
plot(mass_iqr_temp_mid_size_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_mid_size_lmer))
  qqline(resid(mass_iqr_temp_mid_size_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_mid_size_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_temp_mid_size_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_temp_mid_size_lmer, vars=c("Studentized"))

summary(mass_iqr_temp_mid_size_lmer)
confint(mass_iqr_temp_mid_size_lmer) 

# SMALL nestlings
mass_iqr_temp_mid_size_small_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) +
                                        (1|fsite), 
                                      data = subset(late_nestling_parent_care,
                                                    mid_size_order == "min",
                                                    !is.na(x = mass_pre_obs) & 
                                                      !is.na(x = nest_iqr_temp) &
                                                      !is.na(x = mid_size_order)))

summary(mass_iqr_temp_mid_size_small_lmer)
confint(mass_iqr_temp_mid_size_small_lmer) 

# BIG nestlings
mass_iqr_temp_mid_size_big_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) +
                                      (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  mid_size_order == "other",
                                                  !is.na(x = mass_pre_obs) & 
                                                    !is.na(x = nest_iqr_temp) &
                                                    !is.na(x = mid_size_order)))

summary(mass_iqr_temp_mid_size_big_lmer)
confint(mass_iqr_temp_mid_size_big_lmer) 

### Wing unadjusted
wing_iqr_temp_mid_size_lmer <- lmer(rt_wing_length ~ scale(nest_iqr_temp) *
                                  mid_size_order + (1|fsite) +
                                  (1|fnest_id), 
                                data = subset(late_nestling_parent_care,
                                              !is.na(x = rt_wing_length) & 
                                                !is.na(x = mid_size_order)&
                                                !is.na(x = nest_iqr_temp)))

## Check diagnostics for the full model
plot(wing_iqr_temp_mid_size_lmer)
# Normal QQplot
{qqnorm(resid(wing_iqr_temp_mid_size_lmer))
  qqline(resid(wing_iqr_temp_mid_size_lmer))}
# Histogram of residuals
hist(resid(wing_iqr_temp_mid_size_lmer))
# Checking for influential outliers
infIndexPlot(wing_iqr_temp_mid_size_lmer, vars=c("Cook"))
infIndexPlot(wing_iqr_temp_mid_size_lmer, vars=c("Studentized"))

summary(wing_iqr_temp_mid_size_lmer)
confint(wing_iqr_temp_mid_size_lmer)

# SMALL nestlings 
wing_iqr_temp_mid_size_small_lmer <- lmer(rt_wing_length ~ scale(nest_iqr_temp) +
                                        (1|fsite), 
                                      data = subset(late_nestling_parent_care,
                                                    mid_size_order == "min",
                                                    !is.na(x = rt_wing_length) & 
                                                      !is.na(x = mid_size_order)&
                                                      !is.na(x = nest_iqr_temp)))

summary(wing_iqr_temp_mid_size_small_lmer)
confint(wing_iqr_temp_mid_size_small_lmer)

# BIG nestlings 
wing_iqr_temp_mid_size_big_lmer <- lmer(rt_wing_length ~ scale(nest_iqr_temp) +
                                      (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  mid_size_order == "other",
                                                  !is.na(x = rt_wing_length) & 
                                                    !is.na(x = mid_size_order)&
                                                    !is.na(x = nest_iqr_temp)))

summary(wing_iqr_temp_mid_size_big_lmer)
confint(wing_iqr_temp_mid_size_big_lmer)


# Mass adjusted
mass_iqr_temp_mid_size_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) *
                                      mid_size_order + scale(nestling_number) + 
                                      scale(days_summer) + (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  !is.na(x = mass_pre_obs) & 
                                                    !is.na(x = nest_iqr_temp) &
                                                    !is.na(x = mid_size_order)))

## Check diagnostics for the full model
plot(mass_iqr_temp_mid_size_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_mid_size_adj_lmer))
  qqline(resid(mass_iqr_temp_mid_size_adj_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_mid_size_adj_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_temp_mid_size_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_temp_mid_size_adj_lmer, vars=c("Studentized"))

summary(mass_iqr_temp_mid_size_adj_lmer)
confint(mass_iqr_temp_mid_size_adj_lmer) 

# SMALL nestlings
mass_iqr_temp_mid_size_small_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) + 
                                            scale(nestling_number) + 
                                            scale(days_summer) +
                                            (1|fsite), 
                                          data = subset(late_nestling_parent_care,
                                                        mid_size_order == "min",
                                                        !is.na(x = mass_pre_obs) & 
                                                          !is.na(x = nest_iqr_temp) &
                                                          !is.na(x = mid_size_order)))

summary(mass_iqr_temp_mid_size_small_adj_lmer)
confint(mass_iqr_temp_mid_size_small_adj_lmer) 

# BIG nestlings
mass_iqr_temp_mid_size_big_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) + 
                                          scale(nestling_number) + 
                                          scale(days_summer) +
                                          (1|fsite) +
                                          (1|fnest_id), 
                                        data = subset(late_nestling_parent_care,
                                                      mid_size_order == "other",
                                                      !is.na(x = mass_pre_obs) & 
                                                        !is.na(x = nest_iqr_temp) &
                                                        !is.na(x = mid_size_order)))

summary(mass_iqr_temp_mid_size_big_adj_lmer)
confint(mass_iqr_temp_mid_size_big_adj_lmer) 

# Mass adjusted outliers removed
mass_iqr_temp_mid_size_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) *
                                            mid_size_order + scale(nestling_number) + 
                                            scale(days_summer) + (1|fsite) +
                                            (1|fnest_id), 
                                          data = subset(size_outliers_removed,
                                                        !is.na(x = mass_pre_obs) & 
                                                          !is.na(x = nest_iqr_temp) &
                                                          !is.na(x = mid_mid_size_order)))

## Check diagnostics for the full model
plot(mass_iqr_temp_mid_size_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_mid_size_adj_noout_lmer))
  qqline(resid(mass_iqr_temp_mid_size_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_mid_size_adj_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_temp_mid_size_adj_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_temp_mid_size_adj_noout_lmer, vars=c("Studentized"))

summary(mass_iqr_temp_mid_size_adj_noout_lmer)
confint(mass_iqr_temp_mid_size_adj_noout_lmer) 

# SMALL nestlings
mass_iqr_temp_mid_size_small_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) + 
                                                  scale(nestling_number) + 
                                                  scale(days_summer) +
                                                  (1|fsite), 
                                                data = subset(mid_size_outliers_removed,
                                                              mid_size_order == "min",
                                                              !is.na(x = mass_pre_obs) & 
                                                                !is.na(x = nest_iqr_temp) &
                                                                !is.na(x = mid_size_order)))

summary(mass_iqr_temp_mid_size_small_adj_noout_lmer)
confint(mass_iqr_temp_mid_size_small_adj_noout_lmer) 

# BIG nestlings
mass_iqr_temp_mid_size_big_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) + 
                                                scale(nestling_number) + 
                                                scale(days_summer) +
                                                (1|fsite) +
                                                (1|fnest_id), 
                                              data = subset(mid_size_outliers_removed,
                                                            mid_size_order == "other",
                                                            !is.na(x = mass_pre_obs) & 
                                                              !is.na(x = nest_iqr_temp) &
                                                              !is.na(x = mid_size_order)))

summary(mass_iqr_temp_mid_size_big_adj_noout_lmer)
confint(mass_iqr_temp_mid_size_big_adj_noout_lmer) 



########################## LATE SIZE #######################################

#### Minimum temp
### Mass unadjusted
mass_min_temp_size_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) *
                                   size_order + (1|fsite) +
                                   (1|fnest_id), 
                                 data = subset(late_nestling_parent_care,
                                               !is.na(x = size_order) & 
                                                 !is.na(x = nest_min_temp)&
                                                 !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_temp_size_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_size_lmer))
  qqline(resid(mass_min_temp_size_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_size_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_temp_size_lmer, vars=c("Cook"))
infIndexPlot(mass_min_temp_size_lmer, vars=c("Studentized"))

summary(mass_min_temp_size_lmer)
confint(mass_min_temp_size_lmer)


### Wing unadjusted
wing_min_temp_size_lmer <- lmer(rt_wing_length ~ scale(nest_min_temp) *
                                   size_order + (1|fsite) +
                                   (1|fnest_id), 
                                 data = subset(late_nestling_parent_care,
                                               !is.na(x = rt_wing_length) & 
                                                 !is.na(x = size_order)&
                                                 !is.na(x = nest_min_temp)))

## Check diagnostics for the full model
plot(wing_min_temp_size_lmer)
# Normal QQplot
{qqnorm(resid(wing_min_temp_size_lmer))
  qqline(resid(wing_min_temp_size_lmer))}
# Histogram of residuals
hist(resid(wing_min_temp_size_lmer))
# Checking for influential outliers
infIndexPlot(wing_min_temp_size_lmer, vars=c("Cook"))
infIndexPlot(wing_min_temp_size_lmer, vars=c("Studentized"))

summary(wing_min_temp_size_lmer)
confint(wing_min_temp_size_lmer)


### Mass adjusted
mass_min_temp_size_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) *
                                  size_order + scale(nestling_number) +
                                    scale(days_summer) + (1|fsite) +
                                  (1|fnest_id), 
                                data = subset(late_nestling_parent_care,
                                              !is.na(x = size_order) & 
                                                !is.na(x = nest_min_temp)&
                                                !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_temp_size_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_size_adj_lmer))
  qqline(resid(mass_min_temp_size_adj_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_size_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_temp_size_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_min_temp_size_adj_lmer, vars=c("Studentized"))

summary(mass_min_temp_size_adj_lmer)
confint(mass_min_temp_size_adj_lmer)


#### Maximum temp
### Mass unadjusted
mass_max_temp_size_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) *
                                  size_order + (1|fsite) +
                                  (1|fnest_id), 
                                data = subset(late_nestling_parent_care,
                                              !is.na(x = mass_pre_obs) & 
                                                !is.na(x = nest_max_temp) &
                                                !is.na(x = size_order)))

## Check diagnostics for the full model
plot(mass_max_temp_size_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_size_lmer))
  qqline(resid(mass_max_temp_size_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_size_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_temp_size_lmer, vars=c("Cook"))
infIndexPlot(mass_max_temp_size_lmer, vars=c("Studentized"))

summary(mass_max_temp_size_lmer)
confint(mass_max_temp_size_lmer) 

# SMALL nestlings
mass_max_temp_size_small_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                  (1|fsite), 
                                data = subset(late_nestling_parent_care,
                                              size_order == "min",
                                              !is.na(x = mass_pre_obs) & 
                                                !is.na(x = nest_max_temp) &
                                                !is.na(x = size_order)))

summary(mass_max_temp_size_small_lmer)
confint(mass_max_temp_size_small_lmer)

# BIG nestlings
mass_max_temp_size_big_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                        (1|fsite) +
                                        (1|fnest_id), 
                                      data = subset(late_nestling_parent_care,
                                                    size_order == "other",
                                                    !is.na(x = mass_pre_obs) & 
                                                      !is.na(x = nest_max_temp) &
                                                      !is.na(x = size_order)))

summary(mass_max_temp_size_big_lmer)
confint(mass_max_temp_size_big_lmer)


### Wing unadjusted
wing_max_temp_size_lmer <- lmer(rt_wing_length ~ scale(nest_max_temp) *
                                  size_order + (1|fsite) +
                                  (1|fnest_id), 
                                data = subset(late_nestling_parent_care,
                                              !is.na(x = rt_wing_length) & 
                                                !is.na(x = size_order)&
                                                !is.na(x = nest_max_temp)))

## Check diagnostics for the full model
plot(wing_max_temp_size_lmer)
# Normal QQplot
{qqnorm(resid(wing_max_temp_size_lmer))
  qqline(resid(wing_max_temp_size_lmer))}
# Histogram of residuals
hist(resid(wing_max_temp_size_lmer))
# Checking for influential outliers
infIndexPlot(wing_max_temp_size_lmer, vars=c("Cook"))
infIndexPlot(wing_max_temp_size_lmer, vars=c("Studentized"))

summary(wing_max_temp_size_lmer)
confint(wing_max_temp_size_lmer)

# SMALL nestlings 
wing_max_temp_size_small_lmer <- lmer(rt_wing_length ~ scale(nest_max_temp) +
                                  (1|fsite), 
                                data = subset(late_nestling_parent_care,
                                              size_order == "min",
                                              !is.na(x = rt_wing_length) & 
                                                !is.na(x = size_order)&
                                                !is.na(x = nest_max_temp)))

summary(wing_max_temp_size_small_lmer)
confint(wing_max_temp_size_small_lmer)

# LARGE nestlings 
wing_max_temp_size_big_lmer <- lmer(rt_wing_length ~ scale(nest_max_temp) +
                                        (1|fsite) +
                                        (1|fnest_id), 
                                      data = subset(late_nestling_parent_care,
                                                    size_order == "other",
                                                    !is.na(x = rt_wing_length) & 
                                                      !is.na(x = size_order)&
                                                      !is.na(x = nest_max_temp)))

summary(wing_max_temp_size_big_lmer)
confint(wing_max_temp_size_big_lmer)


### Mass adjusted
mass_max_temp_size_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) *
                                  size_order + scale(nestling_number) +
                                  scale(days_summer) + (1|fsite) +
                                  (1|fnest_id), 
                                data = subset(late_nestling_parent_care,
                                              !is.na(x = mass_pre_obs) & 
                                                !is.na(x = nest_max_temp) &
                                                !is.na(x = size_order)))

## Check diagnostics for the full model
plot(mass_max_temp_size_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_size_adj_lmer))
  qqline(resid(mass_max_temp_size_adj_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_size_adj_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_temp_size_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_max_temp_size_adj_lmer, vars=c("Studentized"))

summary(mass_max_temp_size_adj_lmer)
confint(mass_max_temp_size_adj_lmer) 

# SMALL nestlings
mass_max_temp_size_small_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                        scale(nestling_number) + scale(days_summer) +
                                          (1|fsite), 
                                      data = subset(late_nestling_parent_care,
                                                    size_order == "min",
                                                    !is.na(x = mass_pre_obs) & 
                                                      !is.na(x = nest_max_temp) &
                                                      !is.na(x = size_order)))

summary(mass_max_temp_size_small_adj_lmer)
confint(mass_max_temp_size_small_adj_lmer)

# BIG nestlings
mass_max_temp_size_big_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                          scale(nestling_number) + scale(days_summer) +
                                      (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  size_order == "other",
                                                  !is.na(x = mass_pre_obs) & 
                                                    !is.na(x = nest_max_temp) &
                                                    !is.na(x = size_order)))

summary(mass_max_temp_size_big_adj_lmer)
confint(mass_max_temp_size_big_adj_lmer)


size_outliers_removed <- late_nestling_parent_care[-c(44), ]


### Mass adjusted outliers removed
mass_max_temp_size_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) *
                                      size_order + scale(nestling_number) +
                                      scale(days_summer) + (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(size_outliers_removed,
                                                  !is.na(x = mass_pre_obs) & 
                                                    !is.na(x = nest_max_temp) &
                                                    !is.na(x = size_order)))

## Check diagnostics for the full model
plot(mass_max_temp_size_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_size_adj_noout_lmer))
  qqline(resid(mass_max_temp_size_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_size_adj_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_temp_size_adj_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_max_temp_size_adj_noout_lmer, vars=c("Studentized"))

summary(mass_max_temp_size_adj_noout_lmer)
confint(mass_max_temp_size_adj_noout_lmer) 

# SMALL nestlings
mass_max_temp_size_small_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                            scale(nestling_number) + scale(days_summer) +
                                            (1|fsite), 
                                          data = subset(size_outliers_removed,
                                                        size_order == "min",
                                                        !is.na(x = mass_pre_obs) & 
                                                          !is.na(x = nest_max_temp) &
                                                          !is.na(x = size_order)))

summary(mass_max_temp_size_small_adj_noout_lmer)
confint(mass_max_temp_size_small_adj_noout_lmer)

# BIG nestlings
mass_max_temp_size_big_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                          scale(nestling_number) + scale(days_summer) +
                                          (1|fsite) +
                                          (1|fnest_id), 
                                        data = subset(size_outliers_removed,
                                                      size_order == "other",
                                                      !is.na(x = mass_pre_obs) & 
                                                        !is.na(x = nest_max_temp) &
                                                        !is.na(x = size_order)))

summary(mass_max_temp_size_big_adj_noout_lmer)
confint(mass_max_temp_size_big_adj_noout_lmer)


#### IQR temp
### Mass unadjusted
mass_iqr_temp_size_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) *
                                  size_order + (1|fsite) +
                                  (1|fnest_id), 
                                data = subset(late_nestling_parent_care,
                                              !is.na(x = mass_pre_obs) & 
                                                !is.na(x = nest_iqr_temp) &
                                                !is.na(x = size_order)))

## Check diagnostics for the full model
plot(mass_iqr_temp_size_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_size_lmer))
  qqline(resid(mass_iqr_temp_size_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_size_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_temp_size_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_temp_size_lmer, vars=c("Studentized"))

summary(mass_iqr_temp_size_lmer)
confint(mass_iqr_temp_size_lmer) 

# SMALL nestlings
mass_iqr_temp_size_small_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) +
                                  (1|fsite), 
                                data = subset(late_nestling_parent_care,
                                              size_order == "min",
                                              !is.na(x = mass_pre_obs) & 
                                                !is.na(x = nest_iqr_temp) &
                                                !is.na(x = size_order)))

summary(mass_iqr_temp_size_small_lmer)
confint(mass_iqr_temp_size_small_lmer) 

# BIG nestlings
mass_iqr_temp_size_big_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) +
                                        (1|fsite) +
                                        (1|fnest_id), 
                                      data = subset(late_nestling_parent_care,
                                                    size_order == "other",
                                                    !is.na(x = mass_pre_obs) & 
                                                      !is.na(x = nest_iqr_temp) &
                                                      !is.na(x = size_order)))

summary(mass_iqr_temp_size_big_lmer)
confint(mass_iqr_temp_size_big_lmer) 

### Wing unadjusted
wing_iqr_temp_size_lmer <- lmer(rt_wing_length ~ scale(nest_iqr_temp) *
                                  size_order + (1|fsite) +
                                  (1|fnest_id), 
                                data = subset(late_nestling_parent_care,
                                              !is.na(x = rt_wing_length) & 
                                                !is.na(x = size_order)&
                                                !is.na(x = nest_iqr_temp)))

## Check diagnostics for the full model
plot(wing_iqr_temp_size_lmer)
# Normal QQplot
{qqnorm(resid(wing_iqr_temp_size_lmer))
  qqline(resid(wing_iqr_temp_size_lmer))}
# Histogram of residuals
hist(resid(wing_iqr_temp_size_lmer))
# Checking for influential outliers
infIndexPlot(wing_iqr_temp_size_lmer, vars=c("Cook"))
infIndexPlot(wing_iqr_temp_size_lmer, vars=c("Studentized"))

summary(wing_iqr_temp_size_lmer)
confint(wing_iqr_temp_size_lmer)

# SMALL nestlings 
wing_iqr_temp_size_small_lmer <- lmer(rt_wing_length ~ scale(nest_iqr_temp) +
                                  (1|fsite), 
                                data = subset(late_nestling_parent_care,
                                              size_order == "min",
                                              !is.na(x = rt_wing_length) & 
                                                !is.na(x = size_order)&
                                                !is.na(x = nest_iqr_temp)))

summary(wing_iqr_temp_size_small_lmer)
confint(wing_iqr_temp_size_small_lmer)

# BIG nestlings 
wing_iqr_temp_size_big_lmer <- lmer(rt_wing_length ~ scale(nest_iqr_temp) +
                                        (1|fsite) +
                                        (1|fnest_id), 
                                      data = subset(late_nestling_parent_care,
                                                    size_order == "other",
                                                    !is.na(x = rt_wing_length) & 
                                                      !is.na(x = size_order)&
                                                      !is.na(x = nest_iqr_temp)))

summary(wing_iqr_temp_size_big_lmer)
confint(wing_iqr_temp_size_big_lmer)


### Mass adjusted
mass_iqr_temp_size_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) *
                                  size_order + scale(nestling_number) + 
                                    scale(days_summer) + (1|fsite) +
                                  (1|fnest_id), 
                                data = subset(late_nestling_parent_care,
                                              !is.na(x = mass_pre_obs) & 
                                                !is.na(x = nest_iqr_temp) &
                                                !is.na(x = size_order)))

## Check diagnostics for the full model
plot(mass_iqr_temp_size_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_size_adj_lmer))
  qqline(resid(mass_iqr_temp_size_adj_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_size_adj_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_temp_size_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_temp_size_adj_lmer, vars=c("Studentized"))

summary(mass_iqr_temp_size_adj_lmer)
confint(mass_iqr_temp_size_adj_lmer) 

# SMALL nestlings
mass_iqr_temp_size_small_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) + 
                                            scale(nestling_number) + 
                                            scale(days_summer) +
                                        (1|fsite), 
                                      data = subset(late_nestling_parent_care,
                                                    size_order == "min",
                                                    !is.na(x = mass_pre_obs) & 
                                                      !is.na(x = nest_iqr_temp) &
                                                      !is.na(x = size_order)))

summary(mass_iqr_temp_size_small_adj_lmer)
confint(mass_iqr_temp_size_small_adj_lmer) 

# BIG nestlings
mass_iqr_temp_size_big_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) + 
                                      scale(nestling_number) + 
                                      scale(days_summer) +
                                      (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  size_order == "other",
                                                  !is.na(x = mass_pre_obs) & 
                                                    !is.na(x = nest_iqr_temp) &
                                                    !is.na(x = size_order)))

summary(mass_iqr_temp_size_big_adj_lmer)
confint(mass_iqr_temp_size_big_adj_lmer) 

### Mass adjusted outliers removed
mass_iqr_temp_size_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) *
                                      size_order + scale(nestling_number) + 
                                      scale(days_summer) + (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(size_outliers_removed,
                                                  !is.na(x = mass_pre_obs) & 
                                                    !is.na(x = nest_iqr_temp) &
                                                    !is.na(x = size_order)))

## Check diagnostics for the full model
plot(mass_iqr_temp_size_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_size_adj_noout_lmer))
  qqline(resid(mass_iqr_temp_size_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_size_adj_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_temp_size_adj_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_temp_size_adj_noout_lmer, vars=c("Studentized"))

summary(mass_iqr_temp_size_adj_noout_lmer)
confint(mass_iqr_temp_size_adj_noout_lmer) 

# SMALL nestlings
mass_iqr_temp_size_small_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) + 
                                            scale(nestling_number) + 
                                            scale(days_summer) +
                                            (1|fsite), 
                                          data = subset(size_outliers_removed,
                                                        size_order == "min",
                                                        !is.na(x = mass_pre_obs) & 
                                                          !is.na(x = nest_iqr_temp) &
                                                          !is.na(x = size_order)))

summary(mass_iqr_temp_size_small_adj_noout_lmer)
confint(mass_iqr_temp_size_small_adj_noout_lmer) 

# BIG nestlings
mass_iqr_temp_size_big_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) + 
                                          scale(nestling_number) + 
                                          scale(days_summer) +
                                          (1|fsite) +
                                          (1|fnest_id), 
                                        data = subset(size_outliers_removed,
                                                      size_order == "other",
                                                      !is.na(x = mass_pre_obs) & 
                                                        !is.na(x = nest_iqr_temp) &
                                                        !is.na(x = size_order)))

summary(mass_iqr_temp_size_big_adj_noout_lmer)
confint(mass_iqr_temp_size_big_adj_noout_lmer) 



###############################################################################
##############      3. Sensitive development periods      ##############
###############################################################################

#### Minimum temp
### Before thermo
mass_min_before_thermo_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp) + (1|fsite) +
                                  (1|fnest_id), 
                                data = subset(late_nestling_parent_care,
                                                !is.na(x = thermo_bef_min_temp)&
                                                !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_before_thermo_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_before_thermo_lmer))
  qqline(resid(mass_min_before_thermo_lmer))}
# Histogram of residuals
hist(resid(mass_min_before_thermo_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_before_thermo_lmer, vars=c("Cook"))
infIndexPlot(mass_min_before_thermo_lmer, vars=c("Studentized"))

summary(mass_min_before_thermo_lmer)
confint(mass_min_before_thermo_lmer)

### Before thermo adjusted
mass_min_before_thermo_adj_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp) + scale(nestling_number) +
                                      scale(days_summer) + (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  !is.na(x = thermo_bef_min_temp)&
                                                    !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_before_thermo_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_before_thermo_adj_lmer))
  qqline(resid(mass_min_before_thermo_adj_lmer))}
# Histogram of residuals
hist(resid(mass_min_before_thermo_adj_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_before_thermo_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_min_before_thermo_adj_lmer, vars=c("Studentized"))

summary(mass_min_before_thermo_adj_lmer)
confint(mass_min_before_thermo_adj_lmer)

### Before thermo adjusted no outliers
mass_min_before_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp) + scale(nestling_number) +
                                          scale(days_summer) + (1|fsite) +
                                          (1|fnest_id), 
                                        data = subset(mass_outliers_removed,
                                                      !is.na(x = thermo_bef_min_temp)&
                                                        !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_before_thermo_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_before_thermo_adj_noout_lmer))
  qqline(resid(mass_min_before_thermo_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_min_before_thermo_adj_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_before_thermo_adj_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_min_before_thermo_adj_noout_lmer, vars=c("Studentized"))

summary(mass_min_before_thermo_adj_noout_lmer)
confint(mass_min_before_thermo_adj_noout_lmer)

### After thermo
mass_min_after_thermo_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_min_temp) + (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  !is.na(x = thermo_aft_min_temp)&
                                                    !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_after_thermo_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_after_thermo_lmer))
  qqline(resid(mass_min_after_thermo_lmer))}
# Histogram of residuals
hist(resid(mass_min_after_thermo_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_after_thermo_lmer, vars=c("Cook"))
infIndexPlot(mass_min_after_thermo_lmer, vars=c("Studentized"))

summary(mass_min_after_thermo_lmer)
confint(mass_min_after_thermo_lmer)

### After thermo adjusted
mass_min_after_thermo_adj_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_min_temp) + scale(nestling_number) +
                                         scale(days_summer) + (1|fsite) +
                                         (1|fnest_id), 
                                       data = subset(late_nestling_parent_care,
                                                     !is.na(x = thermo_aft_min_temp)&
                                                       !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_after_thermo_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_after_thermo_adj_lmer))
  qqline(resid(mass_min_after_thermo_adj_lmer))}
# Histogram of residuals
hist(resid(mass_min_after_thermo_adj_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_after_thermo_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_min_after_thermo_adj_lmer, vars=c("Studentized"))

summary(mass_min_after_thermo_adj_lmer)
confint(mass_min_after_thermo_adj_lmer)


### After thermo adjusted outliers removed
mass_min_after_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_min_temp) + scale(nestling_number) +
                                               scale(days_summer) + (1|fsite) +
                                               (1|fnest_id), 
                                             data = subset(mass_outliers_removed,
                                                           !is.na(x = thermo_aft_min_temp)&
                                                             !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_after_thermo_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_after_thermo_adj_noout_lmer))
  qqline(resid(mass_min_after_thermo_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_min_after_thermo_adj_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_after_thermo_adj_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_min_after_thermo_adj_noout_lmer, vars=c("Studentized"))

summary(mass_min_after_thermo_adj_noout_lmer)
confint(mass_min_after_thermo_adj_noout_lmer)


### Before and after thermo
mass_min_both_thermo_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp) + 
                                    scale(thermo_aft_min_temp) + (1|fsite) +
                                    (1|fnest_id), 
                                  data = subset(late_nestling_parent_care,
                                                !is.na(x = thermo_bef_min_temp) &
                                                  !is.na(x = thermo_aft_min_temp)&
                                                  !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_both_thermo_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_both_thermo_lmer))
  qqline(resid(mass_min_both_thermo_lmer))}
# Histogram of residuals
hist(resid(mass_min_both_thermo_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_both_thermo_lmer, vars=c("Cook"))
infIndexPlot(mass_min_both_thermo_lmer, vars=c("Studentized"))

summary(mass_min_both_thermo_lmer)
confint(mass_min_both_thermo_lmer)
vif(mass_min_both_thermo_lmer)

### Before and after thermo adjusted
mass_min_both_thermo_adj_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp) + 
                                        scale(thermo_aft_min_temp) + scale(nestling_number) +
                                        scale(days_summer) + (1|fsite) +
                                        (1|fnest_id), 
                                      data = subset(late_nestling_parent_care,
                                                    !is.na(x = thermo_bef_min_temp) &
                                                      !is.na(x = thermo_aft_min_temp)&
                                                      !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_both_thermo_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_both_thermo_adj_lmer))
  qqline(resid(mass_min_both_thermo_adj_lmer))}
# Histogram of residuals
hist(resid(mass_min_both_thermo_adj_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_both_thermo_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_min_both_thermo_adj_lmer, vars=c("Studentized"))

summary(mass_min_both_thermo_adj_lmer)
confint(mass_min_both_thermo_adj_lmer)
vif(mass_min_both_thermo_adj_lmer)

### Before and after thermo outliers removed
thermo_outliers_removed <- late_nestling_parent_care[-c(43, 44), ]
mass_min_both_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp) + 
                                        scale(thermo_aft_min_temp) + scale(nestling_number) +
                                        scale(days_summer) + (1|fsite) +
                                        (1|fnest_id), 
                                      data = subset(thermo_outliers_removed,
                                                    !is.na(x = thermo_bef_min_temp) &
                                                      !is.na(x = thermo_aft_min_temp)&
                                                      !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_both_thermo_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_both_thermo_adj_noout_lmer))
  qqline(resid(mass_min_both_thermo_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_min_both_thermo_adj_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_both_thermo_adj_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_min_both_thermo_adj_noout_lmer, vars=c("Studentized"))

summary(mass_min_both_thermo_adj_noout_lmer)
confint(mass_min_both_thermo_adj_noout_lmer)
vif(mass_min_both_thermo_adj_noout_lmer)

#### Maximum temp
### Before thermo
mass_max_before_thermo_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp) + (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  !is.na(x = thermo_bef_max_temp)&
                                                    !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_before_thermo_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_before_thermo_lmer))
  qqline(resid(mass_max_before_thermo_lmer))}
# Histogram of residuals
hist(resid(mass_max_before_thermo_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_before_thermo_lmer, vars=c("Cook"))
infIndexPlot(mass_max_before_thermo_lmer, vars=c("Studentized"))

summary(mass_max_before_thermo_lmer)
confint(mass_max_before_thermo_lmer)

### Before thermo adjusted
mass_max_before_thermo_adj_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp) + scale(nestling_number) +
                                          scale(days_summer) + (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  !is.na(x = thermo_bef_max_temp)&
                                                    !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_before_thermo_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_before_thermo_adj_lmer))
  qqline(resid(mass_max_before_thermo_adj_lmer))}
# Histogram of residuals
hist(resid(mass_max_before_thermo_adj_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_before_thermo_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_max_before_thermo_adj_lmer, vars=c("Studentized"))

summary(mass_max_before_thermo_adj_lmer)
confint(mass_max_before_thermo_adj_lmer)

### Before thermo adjusted no outliers
mass_max_before_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp) + scale(nestling_number) +
                                          scale(days_summer) + (1|fsite) +
                                          (1|fnest_id), 
                                        data = subset(mass_outliers_removed,
                                                      !is.na(x = thermo_bef_max_temp)&
                                                        !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_before_thermo_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_before_thermo_adj_noout_lmer))
  qqline(resid(mass_max_before_thermo_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_max_before_thermo_adj_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_before_thermo_adj_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_max_before_thermo_adj_noout_lmer, vars=c("Studentized"))

summary(mass_max_before_thermo_adj_noout_lmer)
confint(mass_max_before_thermo_adj_noout_lmer)

### After thermo
mass_max_after_thermo_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_max_temp) + (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  !is.na(x = thermo_aft_max_temp)&
                                                    !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_after_thermo_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_after_thermo_lmer))
  qqline(resid(mass_max_after_thermo_lmer))}
# Histogram of residuals
hist(resid(mass_max_after_thermo_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_after_thermo_lmer, vars=c("Cook"))
infIndexPlot(mass_max_after_thermo_lmer, vars=c("Studentized"))

summary(mass_max_after_thermo_lmer)
confint(mass_max_after_thermo_lmer)

### After thermo adjusted
mass_max_after_thermo_adj_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_max_temp) + scale(nestling_number) +
                                         scale(days_summer) + (1|fsite) +
                                     (1|fnest_id), 
                                   data = subset(late_nestling_parent_care,
                                                 !is.na(x = thermo_aft_max_temp)&
                                                   !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_after_thermo_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_after_thermo_adj_lmer))
  qqline(resid(mass_max_after_thermo_adj_lmer))}
# Histogram of residuals
hist(resid(mass_max_after_thermo_adj_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_after_thermo_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_max_after_thermo_adj_lmer, vars=c("Studentized"))

summary(mass_max_after_thermo_adj_lmer)
confint(mass_max_after_thermo_adj_lmer)

### After thermo adjusted outliers removed
mass_max_after_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_max_temp) + scale(nestling_number) +
                                         scale(days_summer) + (1|fsite) +
                                         (1|fnest_id), 
                                       data = subset(thermo_outliers_removed,
                                                     !is.na(x = thermo_aft_max_temp)&
                                                       !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_after_thermo_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_after_thermo_adj_noout_lmer))
  qqline(resid(mass_max_after_thermo_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_max_after_thermo_adj_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_after_thermo_adj_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_max_after_thermo_adj_noout_lmer, vars=c("Studentized"))

summary(mass_max_after_thermo_adj_noout_lmer)
confint(mass_max_after_thermo_adj_noout_lmer)

### Before and after thermo
mass_max_both_thermo_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp) + 
                                    scale(thermo_aft_max_temp) + (1|fsite) +
                                    (1|fnest_id), 
                                  data = subset(late_nestling_parent_care,
                                                !is.na(x = thermo_bef_max_temp) &
                                                  !is.na(x = thermo_aft_max_temp)&
                                                  !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_both_thermo_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_both_thermo_lmer))
  qqline(resid(mass_max_both_thermo_lmer))}
# Histogram of residuals
hist(resid(mass_max_both_thermo_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_both_thermo_lmer, vars=c("Cook"))
infIndexPlot(mass_max_both_thermo_lmer, vars=c("Studentized"))

summary(mass_max_both_thermo_lmer)
confint(mass_max_both_thermo_lmer)
vif(mass_max_both_thermo_lmer)

### Before and after thermo adjusted
mass_max_both_thermo_adj_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp) + 
                                    scale(thermo_aft_max_temp) + scale(nestling_number) +
                                      scale(days_summer) +  (1|fsite) +
                                    (1|fnest_id), 
                                  data = subset(late_nestling_parent_care,
                                                !is.na(x = thermo_bef_max_temp) &
                                                  !is.na(x = thermo_aft_max_temp)&
                                                  !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_both_thermo_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_both_thermo_adj_lmer))
  qqline(resid(mass_max_both_thermo_adj_lmer))}
# Histogram of residuals
hist(resid(mass_max_both_thermo_adj_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_both_thermo_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_max_both_thermo_adj_lmer, vars=c("Studentized"))

summary(mass_max_both_thermo_adj_lmer)
confint(mass_max_both_thermo_adj_lmer)


### Before and after thermo adjusted no outliers
mass_max_both_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp) + 
                                        scale(thermo_aft_max_temp) + scale(nestling_number) +
                                        scale(days_summer) +  (1|fsite) +
                                        (1|fnest_id), 
                                      data = subset(mass_outliers_removed,
                                                    !is.na(x = thermo_bef_max_temp) &
                                                      !is.na(x = thermo_aft_max_temp)&
                                                      !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_both_thermo_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_both_thermo_adj_noout_lmer))
  qqline(resid(mass_max_both_thermo_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_max_both_thermo_adj_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_both_thermo_adj_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_max_both_thermo_adj_noout_lmer, vars=c("Studentized"))

summary(mass_max_both_thermo_adj_noout_lmer)
confint(mass_max_both_thermo_adj_noout_lmer)

#### IQR of temp
### Before thermo
mass_iqr_before_thermo_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp) + 
                                      (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  !is.na(x = thermo_bef_iqr_temp)&
                                                    !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_before_thermo_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_before_thermo_lmer))
  qqline(resid(mass_max_before_thermo_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_before_thermo_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_before_thermo_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_before_thermo_lmer, vars=c("Studentized"))

summary(mass_iqr_before_thermo_lmer)
confint(mass_iqr_before_thermo_lmer)

### Before thermo adjusted
mass_iqr_before_thermo_adj_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp) + scale(nestling_number) +
                                          scale(days_summer) +  
                                      (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  !is.na(x = thermo_bef_iqr_temp)&
                                                    !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_before_thermo_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_before_thermo_adj_lmer))
  qqline(resid(mass_iqr_before_thermo_adj_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_before_thermo_adj_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_before_thermo_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_before_thermo_adj_lmer, vars=c("Studentized"))

summary(mass_iqr_before_thermo_adj_lmer)
confint(mass_iqr_before_thermo_adj_lmer)


### After thermo
mass_iqr_after_thermo_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_iqr_temp) + 
                                      (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  !is.na(x = thermo_aft_iqr_temp)&
                                                    !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_after_thermo_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_after_thermo_lmer))
  qqline(resid(mass_iqr_after_thermo_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_after_thermo_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_after_thermo_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_after_thermo_lmer, vars=c("Studentized"))

summary(mass_iqr_after_thermo_lmer)
confint(mass_iqr_after_thermo_lmer)

### After thermo adjusted
mass_iqr_after_thermo_adj_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_iqr_temp) + scale(nestling_number) +
                                         scale(days_summer) + 
                                     (1|fsite) +
                                     (1|fnest_id), 
                                   data = subset(late_nestling_parent_care,
                                                 !is.na(x = thermo_aft_iqr_temp)&
                                                   !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_after_thermo_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_after_thermo_adj_lmer))
  qqline(resid(mass_iqr_after_thermo_adj_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_after_thermo_adj_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_after_thermo_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_after_thermo_adj_lmer, vars=c("Studentized"))

summary(mass_iqr_after_thermo_adj_lmer)
confint(mass_iqr_after_thermo_adj_lmer)

### After thermo adjusted no outliers
mass_iqr_after_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_iqr_temp) + scale(nestling_number) +
                                         scale(days_summer) + 
                                         (1|fsite) +
                                         (1|fnest_id), 
                                       data = subset(size_outliers_removed,
                                                     !is.na(x = thermo_aft_iqr_temp)&
                                                       !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_after_thermo_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_after_thermo_adj_noout_lmer))
  qqline(resid(mass_iqr_after_thermo_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_after_thermo_adj_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_after_thermo_adj_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_after_thermo_adj_noout_lmer, vars=c("Studentized"))

summary(mass_iqr_after_thermo_adj_noout_lmer)
confint(mass_iqr_after_thermo_adj_noout_lmer)

### Before and after thermo
mass_iqr_both_thermo_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp) + 
                                    scale(thermo_aft_iqr_temp) + (1|fsite) +
                                    (1|fnest_id), 
                                  data = subset(late_nestling_parent_care,
                                                !is.na(x = thermo_bef_iqr_temp) &
                                                  !is.na(x = thermo_aft_iqr_temp)&
                                                  !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_both_thermo_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_both_thermo_lmer))
  qqline(resid(mass_iqr_both_thermo_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_both_thermo_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_both_thermo_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_both_thermo_lmer, vars=c("Studentized"))

summary(mass_iqr_both_thermo_lmer)
confint(mass_iqr_both_thermo_lmer)
vif(mass_iqr_both_thermo_lmer)

### Before and after thermo adjusted
mass_iqr_both_thermo_adj_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp) + 
                                    scale(thermo_aft_iqr_temp) + scale(nestling_number) +
                                    scale(days_summer) +  (1|fsite) +
                                    (1|fnest_id), 
                                  data = subset(late_nestling_parent_care,
                                                !is.na(x = thermo_bef_iqr_temp) &
                                                  !is.na(x = thermo_aft_iqr_temp)&
                                                  !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_both_thermo_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_both_thermo_adj_lmer))
  qqline(resid(mass_iqr_both_thermo_adj_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_both_thermo_adj_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_both_thermo_adj_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_both_thermo_adj_lmer, vars=c("Studentized"))

summary(mass_iqr_both_thermo_adj_lmer)
confint(mass_iqr_both_thermo_adj_lmer)
vif(mass_iqr_both_thermo_adj_lmer)

### Before and after thermo adjusted
mass_iqr_both_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp) + 
                                        scale(thermo_aft_iqr_temp) + scale(nestling_number) +
                                        scale(days_summer) +  (1|fsite) +
                                        (1|fnest_id), 
                                      data = subset(size_outliers_removed,
                                                    !is.na(x = thermo_bef_iqr_temp) &
                                                      !is.na(x = thermo_aft_iqr_temp)&
                                                      !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_both_thermo_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_both_thermo_adj_noout_lmer))
  qqline(resid(mass_iqr_both_thermo_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_both_thermo_adj_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_both_thermo_adj_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_both_thermo_adj_noout_lmer, vars=c("Studentized"))

summary(mass_iqr_both_thermo_adj_noout_lmer)
confint(mass_iqr_both_thermo_adj_noout_lmer)
vif(mass_iqr_both_thermo_adj_noout_lmer)
