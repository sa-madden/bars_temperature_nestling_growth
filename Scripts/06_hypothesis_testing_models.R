####### Purpose: run models to test hypotheses for barn swallow
####### nest mircoclimate and nestling growth dataset
####### By: XXX
####### Created: 12/19/2022
####### Last modified: 6/10/2025

# Code Blocks
# 1: Configure work space
# 2: Load data
# 3: Growth, temp, and parental care (question 3)
# 4: Growth, temp, and relative size (question 2)
# 5: Sensitive development periods (question 1)
# 6: Model visualizations

### Sections 3 to 5 include the code to run models, check model diagnostics,
### and create results tables for export for each of three questions. 
### For each question, we ran unadjusted and adjusted (including covariates) models
### for each of three temperature variables (maximum, minimum, interquartile range)

### Section 6 includes code to create plots of model predictions and raw data.


###############################################################################
##############                 Configure work space              ##############
###############################################################################

### Global options
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
library('ggpubr')
library('ggeffects')
library('ggnewscale')
library('gt')


## Modelling Packages
library('lme4')
library ('emmeans')
library('MuMIn')
library('car')
library('boot')

# load pbkrtest and lmertest (emmeans dependency)
library('pbkrtest')
library('lmerTest')

## Broom packages
library('broom')
library('broom.mixed')

## Model checking packages
library('DHARMa')


### Get Version and Session Info
R.Version()
sessionInfo()

###############################################################################
##############                      Load Data                    ##############
###############################################################################  

### Load RData
## Load RData tidy barn swallow data
load('Data/Tidy/tidy_parent_nestl_weather_data_10-4_with_BLUPs.RData')

# Make site a factor
late_nestling_parent_care$fsite <- as.factor(late_nestling_parent_care$site)
late_nestling_parent_care$fnest_id <- as.factor(late_nestling_parent_care$nest_id)

## Create three categories for BLUPs to allow stratification 
late_nestling_parent_care <- late_nestling_parent_care %>%
  mutate(feeding_expontd_blups_strat =  as.integer(Hmisc::cut2(feeding_expontd_blups, g=3)))

### Remove outliers in minimum temp x parental feeding analyses (identified based on 
### diagnostic plots later in script)
mass_outliers_removed <- late_nestling_parent_care[-c(4, 40), ]

### Remove outliers in maximum temp x parental feeding analyses (identified based on 
### diagnostic plots later in script)
mass_max_outliers_removed <- late_nestling_parent_care[-c(4, 39, 40), ]

### Remove outliers in sensitive periods before model - min temp (identified based on 
### diagnostic plots later in script)
before_thermo_outliers_removed <- late_nestling_parent_care[-c(4, 39, 40), ]


###############################################################################
##############             Growth, temp, and parental care       ##############
###############################################################################

######################### Minimum temperature #################################

### Inteaction model - unadjusted
mass_min_temp_blups_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) *
                                   scale(feeding_expontd_blups) +
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

# Calculate R squared 
r.squaredGLMM(mass_min_temp_blups_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_blups_lmer<- bootMer(x = mass_min_temp_blups_lmer,
                                FUN = fixef, nsim = 2000,
                                seed = 632760,
                                use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_blups_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_blups_lmer <- boot.ci(boot_mass_min_temp_blups_lmer,
                                 type = c('perc', 'norm', 'basic'),
                                 index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_blups_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_min_temp_blups_lmer_2 <- boot.ci(boot_mass_min_temp_blups_lmer,
                                          type = c('perc', 'norm', 'basic'),
                                          index = 3) # CI for 2nd betas
print(bt_ci_mass_min_temp_blups_lmer_2)

# use 'boot' package to generate 95% CI for 3rd beta
bt_ci_mass_min_temp_blups_lmer_3 <- boot.ci(boot_mass_min_temp_blups_lmer,
                                          type = c('perc', 'norm', 'basic'),
                                          index = 4) # CI for 3rd betas
print(bt_ci_mass_min_temp_blups_lmer_3)



### Stratified models - unadjusted
## LOW parental care model - unadjusted
mass_min_temp_blups_low_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) + 
                                           (1|fnest_id), 
                                         data = subset(late_nestling_parent_care,
                                                       feeding_expontd_blups_strat == 1,
                                                       !is.na(x = mass_pre_obs) & 
                                                         !is.na(x = nest_min_temp) &
                                                         !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_min_temp_blups_low_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_blups_low_lmer))
  qqline(resid(mass_min_temp_blups_low_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_blups_low_lmer))


summary(mass_min_temp_blups_low_lmer)
confint(mass_min_temp_blups_low_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_blups_low_lmer <- bootMer(x = mass_min_temp_blups_low_lmer,
                                                 FUN = fixef, nsim = 2000,
                                                 seed = 632760,
                                                 use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_blups_low_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_blups_low_lmer <- boot.ci(boot_mass_min_temp_blups_low_lmer,
                                                  type = c('perc', 'norm', 'basic'),
                                                  index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_blups_low_lmer)


## MED parental care- unadjusted
mass_min_temp_blups_med_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) + 
                                           (1|fnest_id), 
                                         data = subset(late_nestling_parent_care,
                                                       feeding_expontd_blups_strat == 2,
                                                       !is.na(x = mass_pre_obs) & 
                                                         !is.na(x = nest_min_temp) &
                                                         !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_min_temp_blups_med_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_blups_med_lmer))
  qqline(resid(mass_min_temp_blups_med_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_blups_med_lmer))


summary(mass_min_temp_blups_med_lmer)
confint(mass_min_temp_blups_med_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_blups_med_lmer <- bootMer(x = mass_min_temp_blups_med_lmer,
                                                 FUN = fixef, nsim = 2000,
                                                 seed = 632760,
                                                 use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_blups_med_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_blups_med_lmer <- boot.ci(boot_mass_min_temp_blups_med_lmer,
                                                  type = c('perc', 'norm', 'basic'),
                                                  index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_blups_med_lmer)


## HIGH parental care - unadjusted
mass_min_temp_blups_high_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) +
                                            (1|fnest_id), 
                                          data = subset(late_nestling_parent_care,
                                                        feeding_expontd_blups_strat == 3,
                                                        !is.na(x = mass_pre_obs) & 
                                                          !is.na(x = nest_min_temp) &
                                                          !is.na(x = feeding_expontd_blups)))


## Check diagnostics for the full model
plot(mass_min_temp_blups_high_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_blups_high_lmer))
  qqline(resid(mass_min_temp_blups_high_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_blups_high_lmer))


summary(mass_min_temp_blups_high_lmer)
confint(mass_min_temp_blups_high_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_blups_high_lmer <- bootMer(x = mass_min_temp_blups_high_lmer,
                                                  FUN = fixef, nsim = 2000,
                                                  seed = 632760,
                                                  use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_blups_high_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_blups_high_lmer <- boot.ci(boot_mass_min_temp_blups_high_lmer,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_blups_high_lmer)



### Interaction model - adjusted
mass_min_temp_blups_adj_lmer <- lmerTest::lmer(mass_pre_obs ~ scale(nest_min_temp) *
                                   scale(feeding_expontd_blups) + scale(nestling_number) + 
                                     scale(days_summer) +
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

# Calculate R squared 
r.squaredGLMM(mass_min_temp_blups_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_blups_adj_lmer<- bootMer(x = mass_min_temp_blups_adj_lmer,
                                        FUN = fixef, nsim = 2000,
                                        seed = 632760,
                                        use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_blups_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_blups_adj_lmer <- boot.ci(boot_mass_min_temp_blups_adj_lmer,
                                          type = c('perc', 'norm', 'basic'),
                                          index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_blups_adj_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_min_temp_blups_adj_lmer_2 <- boot.ci(boot_mass_min_temp_blups_adj_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 3) # CI for 2nd betas
print(bt_ci_mass_min_temp_blups_adj_lmer_2)

# use 'boot' package to generate 95% CI for 5th beta
bt_ci_mass_min_temp_blups_adj_lmer_5 <- boot.ci(boot_mass_min_temp_blups_adj_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 6) # CI for 3rd betas
print(bt_ci_mass_min_temp_blups_adj_lmer_5)


### Straified models - adjusted
## LOW parental care model - adjusted
mass_min_temp_blups_low_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                           (1|fnest_id), 
                                         data = subset(late_nestling_parent_care,
                                                       feeding_expontd_blups_strat == 1,
                                                       !is.na(x = mass_pre_obs) & 
                                                         !is.na(x = nest_min_temp) &
                                                         !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_min_temp_blups_low_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_blups_low_adj_lmer))
  qqline(resid(mass_min_temp_blups_low_adj_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_blups_low_adj_lmer))


summary(mass_min_temp_blups_low_adj_lmer)
confint(mass_min_temp_blups_low_adj_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_blups_low_adj_lmer <- bootMer(x = mass_min_temp_blups_low_adj_lmer,
                                            FUN = fixef, nsim = 2000,
                                            seed = 632760,
                                            use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_blups_low_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_blups_low_adj_lmer <- boot.ci(boot_mass_min_temp_blups_low_adj_lmer,
                                              type = c('perc', 'norm', 'basic'),
                                              index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_blups_low_adj_lmer)


## MED parental care - adjusted
mass_min_temp_blups_med_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                       (1|fnest_id), 
                                     data = subset(late_nestling_parent_care,
                                                   feeding_expontd_blups_strat == 2,
                                                   !is.na(x = mass_pre_obs) & 
                                                     !is.na(x = nest_min_temp) &
                                                     !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_min_temp_blups_med_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_blups_med_adj_lmer))
  qqline(resid(mass_min_temp_blups_med_adj_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_blups_med_adj_lmer))


summary(mass_min_temp_blups_med_adj_lmer)
confint(mass_min_temp_blups_med_adj_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_blups_med_adj_lmer <- bootMer(x = mass_min_temp_blups_med_adj_lmer,
                                                 FUN = fixef, nsim = 2000,
                                                 seed = 632760,
                                                 use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_blups_med_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_blups_med_adj_lmer <- boot.ci(boot_mass_min_temp_blups_med_adj_lmer,
                                                  type = c('perc', 'norm', 'basic'),
                                                  index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_blups_med_adj_lmer)


## HIGH parental care - adjusted 
mass_min_temp_blups_high_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) +
                                            scale(nestling_number) + scale(days_summer) + 
                                        (1|fnest_id), 
                                      data = subset(late_nestling_parent_care,
                                                    feeding_expontd_blups_strat == 3,
                                                    !is.na(x = mass_pre_obs) & 
                                                      !is.na(x = nest_min_temp) &
                                                      !is.na(x = feeding_expontd_blups)))


## Check diagnostics for the full model
plot(mass_min_temp_blups_high_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_blups_high_adj_lmer))
  qqline(resid(mass_min_temp_blups_high_adj_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_blups_high_adj_lmer))


summary(mass_min_temp_blups_high_adj_lmer)
confint(mass_min_temp_blups_high_adj_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_blups_high_adj_lmer <- bootMer(x = mass_min_temp_blups_high_adj_lmer,
                                                 FUN = fixef, nsim = 2000,
                                                 seed = 632760,
                                                 use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_blups_high_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_blups_high_adj_lmer <- boot.ci(boot_mass_min_temp_blups_high_adj_lmer,
                                                  type = c('perc', 'norm', 'basic'),
                                                  index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_blups_high_adj_lmer)



### Interaction model - unadjusted, outliers removed
mass_min_temp_blups_noout_lmer <- lmerTest::lmer(mass_pre_obs ~ scale(nest_min_temp) *
                                   scale(feeding_expontd_blups) +
                                   (1|fnest_id), 
                                 data = subset(mass_outliers_removed,
                                               !is.na(x = mass_pre_obs) & 
                                                 !is.na(x = nest_min_temp) &
                                                 !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_min_temp_blups_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_blups_noout_lmer))
  qqline(resid(mass_min_temp_blups_noout_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_blups_noout_lmer))

summary(mass_min_temp_blups_noout_lmer)
confint(mass_min_temp_blups_noout_lmer)  

# Calculate R squared 
r.squaredGLMM(mass_min_temp_blups_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_blups_noout_lmer <- bootMer(x = mass_min_temp_blups_noout_lmer,
                                        FUN = fixef, nsim = 2000,
                                        seed = 632760,
                                        use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_blups_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_blups_noout_lmer <- boot.ci(boot_mass_min_temp_blups_noout_lmer,
                                          type = c('perc', 'norm', 'basic'),
                                          index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_blups_noout_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_min_temp_blups_noout_lmer_2 <- boot.ci(boot_mass_min_temp_blups_noout_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 3) # CI for 2nd betas
print(bt_ci_mass_min_temp_blups_noout_lmer_2)

# use 'boot' package to generate 95% CI for 3rd beta
bt_ci_mass_min_temp_blups_noout_lmer_3 <- boot.ci(boot_mass_min_temp_blups_noout_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 4) # CI for 3rd betas
print(bt_ci_mass_min_temp_blups_noout_lmer_3)


### Interaction model - adjusted, outliers removed
mass_min_temp_blups_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) *
                                       scale(feeding_expontd_blups) + scale(nestling_number) + 
                                       scale(days_summer) +
                                       (1|fnest_id), 
                                     data = subset(mass_outliers_removed,
                                                   !is.na(x = mass_pre_obs) & 
                                                     !is.na(x = nest_min_temp) &
                                                     !is.na(x = feeding_expontd_blups)))

# Check diagnostics for the full model
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

# Calculate R squared 
r.squaredGLMM(mass_min_temp_blups_adj_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_blups_adj_noout_lmer <- bootMer(x = mass_min_temp_blups_adj_noout_lmer,
                                        FUN = fixef, nsim = 2000,
                                        seed = 632760,
                                        use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_blups_adj_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_blups_adj_noout_lmer <- boot.ci(boot_mass_min_temp_blups_adj_noout_lmer,
                                          type = c('perc', 'norm', 'basic'),
                                          index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_blups_adj_noout_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_min_temp_blups_adj_noout_lmer_2 <- boot.ci(boot_mass_min_temp_blups_adj_noout_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 3) # CI for 2nd betas
print(bt_ci_mass_min_temp_blups_adj_noout_lmer_2)

# use 'boot' package to generate 95% CI for 5th beta
bt_ci_mass_min_temp_blups_adj_noout_lmer_5 <- boot.ci(boot_mass_min_temp_blups_adj_noout_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 6) # CI for 3rd betas
print(bt_ci_mass_min_temp_blups_adj_noout_lmer_5)


### Stratified models - adjusted, outliers removed
## LOW parental care model - adjusted, outliers removed
mass_min_temp_blups_low_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                           (1|fnest_id), 
                                         data = subset(mass_outliers_removed,
                                                       feeding_expontd_blups_strat == 1,
                                                       !is.na(x = mass_pre_obs) & 
                                                         !is.na(x = nest_min_temp) &
                                                         !is.na(x = feeding_expontd_blups)))

# Check diagnostics for the full model
plot(mass_min_temp_blups_low_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_blups_low_adj_noout_lmer))
  qqline(resid(mass_min_temp_blups_low_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_blups_low_adj_noout_lmer))

summary(mass_min_temp_blups_low_adj_noout_lmer)
confint(mass_min_temp_blups_low_adj_noout_lmer) 

# Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_blups_low_adj_noout_lmer <- bootMer(x = mass_min_temp_blups_low_adj_noout_lmer,
                                                   FUN = fixef, nsim = 2000,
                                                   seed = 632760,
                                                   use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_blups_low_adj_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_blups_low_adj_noout_lmer <- boot.ci(boot_mass_min_temp_blups_low_adj_noout_lmer,
                                                    type = c('perc', 'norm', 'basic'),
                                                    index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_blups_low_adj_noout_lmer)


## MED parental care - adjusted, outliers removed
mass_min_temp_blups_med_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                           (1|fnest_id), 
                                         data = subset(mass_outliers_removed,
                                                       feeding_expontd_blups_strat == 2,
                                                       !is.na(x = mass_pre_obs) & 
                                                         !is.na(x = nest_min_temp) &
                                                         !is.na(x = feeding_expontd_blups)))

# Check diagnostics for the full model
plot(mass_min_temp_blups_med_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_blups_med_adj_noout_lmer))
  qqline(resid(mass_min_temp_blups_med_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_blups_med_adj_noout_lmer))

summary(mass_min_temp_blups_med_adj_noout_lmer)
confint(mass_min_temp_blups_med_adj_noout_lmer) 

# Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_blups_med_adj_noout_lmer <- bootMer(x = mass_min_temp_blups_med_adj_noout_lmer,
                                                       FUN = fixef, nsim = 2000,
                                                       seed = 632760,
                                                       use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_blups_med_adj_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_blups_med_adj_noout_lmer <- boot.ci(boot_mass_min_temp_blups_med_adj_noout_lmer,
                                                        type = c('perc', 'norm', 'basic'),
                                                        index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_blups_med_adj_noout_lmer)

## HIGH parental care - adjusted, outliers removed
mass_min_temp_blups_high_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) +
                                            scale(nestling_number) + scale(days_summer) + 
                                            (1|fnest_id), 
                                          data = subset(mass_outliers_removed,
                                                        feeding_expontd_blups_strat == 3,
                                                        !is.na(x = mass_pre_obs) & 
                                                          !is.na(x = nest_min_temp) &
                                                          !is.na(x = feeding_expontd_blups)))

# Check diagnostics for the full model
plot(mass_min_temp_blups_high_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_blups_high_adj_noout_lmer))
  qqline(resid(mass_min_temp_blups_high_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_blups_high_adj_noout_lmer))

summary(mass_min_temp_blups_high_adj_noout_lmer)
confint(mass_min_temp_blups_high_adj_noout_lmer) 

# Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_blups_high_adj_noout_lmer <- bootMer(x = mass_min_temp_blups_high_adj_noout_lmer,
                                                       FUN = fixef, nsim = 2000,
                                                       seed = 632760,
                                                       use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_blups_high_adj_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_blups_high_adj_noout_lmer <- boot.ci(boot_mass_min_temp_blups_high_adj_noout_lmer,
                                                        type = c('perc', 'norm', 'basic'),
                                                        index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_blups_high_adj_noout_lmer)


############################## Maximum temperature ############################

## Interaction model - unadjusted 
mass_max_temp_blups_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) *
                                   scale(feeding_expontd_blups) +
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

# Calculate R squared 
r.squaredGLMM(mass_max_temp_blups_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_blups_lmer <- bootMer(x = mass_max_temp_blups_lmer,
                                                   FUN = fixef, nsim = 2000,
                                                   seed = 632760,
                                                   use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_blups_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_blups_lmer <- boot.ci(boot_mass_max_temp_blups_lmer,
                                                    type = c('perc', 'norm', 'basic'),
                                                    index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_blups_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_max_temp_blups_lmer_2 <- boot.ci(boot_mass_max_temp_blups_lmer,
                                                      type = c('perc', 'norm', 'basic'),
                                                      index = 3) # CI for 2nd betas
print(bt_ci_mass_max_temp_blups_lmer_2)

# use 'boot' package to generate 95% CI for 3rd beta
bt_ci_mass_max_temp_blups_lmer_3 <- boot.ci(boot_mass_max_temp_blups_lmer,
                                                      type = c('perc', 'norm', 'basic'),
                                                      index = 4) # CI for 3rd betas
print(bt_ci_mass_max_temp_blups_lmer_3)


### Stratified models - unadjusted
## LOW parental care model - unadjusted
mass_max_temp_blups_low_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) + 
                                       (1|fnest_id), 
                                     data = subset(late_nestling_parent_care,
                                                   feeding_expontd_blups_strat == 1,
                                                   !is.na(x = mass_pre_obs) & 
                                                     !is.na(x = nest_max_temp) &
                                                     !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_max_temp_blups_low_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_blups_low_lmer))
  qqline(resid(mass_max_temp_blups_low_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_blups_low_lmer))


summary(mass_max_temp_blups_low_lmer)
confint(mass_max_temp_blups_low_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_blups_low_lmer <- bootMer(x = mass_max_temp_blups_low_lmer,
                                             FUN = fixef, nsim = 2000,
                                             seed = 632760,
                                             use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_blups_low_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_blups_low_lmer <- boot.ci(boot_mass_max_temp_blups_low_lmer,
                                              type = c('perc', 'norm', 'basic'),
                                              index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_blups_low_lmer)


## MED parental care - unadjusted
mass_max_temp_blups_med_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) + 
                                       (1|fnest_id), 
                                     data = subset(late_nestling_parent_care,
                                                   feeding_expontd_blups_strat == 2,
                                                   !is.na(x = mass_pre_obs) & 
                                                     !is.na(x = nest_max_temp) &
                                                     !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_max_temp_blups_med_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_blups_med_lmer))
  qqline(resid(mass_max_temp_blups_med_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_blups_med_lmer))


summary(mass_max_temp_blups_med_lmer)
confint(mass_max_temp_blups_med_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_blups_med_lmer <- bootMer(x = mass_max_temp_blups_med_lmer,
                                             FUN = fixef, nsim = 2000,
                                             seed = 632760,
                                             use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_blups_med_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_blups_med_lmer <- boot.ci(boot_mass_max_temp_blups_med_lmer,
                                              type = c('perc', 'norm', 'basic'),
                                              index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_blups_med_lmer)


## HIGH parental care - unadjusted 
mass_max_temp_blups_high_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                        (1|fnest_id), 
                                      data = subset(late_nestling_parent_care,
                                                    feeding_expontd_blups_strat == 3,
                                                    !is.na(x = mass_pre_obs) & 
                                                      !is.na(x = nest_max_temp) &
                                                      !is.na(x = feeding_expontd_blups)))


## Check diagnostics for the full model
plot(mass_max_temp_blups_high_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_blups_high_lmer))
  qqline(resid(mass_max_temp_blups_high_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_blups_high_lmer))


summary(mass_max_temp_blups_high_lmer)
confint(mass_max_temp_blups_high_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_blups_high_lmer <- bootMer(x = mass_max_temp_blups_high_lmer,
                                              FUN = fixef, nsim = 2000,
                                              seed = 632760,
                                              use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_blups_high_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_blups_high_lmer <- boot.ci(boot_mass_max_temp_blups_high_lmer,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_blups_high_lmer)


### Interaction model - adjusted
mass_max_temp_blups_adj_lmer <- lmerTest::lmer(mass_pre_obs ~ scale(nest_max_temp) *
                                   scale(feeding_expontd_blups) + 
                                     scale(nestling_number) + scale(days_summer) +
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

# Calculate R squared 
r.squaredGLMM(mass_max_temp_blups_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_blups_adj_lmer <- bootMer(x = mass_max_temp_blups_adj_lmer,
                                         FUN = fixef, nsim = 2000,
                                         seed = 632760,
                                         use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_blups_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_blups_adj_lmer <- boot.ci(boot_mass_max_temp_blups_adj_lmer,
                                          type = c('perc', 'norm', 'basic'),
                                          index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_blups_adj_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_max_temp_blups_adj_lmer_2 <- boot.ci(boot_mass_max_temp_blups_adj_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 3) # CI for 2nd betas
print(bt_ci_mass_max_temp_blups_adj_lmer_2)

# use 'boot' package to generate 95% CI for 5th beta
bt_ci_mass_max_temp_blups_adj_lmer_5 <- boot.ci(boot_mass_max_temp_blups_adj_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 6) # CI for 3rd betas
print(bt_ci_mass_max_temp_blups_adj_lmer_5)


### Stratified models - adjusted
## LOW parental care model - adjusted
mass_max_temp_blups_low_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                           (1|fnest_id), 
                                         data = subset(late_nestling_parent_care,
                                                       feeding_expontd_blups_strat == 1,
                                                       !is.na(x = mass_pre_obs) & 
                                                         !is.na(x = nest_max_temp) &
                                                         !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_max_temp_blups_low_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_blups_low_adj_lmer))
  qqline(resid(mass_max_temp_blups_low_adj_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_blups_low_adj_lmer))


summary(mass_max_temp_blups_low_adj_lmer)
confint(mass_max_temp_blups_low_adj_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_blups_low_adj_lmer <- bootMer(x = mass_max_temp_blups_low_adj_lmer,
                                                 FUN = fixef, nsim = 2000,
                                                 seed = 632760,
                                                 use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_blups_low_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_blups_low_adj_lmer <- boot.ci(boot_mass_max_temp_blups_low_adj_lmer,
                                                  type = c('perc', 'norm', 'basic'),
                                                  index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_blups_low_adj_lmer)


## MED parental care - adjusted
mass_max_temp_blups_med_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                           (1|fnest_id), 
                                         data = subset(late_nestling_parent_care,
                                                       feeding_expontd_blups_strat == 2,
                                                       !is.na(x = mass_pre_obs) & 
                                                         !is.na(x = nest_max_temp) &
                                                         !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_max_temp_blups_med_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_blups_med_adj_lmer))
  qqline(resid(mass_max_temp_blups_med_adj_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_blups_med_adj_lmer))


summary(mass_max_temp_blups_med_adj_lmer)
confint(mass_max_temp_blups_med_adj_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_blups_med_adj_lmer <- bootMer(x = mass_max_temp_blups_med_adj_lmer,
                                                 FUN = fixef, nsim = 2000,
                                                 seed = 632760,
                                                 use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_blups_med_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_blups_med_adj_lmer <- boot.ci(boot_mass_max_temp_blups_med_adj_lmer,
                                                  type = c('perc', 'norm', 'basic'),
                                                  index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_blups_med_adj_lmer)


## HIGH parental care - adjusted 
mass_max_temp_blups_high_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                            scale(nestling_number) + scale(days_summer) + 
                                            (1|fnest_id), 
                                          data = subset(late_nestling_parent_care,
                                                        feeding_expontd_blups_strat == 3,
                                                        !is.na(x = mass_pre_obs) & 
                                                          !is.na(x = nest_max_temp) &
                                                          !is.na(x = feeding_expontd_blups)))


## Check diagnostics for the full model
plot(mass_max_temp_blups_high_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_blups_high_adj_lmer))
  qqline(resid(mass_max_temp_blups_high_adj_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_blups_high_adj_lmer))

summary(mass_max_temp_blups_high_adj_lmer)
confint(mass_max_temp_blups_high_adj_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_blups_high_adj_lmer <- bootMer(x = mass_max_temp_blups_high_adj_lmer,
                                                  FUN = fixef, nsim = 2000,
                                                  seed = 632760,
                                                  use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_blups_high_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_blups_high_adj_lmer <- boot.ci(boot_mass_max_temp_blups_high_adj_lmer,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_blups_high_adj_lmer)



## Interaction model - unadjusted, outliers removed
mass_max_temp_blups_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) *
                                   scale(feeding_expontd_blups) +
                                   (1|fnest_id), 
                                 data = subset(mass_max_outliers_removed,
                                               !is.na(x = mass_pre_obs) & 
                                                 !is.na(x = nest_max_temp) &
                                                 !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_max_temp_blups_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_blups_noout_lmer))
  qqline(resid(mass_max_temp_blups_noout_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_blups_noout_lmer))

summary(mass_max_temp_blups_noout_lmer)
confint(mass_max_temp_blups_noout_lmer) 

# Calculate R squared 
r.squaredGLMM(mass_max_temp_blups_noout_lmer)


## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_blups_noout_lmer <- bootMer(x = mass_max_temp_blups_noout_lmer,
                                         FUN = fixef, nsim = 2000,
                                         seed = 632760,
                                         use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_blups_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_blups_noout_lmer <- boot.ci(boot_mass_max_temp_blups_noout_lmer,
                                          type = c('perc', 'norm', 'basic'),
                                          index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_blups_noout_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_max_temp_blups_noout_lmer_2 <- boot.ci(boot_mass_max_temp_blups_noout_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 3) # CI for 2nd betas
print(bt_ci_mass_max_temp_blups_lmer_2)

# use 'boot' package to generate 95% CI for 3rd beta
bt_ci_mass_max_temp_blups_noout_lmer_3 <- boot.ci(boot_mass_max_temp_blups_noout_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 4) # CI for 3rd betas
print(bt_ci_mass_max_temp_blups_noout_lmer_3)


### Interaction model - adjusted, outliers removed
mass_max_temp_blups_adj_noout_lmer <- lmerTest::lmer(mass_pre_obs ~ scale(nest_max_temp) *
                                       scale(feeding_expontd_blups) + 
                                       scale(nestling_number) + scale(days_summer) +
                                       (1|fnest_id), 
                                     data = subset(mass_max_outliers_removed,
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

# Calculate R squared 
r.squaredGLMM(mass_max_temp_blups_adj_noout_lmer)


## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_blups_adj_noout_lmer <- bootMer(x = mass_max_temp_blups_adj_noout_lmer,
                                             FUN = fixef, nsim = 2000,
                                             seed = 632760,
                                             use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_blups_adj_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_blups_adj_noout_lmer <- boot.ci(boot_mass_max_temp_blups_adj_noout_lmer,
                                              type = c('perc', 'norm', 'basic'),
                                              index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_blups_adj_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_max_temp_blups_adj_noout_lmer_2 <- boot.ci(boot_mass_max_temp_blups_adj_noout_lmer,
                                                type = c('perc', 'norm', 'basic'),
                                                index = 3) # CI for 2nd betas
print(bt_ci_mass_max_temp_blups_adj_noout_lmer_2)

# use 'boot' package to generate 95% CI for 5th beta
bt_ci_mass_max_temp_blups_adj_noout_lmer_5 <- boot.ci(boot_mass_max_temp_blups_adj_noout_lmer,
                                                type = c('perc', 'norm', 'basic'),
                                                index = 6) # CI for 3rd betas
print(bt_ci_mass_max_temp_blups_adj_noout_lmer_5)


### Stratified models - adjusted, outliers removed
## LOW parental care model - adjusted, outliers removed
mass_max_temp_blups_low_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                           (1|fnest_id), 
                                         data = subset(mass_max_outliers_removed,
                                                       feeding_expontd_blups_strat == 1,
                                                       !is.na(x = mass_pre_obs) & 
                                                         !is.na(x = nest_max_temp) &
                                                         !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_max_temp_blups_low_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_blups_low_adj_noout_lmer))
  qqline(resid(mass_max_temp_blups_low_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_blups_low_adj_noout_lmer))


summary(mass_max_temp_blups_low_adj_noout_lmer)
confint(mass_max_temp_blups_low_adj_noout_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_blups_low_adj_noout_lmer <- bootMer(x = mass_max_temp_blups_low_adj_noout_lmer,
                                                 FUN = fixef, nsim = 2000,
                                                 seed = 632760,
                                                 use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_blups_low_adj_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_blups_low_adj_noout_lmer <- boot.ci(boot_mass_max_temp_blups_low_adj_noout_lmer,
                                                  type = c('perc', 'norm', 'basic'),
                                                  index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_blups_low_adj_noout_lmer)


## MED parental care - adjusted, outliers removed
mass_max_temp_blups_med_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                           (1|fnest_id), 
                                         data = subset(mass_max_outliers_removed,
                                                       feeding_expontd_blups_strat == 2,
                                                       !is.na(x = mass_pre_obs) & 
                                                         !is.na(x = nest_max_temp) &
                                                         !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_max_temp_blups_med_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_blups_med_adj_noout_lmer))
  qqline(resid(mass_max_temp_blups_med_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_blups_med_adj_noout_lmer))


summary(mass_max_temp_blups_med_adj_noout_lmer)
confint(mass_max_temp_blups_med_adj_noout_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_blups_med_adj_noout_lmer <- bootMer(x = mass_max_temp_blups_med_adj_noout_lmer,
                                                 FUN = fixef, nsim = 2000,
                                                 seed = 632760,
                                                 use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_blups_med_adj_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_blups_med_adj_noout_lmer <- boot.ci(boot_mass_max_temp_blups_med_adj_noout_lmer,
                                                  type = c('perc', 'norm', 'basic'),
                                                  index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_blups_med_adj_noout_lmer)


## HIGH parental care - adjusted, outliers removed
mass_max_temp_blups_high_noout_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                            scale(nestling_number) + scale(days_summer) + 
                                            (1|fnest_id), 
                                          data = subset(mass_max_outliers_removed,
                                                        feeding_expontd_blups_strat == 3,
                                                        !is.na(x = mass_pre_obs) & 
                                                          !is.na(x = nest_max_temp) &
                                                          !is.na(x = feeding_expontd_blups)))


## Check diagnostics for the full model
plot(mass_max_temp_blups_high_noout_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_blups_high_noout_adj_lmer))
  qqline(resid(mass_max_temp_blups_high_noout_adj_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_blups_high_noout_adj_lmer))


summary(mass_max_temp_blups_high_noout_adj_lmer)
confint(mass_max_temp_blups_high_noout_adj_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_blups_high_noout_adj_lmer <- bootMer(x = mass_max_temp_blups_high_noout_adj_lmer,
                                                  FUN = fixef, nsim = 2000,
                                                  seed = 632760,
                                                  use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_blups_high_noout_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_blups_high_noout_adj_lmer <- boot.ci(boot_mass_max_temp_blups_high_noout_adj_lmer,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_blups_high_noout_adj_lmer)



##################### Temperature interquartile range #########################

### Interaction model - unadjusted
mass_iqr_temp_blups_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) *
                                   scale(feeding_expontd_blups) +
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

# Calculate R squared 
r.squaredGLMM(mass_iqr_temp_blups_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_blups_lmer <- bootMer(x = mass_iqr_temp_blups_lmer,
                                         FUN = fixef, nsim = 2000,
                                         seed = 632760,
                                         use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_blups_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_blups_lmer <- boot.ci(boot_mass_iqr_temp_blups_lmer,
                                          type = c('perc', 'norm', 'basic'),
                                          index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_blups_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_iqr_temp_blups_lmer_2 <- boot.ci(boot_mass_iqr_temp_blups_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 3) # CI for 2nd betas
print(bt_ci_mass_iqr_temp_blups_lmer_2)

# use 'boot' package to generate 95% CI for 3rd beta
bt_ci_mass_iqr_temp_blups_lmer_3 <- boot.ci(boot_mass_iqr_temp_blups_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 4) # CI for 3rd betas
print(bt_ci_mass_iqr_temp_blups_lmer_3)


### Stratified models - unadjusted
## LOW parental care model - unadjusted
mass_iqr_temp_blups_low_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) + 
                                       (1|fnest_id), 
                                     data = subset(late_nestling_parent_care,
                                                   feeding_expontd_blups_strat == 1,
                                                   !is.na(x = mass_pre_obs) & 
                                                     !is.na(x = nest_iqr_temp) &
                                                     !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_iqr_temp_blups_low_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_blups_low_lmer))
  qqline(resid(mass_iqr_temp_blups_low_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_blups_low_lmer))


summary(mass_iqr_temp_blups_low_lmer)
confint(mass_iqr_temp_blups_low_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_blups_low_lmer <- bootMer(x = mass_iqr_temp_blups_low_lmer,
                                             FUN = fixef, nsim = 2000,
                                             seed = 632760,
                                             use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_blups_low_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_blups_low_lmer <- boot.ci(boot_mass_iqr_temp_blups_low_lmer,
                                              type = c('perc', 'norm', 'basic'),
                                              index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_blups_low_lmer)


## MED parental care - unadjusted
mass_iqr_temp_blups_med_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) + 
                                       (1|fnest_id), 
                                     data = subset(late_nestling_parent_care,
                                                   feeding_expontd_blups_strat == 2,
                                                   !is.na(x = mass_pre_obs) & 
                                                     !is.na(x = nest_iqr_temp) &
                                                     !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_iqr_temp_blups_med_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_blups_med_lmer))
  qqline(resid(mass_iqr_temp_blups_med_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_blups_med_lmer))


summary(mass_iqr_temp_blups_med_lmer)
confint(mass_iqr_temp_blups_med_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_blups_med_lmer <- bootMer(x = mass_iqr_temp_blups_med_lmer,
                                             FUN = fixef, nsim = 2000,
                                             seed = 632760,
                                             use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_blups_med_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_blups_med_lmer <- boot.ci(boot_mass_iqr_temp_blups_med_lmer,
                                              type = c('perc', 'norm', 'basic'),
                                              index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_blups_med_lmer)


## HIGH parental care - unadjusted 
mass_iqr_temp_blups_high_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) +
                                        (1|fnest_id), 
                                      data = subset(late_nestling_parent_care,
                                                    feeding_expontd_blups_strat == 3,
                                                    !is.na(x = mass_pre_obs) & 
                                                      !is.na(x = nest_iqr_temp) &
                                                      !is.na(x = feeding_expontd_blups)))


## Check diagnostics for the full model
plot(mass_iqr_temp_blups_high_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_blups_high_lmer))
  qqline(resid(mass_iqr_temp_blups_high_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_blups_high_lmer))


summary(mass_iqr_temp_blups_high_lmer)
confint(mass_iqr_temp_blups_high_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_blups_high_lmer <- bootMer(x = mass_iqr_temp_blups_high_lmer,
                                              FUN = fixef, nsim = 2000,
                                              seed = 632760,
                                              use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_blups_high_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_blups_high_lmer <- boot.ci(boot_mass_iqr_temp_blups_high_lmer,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_blups_high_lmer)


### Interaction model - adjusted
mass_iqr_temp_blups_adj_lmer <- lmerTest::lmer(mass_pre_obs ~ scale(nest_iqr_temp) *
                                   scale(feeding_expontd_blups) + 
                                     scale(nestling_number) + scale(days_summer) +
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

# Calculate R squared 
r.squaredGLMM(mass_iqr_temp_blups_adj_lmer)


## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_blups_adj_lmer <- bootMer(x = mass_iqr_temp_blups_adj_lmer,
                                             FUN = fixef, nsim = 2000,
                                             seed = 632760,
                                             use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_blups_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_blups_adj_lmer <- boot.ci(boot_mass_iqr_temp_blups_adj_lmer,
                                              type = c('perc', 'norm', 'basic'),
                                              index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_blups_adj_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_iqr_temp_blups_adj_lmer_2 <- boot.ci(boot_mass_iqr_temp_blups_adj_lmer,
                                                type = c('perc', 'norm', 'basic'),
                                                index = 3) # CI for 2nd betas
print(bt_ci_mass_iqr_temp_blups_adj_lmer_2)

# use 'boot' package to generate 95% CI for 5th beta
bt_ci_mass_iqr_temp_blups_adj_lmer_5 <- boot.ci(boot_mass_iqr_temp_blups_adj_lmer,
                                                type = c('perc', 'norm', 'basic'),
                                                index = 6) # CI for 3rd betas
print(bt_ci_mass_iqr_temp_blups_adj_lmer_5)


### Stratified models - adjusted
## LOW parental care model - adjusted
mass_iqr_temp_blups_low_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                           (1|fnest_id), 
                                         data = subset(late_nestling_parent_care,
                                                       feeding_expontd_blups_strat == 1,
                                                       !is.na(x = mass_pre_obs) & 
                                                         !is.na(x = nest_iqr_temp) &
                                                         !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_iqr_temp_blups_low_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_blups_low_adj_lmer))
  qqline(resid(mass_iqr_temp_blups_low_adj_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_blups_low_adj_lmer))


summary(mass_iqr_temp_blups_low_adj_lmer)
confint(mass_iqr_temp_blups_low_adj_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_blups_low_adj_lmer <- bootMer(x = mass_iqr_temp_blups_low_adj_lmer,
                                                 FUN = fixef, nsim = 2000,
                                                 seed = 632760,
                                                 use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_blups_low_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_blups_low_adj_lmer <- boot.ci(boot_mass_iqr_temp_blups_low_adj_lmer,
                                                  type = c('perc', 'norm', 'basic'),
                                                  index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_blups_low_adj_lmer)


## MED parental care - adjusted
mass_iqr_temp_blups_med_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                           (1|fnest_id), 
                                         data = subset(late_nestling_parent_care,
                                                       feeding_expontd_blups_strat == 2,
                                                       !is.na(x = mass_pre_obs) & 
                                                         !is.na(x = nest_iqr_temp) &
                                                         !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_iqr_temp_blups_med_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_blups_med_adj_lmer))
  qqline(resid(mass_iqr_temp_blups_med_adj_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_blups_med_adj_lmer))


summary(mass_iqr_temp_blups_med_adj_lmer)
confint(mass_iqr_temp_blups_med_adj_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_blups_med_adj_lmer <- bootMer(x = mass_iqr_temp_blups_med_adj_lmer,
                                                 FUN = fixef, nsim = 2000,
                                                 seed = 632760,
                                                 use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_blups_med_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_blups_med_adj_lmer <- boot.ci(boot_mass_iqr_temp_blups_med_adj_lmer,
                                                  type = c('perc', 'norm', 'basic'),
                                                  index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_blups_med_adj_lmer)


## HIGH parental care - adjusted 
mass_iqr_temp_blups_high_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) +
                                            scale(nestling_number) + scale(days_summer) + 
                                            (1|fnest_id), 
                                          data = subset(late_nestling_parent_care,
                                                        feeding_expontd_blups_strat == 3,
                                                        !is.na(x = mass_pre_obs) & 
                                                          !is.na(x = nest_iqr_temp) &
                                                          !is.na(x = feeding_expontd_blups)))


## Check diagnostics for the full model
plot(mass_iqr_temp_blups_high_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_blups_high_adj_lmer))
  qqline(resid(mass_iqr_temp_blups_high_adj_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_blups_high_adj_lmer))

summary(mass_iqr_temp_blups_high_adj_lmer)
confint(mass_iqr_temp_blups_high_adj_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_blups_high_adj_lmer <- bootMer(x = mass_iqr_temp_blups_high_adj_lmer,
                                                  FUN = fixef, nsim = 2000,
                                                  seed = 632760,
                                                  use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_blups_high_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_blups_high_adj_lmer <- boot.ci(boot_mass_iqr_temp_blups_high_adj_lmer,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_blups_high_adj_lmer)


### Interaction model - unadjusted, outliers removed 
mass_iqr_temp_blups_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) *
                                   scale(feeding_expontd_blups) +
                                   (1|fnest_id), 
                                 data = subset(mass_max_outliers_removed,
                                               !is.na(x = mass_pre_obs) & 
                                                 !is.na(x = nest_iqr_temp) &
                                                 !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_iqr_temp_blups_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_blups_noout_lmer))
  qqline(resid(mass_iqr_temp_blups_noout_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_blups_noout_lmer))

summary(mass_iqr_temp_blups_noout_lmer)
confint(mass_iqr_temp_blups_noout_lmer)  

# Calculate R squared 
r.squaredGLMM(mass_iqr_temp_blups_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_blups_noout_lmer <- bootMer(x = mass_iqr_temp_blups_noout_lmer,
                                         FUN = fixef, nsim = 2000,
                                         seed = 632760,
                                         use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_blups_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_blups_noout_lmer <- boot.ci(boot_mass_iqr_temp_blups_noout_lmer,
                                          type = c('perc', 'norm', 'basic'),
                                          index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_blups_noout_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_iqr_temp_blups_noout_lmer_2 <- boot.ci(boot_mass_iqr_temp_blups_noout_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 3) # CI for 2nd betas
print(bt_ci_mass_iqr_temp_blups_noout_lmer_2)

# use 'boot' package to generate 95% CI for 3rd beta
bt_ci_mass_iqr_temp_blups_noout_lmer_3 <- boot.ci(boot_mass_iqr_temp_blups_noout_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 4) # CI for 3rd betas
print(bt_ci_mass_iqr_temp_blups_noout_lmer_3)


### Interaction model - adjusted, outliers removed 
mass_iqr_temp_blups_adj_noout_lmer <- lmerTest::lmer(mass_pre_obs ~ scale(nest_iqr_temp) *
                                       scale(feeding_expontd_blups) + 
                                       scale(nestling_number) + scale(days_summer) +
                                       (1|fnest_id), 
                                     data = subset(mass_max_outliers_removed,
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

# Calculate R squared 
r.squaredGLMM(mass_iqr_temp_blups_adj_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_blups_adj_noout_lmer <- bootMer(x = mass_iqr_temp_blups_adj_noout_lmer,
                                             FUN = fixef, nsim = 2000,
                                             seed = 632760,
                                             use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_blups_adj_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_blups_adj_noout_lmer <- boot.ci(boot_mass_iqr_temp_blups_adj_noout_lmer,
                                              type = c('perc', 'norm', 'basic'),
                                              index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_blups_adj_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_iqr_temp_blups_adj_noout_lmer_2 <- boot.ci(boot_mass_iqr_temp_blups_adj_noout_lmer,
                                                type = c('perc', 'norm', 'basic'),
                                                index = 3) # CI for 2nd betas
print(bt_ci_mass_iqr_temp_blups_adj_noout_lmer_2)

# use 'boot' package to generate 95% CI for 5th beta
bt_ci_mass_iqr_temp_blups_adj_noout_lmer_5 <- boot.ci(boot_mass_iqr_temp_blups_adj_noout_lmer,
                                                type = c('perc', 'norm', 'basic'),
                                                index = 6) # CI for 3rd betas
print(bt_ci_mass_iqr_temp_blups_adj_noout_lmer_5)


### Stratified models - adjusted, outliers removed
## LOW parental care model - adjusted, outliers removed
mass_iqr_temp_blups_low_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                           (1|fnest_id), 
                                         data = subset(mass_max_outliers_removed,
                                                       feeding_expontd_blups_strat == 1,
                                                       !is.na(x = mass_pre_obs) & 
                                                         !is.na(x = nest_iqr_temp) &
                                                         !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_iqr_temp_blups_low_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_blups_low_adj_noout_lmer))
  qqline(resid(mass_iqr_temp_blups_low_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_blups_low_adj_noout_lmer))


summary(mass_iqr_temp_blups_low_adj_noout_lmer)
confint(mass_iqr_temp_blups_low_adj_noout_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_blups_low_adj_noout_lmer <- bootMer(x = mass_iqr_temp_blups_low_adj_noout_lmer,
                                                 FUN = fixef, nsim = 2000,
                                                 seed = 632760,
                                                 use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_blups_low_adj_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_blups_low_adj_noout_lmer <- boot.ci(boot_mass_iqr_temp_blups_low_adj_noout_lmer,
                                                  type = c('perc', 'norm', 'basic'),
                                                  index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_blups_low_adj_noout_lmer)


## MED parental care - adjusted, outliers removed
mass_iqr_temp_blups_med_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                           (1|fnest_id), 
                                         data = subset(mass_max_outliers_removed,
                                                       feeding_expontd_blups_strat == 2,
                                                       !is.na(x = mass_pre_obs) & 
                                                         !is.na(x = nest_iqr_temp) &
                                                         !is.na(x = feeding_expontd_blups)))

## Check diagnostics for the full model
plot(mass_iqr_temp_blups_med_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_blups_med_adj_noout_lmer))
  qqline(resid(mass_iqr_temp_blups_med_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_blups_med_adj_noout_lmer))


summary(mass_iqr_temp_blups_med_adj_noout_lmer)
confint(mass_iqr_temp_blups_med_adj_noout_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_blups_med_adj_noout_lmer <- bootMer(x = mass_iqr_temp_blups_med_adj_noout_lmer,
                                                 FUN = fixef, nsim = 2000,
                                                 seed = 632760,
                                                 use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_blups_med_adj_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_blups_med_adj_noout_lmer <- boot.ci(boot_mass_iqr_temp_blups_med_adj_noout_lmer,
                                                  type = c('perc', 'norm', 'basic'),
                                                  index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_blups_med_adj_noout_lmer)


## HIGH parental care - adjusted, outliers removed
mass_iqr_temp_blups_high_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) +
                                            scale(nestling_number) + scale(days_summer) + 
                                            (1|fnest_id), 
                                          data = subset(mass_max_outliers_removed,
                                                        feeding_expontd_blups_strat == 3,
                                                        !is.na(x = mass_pre_obs) & 
                                                          !is.na(x = nest_iqr_temp) &
                                                          !is.na(x = feeding_expontd_blups)))


## Check diagnostics for the full model
plot(mass_iqr_temp_blups_high_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_blups_high_adj_noout_lmer))
  qqline(resid(mass_iqr_temp_blups_high_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_blups_high_adj_noout_lmer))

summary(mass_iqr_temp_blups_high_adj_noout_lmer)
confint(mass_iqr_temp_blups_high_adj_noout_lmer) 

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_blups_high_adj_noout_lmer <- bootMer(x = mass_iqr_temp_blups_high_adj_noout_lmer,
                                                  FUN = fixef, nsim = 2000,
                                                  seed = 632760,
                                                  use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_blups_high_adj_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_blups_high_adj_noout_lmer <- boot.ci(boot_mass_iqr_temp_blups_high_adj_noout_lmer,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_blups_high_adj_noout_lmer)



######## Create model results tables for straified parental care mdoels #######

# Extract number of observations for unadjusted models
low_un_n <- c(nobs(mass_min_temp_blups_low_lmer), 
              nobs(mass_max_temp_blups_low_lmer),
              nobs(mass_iqr_temp_blups_low_lmer))

med_un_n <- c(nobs(mass_min_temp_blups_med_lmer), 
              nobs(mass_max_temp_blups_med_lmer),
              nobs(mass_iqr_temp_blups_med_lmer))

high_un_n <- c(nobs(mass_min_temp_blups_high_lmer), 
               nobs(mass_max_temp_blups_high_lmer),
               nobs(mass_iqr_temp_blups_high_lmer))

# Extract Betas for unadjusted models
low_un_b <- c(paste(round(fixef(mass_min_temp_blups_low_lmer)[2], 2), " (", 
                    round(bt_ci_mass_min_temp_blups_low_lmer[[6]][4], 2), ", ", 
                    round(bt_ci_mass_min_temp_blups_low_lmer[[6]][5], 2), ")",
                    sep = ""),
              paste(round(fixef(mass_max_temp_blups_low_lmer)[2], 2), " (", 
                    round(bt_ci_mass_max_temp_blups_low_lmer[[6]][4], 2), ", ", 
                    round(bt_ci_mass_max_temp_blups_low_lmer[[6]][5], 2), ")",
                    sep = ""),
              paste(round(fixef(mass_iqr_temp_blups_low_lmer)[2], 2), " (", 
                    round(bt_ci_mass_iqr_temp_blups_low_lmer[[6]][4], 2), ", ", 
                    round(bt_ci_mass_iqr_temp_blups_low_lmer[[6]][5], 2), ")",
                    sep = ""))


med_un_b <- c(paste(round(fixef(mass_min_temp_blups_med_lmer)[2], 2), " (", 
                    round(bt_ci_mass_min_temp_blups_med_lmer[[6]][4], 2), ", ", 
                    round(bt_ci_mass_min_temp_blups_med_lmer[[6]][5], 2), ")",
                    sep = ""),
              paste(round(fixef(mass_max_temp_blups_med_lmer)[2], 2), " (", 
                    round(bt_ci_mass_max_temp_blups_med_lmer[[6]][4], 2), ", ", 
                    round(bt_ci_mass_max_temp_blups_med_lmer[[6]][5], 2), ")",
                    sep = ""),
              paste(round(fixef(mass_iqr_temp_blups_med_lmer)[2], 2), " (", 
                    round(bt_ci_mass_iqr_temp_blups_med_lmer[[6]][4], 2), ", ", 
                    round(bt_ci_mass_iqr_temp_blups_med_lmer[[6]][5], 2), ")",
                    sep = ""))


high_un_b <- c(paste(round(fixef(mass_min_temp_blups_high_lmer)[2], 2), " (", 
                     round(bt_ci_mass_min_temp_blups_high_lmer[[6]][4], 2), ", ", 
                     round(bt_ci_mass_min_temp_blups_high_lmer[[6]][5], 2), ")",
                     sep = ""),
               paste(round(fixef(mass_max_temp_blups_high_lmer)[2], 2), " (", 
                     round(bt_ci_mass_max_temp_blups_high_lmer[[6]][4], 2), ", ", 
                     round(bt_ci_mass_max_temp_blups_high_lmer[[6]][5], 2), ")",
                     sep = ""),
               paste(round(fixef(mass_iqr_temp_blups_high_lmer)[2], 2), " (", 
                     round(bt_ci_mass_iqr_temp_blups_high_lmer[[6]][4], 2), ", ", 
                     round(bt_ci_mass_iqr_temp_blups_high_lmer[[6]][5], 2), ")",
                     sep = ""))


# Extract number of observations for adjusted models
low_ad_n <- c(nobs(mass_min_temp_blups_low_adj_lmer), 
              nobs(mass_max_temp_blups_low_adj_lmer),
              nobs(mass_iqr_temp_blups_low_adj_lmer))

med_ad_n <- c(nobs(mass_min_temp_blups_med_adj_lmer), 
              nobs(mass_max_temp_blups_med_adj_lmer),
              nobs(mass_iqr_temp_blups_med_adj_lmer))

high_ad_n <- c(nobs(mass_min_temp_blups_high_adj_lmer), 
               nobs(mass_max_temp_blups_high_adj_lmer),
               nobs(mass_iqr_temp_blups_high_adj_lmer))

# Extract Betas for adjusted models
low_ad_b <- c(paste(round(fixef(mass_min_temp_blups_low_adj_lmer)[2], 2), " (", 
                    round(bt_ci_mass_min_temp_blups_low_adj_lmer[[6]][4], 2), ", ", 
                    round(bt_ci_mass_min_temp_blups_low_adj_lmer[[6]][5], 2), ")",
                    sep = ""),
              paste(round(fixef(mass_max_temp_blups_low_adj_lmer)[2], 2), " (", 
                    round(bt_ci_mass_max_temp_blups_low_adj_lmer[[6]][4], 2), ", ", 
                    round(bt_ci_mass_max_temp_blups_low_adj_lmer[[6]][5], 2), ")",
                    sep = ""),
              paste(round(fixef(mass_iqr_temp_blups_low_adj_lmer)[2], 2), " (", 
                    round(bt_ci_mass_iqr_temp_blups_low_adj_lmer[[6]][4], 2), ", ", 
                    round(bt_ci_mass_iqr_temp_blups_low_adj_lmer[[6]][5], 2), ")",
                    sep = ""))

med_ad_b <- c(paste(round(fixef(mass_min_temp_blups_med_adj_lmer)[2], 2), " (", 
                    round(bt_ci_mass_min_temp_blups_med_adj_lmer[[6]][4], 2), ", ", 
                    round(bt_ci_mass_min_temp_blups_med_adj_lmer[[6]][5], 2), ")",
                    sep = ""),
              paste(round(fixef(mass_max_temp_blups_med_adj_lmer)[2], 2), " (", 
                    round(bt_ci_mass_max_temp_blups_med_adj_lmer[[6]][4], 2), ", ", 
                    round(bt_ci_mass_max_temp_blups_med_adj_lmer[[6]][5], 2), ")",
                    sep = ""),
              paste(round(fixef(mass_iqr_temp_blups_med_adj_lmer)[2], 2), " (", 
                    round(bt_ci_mass_iqr_temp_blups_med_adj_lmer[[6]][4], 2), ", ", 
                    round(bt_ci_mass_iqr_temp_blups_med_adj_lmer[[6]][5], 2), ")",
                    sep = ""))


high_ad_b <- c(paste(round(fixef(mass_min_temp_blups_high_adj_lmer)[2], 2), " (", 
                     round(bt_ci_mass_min_temp_blups_high_adj_lmer[[6]][4], 2), ", ", 
                     round(bt_ci_mass_min_temp_blups_high_adj_lmer[[6]][5], 2), ")",
                     sep = ""),
               paste(round(fixef(mass_max_temp_blups_high_adj_lmer)[2], 2), " (", 
                     round(bt_ci_mass_max_temp_blups_high_adj_lmer[[6]][4], 2), ", ", 
                     round(bt_ci_mass_max_temp_blups_high_adj_lmer[[6]][5], 2), ")",
                     sep = ""),
               paste(round(fixef(mass_iqr_temp_blups_high_adj_lmer)[2], 2), " (", 
                     round(bt_ci_mass_iqr_temp_blups_high_adj_lmer[[6]][4], 2), ", ", 
                     round(bt_ci_mass_iqr_temp_blups_high_adj_lmer[[6]][5], 2), ")",
                     sep = ""))



# Bind the values together into a single dataframe
parental_care_strat_ad <- data.frame(
  variable <- c("Minimum temperature", "Maximum temperature", "Temperature variability"),
  type <- c("Adjusted", "Adjusted", "Adjusted"),
  low_ad_n, low_ad_b, 
  med_ad_n, med_ad_b, 
  high_ad_n, high_ad_b
)

parental_care_strat_un <- data.frame(
  variable <- c("Minimum temperature", "Maximum temperature", "Temperature variability"),
  type <- c("Unadjusted", "Unadjusted", "Unadjusted"),
  low_un_n, low_un_b, 
  med_un_n, med_un_b, 
  high_un_n, high_un_b
)

colnames(parental_care_strat_ad) <- c("variable", "type",
                                      "low_n", "low_b", 
                                      "med_n", "med_b", 
                                      "high_n", "high_b")

colnames(parental_care_strat_un) <- c("variable", "type",
                                      "low_n", "low_b", 
                                      "med_n", "med_b", 
                                      "high_n", "high_b")

parental_care_strat <- rbind(parental_care_strat_un, parental_care_strat_ad)

parental_care_strat <- arrange(parental_care_strat,
                               factor(variable, levels = c("Minimum temperature", "Maximum temperature", "Temperature variability")))


str(parental_care_strat)

# Create display table
parental_care_strat_mod_table <- gt(parental_care_strat, rowname_col = "type") %>%
  tab_header(
    title = md("**Supplemental Table 5.** Associations of nestling mass and temperature, assessed in separate models stratified by levels of parental feeding, when influential outliers are included in the dataset. Temperature variability is defined as the interquartile range.")
  ) %>%
  tab_footnote(
    footnote = "Estimated β (95% CI) from stratified linear mixed models in which temperature is the explanatory variable of interest, nestling mass is the outcome of interest, and nest ID was included as a random intercept. Adjusted models include hatch date and number of nestlings in the nest. Continuous predictors are z-score standardized.",
    locations = cells_column_labels(columns = c(low_b, med_b, high_b))
  ) %>%
  tab_footnote(
    footnote = md(paste("R-squared for adjusted minimum temperature models. ",
                        "Low parental feeding model: Marginal R-squared = ", round(r.squaredGLMM(mass_min_temp_blups_low_adj_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_min_temp_blups_low_adj_lmer)[2], 2),
                        "; Medium parental feeding model: Marginal R-squared = ", round(r.squaredGLMM(mass_min_temp_blups_med_adj_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_min_temp_blups_med_adj_lmer)[2], 2),
                        "; High parental feeding model: Marginal R-squared = ", round(r.squaredGLMM(mass_min_temp_blups_high_adj_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_min_temp_blups_high_adj_lmer)[2], 2),
                 sep = "")),
    locations = cells_stub(rows = c(2))
    ) %>%
  tab_footnote(
    footnote = md(paste("R-squared for adjusted maximum temperature models. ",
                        "Low parental feeding model: Marginal R-squared = ", round(r.squaredGLMM(mass_max_temp_blups_low_adj_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_max_temp_blups_low_adj_lmer)[2], 2),
                        "; Medium parental feeding model: Marginal R-squared = ", round(r.squaredGLMM(mass_max_temp_blups_med_adj_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_max_temp_blups_med_adj_lmer)[2], 2),
                        "; High parental feeding model: Marginal R-squared = ", round(r.squaredGLMM(mass_max_temp_blups_high_adj_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_max_temp_blups_high_adj_lmer)[2], 2),
                        sep = "")),
    locations = cells_stub(rows = c(4))
  ) %>%
  tab_footnote(
    footnote = md(paste("R-squared for adjusted temperature variability models. ",
                        "Low parental feeding model: Marginal R-squared = ", round(r.squaredGLMM(mass_iqr_temp_blups_low_adj_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_iqr_temp_blups_low_adj_lmer)[2], 2),
                        "; Medium parental feeding model: Marginal R-squared = ", round(r.squaredGLMM(mass_iqr_temp_blups_med_adj_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_iqr_temp_blups_med_adj_lmer)[2], 2),
                        "; High parental feeding model: Marginal R-squared = ", round(r.squaredGLMM(mass_iqr_temp_blups_high_adj_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_iqr_temp_blups_high_adj_lmer)[2], 2),
                        sep = "")),
    locations = cells_stub(rows = c(6))
  ) %>%
  tab_stubhead(
    label = md("Type")
  ) %>%
  tab_row_group(
    label = md("**Effect of temperature variability**"), 
    rows = c(5:6)
  ) %>%
  tab_row_group(
    label = md("**Effect of maximum temperature**"), 
    rows = c(3:4)
  ) %>%
  tab_row_group(
    label = md("**Effect of minimum temperature**"), 
    rows = c(1:2)
  ) %>%
  tab_spanner(
    label = "Low parental feeding models",
    columns = c(low_n, low_b)
  ) %>%
  tab_spanner(
    label = "Medium parental feeding models",
    columns = c(med_n, med_b)
  ) %>%
  tab_spanner(
    label = "High parental feeding models",
    columns = c(high_n, high_b)
  ) %>% 
  cols_label(
    ends_with("n") ~ "N", 
    ends_with("b") ~ "β (95% CI)",
    type = "Type"
  ) %>%
  cols_hide(variable)  %>%
  opt_table_font(font = "Arial", size = 12)  %>%
  tab_options(footnotes.font.size = 10)


parental_care_strat_mod_table
gtsave(parental_care_strat_mod_table, filename = "Output/parental_care_adjusted_strat.docx") 
gtsave(parental_care_strat_mod_table, filename = "Output/parental_care_adjusted_strat.html")



###############################################################################
##############          Growth, temp, and relative size          ##############
###############################################################################

########## Relative size calculated at mid development ########################

########################### Minimum temp ######################################

### Interaction model - unadjusted
mass_min_temp_mid_size_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) *
                                  mid_size_order +
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

# Calculate R squared 
r.squaredGLMM(mass_min_temp_mid_size_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_mid_size_lmer <- bootMer(x = mass_min_temp_mid_size_lmer,
                                         FUN = fixef, nsim = 2000,
                                         seed = 632760,
                                         use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_mid_size_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_mid_size_lmer <- boot.ci(boot_mass_min_temp_mid_size_lmer,
                                          type = c('perc', 'norm', 'basic'),
                                          index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_mid_size_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_min_temp_mid_size_lmer_2 <- boot.ci(boot_mass_min_temp_mid_size_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 3) # CI for 2nd betas
print(bt_ci_mass_min_temp_mid_size_lmer_2)

# use 'boot' package to generate 95% CI for 3rd beta
bt_ci_mass_min_temp_mid_size_lmer_3 <- boot.ci(boot_mass_min_temp_mid_size_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 4) # CI for 3rd betas
print(bt_ci_mass_min_temp_mid_size_lmer_3)


### Stratified models - unadjusted
# SMALL nestlings - unadjusted
mass_min_temp_mid_size_small_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) +
                                                (1|fnest_id), 
                                              data = subset(late_nestling_parent_care,
                                                            mid_size_order == "min",
                                                            !is.na(x = mass_pre_obs) & 
                                                              !is.na(x = nest_min_temp) &
                                                              !is.na(x = mid_size_order)))

## Check diagnostics for the full model
plot(mass_min_temp_mid_size_small_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_mid_size_small_lmer))
  qqline(resid(mass_min_temp_mid_size_small_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_mid_size_small_lmer))

summary(mass_min_temp_mid_size_small_lmer)
confint(mass_min_temp_mid_size_small_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_mid_size_small_lmer <- bootMer(x = mass_min_temp_mid_size_small_lmer,
                                                      FUN = fixef, nsim = 2000,
                                                      seed = 632760,
                                                      use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_mid_size_small_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_mid_size_small_lmer <- boot.ci(boot_mass_min_temp_mid_size_small_lmer,
                                                       type = c('perc', 'norm', 'basic'),
                                                       index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_mid_size_small_lmer)


# BIG nestlings - unadjusted
mass_min_temp_mid_size_big_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) +
                                              (1|fnest_id), 
                                            data = subset(late_nestling_parent_care,
                                                          mid_size_order == "other",
                                                          !is.na(x = mass_pre_obs) & 
                                                            !is.na(x = nest_min_temp) &
                                                            !is.na(x = mid_size_order)))

## Check diagnostics for the full model
plot(mass_min_temp_mid_size_big_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_mid_size_big_lmer))
  qqline(resid(mass_min_temp_mid_size_big_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_mid_size_big_lmer))

summary(mass_min_temp_mid_size_big_lmer)
confint(mass_min_temp_mid_size_big_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_mid_size_big_lmer <- bootMer(x = mass_min_temp_mid_size_big_lmer,
                                                    FUN = fixef, nsim = 2000,
                                                    seed = 632760,
                                                    use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_mid_size_big_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_mid_size_big_lmer <- boot.ci(boot_mass_min_temp_mid_size_big_lmer,
                                                     type = c('perc', 'norm', 'basic'),
                                                     index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_mid_size_big_lmer)


### Interaction model - adjusted
mass_min_temp_mid_size_adj_lmer <- lmerTest::lmer(mass_pre_obs ~ scale(nest_min_temp) *
                                      mid_size_order + scale(nestling_number) +
                                      scale(days_summer) +
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

# Calculate R squared 
r.squaredGLMM(mass_min_temp_mid_size_adj_lmer)


## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_mid_size_adj_lmer <- bootMer(x = mass_min_temp_mid_size_adj_lmer,
                                                   FUN = fixef, nsim = 2000,
                                                   seed = 632760,
                                                   use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_mid_size_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_mid_size_adj_lmer <- boot.ci(boot_mass_min_temp_mid_size_adj_lmer,
                                                    type = c('perc', 'norm', 'basic'),
                                                    index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_mid_size_adj_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_min_temp_mid_size_adj_lmer_2 <- boot.ci(boot_mass_min_temp_mid_size_adj_lmer,
                                                      type = c('perc', 'norm', 'basic'),
                                                      index = 3) # CI for 2nd betas
print(bt_ci_mass_min_temp_mid_size_adj_lmer_2)

# use 'boot' package to generate 95% CI for 5th beta
bt_ci_mass_min_temp_mid_size_adj_lmer_5 <- boot.ci(boot_mass_min_temp_mid_size_adj_lmer,
                                                      type = c('perc', 'norm', 'basic'),
                                                      index = 6) # CI for 3rd betas
print(bt_ci_mass_min_temp_mid_size_adj_lmer_5)


## Stratified models - adjusted
# SMALL nestlings - adjusted
mass_min_temp_mid_size_small_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) +
                                                scale(nestling_number) + scale(days_summer) +
                                                (1|fnest_id), 
                                              data = subset(late_nestling_parent_care,
                                                            mid_size_order == "min",
                                                            !is.na(x = mass_pre_obs) & 
                                                              !is.na(x = nest_min_temp) &
                                                              !is.na(x = mid_size_order)))

## Check diagnostics for the full model
plot(mass_min_temp_mid_size_small_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_mid_size_small_adj_lmer))
  qqline(resid(mass_min_temp_mid_size_small_adj_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_mid_size_small_adj_lmer))

summary(mass_min_temp_mid_size_small_adj_lmer)
confint(mass_min_temp_mid_size_small_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_mid_size_small_adj_lmer <- bootMer(x = mass_min_temp_mid_size_small_adj_lmer,
                                                FUN = fixef, nsim = 2000,
                                                seed = 632760,
                                                use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_mid_size_small_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_mid_size_small_adj_lmer <- boot.ci(boot_mass_min_temp_mid_size_small_adj_lmer,
                                                 type = c('perc', 'norm', 'basic'),
                                                 index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_mid_size_small_adj_lmer)


# BIG nestlings - adjusted
mass_min_temp_mid_size_big_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) +
                                              scale(nestling_number) + scale(days_summer) +
                                              (1|fnest_id), 
                                            data = subset(late_nestling_parent_care,
                                                          mid_size_order == "other",
                                                          !is.na(x = mass_pre_obs) & 
                                                            !is.na(x = nest_min_temp) &
                                                            !is.na(x = mid_size_order)))

## Check diagnostics for the full model
plot(mass_min_temp_mid_size_small_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_mid_size_small_adj_lmer))
  qqline(resid(mass_min_temp_mid_size_small_adj_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_mid_size_small_adj_lmer))

summary(mass_min_temp_mid_size_big_adj_lmer)
confint(mass_min_temp_mid_size_big_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_mid_size_big_adj_lmer <- bootMer(x = mass_min_temp_mid_size_big_adj_lmer,
                                                      FUN = fixef, nsim = 2000,
                                                      seed = 632760,
                                                      use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_mid_size_big_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_mid_size_big_adj_lmer <- boot.ci(boot_mass_min_temp_mid_size_big_adj_lmer,
                                                       type = c('perc', 'norm', 'basic'),
                                                       index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_mid_size_big_adj_lmer)


############################## Maximum temperature ############################

### Interaction model - unadjusted
mass_max_temp_mid_size_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) *
                                  mid_size_order +
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

# Calculate R squared 
r.squaredGLMM(mass_max_temp_mid_size_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_mid_size_lmer<- bootMer(x = mass_max_temp_mid_size_lmer,
                                            FUN = fixef, nsim = 2000,
                                            seed = 632760,
                                            use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_mid_size_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_mid_size_lmer <- boot.ci(boot_mass_max_temp_mid_size_lmer,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_mid_size_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_max_temp_mid_size_lmer_2 <- boot.ci(boot_mass_max_temp_mid_size_lmer,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 3) # CI for 2nd betas
print(bt_ci_mass_max_temp_mid_size_lmer_2)

# use 'boot' package to generate 95% CI for 3rd beta
bt_ci_mass_max_temp_mid_size_lmer_3 <- boot.ci(boot_mass_max_temp_mid_size_lmer,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 4) # CI for 3rd betas
print(bt_ci_mass_max_temp_mid_size_lmer_3)


### Stratified models - unadjusted
# SMALL nestlings - unadjusted
mass_max_temp_mid_size_small_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                                (1|fnest_id), 
                                              data = subset(late_nestling_parent_care,
                                                            mid_size_order == "min",
                                                            !is.na(x = mass_pre_obs) & 
                                                              !is.na(x = nest_max_temp) &
                                                              !is.na(x = mid_size_order)))

## Check diagnostics for the full model
plot(mass_max_temp_mid_size_small_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_mid_size_small_lmer))
  qqline(resid(mass_max_temp_mid_size_small_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_mid_size_small_lmer))

summary(mass_max_temp_mid_size_small_lmer)
confint(mass_max_temp_mid_size_small_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_mid_size_small_lmer <- bootMer(x = mass_max_temp_mid_size_small_lmer,
                                                      FUN = fixef, nsim = 2000,
                                                      seed = 632760,
                                                      use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_mid_size_small_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_mid_size_small_lmer <- boot.ci(boot_mass_max_temp_mid_size_small_lmer,
                                                       type = c('perc', 'norm', 'basic'),
                                                       index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_mid_size_small_lmer)


# BIG nestlings - unadjusted
mass_max_temp_mid_size_big_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                              (1|fnest_id), 
                                            data = subset(late_nestling_parent_care,
                                                          mid_size_order == "other",
                                                          !is.na(x = mass_pre_obs) & 
                                                            !is.na(x = nest_max_temp) &
                                                            !is.na(x = mid_size_order)))

## Check diagnostics for the full model
plot(mass_max_temp_mid_size_big_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_mid_size_big_lmer))
  qqline(resid(mass_max_temp_mid_size_big_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_mid_size_big_lmer))

summary(mass_max_temp_mid_size_big_lmer)
confint(mass_max_temp_mid_size_big_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_mid_size_big_lmer <- bootMer(x = mass_max_temp_mid_size_big_lmer,
                                                    FUN = fixef, nsim = 2000,
                                                    seed = 632760,
                                                    use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_mid_size_big_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_mid_size_big_lmer <- boot.ci(boot_mass_max_temp_mid_size_big_lmer,
                                                     type = c('perc', 'norm', 'basic'),
                                                     index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_mid_size_big_lmer)


### Interaction model - adjusted
mass_max_temp_mid_size_adj_lmer <- lmerTest::lmer(mass_pre_obs ~ scale(nest_max_temp) *
                                      mid_size_order + scale(nestling_number) +
                                      scale(days_summer) +
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

# Calculate R squared 
r.squaredGLMM(mass_max_temp_mid_size_adj_lmer)


## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_mid_size_adj_lmer <- bootMer(x = mass_max_temp_mid_size_adj_lmer,
                                                FUN = fixef, nsim = 2000,
                                                seed = 632760,
                                                use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_mid_size_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_mid_size_adj_lmer <- boot.ci(boot_mass_max_temp_mid_size_adj_lmer,
                                                 type = c('perc', 'norm', 'basic'),
                                                 index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_mid_size_adj_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_max_temp_mid_size_adj_lmer_2 <- boot.ci(boot_mass_max_temp_mid_size_adj_lmer,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 3) # CI for 2nd betas
print(bt_ci_mass_max_temp_mid_size_adj_lmer_2)

# use 'boot' package to generate 95% CI for 5th beta
bt_ci_mass_max_temp_mid_size_adj_lmer_5 <- boot.ci(boot_mass_max_temp_mid_size_adj_lmer,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 6) # CI for 3rd betas
print(bt_ci_mass_max_temp_mid_size_adj_lmer_5)


### Stratified models - adjusted
# SMALL nestlings - adjusted
mass_max_temp_mid_size_small_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                            scale(nestling_number) + scale(days_summer) +
                                              (1|fnest_id), 
                                          data = subset(late_nestling_parent_care,
                                                        mid_size_order == "min",
                                                        !is.na(x = mass_pre_obs) & 
                                                          !is.na(x = nest_max_temp) &
                                                          !is.na(x = mid_size_order)))

## Check diagnostics for the full model
plot(mass_max_temp_mid_size_small_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_mid_size_small_adj_lmer))
  qqline(resid(mass_max_temp_mid_size_small_adj_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_mid_size_small_adj_lmer))

summary(mass_max_temp_mid_size_small_adj_lmer)
confint(mass_max_temp_mid_size_small_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_mid_size_small_adj_lmer <- bootMer(x = mass_max_temp_mid_size_small_adj_lmer,
                                                    FUN = fixef, nsim = 2000,
                                                    seed = 632760,
                                                    use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_mid_size_small_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_mid_size_small_adj_lmer <- boot.ci(boot_mass_max_temp_mid_size_small_adj_lmer,
                                                     type = c('perc', 'norm', 'basic'),
                                                     index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_mid_size_small_adj_lmer)


# BIG nestlings - adjusted
mass_max_temp_mid_size_big_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                          scale(nestling_number) + scale(days_summer) +
                                          (1|fnest_id), 
                                        data = subset(late_nestling_parent_care,
                                                      mid_size_order == "other",
                                                      !is.na(x = mass_pre_obs) & 
                                                        !is.na(x = nest_max_temp) &
                                                        !is.na(x = mid_size_order)))

## Check diagnostics for the full model
plot(mass_max_temp_mid_size_big_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_mid_size_big_adj_lmer))
  qqline(resid(mass_max_temp_mid_size_big_adj_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_mid_size_big_adj_lmer))

summary(mass_max_temp_mid_size_big_adj_lmer)
confint(mass_max_temp_mid_size_big_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_mid_size_big_adj_lmer <- bootMer(x = mass_max_temp_mid_size_big_adj_lmer,
                                                    FUN = fixef, nsim = 2000,
                                                    seed = 632760,
                                                    use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_mid_size_big_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_mid_size_big_adj_lmer <- boot.ci(boot_mass_max_temp_mid_size_big_adj_lmer,
                                                     type = c('perc', 'norm', 'basic'),
                                                     index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_mid_size_big_adj_lmer)


######################### Temperature interquartile range #####################

### Interaction model - unadjusted
mass_iqr_temp_mid_size_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) *
                                  mid_size_order +
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

# Calculate R squared 
r.squaredGLMM(mass_iqr_temp_mid_size_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_mid_size_lmer <- bootMer(x = mass_iqr_temp_mid_size_lmer,
                                           FUN = fixef, nsim = 2000,
                                           seed = 632760,
                                           use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_mid_size_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_mid_size_lmer <- boot.ci(boot_mass_iqr_temp_mid_size_lmer,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_mid_size_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_iqr_temp_mid_size_lmer_2 <- boot.ci(boot_mass_iqr_temp_mid_size_lmer,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 3) # CI for 2nd betas
print(bt_ci_mass_iqr_temp_mid_size_lmer_2)

# use 'boot' package to generate 95% CI for 3rd beta
bt_ci_mass_iqr_temp_mid_size_lmer_3 <- boot.ci(boot_mass_iqr_temp_mid_size_lmer,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 4) # CI for 3rd betas
print(bt_ci_mass_iqr_temp_mid_size_lmer_3)


### Stratified models - unadjusted
# SMALL nestlings - unadjusted
mass_iqr_temp_mid_size_small_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) + 
                                                (1|fnest_id), 
                                              data = subset(late_nestling_parent_care,
                                                            mid_size_order == "min",
                                                            !is.na(x = mass_pre_obs) & 
                                                              !is.na(x = nest_iqr_temp) &
                                                              !is.na(x = mid_size_order)))

## Check diagnostics for the full model
plot(mass_iqr_temp_mid_size_small_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_mid_size_small_lmer))
  qqline(resid(mass_iqr_temp_mid_size_small_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_mid_size_small_lmer))

summary(mass_iqr_temp_mid_size_small_lmer)
confint(mass_iqr_temp_mid_size_small_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_mid_size_small_lmer <- bootMer(x = mass_iqr_temp_mid_size_small_lmer,
                                                      FUN = fixef, nsim = 2000,
                                                      seed = 632760,
                                                      use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_mid_size_small_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_mid_size_small_lmer <- boot.ci(boot_mass_iqr_temp_mid_size_small_lmer,
                                                       type = c('perc', 'norm', 'basic'),
                                                       index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_mid_size_small_lmer)


# BIG nestlings - unadjusted
mass_iqr_temp_mid_size_big_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) + 
                                              (1|fnest_id), 
                                            data = subset(late_nestling_parent_care,
                                                          mid_size_order == "other",
                                                          !is.na(x = mass_pre_obs) & 
                                                            !is.na(x = nest_iqr_temp) &
                                                            !is.na(x = mid_size_order)))


## Check diagnostics for the full model
plot(mass_iqr_temp_mid_size_big_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_mid_size_big_lmer))
  qqline(resid(mass_iqr_temp_mid_size_big_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_mid_size_big_lmer))

summary(mass_iqr_temp_mid_size_big_lmer)
confint(mass_iqr_temp_mid_size_big_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_mid_size_big_lmer <- bootMer(x = mass_iqr_temp_mid_size_big_lmer,
                                                    FUN = fixef, nsim = 2000,
                                                    seed = 632760,
                                                    use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_mid_size_big_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_mid_size_big_lmer <- boot.ci(boot_mass_iqr_temp_mid_size_big_lmer,
                                                     type = c('perc', 'norm', 'basic'),
                                                     index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_mid_size_big_lmer)


### Interaction model - adjusted
mass_iqr_temp_mid_size_adj_lmer <- lmerTest::lmer(mass_pre_obs ~ scale(nest_iqr_temp) *
                                      mid_size_order + scale(nestling_number) + 
                                      scale(days_summer) +
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

# Calculate R squared 
r.squaredGLMM(mass_iqr_temp_mid_size_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_mid_size_adj_lmer <- bootMer(x = mass_iqr_temp_mid_size_adj_lmer,
                                                FUN = fixef, nsim = 2000,
                                                seed = 632760,
                                                use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_mid_size_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_mid_size_adj_lmer <- boot.ci(boot_mass_iqr_temp_mid_size_adj_lmer,
                                                 type = c('perc', 'norm', 'basic'),
                                                 index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_mid_size_adj_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_iqr_temp_mid_size_adj_lmer_2 <- boot.ci(boot_mass_iqr_temp_mid_size_adj_lmer,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 3) # CI for 2nd betas
print(bt_ci_mass_iqr_temp_mid_size_adj_lmer_2)

# use 'boot' package to generate 95% CI for 5th beta
bt_ci_mass_iqr_temp_mid_size_adj_lmer_5 <- boot.ci(boot_mass_iqr_temp_mid_size_adj_lmer,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 6) # CI for 3rd betas
print(bt_ci_mass_iqr_temp_mid_size_adj_lmer_5)


### Stratified models - adjusted
# SMALL nestlings - adjusted
mass_iqr_temp_mid_size_small_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) + 
                                            scale(nestling_number) + 
                                            scale(days_summer) +
                                              (1|fnest_id), 
                                          data = subset(late_nestling_parent_care,
                                                        mid_size_order == "min",
                                                        !is.na(x = mass_pre_obs) & 
                                                          !is.na(x = nest_iqr_temp) &
                                                          !is.na(x = mid_size_order)))

## Check diagnostics for the full model
plot(mass_iqr_temp_mid_size_small_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_mid_size_small_adj_lmer))
  qqline(resid(mass_iqr_temp_mid_size_small_adj_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_mid_size_small_adj_lmer))

summary(mass_iqr_temp_mid_size_small_adj_lmer)
confint(mass_iqr_temp_mid_size_small_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_mid_size_small_adj_lmer <- bootMer(x = mass_iqr_temp_mid_size_small_adj_lmer,
                                                    FUN = fixef, nsim = 2000,
                                                    seed = 632760,
                                                    use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_mid_size_small_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_mid_size_small_adj_lmer <- boot.ci(boot_mass_iqr_temp_mid_size_small_adj_lmer,
                                                     type = c('perc', 'norm', 'basic'),
                                                     index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_mid_size_small_adj_lmer)


# BIG nestlings - adjusted
mass_iqr_temp_mid_size_big_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) + 
                                          scale(nestling_number) + 
                                          scale(days_summer) +
                                          (1|fnest_id), 
                                        data = subset(late_nestling_parent_care,
                                                      mid_size_order == "other",
                                                      !is.na(x = mass_pre_obs) & 
                                                        !is.na(x = nest_iqr_temp) &
                                                        !is.na(x = mid_size_order)))


## Check diagnostics for the full model
plot(mass_iqr_temp_mid_size_big_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_mid_size_big_adj_lmer))
  qqline(resid(mass_iqr_temp_mid_size_big_adj_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_mid_size_big_adj_lmer))

summary(mass_iqr_temp_mid_size_big_adj_lmer)
confint(mass_iqr_temp_mid_size_big_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_mid_size_big_adj_lmer <- bootMer(x = mass_iqr_temp_mid_size_big_adj_lmer,
                                                      FUN = fixef, nsim = 2000,
                                                      seed = 632760,
                                                      use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_mid_size_big_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_mid_size_big_adj_lmer <- boot.ci(boot_mass_iqr_temp_mid_size_big_adj_lmer,
                                                       type = c('perc', 'norm', 'basic'),
                                                       index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_mid_size_big_adj_lmer)


######## Create model results tables for stratified relative size models ######

# Extract number of observations for unadjusted models
small_un_n <- c(nobs(mass_min_temp_mid_size_small_lmer), 
                nobs(mass_max_temp_mid_size_small_lmer),
                nobs(mass_iqr_temp_mid_size_small_lmer))

big_un_n <- c(nobs(mass_min_temp_mid_size_big_lmer), 
              nobs(mass_max_temp_mid_size_big_lmer),
              nobs(mass_iqr_temp_mid_size_big_lmer))


# Extract Betas for unadjusted models
small_un_b <- c(paste(round(fixef(mass_min_temp_mid_size_small_lmer)[2], 2), " (", 
                      round(bt_ci_mass_min_temp_mid_size_small_lmer[[6]][4], 2), ", ", 
                      round(bt_ci_mass_min_temp_mid_size_small_lmer[[6]][5], 2), ")",
                      sep = ""),
                paste(round(fixef(mass_max_temp_mid_size_small_lmer)[2], 2), " (", 
                      round(bt_ci_mass_max_temp_mid_size_small_lmer[[6]][4], 2), ", ", 
                      round(bt_ci_mass_max_temp_mid_size_small_lmer[[6]][5], 2), ")",
                      sep = ""),
                paste(round(fixef(mass_iqr_temp_mid_size_small_lmer)[2], 2), " (", 
                      round(bt_ci_mass_iqr_temp_mid_size_small_lmer[[6]][4], 2), ", ", 
                      round(bt_ci_mass_iqr_temp_mid_size_small_lmer[[6]][5], 2), ")",
                      sep = ""))


big_un_b <- c(paste(round(fixef(mass_min_temp_mid_size_big_lmer)[2], 2), " (", 
                    round(bt_ci_mass_min_temp_mid_size_big_lmer[[6]][4], 2), ", ", 
                    round(bt_ci_mass_min_temp_mid_size_big_lmer[[6]][5], 2), ")",
                    sep = ""),
              paste(round(fixef(mass_max_temp_mid_size_big_lmer)[2], 2), " (", 
                    round(bt_ci_mass_max_temp_mid_size_big_lmer[[6]][4], 2), ", ", 
                    round(bt_ci_mass_max_temp_mid_size_big_lmer[[6]][5], 2), ")",
                    sep = ""),
              paste(round(fixef(mass_iqr_temp_mid_size_big_lmer)[2], 2), " (", 
                    round(bt_ci_mass_iqr_temp_mid_size_big_lmer[[6]][4], 2), ", ", 
                    round(bt_ci_mass_iqr_temp_mid_size_big_lmer[[6]][5], 2), ")",
                    sep = ""))

# Extract number of observations for adjusted models
small_ad_n <- c(nobs(mass_min_temp_mid_size_small_adj_lmer), 
                nobs(mass_max_temp_mid_size_small_adj_lmer),
                nobs(mass_iqr_temp_mid_size_small_adj_lmer))

big_ad_n <- c(nobs(mass_min_temp_mid_size_big_adj_lmer), 
              nobs(mass_max_temp_mid_size_big_adj_lmer),
              nobs(mass_iqr_temp_mid_size_big_adj_lmer))

# Extract Betas for adjusted models
small_ad_b <- c(paste(round(fixef(mass_min_temp_mid_size_small_adj_lmer)[2], 2), " (", 
                      round(bt_ci_mass_min_temp_mid_size_small_adj_lmer[[6]][4], 2), ", ", 
                      round(bt_ci_mass_min_temp_mid_size_small_adj_lmer[[6]][5], 2), ")",
                      sep = ""),
                paste(round(fixef(mass_max_temp_mid_size_small_adj_lmer)[2], 2), " (", 
                      round(bt_ci_mass_max_temp_mid_size_small_adj_lmer[[6]][4], 2), ", ", 
                      round(bt_ci_mass_max_temp_mid_size_small_adj_lmer[[6]][5], 2), ")",
                      sep = ""),
                paste(round(fixef(mass_iqr_temp_mid_size_small_adj_lmer)[2], 2), " (", 
                      round(bt_ci_mass_iqr_temp_mid_size_small_adj_lmer[[6]][4], 2), ", ", 
                      round(bt_ci_mass_iqr_temp_mid_size_small_adj_lmer[[6]][5], 2), ")",
                      sep = ""))


big_ad_b <- c(paste(round(fixef(mass_min_temp_mid_size_big_adj_lmer)[2], 2), " (", 
                    round(bt_ci_mass_min_temp_mid_size_big_adj_lmer[[6]][4], 2), ", ", 
                    round(bt_ci_mass_min_temp_mid_size_big_adj_lmer[[6]][5], 2), ")",
                    sep = ""),
              paste(round(fixef(mass_max_temp_mid_size_big_adj_lmer)[2], 2), " (", 
                    round(bt_ci_mass_max_temp_mid_size_big_adj_lmer[[6]][4], 2), ", ", 
                    round(bt_ci_mass_max_temp_mid_size_big_adj_lmer[[6]][5], 2), ")",
                    sep = ""),
              paste(round(fixef(mass_iqr_temp_mid_size_big_adj_lmer)[2], 2), " (", 
                    round(bt_ci_mass_iqr_temp_mid_size_big_adj_lmer[[6]][4], 2), ", ", 
                    round(bt_ci_mass_iqr_temp_mid_size_big_adj_lmer[[6]][5], 2), ")",
                    sep = ""))

# Bind the values together into a single dataframe
size_strat_ad <- data.frame(
  variable <- c("Minimum temperature", "Maximum temperature", "Temperature variability"),
  type <- c("Adjusted", "Adjusted", "Adjusted"),
  small_ad_n, small_ad_b, 
  big_ad_n, big_ad_b
)

size_strat_un <- data.frame(
  variable <- c("Minimum temperature", "Maximum temperature", "Temperature variability"),
  type <- c("Unadjusted", "Unadjusted", "Unadjusted"),
  small_un_n, small_un_b, 
  big_un_n, big_un_b
)

colnames(size_strat_ad) <- c("variable", "type",
                             "small_n", "small_b", 
                             "big_n", "big_b")

colnames(size_strat_un) <- c("variable", "type",
                             "small_n", "small_b", 
                             "big_n", "big_b")

size_strat <- rbind(size_strat_un, size_strat_ad)

size_strat <- arrange(size_strat, factor(variable, levels = c("Minimum temperature", "Maximum temperature", "Temperature variability")))


str(size_strat)

# Create display table
size_strat_mod_table <- gt(size_strat, rowname_col = "type") %>%
  tab_header(
    title = md("**Supplemental Table 4.** Associations of nestling mass and temperature, assessed in separate models stratified by relative nestling size at mid development measure (smallest vs. other). Temperature variability is defined as the interquartile range. ")
  ) %>%
  tab_footnote(
    footnote = "Estimated β (95% CI) from stratified linear mixed models in which temperature is the explanatory variable of interest, nestling mass is the outcome of interest, and nest ID was included as a random intercept. Adjusted models include hatch date and number of nestlings in the nest. Continuous predictors are z-score standardized.",
    locations = cells_column_labels(columns = c(small_b, big_b))
  ) %>%
  tab_footnote(
    footnote = md(paste("R-squared for adjusted minimum temperature models. ",
                        "Small size model: Marginal R-squared = ", round(r.squaredGLMM(mass_min_temp_mid_size_small_adj_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_min_temp_mid_size_small_adj_lmer)[2], 2),
                        "; Other size model: Marginal R-squared = ", round(r.squaredGLMM(mass_min_temp_mid_size_big_adj_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_min_temp_mid_size_big_adj_lmer)[2], 2),
                        sep = "")),
    locations = cells_stub(rows = c(2))
  ) %>%
  tab_footnote(
    footnote = md(paste("R-squared for adjusted maximum temperature models. ",
                        "Small size model: Marginal R-squared = ", round(r.squaredGLMM(mass_max_temp_mid_size_small_adj_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_max_temp_mid_size_small_adj_lmer)[2], 2),
                        "; Other size model: Marginal R-squared = ", round(r.squaredGLMM(mass_max_temp_mid_size_big_adj_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_max_temp_mid_size_big_adj_lmer)[2], 2),
                        sep = "")),
    locations = cells_stub(rows = c(4))
  ) %>%  
  tab_footnote(
    footnote = md(paste("R-squared for adjusted temperature variability models. ",
                        "Small size model: Marginal R-squared = ", round(r.squaredGLMM(mass_iqr_temp_mid_size_small_adj_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_iqr_temp_mid_size_small_adj_lmer)[2], 2),
                        "; Other size model: Marginal R-squared = ", round(r.squaredGLMM(mass_iqr_temp_mid_size_big_adj_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_iqr_temp_mid_size_big_adj_lmer)[2], 2),
                        sep = "")),
    locations = cells_stub(rows = c(6))
  ) %>%
  tab_stubhead(
    label = md("Type")
  ) %>%
  tab_row_group(
    label = md("**Effect of temperature variability**"), 
    rows = c(5:6)
  ) %>%
  tab_row_group(
    label = md("**Effect of maximum temperature**"), 
    rows = c(3:4)
  ) %>%
  tab_row_group(
    label = md("**Effect of minimum temperature**"), 
    rows = c(1:2)
  ) %>%
  tab_spanner(
    label = "Small size models",
    columns = c(small_n, small_b)
  ) %>%
  tab_spanner(
    label = "Other size models",
    columns = c(big_n, big_b)
  )  %>% 
  cols_label(
    ends_with("n") ~ "N", 
    ends_with("b") ~ "β (95% CI)",
    type = "Type"
  ) %>%
  cols_hide(variable)  %>%
  opt_table_font(font = "Arial", size = 12)  %>%
  tab_options(footnotes.font.size = 10)


size_strat_mod_table
gtsave(size_strat_mod_table, filename = "Output/size_adjusted_strat.docx")
gtsave(size_strat_mod_table, filename = "Output/size_adjusted_strat.html")



###############################################################################
##############           Sensitive development periods           ##############
###############################################################################

## This section includes models with several different thresholds set for 
## development of thermoregulatory independence. The threshold used in the 
## main text is six days post-hatch. We also provide models here for thresholds
## of five and seven days

######################## 6 days threshold #####################################

####################### Minimum temperature ###################################

### Before model - unadjusted
mass_min_before_thermo_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp) +
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

# Calculate R squared 
r.squaredGLMM(mass_min_before_thermo_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_before_thermo_lmer <- bootMer(x = mass_min_before_thermo_lmer,
                                                FUN = fixef, nsim = 2000,
                                                seed = 632760,
                                                use.u = F, type = 'parametric')
tidy(boot_mass_min_before_thermo_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_before_thermo_lmer <- boot.ci(boot_mass_min_before_thermo_lmer,
                                                 type = c('perc', 'norm', 'basic'),
                                                 index = 2) # CI for 1st betas
print(bt_ci_mass_min_before_thermo_lmer)


### Before model -  adjusted
mass_min_before_thermo_adj_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp) + scale(nestling_number) +
                                      scale(days_summer) +
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

# Calculate R squared 
r.squaredGLMM(mass_min_before_thermo_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_before_thermo_adj_lmer <- bootMer(x = mass_min_before_thermo_adj_lmer,
                                            FUN = fixef, nsim = 2000,
                                            seed = 632760,
                                            use.u = F, type = 'parametric')
tidy(boot_mass_min_before_thermo_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_before_thermo_adj_lmer <- boot.ci(boot_mass_min_before_thermo_adj_lmer,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 2) # CI for 1st betas
print(bt_ci_mass_min_before_thermo_adj_lmer)


## Before model - unadjusted, outliers removed
mass_min_before_thermo_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp) +
                                      (1|fnest_id), 
                                    data = subset(before_thermo_outliers_removed,
                                                  !is.na(x = thermo_bef_min_temp)&
                                                    !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_before_thermo_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_before_thermo_noout_lmer))
  qqline(resid(mass_min_before_thermo_noout_lmer))}
# Histogram of residuals
hist(resid(mass_min_before_thermo_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_before_thermo_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_min_before_thermo_noout_lmer, vars=c("Studentized"))

summary(mass_min_before_thermo_noout_lmer)
confint(mass_min_before_thermo_noout_lmer)

# Calculate R squared 
r.squaredGLMM(mass_min_before_thermo_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_before_thermo_noout_lmer <- bootMer(x = mass_min_before_thermo_noout_lmer,
                                            FUN = fixef, nsim = 2000,
                                            seed = 632760,
                                            use.u = F, type = 'parametric')
tidy(boot_mass_min_before_thermo_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_before_thermo_noout_lmer<- boot.ci(boot_mass_min_before_thermo_noout_lmer,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 2) # CI for 1st betas
print(bt_ci_mass_min_before_thermo_noout_lmer)


## Before model - adjusted, outliers removed
mass_min_before_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp) + scale(nestling_number) +
                                          scale(days_summer) +
                                          (1|fnest_id), 
                                        data = subset(before_thermo_outliers_removed,
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

# Calculate R squared 
r.squaredGLMM(mass_min_before_thermo_adj_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_before_thermo_adj_noout_lmer <- bootMer(x = mass_min_before_thermo_adj_noout_lmer,
                                                FUN = fixef, nsim = 2000,
                                                seed = 632760,
                                                use.u = F, type = 'parametric')
tidy(boot_mass_min_before_thermo_adj_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_before_thermo_adj_noout_lmer <- boot.ci(boot_mass_min_before_thermo_adj_noout_lmer,
                                                 type = c('perc', 'norm', 'basic'),
                                                 index = 2) # CI for 1st betas
print(bt_ci_mass_min_before_thermo_adj_noout_lmer)


### After model - unadjusted
mass_min_after_thermo_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_min_temp) +
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

# Calculate R squared 
r.squaredGLMM(mass_min_after_thermo_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_after_thermo_lmer <- bootMer(x = mass_min_after_thermo_lmer,
                                                      FUN = fixef, nsim = 2000,
                                                      seed = 632760,
                                                      use.u = F, type = 'parametric')
tidy(boot_mass_min_after_thermo_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_after_thermo_lmer <- boot.ci(boot_mass_min_after_thermo_lmer,
                                                       type = c('perc', 'norm', 'basic'),
                                                       index = 2) # CI for 1st betas
print(bt_ci_mass_min_after_thermo_lmer)


### After model - adjusted
mass_min_after_thermo_adj_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_min_temp) + scale(nestling_number) +
                                         scale(days_summer) +
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

# Calculate R squared 
r.squaredGLMM(mass_min_after_thermo_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_after_thermo_adj_lmer <- bootMer(x = mass_min_after_thermo_adj_lmer,
                                           FUN = fixef, nsim = 2000,
                                           seed = 632760,
                                           use.u = F, type = 'parametric')
tidy(boot_mass_min_after_thermo_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_after_thermo_adj_lmer <- boot.ci(boot_mass_min_after_thermo_adj_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 2) # CI for 1st betas
print(bt_ci_mass_min_after_thermo_adj_lmer)


### After model - unadjusted, outliers removed
mass_min_after_thermo_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_min_temp) +
                                            (1|fnest_id), 
                                          data = subset(before_thermo_outliers_removed,
                                                        !is.na(x = thermo_aft_min_temp)&
                                                          !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_after_thermo_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_after_thermo_noout_lmer))
  qqline(resid(mass_min_after_thermo_noout_lmer))}
# Histogram of residuals
hist(resid(mass_min_after_thermo_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_after_thermo_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_min_after_thermo_noout_lmer, vars=c("Studentized"))

summary(mass_min_after_thermo_noout_lmer)
confint(mass_min_after_thermo_noout_lmer)

# Calculate R squared 
r.squaredGLMM(mass_min_after_thermo_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_after_thermo_noout_lmer <- bootMer(x = mass_min_after_thermo_noout_lmer,
                                                  FUN = fixef, nsim = 2000,
                                                  seed = 632760,
                                                  use.u = F, type = 'parametric')
tidy(boot_mass_min_after_thermo_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_after_thermo_noout_lmer <- boot.ci(boot_mass_min_after_thermo_noout_lmer,
                                                  type = c('perc', 'norm', 'basic'),
                                                  index = 2) # CI for 1st betas
print(bt_ci_mass_min_after_thermo_noout_lmer)


### After model - adjusted, outliers removed
mass_min_after_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_min_temp) + scale(nestling_number) +
                                               scale(days_summer) +
                                               (1|fnest_id), 
                                             data = subset(before_thermo_outliers_removed,
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

# Calculate R squared 
r.squaredGLMM(mass_min_after_thermo_adj_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_after_thermo_adj_noout_lmer <- bootMer(x = mass_min_after_thermo_adj_noout_lmer,
                                                 FUN = fixef, nsim = 2000,
                                                 seed = 632760,
                                                 use.u = F, type = 'parametric')
tidy(boot_mass_min_after_thermo_adj_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_after_thermo_adj_noout_lmer <- boot.ci(boot_mass_min_after_thermo_adj_noout_lmer,
                                                  type = c('perc', 'norm', 'basic'),
                                                  index = 2) # CI for 1st betas
print(bt_ci_mass_min_after_thermo_adj_noout_lmer)



########################## Maximum temperature ################################

### Before model - unadjusted
mass_max_before_thermo_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp) +
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

# Calculate R squared 
r.squaredGLMM(mass_max_before_thermo_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_before_thermo_lmer <- bootMer(x = mass_max_before_thermo_lmer,
                                                    FUN = fixef, nsim = 2000,
                                                    seed = 632760,
                                                    use.u = F, type = 'parametric')
tidy(boot_mass_max_before_thermo_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_before_thermo_lmer <- boot.ci(boot_mass_max_before_thermo_lmer,
                                                     type = c('perc', 'norm', 'basic'),
                                                     index = 2) # CI for 1st betas
print(bt_ci_mass_max_before_thermo_lmer)


### Before model - adjusted
mass_max_before_thermo_adj_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp) + scale(nestling_number) +
                                          scale(days_summer) +
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

# Calculate R squared 
r.squaredGLMM(mass_max_before_thermo_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_before_thermo_adj_lmer <- bootMer(x = mass_max_before_thermo_adj_lmer,
                                            FUN = fixef, nsim = 2000,
                                            seed = 632760,
                                            use.u = F, type = 'parametric')
tidy(boot_mass_max_before_thermo_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_before_thermo_adj_lmer <- boot.ci(boot_mass_max_before_thermo_adj_lmer,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 2) # CI for 1st betas
print(bt_ci_mass_max_before_thermo_adj_lmer)


### Before model - unadjusted, outliers removed
mass_max_before_thermo_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp) +
                                      (1|fnest_id), 
                                    data = subset(before_thermo_outliers_removed,
                                                  !is.na(x = thermo_bef_max_temp)&
                                                    !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_before_thermo_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_before_thermo_noout_lmer))
  qqline(resid(mass_max_before_thermo_noout_lmer))}
# Histogram of residuals
hist(resid(mass_max_before_thermo_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_before_thermo_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_max_before_thermo_noout_lmer, vars=c("Studentized"))

summary(mass_max_before_thermo_noout_lmer)
confint(mass_max_before_thermo_noout_lmer)

# Calculate R squared 
r.squaredGLMM(mass_max_before_thermo_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_before_thermo_noout_lmer <- bootMer(x = mass_max_before_thermo_noout_lmer,
                                            FUN = fixef, nsim = 2000,
                                            seed = 632760,
                                            use.u = F, type = 'parametric')
tidy(boot_mass_max_before_thermo_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_before_thermo_noout_lmer <- boot.ci(boot_mass_max_before_thermo_noout_lmer,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 2) # CI for 1st betas
print(bt_ci_mass_max_before_thermo_noout_lmer)


### Before model - adjusted, outliers removed
mass_max_before_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp) + scale(nestling_number) +
                                          scale(days_summer) +
                                          (1|fnest_id), 
                                        data = subset(before_thermo_outliers_removed,
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

# Calculate R squared 
r.squaredGLMM(mass_max_before_thermo_adj_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_before_thermo_adj_noout_lmer <- bootMer(x = mass_max_before_thermo_adj_noout_lmer,
                                                  FUN = fixef, nsim = 2000,
                                                  seed = 632760,
                                                  use.u = F, type = 'parametric')
tidy(boot_mass_max_before_thermo_adj_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_before_thermo_adj_noout_lmer <- boot.ci(boot_mass_max_before_thermo_adj_noout_lmer,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 2) # CI for 1st betas
print(bt_ci_mass_max_before_thermo_adj_noout_lmer)


### After model - unadjusted
mass_max_after_thermo_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_max_temp) +
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

# Calculate R squared 
r.squaredGLMM(mass_max_after_thermo_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_after_thermo_lmer <- bootMer(x = mass_max_after_thermo_lmer,
                                                      FUN = fixef, nsim = 2000,
                                                      seed = 632760,
                                                      use.u = F, type = 'parametric')
tidy(boot_mass_max_after_thermo_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_after_thermo_lmer <- boot.ci(boot_mass_max_after_thermo_lmer,
                                                       type = c('perc', 'norm', 'basic'),
                                                       index = 2) # CI for 1st betas
print(bt_ci_mass_max_after_thermo_lmer)


### After model - adjusted
mass_max_after_thermo_adj_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_max_temp) + scale(nestling_number) +
                                         scale(days_summer) +
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

# Calculate R squared 
r.squaredGLMM(mass_max_after_thermo_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_after_thermo_adj_lmer <- bootMer(x = mass_max_after_thermo_adj_lmer,
                                           FUN = fixef, nsim = 2000,
                                           seed = 632760,
                                           use.u = F, type = 'parametric')
tidy(boot_mass_max_after_thermo_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_after_thermo_adj_lmer <- boot.ci(boot_mass_max_after_thermo_adj_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 2) # CI for 1st betas
print(bt_ci_mass_max_after_thermo_adj_lmer)


### After model - unadjusted, outliers removed
mass_max_after_thermo_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_max_temp) +
                                     (1|fnest_id), 
                                   data = subset(before_thermo_outliers_removed,
                                                 !is.na(x = thermo_aft_max_temp)&
                                                   !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_after_thermo_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_after_thermo_noout_lmer))
  qqline(resid(mass_max_after_thermo_noout_lmer))}
# Histogram of residuals
hist(resid(mass_max_after_thermo_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_after_thermo_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_max_after_thermo_noout_lmer, vars=c("Studentized"))

summary(mass_max_after_thermo_noout_lmer)
confint(mass_max_after_thermo_noout_lmer)

# Calculate R squared 
r.squaredGLMM(mass_max_after_thermo_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_after_thermo_noout_lmer <- bootMer(x = mass_max_after_thermo_noout_lmer,
                                           FUN = fixef, nsim = 2000,
                                           seed = 632760,
                                           use.u = F, type = 'parametric')
tidy(boot_mass_max_after_thermo_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_after_thermo_noout_lmer <- boot.ci(boot_mass_max_after_thermo_noout_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 2) # CI for 1st betas
print(bt_ci_mass_max_after_thermo_noout_lmer)

### After model - adjusted, outliers removed
mass_max_after_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_max_temp) + scale(nestling_number) +
                                         scale(days_summer) +
                                         (1|fnest_id), 
                                       data = subset(before_thermo_outliers_removed,
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

# Calculate R squared 
r.squaredGLMM(mass_max_after_thermo_adj_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_after_thermo_adj_noout_lmer <- bootMer(x = mass_max_after_thermo_adj_noout_lmer,
                                                 FUN = fixef, nsim = 2000,
                                                 seed = 632760,
                                                 use.u = F, type = 'parametric')
tidy(boot_mass_max_after_thermo_adj_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_after_thermo_adj_noout_lmer <- boot.ci(boot_mass_max_after_thermo_adj_noout_lmer,
                                                  type = c('perc', 'norm', 'basic'),
                                                  index = 2) # CI for 1st betas
print(bt_ci_mass_max_after_thermo_adj_noout_lmer)


####################### Temperature interquartile range #######################
### Before model - unadjusted
mass_iqr_before_thermo_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp) + 
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

# Calculate R squared 
r.squaredGLMM(mass_iqr_before_thermo_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_before_thermo_lmer <- bootMer(x = mass_iqr_before_thermo_lmer,
                                            FUN = fixef, nsim = 2000,
                                            seed = 632760,
                                            use.u = F, type = 'parametric')
tidy(boot_mass_iqr_before_thermo_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_before_thermo_lmer <- boot.ci(boot_mass_iqr_before_thermo_lmer,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_before_thermo_lmer)


### Before model - adjusted
mass_iqr_before_thermo_adj_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp) + scale(nestling_number) +
                                          scale(days_summer) +
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

# Calculate R squared 
r.squaredGLMM(mass_iqr_before_thermo_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_before_thermo_adj_lmer <- bootMer(x = mass_iqr_before_thermo_adj_lmer,
                                            FUN = fixef, nsim = 2000,
                                            seed = 632760,
                                            use.u = F, type = 'parametric')
tidy(boot_mass_iqr_before_thermo_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_before_thermo_adj_lmer <- boot.ci(boot_mass_iqr_before_thermo_adj_lmer,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_before_thermo_adj_lmer)


### Before model - unadjusted, outliers removed
mass_iqr_before_thermo_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp) + 
                                      (1|fnest_id), 
                                    data = subset(before_thermo_outliers_removed,
                                                  !is.na(x = thermo_bef_iqr_temp)&
                                                    !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_before_thermo_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_before_thermo_noout_lmer))
  qqline(resid(mass_iqr_before_thermo_noout_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_before_thermo_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_before_thermo_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_before_thermo_noout_lmer, vars=c("Studentized"))

summary(mass_iqr_before_thermo_noout_lmer)
confint(mass_iqr_before_thermo_noout_lmer)

# Calculate R squared 
r.squaredGLMM(mass_iqr_before_thermo_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_before_thermo_noout_lmer <- bootMer(x = mass_iqr_before_thermo_noout_lmer,
                                            FUN = fixef, nsim = 2000,
                                            seed = 632760,
                                            use.u = F, type = 'parametric')
tidy(boot_mass_iqr_before_thermo_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_before_thermo_noout_lmer <- boot.ci(boot_mass_iqr_before_thermo_noout_lmer,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_before_thermo_noout_lmer)


### Before model - adjusted, outliers removed
mass_iqr_before_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp) + scale(nestling_number) +
                                          scale(days_summer) +  
                                          (1|fnest_id), 
                                        data = subset(before_thermo_outliers_removed,
                                                      !is.na(x = thermo_bef_iqr_temp)&
                                                        !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_before_thermo_adj_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_before_thermo_adj_noout_lmer))
  qqline(resid(mass_iqr_before_thermo_adj_noout_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_before_thermo_adj_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_before_thermo_adj_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_before_thermo_adj_noout_lmer, vars=c("Studentized"))

summary(mass_iqr_before_thermo_adj_noout_lmer)
confint(mass_iqr_before_thermo_adj_noout_lmer)

# Calculate R squared 
r.squaredGLMM(mass_iqr_before_thermo_adj_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_before_thermo_adj_noout_lmer <- bootMer(x = mass_iqr_before_thermo_adj_noout_lmer,
                                                FUN = fixef, nsim = 2000,
                                                seed = 632760,
                                                use.u = F, type = 'parametric')
tidy(boot_mass_iqr_before_thermo_adj_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_before_thermo_adj_noout_lmer <- boot.ci(boot_mass_iqr_before_thermo_adj_noout_lmer,
                                                 type = c('perc', 'norm', 'basic'),
                                                 index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_before_thermo_adj_noout_lmer)


### After model - unadjusted
mass_iqr_after_thermo_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_iqr_temp) + 
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

# Calculate R squared 
r.squaredGLMM(mass_iqr_after_thermo_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_after_thermo_lmer <- bootMer(x = mass_iqr_after_thermo_lmer,
                                                      FUN = fixef, nsim = 2000,
                                                      seed = 632760,
                                                      use.u = F, type = 'parametric')
tidy(boot_mass_iqr_after_thermo_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_after_thermo_lmer <- boot.ci(boot_mass_iqr_after_thermo_lmer,
                                                       type = c('perc', 'norm', 'basic'),
                                                       index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_after_thermo_lmer)


### After model - adjusted
mass_iqr_after_thermo_adj_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_iqr_temp) + scale(nestling_number) +
                                         scale(days_summer) + 
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

# Calculate R squared 
r.squaredGLMM(mass_iqr_after_thermo_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_after_thermo_adj_lmer <- bootMer(x = mass_iqr_after_thermo_adj_lmer,
                                           FUN = fixef, nsim = 2000,
                                           seed = 632760,
                                           use.u = F, type = 'parametric')
tidy(boot_mass_iqr_after_thermo_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_after_thermo_adj_lmer <- boot.ci(boot_mass_iqr_after_thermo_adj_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_after_thermo_adj_lmer)

### After model - unadjusted, outliers removed
mass_iqr_after_thermo_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_iqr_temp) + 
                                     (1|fnest_id), 
                                   data = subset(before_thermo_outliers_removed,
                                                 !is.na(x = thermo_aft_iqr_temp)&
                                                   !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_after_thermo_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_after_thermo_noout_lmer))
  qqline(resid(mass_iqr_after_thermo_noout_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_after_thermo_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_after_thermo_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_after_thermo_noout_lmer, vars=c("Studentized"))

summary(mass_iqr_after_thermo_noout_lmer)
confint(mass_iqr_after_thermo_noout_lmer)

# Calculate R squared 
r.squaredGLMM(mass_iqr_after_thermo_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_after_thermo_noout_lmer <- bootMer(x = mass_iqr_after_thermo_noout_lmer,
                                           FUN = fixef, nsim = 2000,
                                           seed = 632760,
                                           use.u = F, type = 'parametric')
tidy(boot_mass_iqr_after_thermo_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_after_thermo_noout_lmer <- boot.ci(boot_mass_iqr_after_thermo_noout_lmer,
                                            type = c('perc', 'norm', 'basic'),
                                            index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_after_thermo_noout_lmer)


### After model - adjusted, outliers removed
mass_iqr_after_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_iqr_temp) + scale(nestling_number) +
                                         scale(days_summer) + 
                                         (1|fnest_id), 
                                       data = subset(mass_outliers_removed,
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

# Calculate R squared 
r.squaredGLMM(mass_iqr_after_thermo_adj_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_after_thermo_adj_noout_lmer <- bootMer(x = mass_iqr_after_thermo_adj_noout_lmer,
                                                 FUN = fixef, nsim = 2000,
                                                 seed = 632760,
                                                 use.u = F, type = 'parametric')
tidy(boot_mass_iqr_after_thermo_adj_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_after_thermo_adj_noout_lmer <- boot.ci(boot_mass_iqr_after_thermo_adj_noout_lmer,
                                                  type = c('perc', 'norm', 'basic'),
                                                  index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_after_thermo_adj_noout_lmer)



########## Create model results tables for developmental stages models ########

# Extract number of observations for unadjusted models
before_un_n <- c(nobs(mass_min_before_thermo_lmer), 
                nobs(mass_max_before_thermo_lmer),
                nobs(mass_iqr_before_thermo_lmer))

after_un_n <- c(nobs(mass_min_after_thermo_lmer), 
                 nobs(mass_max_after_thermo_lmer),
                 nobs(mass_iqr_after_thermo_lmer))


# Extract Betas for unadjusted models
before_un_b <- c(paste(round(fixef(mass_min_before_thermo_lmer)[2], 2), " (", 
                      round(bt_ci_mass_min_before_thermo_lmer[[6]][4], 2), ", ", 
                      round(bt_ci_mass_min_before_thermo_lmer[[6]][5], 2), ")",
                      sep = ""),
                paste(round(fixef(mass_max_before_thermo_lmer)[2], 2), " (", 
                      round(bt_ci_mass_max_before_thermo_lmer[[6]][4], 2), ", ", 
                      round(bt_ci_mass_max_before_thermo_lmer[[6]][5], 2), ")",
                      sep = ""),
                paste(round(fixef(mass_iqr_before_thermo_lmer)[2], 2), " (", 
                      round(bt_ci_mass_iqr_before_thermo_lmer[[6]][4], 2), ", ", 
                      round(bt_ci_mass_iqr_before_thermo_lmer[[6]][5], 2), ")",
                      sep = ""))


after_un_b <- c(paste(round(fixef(mass_min_after_thermo_lmer)[2], 2), " (", 
                       round(bt_ci_mass_min_after_thermo_lmer[[6]][4], 2), ", ", 
                       round(bt_ci_mass_min_after_thermo_lmer[[6]][5], 2), ")",
                       sep = ""),
                 paste(round(fixef(mass_max_after_thermo_lmer)[2], 2), " (", 
                       round(bt_ci_mass_max_after_thermo_lmer[[6]][4], 2), ", ", 
                       round(bt_ci_mass_max_after_thermo_lmer[[6]][5], 2), ")",
                       sep = ""),
                 paste(round(fixef(mass_iqr_after_thermo_lmer)[2], 2), " (", 
                       round(bt_ci_mass_iqr_after_thermo_lmer[[6]][4], 2), ", ", 
                       round(bt_ci_mass_iqr_after_thermo_lmer[[6]][5], 2), ")",
                       sep = ""))



# Extract number of observations for adjusted models
before_ad_n <- c(nobs(mass_min_before_thermo_adj_lmer), 
                 nobs(mass_max_before_thermo_adj_lmer),
                 nobs(mass_iqr_before_thermo_adj_lmer))

after_ad_n <- c(nobs(mass_min_after_thermo_adj_lmer), 
                nobs(mass_max_after_thermo_adj_lmer),
                nobs(mass_iqr_after_thermo_adj_lmer))



# Extract Betas for unadjusted models
before_ad_b <- c(paste(round(fixef(mass_min_before_thermo_adj_lmer)[2], 2), " (", 
                       round(bt_ci_mass_min_before_thermo_adj_lmer[[6]][4], 2), ", ", 
                       round(bt_ci_mass_min_before_thermo_adj_lmer[[6]][5], 2), ")",
                       sep = ""),
                 paste(round(fixef(mass_max_before_thermo_adj_lmer)[2], 2), " (", 
                       round(bt_ci_mass_max_before_thermo_adj_lmer[[6]][4], 2), ", ", 
                       round(bt_ci_mass_max_before_thermo_adj_lmer[[6]][5], 2), ")",
                       sep = ""),
                 paste(round(fixef(mass_iqr_before_thermo_adj_lmer)[2], 2), " (", 
                       round(bt_ci_mass_iqr_before_thermo_adj_lmer[[6]][4], 2), ", ", 
                       round(bt_ci_mass_iqr_before_thermo_adj_lmer[[6]][5], 2), ")",
                       sep = ""))


after_ad_b <- c(paste(round(fixef(mass_min_after_thermo_adj_lmer)[2], 2), " (", 
                      round(bt_ci_mass_min_after_thermo_adj_lmer[[6]][4], 2), ", ", 
                      round(bt_ci_mass_min_after_thermo_adj_lmer[[6]][5], 2), ")",
                      sep = ""),
                paste(round(fixef(mass_max_after_thermo_adj_lmer)[2], 2), " (", 
                      round(bt_ci_mass_max_after_thermo_adj_lmer[[6]][4], 2), ", ", 
                      round(bt_ci_mass_max_after_thermo_adj_lmer[[6]][5], 2), ")",
                      sep = ""),
                paste(round(fixef(mass_iqr_after_thermo_adj_lmer)[2], 2), " (", 
                      round(bt_ci_mass_iqr_after_thermo_adj_lmer[[6]][4], 2), ", ", 
                      round(bt_ci_mass_iqr_after_thermo_adj_lmer[[6]][5], 2), ")",
                      sep = ""))



# Bind the values together into a single dataframe
devel_ad <- data.frame(
  variable <- c("Minimum temperature", "Maximum temperature", "Temperature variability"),
  type <- c("Adjusted", "Adjusted", "Adjusted"),
  before_ad_n, before_ad_b, 
  after_ad_n, after_ad_b
)

devel_un <- data.frame(
  variable <- c("Minimum temperature", "Maximum temperature", "Temperature variability"),
  type <- c("Unadjusted", "Unadjusted", "Unadjusted"),
  before_un_n, before_un_b, 
  after_un_n, after_un_b
)

colnames(devel_ad) <- c("variable", "type",
                                      "before_n", "before_b", 
                                      "after_n", "after_b")

colnames(devel_un) <- c("variable", "type",
                                      "before_n", "before_b", 
                                      "after_n", "after_b")

devel <- rbind(devel_un, devel_ad)

devel <- arrange(devel,
                factor(variable, levels = c("Minimum temperature", "Maximum temperature", "Temperature variability")))


str(devel)

# Create display table
devel_mod_table <- gt(devel, rowname_col = "type") %>%
  tab_header(
    title = md("**Supplemental Table 3.** Associations of three temperature variables before or after thermoregulatory independence (defined as occurring six days post-hatch) with nestling mass. Temperature variability is defined as the interquartile range. ")
  ) %>%
  tab_footnote(
    footnote = "Estimated β (95% CI) from stratified linear mixed models in which temperature before or after thermoregulatory independence are the explanatory variables of interest, nestling mass is the outcome of interest, and nest ID was included as a random intercept. Adjusted models include hatch date and number of nestlings in the nest. Continuous predictors as z-score standardized.",
    locations = cells_column_labels(columns = c(before_b, after_b))
  )  %>%
  tab_footnote(
    footnote = md(paste("R-squared for adjusted minimum temperature models. ",
                        "Before model: Marginal R-squared = ", round(r.squaredGLMM(mass_min_before_thermo_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_min_before_thermo_lmer)[2], 2),
                        "; After model: Marginal R-squared = ", round(r.squaredGLMM(mass_min_after_thermo_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_min_after_thermo_lmer)[2], 2),
                        sep = "")),
    locations = cells_stub(rows = c(2))
  ) %>%
  tab_footnote(
    footnote = md(paste("R-squared for adjusted maximum temperature models. ",
                        "Before model: Marginal R-squared = ", round(r.squaredGLMM(mass_max_before_thermo_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_max_before_thermo_lmer)[2], 2),
                        "; After model: Marginal R-squared = ", round(r.squaredGLMM(mass_max_after_thermo_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_max_after_thermo_lmer)[2], 2),
                        sep = "")),
    locations = cells_stub(rows = c(4))
  ) %>%
  tab_footnote(
    footnote = md(paste("R-squared for adjusted temperature variability models. ",
                        "Before model: Marginal R-squared = ", round(r.squaredGLMM(mass_iqr_before_thermo_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_iqr_before_thermo_lmer)[2], 2),
                        "; After model: Marginal R-squared = ", round(r.squaredGLMM(mass_iqr_after_thermo_lmer)[1], 2),
                        ", Conditional R-squared = ", round(r.squaredGLMM(mass_iqr_after_thermo_lmer)[2], 2),
                        sep = "")),
    locations = cells_stub(rows = c(6))
  ) %>%
  tab_stubhead(
    label = md("Type")
  ) %>%
  tab_row_group(
    label = md("**Effect of temperature variability**"), 
    rows = c(5:6)
  ) %>%
  tab_row_group(
    label = md("**Effect of maximum temperature**"), 
    rows = c(3:4)
  ) %>%
  tab_row_group(
    label = md("**Effect of minimum temperature**"), 
    rows = c(1:2)
  ) %>%
  tab_spanner(
    label = "Before models",
    columns = c(before_n, before_b)
  ) %>%
  tab_spanner(
    label = "After models",
    columns = c(after_n, after_b)
  ) %>%
  cols_label(
    ends_with("n") ~ "N", 
    ends_with("b") ~ "β (95% CI)",
    type = "Type"
  ) %>%
  cols_hide(variable)  %>%
  opt_table_font(font = "Arial", size = 12)  %>%
  tab_options(footnotes.font.size = 10)


devel_mod_table
gtsave(devel_mod_table, filename = "Output/devel_mod_table.docx")
gtsave(devel_mod_table, filename = "Output/devel_mod_table.html")


############################ 5 days threshold #################################

########################## Minimum temperature ################################

### Before model - unadjusted
mass_min_before_thermo_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp_5) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  !is.na(x = thermo_bef_min_temp_5)&
                                                    !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_before_thermo_lmer_5)
# Normal QQplot
{qqnorm(resid(mass_min_before_thermo_lmer_5))
  qqline(resid(mass_min_before_thermo_lmer_5))}
# Histogram of residuals
hist(resid(mass_min_before_thermo_lmer_5))
# Checking for influential outliers
infIndexPlot(mass_min_before_thermo_lmer_5, vars=c("Cook"))
infIndexPlot(mass_min_before_thermo_lmer_5, vars=c("Studentized"))

summary(mass_min_before_thermo_lmer_5)
confint(mass_min_before_thermo_lmer_5)

# Calculate R squared 
r.squaredGLMM(mass_min_before_thermo_lmer_5)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_before_thermo_lmer_5 <- bootMer(x = mass_min_before_thermo_lmer_5,
                                            FUN = fixef, nsim = 2000,
                                            seed = 632760,
                                            use.u = F, type = 'parametric')
tidy(boot_mass_min_before_thermo_lmer_5) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_before_thermo_lmer_5 <- boot.ci(boot_mass_min_before_thermo_lmer_5,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 2) # CI for 1st betas
print(bt_ci_mass_min_before_thermo_lmer_5)


### Before model - adjusted
mass_min_before_thermo_adj_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp_5) + scale(nestling_number) +
                                          scale(days_summer) +
                                          (1|fnest_id), 
                                        data = subset(late_nestling_parent_care,
                                                      !is.na(x = thermo_bef_min_temp_5)&
                                                        !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_before_thermo_adj_lmer_5)
# Normal QQplot
{qqnorm(resid(mass_min_before_thermo_adj_lmer_5))
  qqline(resid(mass_min_before_thermo_adj_lmer_5))}
# Histogram of residuals
hist(resid(mass_min_before_thermo_adj_lmer_5))
# Checking for influential outliers
infIndexPlot(mass_min_before_thermo_adj_lmer_5, vars=c("Cook"))
infIndexPlot(mass_min_before_thermo_adj_lmer_5, vars=c("Studentized"))

summary(mass_min_before_thermo_adj_lmer_5)
confint(mass_min_before_thermo_adj_lmer_5)

# Calculate R squared 
r.squaredGLMM(mass_min_before_thermo_adj_lmer_5)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_before_thermo_adj_lmer_5 <- bootMer(x = mass_min_before_thermo_adj_lmer_5,
                                              FUN = fixef, nsim = 2000,
                                              seed = 632760,
                                              use.u = F, type = 'parametric')
tidy(boot_mass_min_before_thermo_adj_lmer_5) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_before_thermo_adj_lmer_5 <- boot.ci(boot_mass_min_before_thermo_adj_lmer_5,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 2) # CI for 1st betas
print(bt_ci_mass_min_before_thermo_adj_lmer_5)


### After model - unadjusted
mass_min_after_thermo_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_aft_min_temp_5) +
                                     (1|fnest_id), 
                                   data = subset(late_nestling_parent_care,
                                                 !is.na(x = thermo_aft_min_temp_5)&
                                                   !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_after_thermo_lmer_5)
# Normal QQplot
{qqnorm(resid(mass_min_after_thermo_lmer_5))
  qqline(resid(mass_min_after_thermo_lmer_5))}
# Histogram of residuals
hist(resid(mass_min_after_thermo_lmer_5))
# Checking for influential outliers
infIndexPlot(mass_min_after_thermo_lmer_5, vars=c("Cook"))
infIndexPlot(mass_min_after_thermo_lmer_5, vars=c("Studentized"))

summary(mass_min_after_thermo_lmer_5)
confint(mass_min_after_thermo_lmer_5)

# Calculate R squared 
r.squaredGLMM(mass_min_after_thermo_lmer_5)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_after_thermo_lmer_5 <- bootMer(x = mass_min_after_thermo_lmer_5,
                                                  FUN = fixef, nsim = 2000,
                                                  seed = 632760,
                                                  use.u = F, type = 'parametric')
tidy(boot_mass_min_after_thermo_lmer_5) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_after_thermo_lmer_5 <- boot.ci(boot_mass_min_after_thermo_lmer_5,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 2) # CI for 1st betas
print(bt_ci_mass_min_after_thermo_lmer_5)


### After model - adjusted
mass_min_after_thermo_adj_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_aft_min_temp_5) + scale(nestling_number) +
                                         scale(days_summer) +
                                         (1|fnest_id), 
                                       data = subset(late_nestling_parent_care,
                                                     !is.na(x = thermo_aft_min_temp_5)&
                                                       !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_after_thermo_adj_lmer_5)
# Normal QQplot
{qqnorm(resid(mass_min_after_thermo_adj_lmer_5))
  qqline(resid(mass_min_after_thermo_adj_lmer_5))}
# Histogram of residuals
hist(resid(mass_min_after_thermo_adj_lmer_5))
# Checking for influential outliers
infIndexPlot(mass_min_after_thermo_adj_lmer_5, vars=c("Cook"))
infIndexPlot(mass_min_after_thermo_adj_lmer_5, vars=c("Studentized"))

summary(mass_min_after_thermo_adj_lmer_5)
confint(mass_min_after_thermo_adj_lmer_5)

# Calculate R squared 
r.squaredGLMM(mass_min_after_thermo_adj_lmer_5)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_after_thermo_adj_lmer_5 <- bootMer(x = mass_min_after_thermo_adj_lmer_5,
                                             FUN = fixef, nsim = 2000,
                                             seed = 632760,
                                             use.u = F, type = 'parametric')
tidy(boot_mass_min_after_thermo_adj_lmer_5) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_after_thermo_adj_lmer_5 <- boot.ci(boot_mass_min_after_thermo_adj_lmer_5,
                                              type = c('perc', 'norm', 'basic'),
                                              index = 2) # CI for 1st betas
print(bt_ci_mass_min_after_thermo_adj_lmer_5)


########################### Maximum temperature ###############################
### Before model - unadjusted
mass_max_before_thermo_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp_5) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  !is.na(x = thermo_bef_max_temp_5)&
                                                    !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_before_thermo_lmer_5)
# Normal QQplot
{qqnorm(resid(mass_max_before_thermo_lmer_5))
  qqline(resid(mass_max_before_thermo_lmer_5))}
# Histogram of residuals
hist(resid(mass_max_before_thermo_lmer_5))
# Checking for influential outliers
infIndexPlot(mass_max_before_thermo_lmer_5, vars=c("Cook"))
infIndexPlot(mass_max_before_thermo_lmer_5, vars=c("Studentized"))

summary(mass_max_before_thermo_lmer_5)
confint(mass_max_before_thermo_lmer_5)

# Calculate R squared 
r.squaredGLMM(mass_max_before_thermo_lmer_5)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_before_thermo_lmer_5 <- bootMer(x = mass_max_before_thermo_lmer_5,
                                             FUN = fixef, nsim = 2000,
                                             seed = 632760,
                                             use.u = F, type = 'parametric')
tidy(boot_mass_max_before_thermo_lmer_5) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_before_thermo_lmer_5 <- boot.ci(boot_mass_max_before_thermo_lmer_5,
                                              type = c('perc', 'norm', 'basic'),
                                              index = 2) # CI for 1st betas
print(bt_ci_mass_max_before_thermo_lmer_5)

### Before model - adjusted
mass_max_before_thermo_adj_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp_5) + scale(nestling_number) +
                                          scale(days_summer) +
                                          (1|fnest_id), 
                                        data = subset(late_nestling_parent_care,
                                                      !is.na(x = thermo_bef_max_temp_5)&
                                                        !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_before_thermo_adj_lmer_5)
# Normal QQplot
{qqnorm(resid(mass_max_before_thermo_adj_lmer_5))
  qqline(resid(mass_max_before_thermo_adj_lmer_5))}
# Histogram of residuals
hist(resid(mass_max_before_thermo_adj_lmer_5))
# Checking for influential outliers
infIndexPlot(mass_max_before_thermo_adj_lmer_5, vars=c("Cook"))
infIndexPlot(mass_max_before_thermo_adj_lmer_5, vars=c("Studentized"))

summary(mass_max_before_thermo_adj_lmer_5)
confint(mass_max_before_thermo_adj_lmer_5)

# Calculate R squared 
r.squaredGLMM(mass_max_before_thermo_adj_lmer_5)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_before_thermo_adj_lmer_5 <- bootMer(x = mass_max_before_thermo_adj_lmer_5,
                                              FUN = fixef, nsim = 2000,
                                              seed = 632760,
                                              use.u = F, type = 'parametric')
tidy(boot_mass_max_before_thermo_adj_lmer_5) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_before_thermo_adj_lmer_5 <- boot.ci(boot_mass_max_before_thermo_adj_lmer_5,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 2) # CI for 1st betas
print(bt_ci_mass_max_before_thermo_adj_lmer_5)


### After model - unadjusted
mass_max_after_thermo_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_aft_max_temp_5) +
                                     (1|fnest_id), 
                                   data = subset(late_nestling_parent_care,
                                                 !is.na(x = thermo_aft_max_temp_5)&
                                                   !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_after_thermo_lmer_5)
# Normal QQplot
{qqnorm(resid(mass_max_after_thermo_lmer_5))
  qqline(resid(mass_max_after_thermo_lmer_5))}
# Histogram of residuals
hist(resid(mass_max_after_thermo_lmer_5))
# Checking for influential outliers
infIndexPlot(mass_max_after_thermo_lmer_5, vars=c("Cook"))
infIndexPlot(mass_max_after_thermo_lmer_5, vars=c("Studentized"))

summary(mass_max_after_thermo_lmer_5)
confint(mass_max_after_thermo_lmer_5)

# Calculate R squared 
r.squaredGLMM(mass_max_after_thermo_lmer_5)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_after_thermo_lmer_5 <- bootMer(x = mass_max_after_thermo_lmer_5,
                                                  FUN = fixef, nsim = 2000,
                                                  seed = 632760,
                                                  use.u = F, type = 'parametric')
tidy(boot_mass_max_after_thermo_lmer_5) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_after_thermo_lmer_5 <- boot.ci(boot_mass_max_after_thermo_lmer_5,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 2) # CI for 1st betas
print(bt_ci_mass_max_after_thermo_lmer_5)


### After model - adjusted
mass_max_after_thermo_adj_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_aft_max_temp_5) + scale(nestling_number) +
                                         scale(days_summer) +
                                         (1|fnest_id), 
                                       data = subset(late_nestling_parent_care,
                                                     !is.na(x = thermo_aft_max_temp_5)&
                                                       !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_after_thermo_adj_lmer_5)
# Normal QQplot
{qqnorm(resid(mass_max_after_thermo_adj_lmer_5))
  qqline(resid(mass_max_after_thermo_adj_lmer_5))}
# Histogram of residuals
hist(resid(mass_max_after_thermo_adj_lmer_5))
# Checking for influential outliers
infIndexPlot(mass_max_after_thermo_adj_lmer_5, vars=c("Cook"))
infIndexPlot(mass_max_after_thermo_adj_lmer_5, vars=c("Studentized"))

summary(mass_max_after_thermo_adj_lmer_5)
confint(mass_max_after_thermo_adj_lmer_5)

# Calculate R squared 
r.squaredGLMM(mass_max_after_thermo_adj_lmer_5)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_after_thermo_adj_lmer_5 <- bootMer(x = mass_max_after_thermo_adj_lmer_5,
                                             FUN = fixef, nsim = 2000,
                                             seed = 632760,
                                             use.u = F, type = 'parametric')
tidy(boot_mass_max_after_thermo_adj_lmer_5) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_after_thermo_adj_lmer_5 <- boot.ci(boot_mass_max_after_thermo_adj_lmer_5,
                                              type = c('perc', 'norm', 'basic'),
                                              index = 2) # CI for 1st betas
print(bt_ci_mass_max_after_thermo_adj_lmer_5)



######################## Temperature interquartile range ######################

### Before model - unadjusted
mass_iqr_before_thermo_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp_5) + 
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  !is.na(x = thermo_bef_iqr_temp_5)&
                                                    !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_before_thermo_lmer_5)
# Normal QQplot
{qqnorm(resid(mass_iqr_before_thermo_lmer_5))
  qqline(resid(mass_iqr_before_thermo_lmer_5))}
# Histogram of residuals
hist(resid(mass_iqr_before_thermo_lmer_5))
# Checking for influential outliers
infIndexPlot(mass_iqr_before_thermo_lmer_5, vars=c("Cook"))
infIndexPlot(mass_iqr_before_thermo_lmer_5, vars=c("Studentized"))

summary(mass_iqr_before_thermo_lmer_5)
confint(mass_iqr_before_thermo_lmer_5)

# Calculate R squared 
r.squaredGLMM(mass_iqr_before_thermo_lmer_5)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_before_thermo_lmer_5 <- bootMer(x = mass_iqr_before_thermo_lmer_5,
                                             FUN = fixef, nsim = 2000,
                                             seed = 632760,
                                             use.u = F, type = 'parametric')
tidy(boot_mass_iqr_before_thermo_lmer_5) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_before_thermo_lmer_5 <- boot.ci(boot_mass_iqr_before_thermo_lmer_5,
                                              type = c('perc', 'norm', 'basic'),
                                              index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_before_thermo_lmer_5)


### Before model - adjusted
mass_iqr_before_thermo_adj_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp_5) + scale(nestling_number) +
                                          scale(days_summer) +  
                                          (1|fnest_id), 
                                        data = subset(late_nestling_parent_care,
                                                      !is.na(x = thermo_bef_iqr_temp_5)&
                                                        !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_before_thermo_adj_lmer_5)
# Normal QQplot
{qqnorm(resid(mass_iqr_before_thermo_adj_lmer_5))
  qqline(resid(mass_iqr_before_thermo_adj_lmer_5))}
# Histogram of residuals
hist(resid(mass_iqr_before_thermo_adj_lmer_5))
# Checking for influential outliers
infIndexPlot(mass_iqr_before_thermo_adj_lmer_5, vars=c("Cook"))
infIndexPlot(mass_iqr_before_thermo_adj_lmer_5, vars=c("Studentized"))

summary(mass_iqr_before_thermo_adj_lmer_5)
confint(mass_iqr_before_thermo_adj_lmer_5)

# Calculate R squared 
r.squaredGLMM(mass_iqr_before_thermo_adj_lmer_5)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_before_thermo_adj_lmer_5 <- bootMer(x = mass_iqr_before_thermo_adj_lmer_5,
                                              FUN = fixef, nsim = 2000,
                                              seed = 632760,
                                              use.u = F, type = 'parametric')
tidy(boot_mass_iqr_before_thermo_adj_lmer_5) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_before_thermo_adj_lmer_5 <- boot.ci(boot_mass_iqr_before_thermo_adj_lmer_5,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_before_thermo_adj_lmer_5)


### After model - unadjusted
mass_iqr_after_thermo_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_aft_iqr_temp_5) + 
                                     (1|fnest_id), 
                                   data = subset(late_nestling_parent_care,
                                                 !is.na(x = thermo_aft_iqr_temp_5)&
                                                   !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_after_thermo_lmer_5)
# Normal QQplot
{qqnorm(resid(mass_iqr_after_thermo_lmer_5))
  qqline(resid(mass_iqr_after_thermo_lmer_5))}
# Histogram of residuals
hist(resid(mass_iqr_after_thermo_lmer_5))
# Checking for influential outliers
infIndexPlot(mass_iqr_after_thermo_lmer_5, vars=c("Cook"))
infIndexPlot(mass_iqr_after_thermo_lmer_5, vars=c("Studentized"))

summary(mass_iqr_after_thermo_lmer_5)
confint(mass_iqr_after_thermo_lmer_5)

# Calculate R squared 
r.squaredGLMM(mass_iqr_after_thermo_lmer_5)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_after_thermo_lmer_5 <- bootMer(x = mass_iqr_after_thermo_lmer_5,
                                                  FUN = fixef, nsim = 2000,
                                                  seed = 632760,
                                                  use.u = F, type = 'parametric')
tidy(boot_mass_iqr_after_thermo_lmer_5) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_after_thermo_lmer_5 <- boot.ci(boot_mass_iqr_after_thermo_lmer_5,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_after_thermo_lmer_5)


### After model - adjusted
mass_iqr_after_thermo_adj_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_aft_iqr_temp_5) + scale(nestling_number) +
                                         scale(days_summer) + 
                                         (1|fnest_id), 
                                       data = subset(late_nestling_parent_care,
                                                     !is.na(x = thermo_aft_iqr_temp_5)&
                                                       !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_after_thermo_adj_lmer_5)
# Normal QQplot
{qqnorm(resid(mass_iqr_after_thermo_adj_lmer_5))
  qqline(resid(mass_iqr_after_thermo_adj_lmer_5))}
# Histogram of residuals
hist(resid(mass_iqr_after_thermo_adj_lmer_5))
# Checking for influential outliers
infIndexPlot(mass_iqr_after_thermo_adj_lmer_5, vars=c("Cook"))
infIndexPlot(mass_iqr_after_thermo_adj_lmer_5, vars=c("Studentized"))

summary(mass_iqr_after_thermo_adj_lmer_5)
confint(mass_iqr_after_thermo_adj_lmer_5)

# Calculate R squared 
r.squaredGLMM(mass_iqr_after_thermo_adj_lmer_5)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_after_thermo_adj_lmer_5 <- bootMer(x = mass_iqr_after_thermo_adj_lmer_5,
                                             FUN = fixef, nsim = 2000,
                                             seed = 632760,
                                             use.u = F, type = 'parametric')
tidy(boot_mass_iqr_after_thermo_adj_lmer_5) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_after_thermo_adj_lmer_5 <- boot.ci(boot_mass_iqr_after_thermo_adj_lmer_5,
                                              type = c('perc', 'norm', 'basic'),
                                              index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_after_thermo_adj_lmer_5)



############################ 7 days threshold #################################

########################### Minimum temperature ###############################

### Before model - unadjusted
mass_min_before_thermo_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp_7) +
                                        (1|fnest_id), 
                                      data = subset(late_nestling_parent_care,
                                                    !is.na(x = thermo_bef_min_temp_7)&
                                                      !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_before_thermo_lmer_7)
# Normal QQplot
{qqnorm(resid(mass_min_before_thermo_lmer_7))
  qqline(resid(mass_min_before_thermo_lmer_7))}
# Histogram of residuals
hist(resid(mass_min_before_thermo_lmer_7))
# Checking for influential outliers
infIndexPlot(mass_min_before_thermo_lmer_7, vars=c("Cook"))
infIndexPlot(mass_min_before_thermo_lmer_7, vars=c("Studentized"))

summary(mass_min_before_thermo_lmer_7)
confint(mass_min_before_thermo_lmer_7)

# Calculate R squared 
r.squaredGLMM(mass_min_before_thermo_lmer_7)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_before_thermo_lmer_7 <- bootMer(x = mass_min_before_thermo_lmer_7,
                                              FUN = fixef, nsim = 2000,
                                              seed = 632760,
                                              use.u = F, type = 'parametric')
tidy(boot_mass_min_before_thermo_lmer_7) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_before_thermo_lmer_7 <- boot.ci(boot_mass_min_before_thermo_lmer_7,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 2) # CI for 1st betas
print(bt_ci_mass_min_before_thermo_lmer_7)


### Before model - adjusted
mass_min_before_thermo_adj_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp_7) + scale(nestling_number) +
                                            scale(days_summer) +
                                            (1|fnest_id), 
                                          data = subset(late_nestling_parent_care,
                                                        !is.na(x = thermo_bef_min_temp_7)&
                                                          !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_before_thermo_adj_lmer_7)
# Normal QQplot
{qqnorm(resid(mass_min_before_thermo_adj_lmer_7))
  qqline(resid(mass_min_before_thermo_adj_lmer_7))}
# Histogram of residuals
hist(resid(mass_min_before_thermo_adj_lmer_7))
# Checking for influential outliers
infIndexPlot(mass_min_before_thermo_adj_lmer_7, vars=c("Cook"))
infIndexPlot(mass_min_before_thermo_adj_lmer_7, vars=c("Studentized"))

summary(mass_min_before_thermo_adj_lmer_7)
confint(mass_min_before_thermo_adj_lmer_7)

# Calculate R squared 
r.squaredGLMM(mass_min_before_thermo_adj_lmer_7)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_before_thermo_adj_lmer_7 <- bootMer(x = mass_min_before_thermo_adj_lmer_7,
                                                  FUN = fixef, nsim = 2000,
                                                  seed = 632760,
                                                  use.u = F, type = 'parametric')
tidy(boot_mass_min_before_thermo_adj_lmer_7) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_before_thermo_adj_lmer_7 <- boot.ci(boot_mass_min_before_thermo_adj_lmer_7,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 2) # CI for 1st betas
print(bt_ci_mass_min_before_thermo_adj_lmer_7)


### After model - unadjusted
mass_min_after_thermo_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_aft_min_temp_7) +
                                       (1|fnest_id), 
                                     data = subset(late_nestling_parent_care,
                                                   !is.na(x = thermo_aft_min_temp_7)&
                                                     !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_after_thermo_lmer_7)
# Normal QQplot
{qqnorm(resid(mass_min_after_thermo_lmer_7))
  qqline(resid(mass_min_after_thermo_lmer_7))}
# Histogram of residuals
hist(resid(mass_min_after_thermo_lmer_7))
# Checking for influential outliers
infIndexPlot(mass_min_after_thermo_lmer_7, vars=c("Cook"))
infIndexPlot(mass_min_after_thermo_lmer_7, vars=c("Studentized"))

summary(mass_min_after_thermo_lmer_7)
confint(mass_min_after_thermo_lmer_7)

# Calculate R squared 
r.squaredGLMM(mass_min_after_thermo_lmer_7)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_after_thermo_lmer_7 <- bootMer(x = mass_min_after_thermo_lmer_7,
                                             FUN = fixef, nsim = 2000,
                                             seed = 632760,
                                             use.u = F, type = 'parametric')
tidy(boot_mass_min_after_thermo_lmer_7) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_after_thermo_lmer_7 <- boot.ci(boot_mass_min_after_thermo_lmer_7,
                                              type = c('perc', 'norm', 'basic'),
                                              index = 2) # CI for 1st betas
print(bt_ci_mass_min_after_thermo_lmer_7)


### After model - adjusted
mass_min_after_thermo_adj_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_aft_min_temp_7) + scale(nestling_number) +
                                           scale(days_summer) +
                                           (1|fnest_id), 
                                         data = subset(late_nestling_parent_care,
                                                       !is.na(x = thermo_aft_min_temp_7)&
                                                         !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_after_thermo_adj_lmer_7)
# Normal QQplot
{qqnorm(resid(mass_min_after_thermo_adj_lmer_7))
  qqline(resid(mass_min_after_thermo_adj_lmer_7))}
# Histogram of residuals
hist(resid(mass_min_after_thermo_adj_lmer_7))
# Checking for influential outliers
infIndexPlot(mass_min_after_thermo_adj_lmer_7, vars=c("Cook"))
infIndexPlot(mass_min_after_thermo_adj_lmer_7, vars=c("Studentized"))

summary(mass_min_after_thermo_adj_lmer_7)
confint(mass_min_after_thermo_adj_lmer_7)

# Calculate R squared 
r.squaredGLMM(mass_min_after_thermo_adj_lmer_7)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_after_thermo_adj_lmer_7 <- bootMer(x = mass_min_after_thermo_adj_lmer_7,
                                                 FUN = fixef, nsim = 2000,
                                                 seed = 632760,
                                                 use.u = F, type = 'parametric')
tidy(boot_mass_min_after_thermo_adj_lmer_7) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_after_thermo_adj_lmer_7 <- boot.ci(boot_mass_min_after_thermo_adj_lmer_7,
                                                  type = c('perc', 'norm', 'basic'),
                                                  index = 2) # CI for 1st betas
print(bt_ci_mass_min_after_thermo_adj_lmer_7)


############################# Maximum temperature #############################

### Before model - unadjusted
mass_max_before_thermo_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp_7) +
                                        (1|fnest_id), 
                                      data = subset(late_nestling_parent_care,
                                                    !is.na(x = thermo_bef_max_temp_7)&
                                                      !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_before_thermo_lmer_7)
# Normal QQplot
{qqnorm(resid(mass_max_before_thermo_lmer_7))
  qqline(resid(mass_max_before_thermo_lmer_7))}
# Histogram of residuals
hist(resid(mass_max_before_thermo_lmer_7))
# Checking for influential outliers
infIndexPlot(mass_max_before_thermo_lmer_7, vars=c("Cook"))
infIndexPlot(mass_max_before_thermo_lmer_7, vars=c("Studentized"))

summary(mass_max_before_thermo_lmer_7)
confint(mass_max_before_thermo_lmer_7)

# Calculate R squared 
r.squaredGLMM(mass_max_before_thermo_lmer_7)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_before_thermo_lmer_7 <- bootMer(x = mass_max_before_thermo_lmer_7,
                                              FUN = fixef, nsim = 2000,
                                              seed = 632760,
                                              use.u = F, type = 'parametric')
tidy(boot_mass_max_before_thermo_lmer_7) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_before_thermo_lmer_7 <- boot.ci(boot_mass_max_before_thermo_lmer_7,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 2) # CI for 1st betas
print(bt_ci_mass_max_before_thermo_lmer_7)

### Before model - adjusted
mass_max_before_thermo_adj_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp_7) + scale(nestling_number) +
                                            scale(days_summer) +
                                            (1|fnest_id), 
                                          data = subset(late_nestling_parent_care,
                                                        !is.na(x = thermo_bef_max_temp_7)&
                                                          !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_before_thermo_adj_lmer_7)
# Normal QQplot
{qqnorm(resid(mass_max_before_thermo_adj_lmer_7))
  qqline(resid(mass_max_before_thermo_adj_lmer_7))}
# Histogram of residuals
hist(resid(mass_max_before_thermo_adj_lmer_7))
# Checking for influential outliers
infIndexPlot(mass_max_before_thermo_adj_lmer_7, vars=c("Cook"))
infIndexPlot(mass_max_before_thermo_adj_lmer_7, vars=c("Studentized"))

summary(mass_max_before_thermo_adj_lmer_7)
confint(mass_max_before_thermo_adj_lmer_7)

# Calculate R squared 
r.squaredGLMM(mass_max_before_thermo_adj_lmer_7)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_before_thermo_adj_lmer_7 <- bootMer(x = mass_max_before_thermo_adj_lmer_7,
                                                  FUN = fixef, nsim = 2000,
                                                  seed = 632760,
                                                  use.u = F, type = 'parametric')
tidy(boot_mass_max_before_thermo_adj_lmer_7) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_before_thermo_adj_lmer_7 <- boot.ci(boot_mass_max_before_thermo_adj_lmer_7,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 2) # CI for 1st betas
print(bt_ci_mass_max_before_thermo_adj_lmer_7)


### After model - unadjusted
mass_max_after_thermo_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_aft_max_temp_7) +
                                       (1|fnest_id), 
                                     data = subset(late_nestling_parent_care,
                                                   !is.na(x = thermo_aft_max_temp_7)&
                                                     !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_after_thermo_lmer_7)
# Normal QQplot
{qqnorm(resid(mass_max_after_thermo_lmer_7))
  qqline(resid(mass_max_after_thermo_lmer_7))}
# Histogram of residuals
hist(resid(mass_max_after_thermo_lmer_7))
# Checking for influential outliers
infIndexPlot(mass_max_after_thermo_lmer_7, vars=c("Cook"))
infIndexPlot(mass_max_after_thermo_lmer_7, vars=c("Studentized"))

summary(mass_max_after_thermo_lmer_7)
confint(mass_max_after_thermo_lmer_7)

# Calculate R squared 
r.squaredGLMM(mass_max_after_thermo_lmer_7)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_after_thermo_lmer_7 <- bootMer(x = mass_max_after_thermo_lmer_7,
                                             FUN = fixef, nsim = 2000,
                                             seed = 632760,
                                             use.u = F, type = 'parametric')
tidy(boot_mass_max_after_thermo_lmer_7) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_after_thermo_lmer_7 <- boot.ci(boot_mass_max_after_thermo_lmer_7,
                                              type = c('perc', 'norm', 'basic'),
                                              index = 2) # CI for 1st betas
print(bt_ci_mass_max_after_thermo_lmer_7)


### After model - adjusted
mass_max_after_thermo_adj_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_aft_max_temp_7) + scale(nestling_number) +
                                           scale(days_summer) +
                                           (1|fnest_id), 
                                         data = subset(late_nestling_parent_care,
                                                       !is.na(x = thermo_aft_max_temp_7)&
                                                         !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_after_thermo_adj_lmer_7)
# Normal QQplot
{qqnorm(resid(mass_max_after_thermo_adj_lmer_7))
  qqline(resid(mass_max_after_thermo_adj_lmer_7))}
# Histogram of residuals
hist(resid(mass_max_after_thermo_adj_lmer_7))
# Checking for influential outliers
infIndexPlot(mass_max_after_thermo_adj_lmer_7, vars=c("Cook"))
infIndexPlot(mass_max_after_thermo_adj_lmer_7, vars=c("Studentized"))

summary(mass_max_after_thermo_adj_lmer_7)
confint(mass_max_after_thermo_adj_lmer_7)

# Calculate R squared 
r.squaredGLMM(mass_max_after_thermo_adj_lmer_7)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_after_thermo_adj_lmer_7 <- bootMer(x = mass_max_after_thermo_adj_lmer_7,
                                                 FUN = fixef, nsim = 2000,
                                                 seed = 632760,
                                                 use.u = F, type = 'parametric')
tidy(boot_mass_max_after_thermo_adj_lmer_7) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_after_thermo_adj_lmer_7 <- boot.ci(boot_mass_max_after_thermo_adj_lmer_7,
                                                  type = c('perc', 'norm', 'basic'),
                                                  index = 2) # CI for 1st betas
print(bt_ci_mass_max_after_thermo_adj_lmer_7)


####################### Temperature interquartile range #######################

### Before model - unadjusted
mass_iqr_before_thermo_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp_7) + 
                                        (1|fnest_id), 
                                      data = subset(late_nestling_parent_care,
                                                    !is.na(x = thermo_bef_iqr_temp_7)&
                                                      !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_before_thermo_lmer_7)
# Normal QQplot
{qqnorm(resid(mass_iqr_before_thermo_lmer_7))
  qqline(resid(mass_iqr_before_thermo_lmer_7))}
# Histogram of residuals
hist(resid(mass_iqr_before_thermo_lmer_7))
# Checking for influential outliers
infIndexPlot(mass_iqr_before_thermo_lmer_7, vars=c("Cook"))
infIndexPlot(mass_iqr_before_thermo_lmer_7, vars=c("Studentized"))

summary(mass_iqr_before_thermo_lmer_7)
confint(mass_iqr_before_thermo_lmer_7)

# Calculate R squared 
r.squaredGLMM(mass_iqr_before_thermo_lmer_7)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_before_thermo_lmer_7 <- bootMer(x = mass_iqr_before_thermo_lmer_7,
                                              FUN = fixef, nsim = 2000,
                                              seed = 632760,
                                              use.u = F, type = 'parametric')
tidy(boot_mass_iqr_before_thermo_lmer_7) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_before_thermo_lmer_7 <- boot.ci(boot_mass_iqr_before_thermo_lmer_7,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_before_thermo_lmer_7)


### Before model - adjusted
mass_iqr_before_thermo_adj_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp_7) + scale(nestling_number) +
                                            scale(days_summer) +  
                                            (1|fnest_id), 
                                          data = subset(late_nestling_parent_care,
                                                        !is.na(x = thermo_bef_iqr_temp_7)&
                                                          !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_before_thermo_adj_lmer_7)
# Normal QQplot
{qqnorm(resid(mass_iqr_before_thermo_adj_lmer_7))
  qqline(resid(mass_iqr_before_thermo_adj_lmer_7))}
# Histogram of residuals
hist(resid(mass_iqr_before_thermo_adj_lmer_7))
# Checking for influential outliers
infIndexPlot(mass_iqr_before_thermo_adj_lmer_7, vars=c("Cook"))
infIndexPlot(mass_iqr_before_thermo_adj_lmer_7, vars=c("Studentized"))

summary(mass_iqr_before_thermo_adj_lmer_7)
confint(mass_iqr_before_thermo_adj_lmer_7)

# Calculate R squared 
r.squaredGLMM(mass_iqr_before_thermo_adj_lmer_7)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_before_thermo_adj_lmer_7 <- bootMer(x = mass_iqr_before_thermo_adj_lmer_7,
                                                  FUN = fixef, nsim = 2000,
                                                  seed = 632760,
                                                  use.u = F, type = 'parametric')
tidy(boot_mass_iqr_before_thermo_adj_lmer_7) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_before_thermo_adj_lmer_7 <- boot.ci(boot_mass_iqr_before_thermo_adj_lmer_7,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_before_thermo_adj_lmer_7)


### After model - unadjusted
mass_iqr_after_thermo_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_aft_iqr_temp_7) + 
                                       (1|fnest_id), 
                                     data = subset(late_nestling_parent_care,
                                                   !is.na(x = thermo_aft_iqr_temp_7)&
                                                     !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_after_thermo_lmer_7)
# Normal QQplot
{qqnorm(resid(mass_iqr_after_thermo_lmer_7))
  qqline(resid(mass_iqr_after_thermo_lmer_7))}
# Histogram of residuals
hist(resid(mass_iqr_after_thermo_lmer_7))
# Checking for influential outliers
infIndexPlot(mass_iqr_after_thermo_lmer_7, vars=c("Cook"))
infIndexPlot(mass_iqr_after_thermo_lmer_7, vars=c("Studentized"))

summary(mass_iqr_after_thermo_lmer_7)
confint(mass_iqr_after_thermo_lmer_7)

# Calculate R squared 
r.squaredGLMM(mass_iqr_after_thermo_lmer_7)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_after_thermo_lmer_7 <- bootMer(x = mass_iqr_after_thermo_lmer_7,
                                             FUN = fixef, nsim = 2000,
                                             seed = 632760,
                                             use.u = F, type = 'parametric')
tidy(boot_mass_iqr_after_thermo_lmer_7) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_after_thermo_lmer_7 <- boot.ci(boot_mass_iqr_after_thermo_lmer_7,
                                              type = c('perc', 'norm', 'basic'),
                                              index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_after_thermo_lmer_7)


### After model - adjusted
mass_iqr_after_thermo_adj_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_aft_iqr_temp_7) + scale(nestling_number) +
                                           scale(days_summer) + 
                                           (1|fnest_id), 
                                         data = subset(late_nestling_parent_care,
                                                       !is.na(x = thermo_aft_iqr_temp_7)&
                                                         !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_after_thermo_adj_lmer_7)
# Normal QQplot
{qqnorm(resid(mass_iqr_after_thermo_adj_lmer_7))
  qqline(resid(mass_iqr_after_thermo_adj_lmer_7))}
# Histogram of residuals
hist(resid(mass_iqr_after_thermo_adj_lmer_7))
# Checking for influential outliers
infIndexPlot(mass_iqr_after_thermo_adj_lmer_7, vars=c("Cook"))
infIndexPlot(mass_iqr_after_thermo_adj_lmer_7, vars=c("Studentized"))

summary(mass_iqr_after_thermo_adj_lmer_7)
confint(mass_iqr_after_thermo_adj_lmer_7)

# Calculate R squared 
r.squaredGLMM(mass_iqr_after_thermo_adj_lmer_7)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_after_thermo_adj_lmer_7 <- bootMer(x = mass_iqr_after_thermo_adj_lmer_7,
                                                 FUN = fixef, nsim = 2000,
                                                 seed = 632760,
                                                 use.u = F, type = 'parametric')
tidy(boot_mass_iqr_after_thermo_adj_lmer_7) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_after_thermo_adj_lmer_7 <- boot.ci(boot_mass_iqr_after_thermo_adj_lmer_7,
                                                  type = c('perc', 'norm', 'basic'),
                                                  index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_after_thermo_adj_lmer_7)



###############################################################################
##############                 Model visualizations              ##############
###############################################################################

################### Parental feeding models (question 3) #####################

### Minimum temperature

## Compute predicted values from stratified models and bind them together
# Low
pred_low <- ggpredict(mass_min_temp_blups_low_adj_lmer, c("nest_min_temp"))
pred_low$feeding <- "1"

# Med
pred_med <- ggpredict(mass_min_temp_blups_med_adj_lmer, c("nest_min_temp"))
pred_med$feeding <- "2"

# High
pred_high <- ggpredict(mass_min_temp_blups_high_adj_lmer, c("nest_min_temp"))
pred_high$feeding <- "3"

pred_feeding <- rbind(pred_low, pred_med, pred_high)

## create a subset of the raw data excluding missing data
sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_min_temp) &
                !is.na(x = feeding_expontd_blups))

## Create vectors of colors for plot
cols <- c("#481567FF", "#20A387FF", "#95D840FF")

## Create plot of model predictions and raw data
temp_min_feeding_blups_predicted <- ggplot() +
  geom_line(data = pred_feeding, aes(x = x, y = predicted, col = feeding), size = 1.5) +
  geom_point(data = sub, aes(x = nest_min_temp, y = mass_pre_obs,
                             col = as.character(feeding_expontd_blups_strat)), 
             show.legend = FALSE, size = 1.5, alpha = 0.5) +
  theme_classic() +
  labs(x = "Minimum temperature (C)", y = "Nestling mass (g)",
       colour = "Parent feeding level", fill = "Parent feeding level") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11),
        legend.text = element_text(size = 9), legend.title = element_text(size = 11)) +
  scale_color_manual(values = cols, labels=c("Low", "Med", "High"),
                     aesthetics = c("colour", "fill"))

# Print plot
print(temp_min_feeding_blups_predicted)



### Maximum temperature

## Compute predicted values from stratified models and bind them together
# Low
pred_low <- ggpredict(mass_max_temp_blups_low_adj_lmer, c("nest_max_temp"))
pred_low$feeding <- "1"

# Med
pred_med <- ggpredict(mass_max_temp_blups_med_adj_lmer, c("nest_max_temp"))
pred_med$feeding <- "2"

# High
pred_high <- ggpredict(mass_max_temp_blups_high_adj_lmer, c("nest_max_temp"))
pred_high$feeding <- "3"

pred_feeding <- rbind(pred_low, pred_med, pred_high)

## Create subset of raw data excluding missing data
sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_max_temp) &
                !is.na(x = feeding_expontd_blups))

## Plot models predictions and raw data
temp_max_feeding_blups_predicted <- ggplot() +
  geom_line(data = pred_feeding, aes(x = x, y = predicted, col = feeding), size = 1.5) +
  geom_point(data = sub, aes(x = nest_max_temp, y = mass_pre_obs,
                             col = as.character(feeding_expontd_blups_strat)), 
             show.legend = FALSE, size = 1.5, alpha = 0.5) +
  theme_classic() +
  labs(x = "Maximum temperature (C)", y = "Nestling mass (g)",
       colour = "Parent feeding level", fill = "Parent feeding level") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11),
        legend.text = element_text(size = 9), legend.title = element_text(size = 11)) +
  scale_color_manual(values = cols, labels=c("Low", "Med", "High"),
                     aesthetics = c("colour", "fill")) 


# Print plot
print(temp_max_feeding_blups_predicted)


### Temperature interquartile range

## Compute predicted values from stratified models and bind them together
# Low
pred_low <- ggpredict(mass_iqr_temp_blups_low_adj_lmer, c("nest_iqr_temp"))
pred_low$feeding <- "1"

# Med
pred_med <- ggpredict(mass_iqr_temp_blups_med_adj_lmer, c("nest_iqr_temp"))
pred_med$feeding <- "2"

# High
pred_high <- ggpredict(mass_iqr_temp_blups_high_adj_lmer, c("nest_iqr_temp"))
pred_high$feeding <- "3"

pred_feeding <- rbind(pred_low, pred_med, pred_high)

# Create subset of raw data excluding missing data
sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_iqr_temp) &
                !is.na(x = feeding_expontd_blups))


# Plot model predictions and raw data
temp_iqr_feeding_blups_predicted <- ggplot() +
  geom_line(data = pred_feeding, aes(x = x, y = predicted, col = feeding), size = 1.5) +
  geom_point(data = sub, aes(x = nest_iqr_temp, y = mass_pre_obs,
                             col = as.character(feeding_expontd_blups_strat)), 
             show.legend = FALSE, size = 1.5, alpha = 0.5) +
  theme_classic() +
  labs(x = "Temperature variability (C)", y = "Nestling mass (g)",
       colour = "Parent feeding level", fill = "Parent feeding level") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11),
        legend.text = element_text(size = 9), legend.title = element_text(size = 11)) +
  scale_color_manual(values = cols, labels=c("Low", "Med", "High"),
                     aesthetics = c("colour", "fill"))

# Print plot
print(temp_iqr_feeding_blups_predicted)


## Combine plots for the three temperature variables
combined_temp_feeding_predictions <- 
  ggarrange(temp_min_feeding_blups_predicted, temp_max_feeding_blups_predicted,
            temp_iqr_feeding_blups_predicted,
            ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom", 
            labels = c("(a)", "(b)", "(c)"), 
            font.label = list(size = 11, face = "bold", color = "red"), hjust = -0.1)

# Print clombined plot
print(combined_temp_feeding_predictions)

# Save combined plot
ggsave('combined_temp_feeding_predictions_solid.png', plot = combined_temp_feeding_predictions, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 4, 
       height = 10, 
       units = c('in'), dpi = 300, limitsize = TRUE)



####################### Relative size models (question 2) #####################

### Minimum temperature

## Compute predicted values from stratified models and bind them together
# Small
pred_small <- ggpredict(mass_min_temp_mid_size_small_adj_lmer, c("nest_min_temp"))
pred_small$size <- "min"

# Other
pred_big <- ggpredict(mass_min_temp_mid_size_big_adj_lmer, c("nest_min_temp"))
pred_big$size <- "other"

pred_size <- rbind(pred_small, pred_big)

# Create subset of raw data excluding missing data
sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_min_temp) &
                !is.na(x = mid_size_order))

# Create vector of colors for plot
cols <- c("#481567FF","#95D840FF")

# Plot model predictions and raw data
temp_min_size_predicted <- ggplot() +
  geom_line(data = pred_size, aes(x = x, y = predicted, col = size), size = 1.5) +
  geom_point(data = sub, aes(x = nest_min_temp, y = mass_pre_obs,
                             col = as.character(mid_size_order)), 
             show.legend = FALSE, size = 1.5, alpha = 0.5) +
  theme_classic() +
  labs(x = "Minimum temperature (C)", y = "Nestling mass (g)",
       colour = "Relative nestling size", fill = "Relative nestling size") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11),
        legend.text = element_text(size = 9), legend.title = element_text(size = 11)) +
  scale_color_manual(values = cols, labels=c("Smallest", "Other"),
                     aesthetics = c("colour", "fill")) 

# Print plot
print(temp_min_size_predicted)


### Maximum temperature

## Compute predicted values from stratified models and bind them together
# Small
pred_small <- ggpredict(mass_max_temp_mid_size_small_adj_lmer, c("nest_max_temp"))
pred_small$size <- "min"

# Other
pred_big <- ggpredict(mass_max_temp_mid_size_big_adj_lmer, c("nest_max_temp"))
pred_big$size <- "other"

pred_size <- rbind(pred_small, pred_big)

# Create subset of raw data excluding missing data
sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_max_temp) &
                !is.na(x = mid_size_order))

# Plot model predictions and raw data
temp_max_size_predicted <- ggplot() +
  geom_line(data = pred_size, aes(x = x, y = predicted, col = size), size = 1.5) +
  geom_point(data = sub, aes(x = nest_max_temp, y = mass_pre_obs,
                             col = as.character(mid_size_order)), 
             show.legend = FALSE, size = 1.5, alpha = 0.5) +
  theme_classic() +
  labs(x = "Maximum temperature (C)", y = "Nestling mass (g)",
       colour = "Relative nestling size", fill = "Relative nestling size") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11),
        legend.text = element_text(size = 9), legend.title = element_text(size = 11)) +
  scale_color_manual(values = cols, labels=c("Smallest", "Other"),
                     aesthetics = c("colour", "fill"))

# Print plot
print(temp_max_size_predicted)


### Temperature interquartile range

## Compute predicted values from stratified models and bind them together
# Small
pred_small <- ggpredict(mass_iqr_temp_mid_size_small_adj_lmer, c("nest_iqr_temp"))
pred_small$size <- "min"

# Other
pred_big <- ggpredict(mass_iqr_temp_mid_size_big_adj_lmer, c("nest_iqr_temp"))
pred_big$size <- "other"

pred_size <- rbind(pred_small, pred_big)

# Create subset of raw data excluding missing data
sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_iqr_temp) &
                !is.na(x = mid_size_order))

# Plot model predictions and raw data
temp_iqr_size_predicted <- ggplot() +
  geom_line(data = pred_size, aes(x = x, y = predicted, col = size), size = 1.5) +
  geom_point(data = sub, aes(x = nest_iqr_temp, y = mass_pre_obs,
                             col = as.character(mid_size_order)), 
             show.legend = FALSE, size = 1.5, alpha = 0.5) +
  theme_classic() +
  labs(x = "Temperature variability (C)", y = "Nestling mass (g)",
       colour = "Relative nestling size", fill = "Relative nestling size") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11),
        legend.text = element_text(size = 9), legend.title = element_text(size = 11)) +
  scale_color_manual(values = cols, labels=c("Smallest", "Other"),
                     aesthetics = c("colour", "fill")) 

# Print plot
print(temp_iqr_size_predicted)

## Combine plots for three temeprature variables
combined_temp_size_predictions <- 
  ggarrange(temp_min_size_predicted, temp_max_size_predicted,
            temp_iqr_size_predicted,
            ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom", 
            labels = c("(a)", "(b)", "(c)"), 
            font.label = list(size = 11, face = "bold", color = "red"), hjust = -0.1)

# Print combined plot
print(combined_temp_size_predictions)

# Save combined plot
ggsave('combined_temp_size_predictions_solid.png', plot = combined_temp_size_predictions, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 4, 
       height = 10, 
       units = c('in'), dpi = 300, limitsize = TRUE)


################## Sensitive periods models (question 1) ######################

### Minimum temperature 

## Compute predicted values from stratified models and bind them together
# Before
pred_bef <- ggpredict(mass_min_before_thermo_adj_lmer, c("thermo_bef_min_temp"))

# After
pred_aft <- ggpredict(mass_min_after_thermo_adj_lmer, c("thermo_aft_min_temp"))

# Create subset of raw data exlcuding missing data
sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = thermo_bef_min_temp) &
                !is.na(x = thermo_aft_min_temp))

# Create vector of colors with labels for plot
colors <- c("Before day 6" = "#481567FF", "After day 6" = "#95D840FF")

# Plot model predictions and raw data
temp_thermo_both_min_predicted <- ggplot() + 
  theme_classic() +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11),
        legend.text = element_text(size = 9), legend.title = element_text(size = 11)) +
  geom_line(data= pred_bef, size = 1.5, aes(x = x, y = predicted, col = "Before day 6")) +
  geom_line(data= pred_aft, size = 1.5, aes(x = x, y = predicted, col = "After day 6")) +
  geom_point(data = sub, aes(x = thermo_bef_min_temp, y = mass_pre_obs, col = "Before day 6"), 
             size = 1, alpha = 0.5, show.legend = FALSE) +
  geom_point(data = sub, aes(x = thermo_aft_min_temp, y = mass_pre_obs, col = "After day 6"), 
             size = 1, alpha = 0.5, show.legend = FALSE) +
  labs(x = "Minimum temperature (C)", y = "Nestling mass (g)", 
       color = "Before vs. after day 6", fill = "Before vs. after day 6") +
  scale_color_manual(values = colors, breaks = c("Before day 6", "After day 6"),
                     aesthetics = c("color", "fill"))

# Print plot
print(temp_thermo_both_min_predicted)


### Maximum temperature

## Compute predicted values from stratified models and bind them together
# Before
pred_bef <- ggpredict(mass_max_before_thermo_adj_lmer, c("thermo_bef_max_temp"))

# After
pred_aft <- ggpredict(mass_max_after_thermo_adj_lmer, c("thermo_aft_max_temp"))

# Create subset of raw data excluding missing data
sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = thermo_bef_max_temp) &
                !is.na(x = thermo_aft_max_temp))

# Plot models predictions and raw data
temp_thermo_both_max_predicted <- ggplot() + 
  theme_classic() +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11),
        legend.text = element_text(size = 9), legend.title = element_text(size = 11)) +
  geom_line(data= pred_bef, size = 1.5, aes(x = x, y = predicted, col = "Before day 6")) +
  geom_line(data= pred_aft, size = 1.5, aes(x = x, y = predicted, col = "After day 6")) +
  geom_point(data = sub, aes(x = thermo_bef_max_temp, y = mass_pre_obs, col = "Before day 6"), 
             size = 1, alpha = 0.5, show.legend = FALSE) +
  geom_point(data = sub, aes(x = thermo_aft_max_temp, y = mass_pre_obs, col = "After day 6"), 
             size = 1, alpha = 0.5, show.legend = FALSE) +
  labs(x = "Maximum temperature (C)", y = "Nestling mass (g)", 
       color = "Before vs. after day 6", fill = "Before vs. after day 6") +
  scale_color_manual(values = colors, breaks = c("Before day 6", "After day 6"),
                     aesthetics = c("color", "fill"))

# Print plot
print(temp_thermo_both_max_predicted)


### Temperature interquartile range
# Before
pred_bef <- ggpredict(mass_iqr_before_thermo_adj_lmer, c("thermo_bef_iqr_temp"))

# After
pred_aft <- ggpredict(mass_iqr_after_thermo_adj_lmer, c("thermo_aft_iqr_temp"))

# Create subset of raw data excluding missing data
sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = thermo_bef_iqr_temp) &
                !is.na(x = thermo_aft_iqr_temp))

# Plot model predictions and raw data
temp_thermo_both_iqr_predicted <- ggplot() + 
  theme_classic() +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11),
        legend.text = element_text(size = 9), legend.title = element_text(size = 11)) +
  geom_line(data= pred_bef, size = 1.5, aes(x = x, y = predicted, col = "Before day 6")) +
  geom_line(data= pred_aft, size = 1.5, aes(x = x, y = predicted, col = "After day 6")) +
  geom_point(data = sub, aes(x = thermo_bef_iqr_temp, y = mass_pre_obs, col = "Before day 6"), 
             size = 1, alpha = 0.5, show.legend = FALSE) +
  geom_point(data = sub, aes(x = thermo_aft_iqr_temp, y = mass_pre_obs, col = "After day 6"), 
             size = 1, alpha = 0.5, show.legend = FALSE) +
  labs(x = "Temperature variability (C)", y = "Nestling mass (g)", 
       color = "Before vs. after day 6", fill = "Before vs. after day 6") +
  scale_color_manual(values = colors, breaks = c("Before day 6", "After day 6"),
                     aesthetics = c("color", "fill"))

# Print plot
print(temp_thermo_both_iqr_predicted)


### Combine plots for three temperature variables
combined_thermo_both_predictions_long <- 
  ggarrange(temp_thermo_both_min_predicted, 
            temp_thermo_both_max_predicted, temp_thermo_both_iqr_predicted,
            ncol = 1, nrow = 3, common.legend = TRUE,
            legend = "bottom", labels = c("(a)", "(b)", "(c)"), 
            font.label = list(size = 11, face = "bold", color = "red"), hjust = -0.1)

# Print combined plot
print(combined_thermo_both_predictions_long)

# Save combined plot
ggsave('combined_thermo_both_predictions_long.png', plot = combined_thermo_both_predictions_long, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 4, 
       height = 10, 
       units = c('in'), dpi = 300, limitsize = TRUE)

