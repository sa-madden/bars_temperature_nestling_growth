####### Purpose: run covariate models for barn swallow
####### nest microclimate and nestling growth dataset
####### By: Sage Madden
####### Created: 1/16/2023
####### Last modified: 10/21/2024

# Code Blocks
# 1: Configure work space
# 2: Load data
# 3: Nestling growth vs. covariates
# 4: Temperature vs. covariates
# 5: Parental care vs. covariates
# 6: Growth vs. parental care


###############################################################################
##############                 Configure work space              ##############
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
library('boot')

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
##############                     Load Data                     ##############
###############################################################################  

### 2.1 Load RData
## Load RData tidy barn swallow data
load('Data/Tidy/tidy_parent_nestl_weather_data_10-4_with_BLUPs.RData')

# Make site a factor
late_nestling_parent_care$fsite <- as.factor(late_nestling_parent_care$site)
late_nestling_parent_care$fnest_id <- as.factor(late_nestling_parent_care$nest_id)

late_nestling_parent_care$num_pairs <- as.numeric(late_nestling_parent_care$num_pairs)


###############################################################################
##############           Nestling growth vs. covariates          ##############
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

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_hatch_lmer <- bootMer(x = mass_hatch_lmer,
                                      FUN = fixef, nsim = 2000,
                                      seed = 632760,
                                      use.u = F, type = 'parametric')
tidy(boot_mass_hatch_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI
bt_ci_mass_hatch_lmer <- boot.ci(boot_mass_hatch_lmer,
                                       type = c('perc', 'norm', 'basic'),
                                       index = 2) # CI for 1st betas
print(bt_ci_mass_hatch_lmer)


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

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_colony_lmer <- bootMer(x = mass_colony_lmer,
                                FUN = fixef, nsim = 2000,
                                seed = 632760,
                                use.u = F, type = 'parametric')
tidy(boot_mass_colony_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI
bt_ci_mass_colony_lmer <- boot.ci(boot_mass_colony_lmer,
                                 type = c('perc', 'norm', 'basic'),
                                 index = 2) # CI for 1st betas
print(bt_ci_mass_colony_lmer)


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

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_brood_lmer <- bootMer(x = mass_brood_lmer,
                                 FUN = fixef, nsim = 2000,
                                 seed = 632760,
                                 use.u = F, type = 'parametric')
tidy(boot_mass_brood_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI
bt_ci_mass_brood_lmer <- boot.ci(boot_mass_brood_lmer,
                                  type = c('perc', 'norm', 'basic'),
                                  index = 2) # CI for 1st betas
print(bt_ci_mass_brood_lmer)

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

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_age_lmer <- bootMer(x = mass_age_lmer,
                                FUN = fixef, nsim = 2000,
                                seed = 632760,
                                use.u = F, type = 'parametric')
tidy(boot_mass_age_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI
bt_ci_mass_age_lmer <- boot.ci(boot_mass_age_lmer,
                                 type = c('perc', 'norm', 'basic'),
                                 index = 2) # CI for 1st betas
print(bt_ci_mass_age_lmer)


###############################################################################
##############        Temperature vs. covariates                 ##############
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

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_max_temp_hatch_lmer <- bootMer(x = max_temp_hatch_lmer,
                                FUN = fixef, nsim = 2000,
                                seed = 632760,
                                use.u = F, type = 'parametric')
tidy(boot_max_temp_hatch_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI
bt_ci_max_temp_hatch_lmer <- boot.ci(boot_max_temp_hatch_lmer,
                                 type = c('perc', 'norm', 'basic'),
                                 index = 2) # CI for 1st betas
print(bt_ci_max_temp_hatch_lmer)


## min_temp and hatch date -- significant positive relationship 
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

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_min_temp_hatch_lmer <- bootMer(x = min_temp_hatch_lmer,
                                    FUN = fixef, nsim = 2000,
                                    seed = 632760,
                                    use.u = F, type = 'parametric')
tidy(boot_min_temp_hatch_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI
bt_ci_min_temp_hatch_lmer <- boot.ci(boot_min_temp_hatch_lmer,
                                     type = c('perc', 'norm', 'basic'),
                                     index = 2) # CI for 1st betas
print(bt_ci_min_temp_hatch_lmer)


## iqr_temp and hatch date -- significant negative relationship
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

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_iqr_temp_hatch_lmer <- bootMer(x = iqr_temp_hatch_lmer,
                                    FUN = fixef, nsim = 2000,
                                    seed = 632760,
                                    use.u = F, type = 'parametric')
tidy(boot_iqr_temp_hatch_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI
bt_ci_iqr_temp_hatch_lmer <- boot.ci(boot_iqr_temp_hatch_lmer,
                                     type = c('perc', 'norm', 'basic'),
                                     index = 2) # CI for 1st betas
print(bt_ci_iqr_temp_hatch_lmer)


###############################################################################
##############         Parental care  vs. covariates             ##############
###############################################################################
## feeding_blups and hatch date -- not significant
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

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_feeding_blups_hatch_lmer <- bootMer(x = feeding_blups_hatch_lmer,
                                FUN = fixef, nsim = 2000,
                                seed = 632760,
                                use.u = F, type = 'parametric')
tidy(boot_feeding_blups_hatch_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI
bt_ci_feeding_blups_hatch_lmer <- boot.ci(boot_feeding_blups_hatch_lmer,
                                 type = c('perc', 'norm', 'basic'),
                                 index = 2) # CI for 1st betas
print(bt_ci_feeding_blups_hatch_lmer)


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


## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_feeding_blups_hatch_lmer <- bootMer(x = feeding_blups_hatch_lmer,
                                         FUN = fixef, nsim = 2000,
                                         seed = 632760,
                                         use.u = F, type = 'parametric')
tidy(boot_feeding_blups_hatch_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI
bt_ci_feeding_blups_hatch_lmer <- boot.ci(boot_feeding_blups_hatch_lmer,
                                          type = c('perc', 'norm', 'basic'),
                                          index = 2) # CI for 1st betas
print(bt_ci_feeding_blups_hatch_lmer)


###############################################################################
##############             Growth vs. parental care              ##############
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

# Calculate avg feeding rates as an alternative to BLUPs
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
