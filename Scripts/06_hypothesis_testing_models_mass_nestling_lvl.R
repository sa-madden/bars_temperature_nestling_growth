####### Purpose: run models to test hypotheses for barn swallow
####### nest mircoclimate and nestling growth dataset
####### By: Sage Madden
####### Created: 12/19/2022
####### Last modified: 12/9/2024

# Code Blocks
# 1: Configure work space
# 2: Load data
# 3: Growth, temp, and parental care
# 4: Growth, temp, and relative size
# 5: Sensitive development periods
# 6: Model visualizations


###############################################################################
##############                 Configure work space              ##############
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
library('ggpubr')
library('ggeffects')
library('ggnewscale')


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

late_nestling_parent_care$num_pairs <- as.numeric(late_nestling_parent_care$num_pairs)


###############################################################################
##############             Growth, temp, and parental care       ##############
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
                                         (1|fsite) +
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
                                        (1|fsite) +
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



### Outliers removed
mass_outliers_removed <- late_nestling_parent_care[-c(4, 40), ]

### Mass unadjusted with outliers removed
mass_min_temp_blups_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) *
                                   scale(feeding_expontd_blups) + (1|fsite) +
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

## Mass adjusted with outliers removed
mass_min_temp_blups_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) *
                                       scale(feeding_expontd_blups) + scale(nestling_number) + 
                                       scale(days_summer) + (1|fsite) +
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


## LOW parental care model - adjusted
mass_max_temp_blups_low_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                           (1|fsite) +
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
                                           (1|fsite) +
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
                                            (1|fsite) +
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

### Outliers removed
mass_max_outliers_removed <- late_nestling_parent_care[-c(4, 39, 40), ]

## Mass unadjusted, outliers removed
mass_max_temp_blups_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) *
                                   scale(feeding_expontd_blups) + (1|fsite) +
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


### Mass adjusted -- outliers removed
mass_max_temp_blups_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) *
                                       scale(feeding_expontd_blups) + 
                                       scale(nestling_number) + scale(days_summer) +
                                       (1|fsite) +
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


## LOW parental care model - adjusted
mass_max_temp_blups_low_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                           (1|fsite) +
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


## MED parental care - adjusted
mass_max_temp_blups_med_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                           (1|fsite) +
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


## HIGH parental care - adjusted 
mass_max_temp_blups_high_noout_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                            scale(nestling_number) + scale(days_summer) + 
                                            (1|fsite) +
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


## LOW parental care model - adjusted
mass_iqr_temp_blups_low_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                           (1|fsite) +
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
                                           (1|fsite) +
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
                                            (1|fsite) +
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



### Outliers removed

### Mass unadjusted -- outliers removed 
mass_iqr_temp_blups_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) *
                                   scale(feeding_expontd_blups) + (1|fsite) +
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

### Mass adjusted -- outliers removed 
mass_iqr_temp_blups_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) *
                                       scale(feeding_expontd_blups) + 
                                       scale(nestling_number) + scale(days_summer) +
                                       (1|fsite) +
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


## LOW parental care model - adjusted, no outliers
mass_iqr_temp_blups_low_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(nest_iqr_temp) + 
                                           scale(nestling_number) + scale(days_summer) + 
                                           (1|fsite) +
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
                                           (1|fsite) +
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
                                            (1|fsite) +
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


###############################################################################
##############          Growth, temp, and relative size          ##############
###############################################################################

########## Relative size calculated at mid development #########################

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

## Stratified by relative size

# SMALL nestlings
mass_min_temp_mid_size_small_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) +
                                                scale(nestling_number) + scale(days_summer) +
                                                (1|fsite), 
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


# BIG nestlings
mass_min_temp_mid_size_big_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) +
                                              scale(nestling_number) + scale(days_summer) +
                                              (1|fsite) +
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

summary(mass_min_temp_mid_size_small_adj_lmer)
confint(mass_min_temp_mid_size_small_adj_lmer)

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


## Stratified by relative size

# SMALL nestlings
mass_max_temp_mid_size_small_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                            scale(nestling_number) + scale(days_summer) +
                                            (1|fsite), 
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


### Mass adjusted
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

# Calculate R squared 
r.squaredGLMM(mass_min_temp_size_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_size_lmer <- bootMer(x = mass_min_temp_size_lmer,
                                           FUN = fixef, nsim = 2000,
                                           seed = 632760,
                                           use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_size_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_size_lmer <- boot.ci(boot_mass_min_temp_size_lmer,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_size_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_min_temp_size_lmer_2 <- boot.ci(boot_mass_min_temp_size_lmer,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 3) # CI for 2nd betas
print(bt_ci_mass_min_temp_size_lmer_2)

# use 'boot' package to generate 95% CI for 3rd beta
bt_ci_mass_min_temp_size_lmer_3 <- boot.ci(boot_mass_min_temp_size_lmer,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 4) # CI for 3rd betas
print(bt_ci_mass_min_temp_size_lmer_3)


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

# Calculate R squared 
r.squaredGLMM(mass_min_temp_size_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_size_adj_lmer <- bootMer(x = mass_min_temp_size_adj_lmer,
                                                FUN = fixef, nsim = 2000,
                                                seed = 632760,
                                                use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_size_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_size_adj_lmer <- boot.ci(boot_mass_min_temp_size_adj_lmer,
                                                 type = c('perc', 'norm', 'basic'),
                                                 index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_size_adj_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_min_temp_size_adj_lmer_2 <- boot.ci(boot_mass_min_temp_size_adj_lmer,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 3) # CI for 2nd betas
print(bt_ci_mass_min_temp_size_adj_lmer_2)

# use 'boot' package to generate 95% CI for 5th beta
bt_ci_mass_min_temp_size_adj_lmer_5 <- boot.ci(boot_mass_min_temp_size_adj_lmer,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 6) # CI for 3rd betas
print(bt_ci_mass_min_temp_size_adj_lmer_5)


# SMALL nestlings
mass_min_temp_size_small_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) + 
                                                scale(nestling_number) + 
                                                scale(days_summer) +
                                                (1|fsite), 
                                              data = subset(late_nestling_parent_care,
                                                            size_order == "min",
                                                            !is.na(x = mass_pre_obs) & 
                                                              !is.na(x = nest_min_temp) &
                                                              !is.na(x = size_order)))

## Check diagnostics for the full model
plot(mass_min_temp_size_small_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_size_small_adj_lmer))
  qqline(resid(mass_min_temp_size_small_adj_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_size_small_adj_lmer))

summary(mass_min_temp_size_small_adj_lmer)
confint(mass_min_temp_size_small_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_size_small_adj_lmer <- bootMer(x = mass_min_temp_size_small_adj_lmer,
                                                      FUN = fixef, nsim = 2000,
                                                      seed = 632760,
                                                      use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_size_small_adj_lmer_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_size_small_adj_lmer <- boot.ci(boot_mass_min_temp_size_small_adj_lmer,
                                                       type = c('perc', 'norm', 'basic'),
                                                       index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_size_small_adj_lmer)


# BIG nestlings
mass_min_temp_size_big_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_min_temp) + 
                                              scale(nestling_number) + 
                                              scale(days_summer) +
                                              (1|fsite) +
                                              (1|fnest_id), 
                                            data = subset(late_nestling_parent_care,
                                                          size_order == "other",
                                                          !is.na(x = mass_pre_obs) & 
                                                            !is.na(x = nest_min_temp) &
                                                            !is.na(x = size_order)))


## Check diagnostics for the full model
plot(mass_min_temp_size_big_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_temp_size_big_adj_lmer))
  qqline(resid(mass_min_temp_size_big_adj_lmer))}
# Histogram of residuals
hist(resid(mass_min_temp_size_big_adj_lmer))

summary(mass_min_temp_size_big_adj_lmer)
confint(mass_min_temp_size_big_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_temp_size_big_adj_lmer <- bootMer(x = mass_min_temp_size_big_adj_lmer,
                                                    FUN = fixef, nsim = 2000,
                                                    seed = 632760,
                                                    use.u = F, type = 'parametric')
tidy(boot_mass_min_temp_size_big_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_temp_size_big_adj_lmer <- boot.ci(boot_mass_min_temp_size_big_adj_lmer,
                                                     type = c('perc', 'norm', 'basic'),
                                                     index = 2) # CI for 1st betas
print(bt_ci_mass_min_temp_size_big_adj_lmer)


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

# Calculate R squared 
r.squaredGLMM(mass_max_temp_size_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_size_lmer <- bootMer(x = mass_max_temp_size_lmer,
                                        FUN = fixef, nsim = 2000,
                                        seed = 632760,
                                        use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_size_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_size_lmer <- boot.ci(boot_mass_max_temp_size_lmer,
                                         type = c('perc', 'norm', 'basic'),
                                         index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_size_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_max_temp_size_lmer_2 <- boot.ci(boot_mass_max_temp_size_lmer,
                                           type = c('perc', 'norm', 'basic'),
                                           index = 3) # CI for 2nd betas
print(bt_ci_mass_max_temp_size_lmer_2)

# use 'boot' package to generate 95% CI for 3rd beta
bt_ci_mass_max_temp_size_lmer_3 <- boot.ci(boot_mass_max_temp_size_lmer,
                                           type = c('perc', 'norm', 'basic'),
                                           index = 4) # CI for 3rd betas
print(bt_ci_mass_max_temp_size_lmer_3)


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

# Calculate R squared 
r.squaredGLMM(mass_max_temp_size_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_size_adj_lmer <- bootMer(x = mass_max_temp_size_adj_lmer,
                                            FUN = fixef, nsim = 2000,
                                            seed = 632760,
                                            use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_size_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_size_adj_lmer <- boot.ci(boot_mass_max_temp_size_adj_lmer,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_size_adj_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_max_temp_size_adj_lmer_2 <- boot.ci(boot_mass_max_temp_size_adj_lmer,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 3) # CI for 2nd betas
print(bt_ci_mass_max_temp_size_adj_lmer_2)

# use 'boot' package to generate 95% CI for 5th beta
bt_ci_mass_max_temp_size_adj_lmer_5 <- boot.ci(boot_mass_max_temp_size_adj_lmer,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 6) # CI for 3rd betas
print(bt_ci_mass_max_temp_size_adj_lmer_5)


# SMALL nestlings
mass_max_temp_size_small_adj_lmer <- lmer(mass_pre_obs ~ scale(nest_max_temp) +
                                        scale(nestling_number) + scale(days_summer) +
                                          (1|fsite), 
                                      data = subset(late_nestling_parent_care,
                                                    size_order == "min",
                                                    !is.na(x = mass_pre_obs) & 
                                                      !is.na(x = nest_max_temp) &
                                                      !is.na(x = size_order)))


## Check diagnostics for the full model
plot(mass_max_temp_size_small_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_size_small_adj_lmer))
  qqline(resid(mass_max_temp_size_small_adj_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_size_small_adj_lmer))

summary(mass_max_temp_size_small_adj_lmer)
confint(mass_max_temp_size_small_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_size_small_adj_lmer <- bootMer(x = mass_max_temp_size_small_adj_lmer,
                                                  FUN = fixef, nsim = 2000,
                                                  seed = 632760,
                                                  use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_size_small_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_size_small_adj_lmer <- boot.ci(boot_mass_max_temp_size_small_adj_lmer,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_size_small_adj_lmer)


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


## Check diagnostics for the full model
plot(mass_max_temp_size_big_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_temp_size_big_adj_lmer))
  qqline(resid(mass_max_temp_size_big_adj_lmer))}
# Histogram of residuals
hist(resid(mass_max_temp_size_big_adj_lmer))

summary(mass_max_temp_size_big_adj_lmer)
confint(mass_max_temp_size_big_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_temp_size_big_adj_lmer <- bootMer(x = mass_max_temp_size_big_adj_lmer,
                                                  FUN = fixef, nsim = 2000,
                                                  seed = 632760,
                                                  use.u = F, type = 'parametric')
tidy(boot_mass_max_temp_size_big_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_temp_size_big_adj_lmer <- boot.ci(boot_mass_max_temp_size_big_adj_lmer,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 2) # CI for 1st betas
print(bt_ci_mass_max_temp_size_big_adj_lmer)


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

# Calculate R squared 
r.squaredGLMM(mass_iqr_temp_size_lmer)


## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_size_lmer <- bootMer(x = mass_iqr_temp_size_lmer,
                                        FUN = fixef, nsim = 2000,
                                        seed = 632760,
                                        use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_size_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_size_lmer <- boot.ci(boot_mass_iqr_temp_size_lmer,
                                         type = c('perc', 'norm', 'basic'),
                                         index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_size_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_iqr_temp_size_lmer_2 <- boot.ci(boot_mass_iqr_temp_size_lmer,
                                           type = c('perc', 'norm', 'basic'),
                                           index = 3) # CI for 2nd betas
print(bt_ci_mass_iqr_temp_size_lmer_2)

# use 'boot' package to generate 95% CI for 3rd beta
bt_ci_mass_iqr_temp_size_lmer_3 <- boot.ci(boot_mass_iqr_temp_size_lmer,
                                           type = c('perc', 'norm', 'basic'),
                                           index = 4) # CI for 3rd betas
print(bt_ci_mass_iqr_temp_size_lmer_3)


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

# Calculate R squared 
r.squaredGLMM(mass_iqr_temp_size_adj_lmer)


## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_size_adj_lmer <- bootMer(x = mass_iqr_temp_size_adj_lmer,
                                            FUN = fixef, nsim = 2000,
                                            seed = 632760,
                                            use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_size_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_size_adj_lmer <- boot.ci(boot_mass_iqr_temp_size_adj_lmer,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_size_adj_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_iqr_temp_size_adj_lmer_2 <- boot.ci(boot_mass_iqr_temp_size_adj_lmer,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 3) # CI for 2nd betas
print(bt_ci_mass_iqr_temp_size_adj_lmer_2)

# use 'boot' package to generate 95% CI for 5th beta
bt_ci_mass_iqr_temp_size_adj_lmer_5 <- boot.ci(boot_mass_iqr_temp_size_adj_lmer,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 6) # CI for 3rd betas
print(bt_ci_mass_iqr_temp_size_adj_lmer_5)


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

## Check diagnostics for the full model
plot(mass_iqr_temp_size_small_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_size_small_adj_lmer))
  qqline(resid(mass_iqr_temp_size_small_adj_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_size_small_adj_lmer))

summary(mass_iqr_temp_size_small_adj_lmer)
confint(mass_iqr_temp_size_small_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_size_small_adj_lmer <- bootMer(x = mass_iqr_temp_size_small_adj_lmer,
                                                FUN = fixef, nsim = 2000,
                                                seed = 632760,
                                                use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_size_small_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_size_small_adj_lmer <- boot.ci(boot_mass_iqr_temp_size_small_adj_lmer,
                                                 type = c('perc', 'norm', 'basic'),
                                                 index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_size_small_adj_lmer)


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

## Check diagnostics for the full model
plot(mass_iqr_temp_size_big_adj_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_temp_size_big_adj_lmer))
  qqline(resid(mass_iqr_temp_size_big_adj_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_temp_size_big_adj_lmer))

summary(mass_iqr_temp_size_big_adj_lmer)
confint(mass_iqr_temp_size_big_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_temp_size_big_adj_lmer <- bootMer(x = mass_iqr_temp_size_big_adj_lmer,
                                                FUN = fixef, nsim = 2000,
                                                seed = 632760,
                                                use.u = F, type = 'parametric')
tidy(boot_mass_iqr_temp_size_big_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_temp_size_big_adj_lmer <- boot.ci(boot_mass_iqr_temp_size_big_adj_lmer,
                                                 type = c('perc', 'norm', 'basic'),
                                                 index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_temp_size_big_adj_lmer)


###############################################################################
##############           Sensitive development periods           ##############
###############################################################################

######################## 6 days cut off #######################################

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

### Before thermo no outliers
before_thermo_outliers_removed <- late_nestling_parent_care[-c(4, 39, 40), ]

## Unadjusted
mass_min_before_thermo_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp) + (1|fsite) +
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


## Adjusted
mass_min_before_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp) + scale(nestling_number) +
                                          scale(days_summer) + (1|fsite) +
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


### After thermo unadjusted outliers removed
mass_min_after_thermo_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_min_temp) + (1|fsite) +
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


### After thermo adjusted outliers removed
mass_min_after_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_min_temp) + scale(nestling_number) +
                                               scale(days_summer) + (1|fsite) +
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

# Calculate R squared 
r.squaredGLMM(mass_min_both_thermo_lmer)

# Calculate VIF
vif(mass_min_both_thermo_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_both_thermo_lmer <- bootMer(x = mass_min_both_thermo_lmer,
                                                     FUN = fixef, nsim = 2000,
                                                     seed = 632760,
                                                     use.u = F, type = 'parametric')
tidy(boot_mass_min_both_thermo_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_both_thermo_lmer <- boot.ci(boot_mass_min_both_thermo_lmer,
                                                      type = c('perc', 'norm', 'basic'),
                                                      index = 2) # CI for 1st betas
print(bt_ci_mass_min_both_thermo_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_min_both_thermo_lmer_2 <- boot.ci(boot_mass_min_both_thermo_lmer,
                                           type = c('perc', 'norm', 'basic'),
                                           index = 3) # CI for 2nd beta
print(bt_ci_mass_min_both_thermo_lmer_2)


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

# Calculate R squared 
r.squaredGLMM(mass_min_both_thermo_adj_lmer)

# Calculate VIF
vif(mass_min_both_thermo_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_both_thermo_adj_lmer <- bootMer(x = mass_min_both_thermo_adj_lmer,
                                          FUN = fixef, nsim = 2000,
                                          seed = 632760,
                                          use.u = F, type = 'parametric')
tidy(boot_mass_min_both_thermo_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_both_thermo_adj_lmer <- boot.ci(boot_mass_min_both_thermo_adj_lmer,
                                           type = c('perc', 'norm', 'basic'),
                                           index = 2) # CI for 1st betas
print(bt_ci_mass_min_both_thermo_adj_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_min_both_thermo_adj_lmer_2 <- boot.ci(boot_mass_min_both_thermo_adj_lmer,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 3) # CI for 2nd beta
print(bt_ci_mass_min_both_thermo_adj_lmer_2)


### Before and after thermo outliers removed
## Unadjusted
mass_min_both_thermo_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp) + 
                                    scale(thermo_aft_min_temp) + (1|fsite) +
                                    (1|fnest_id), 
                                  data = subset(before_thermo_outliers_removed,
                                                !is.na(x = thermo_bef_min_temp) &
                                                  !is.na(x = thermo_aft_min_temp)&
                                                  !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_both_thermo_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_min_both_thermo_noout_lmer))
  qqline(resid(mass_min_both_thermo_noout_lmer))}
# Histogram of residuals
hist(resid(mass_min_both_thermo_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_min_both_thermo_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_min_both_thermo_noout_lmer, vars=c("Studentized"))

summary(mass_min_both_thermo_noout_lmer)
confint(mass_min_both_thermo_noout_lmer)

# Calculate R squared 
r.squaredGLMM(mass_min_both_thermo_noout_lmer)

# Calculate VIF
vif(mass_min_both_thermo_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_both_thermo_noout_lmer <- bootMer(x = mass_min_both_thermo_noout_lmer,
                                          FUN = fixef, nsim = 2000,
                                          seed = 632760,
                                          use.u = F, type = 'parametric')
tidy(boot_mass_min_both_thermo_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_both_thermo_noout_lmer <- boot.ci(boot_mass_min_both_thermo_noout_lmer,
                                           type = c('perc', 'norm', 'basic'),
                                           index = 2) # CI for 1st betas
print(bt_ci_mass_min_both_thermo_noout_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_min_both_thermo_noout_lmer_2 <- boot.ci(boot_mass_min_both_thermo_noout_lmer,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 3) # CI for 2nd beta
print(bt_ci_mass_min_both_thermo_noout_lmer_2)


## Adjusted
mass_min_both_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp) + 
                                        scale(thermo_aft_min_temp) + scale(nestling_number) +
                                        scale(days_summer) + (1|fsite) +
                                        (1|fnest_id), 
                                      data = subset(before_thermo_outliers_removed,
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

# Calculate R squared 
r.squaredGLMM(mass_min_both_thermo_adj_noout_lmer)

# Calculate VIF
vif(mass_min_both_thermo_adj_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_both_thermo_adj_noout_lmer <- bootMer(x = mass_min_both_thermo_adj_noout_lmer,
                                                FUN = fixef, nsim = 2000,
                                                seed = 632760,
                                                use.u = F, type = 'parametric')
tidy(boot_mass_min_both_thermo_adj_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_both_thermo_adj_noout_lmer <- boot.ci(boot_mass_min_both_thermo_adj_noout_lmer,
                                                 type = c('perc', 'norm', 'basic'),
                                                 index = 2) # CI for 1st betas
print(bt_ci_mass_min_both_thermo_adj_noout_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_min_both_thermo_adj_noout_lmer_2 <- boot.ci(boot_mass_min_both_thermo_adj_noout_lmer,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 3) # CI for 2nd beta
print(bt_ci_mass_min_both_thermo_adj_noout_lmer_2)

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


### Before thermo no outliers
## Unadjusted
mass_max_before_thermo_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp) + (1|fsite) +
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


# Adjusted
mass_max_before_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp) + scale(nestling_number) +
                                          scale(days_summer) + (1|fsite) +
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


### After thermo outliers removed
## Unadjusted
mass_max_after_thermo_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_max_temp) + (1|fsite) +
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

# Adjusted
mass_max_after_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_max_temp) + scale(nestling_number) +
                                         scale(days_summer) + (1|fsite) +
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

# Calculate R squared 
r.squaredGLMM(mass_max_both_thermo_lmer)

# Calcilate VIF
vif(mass_max_both_thermo_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_both_thermo_lmer <- bootMer(x = mass_max_both_thermo_lmer,
                                                    FUN = fixef, nsim = 2000,
                                                    seed = 632760,
                                                    use.u = F, type = 'parametric')
tidy(boot_mass_max_both_thermo_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_both_thermo_lmer <- boot.ci(boot_mass_max_both_thermo_lmer,
                                                     type = c('perc', 'norm', 'basic'),
                                                     index = 2) # CI for 1st betas
print(bt_ci_mass_max_both_thermo_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_max_both_thermo_lmer_2 <- boot.ci(boot_mass_max_both_thermo_lmer,
                                                       type = c('perc', 'norm', 'basic'),
                                                       index = 3) # CI for 2nd beta
print(bt_ci_mass_max_both_thermo_lmer_2)


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
vif(mass_max_both_thermo_adj_lmer)

# Calculate R squared 
r.squaredGLMM(mass_max_both_thermo_adj_lmer)

# Calcilate VIF
vif(mass_max_both_thermo_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_both_thermo_adj_lmer <- bootMer(x = mass_max_both_thermo_adj_lmer,
                                          FUN = fixef, nsim = 2000,
                                          seed = 632760,
                                          use.u = F, type = 'parametric')
tidy(boot_mass_max_both_thermo_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_both_thermo_adj_lmer <- boot.ci(boot_mass_max_both_thermo_adj_lmer,
                                           type = c('perc', 'norm', 'basic'),
                                           index = 2) # CI for 1st betas
print(bt_ci_mass_max_both_thermo_adj_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_max_both_thermo_adj_lmer_2 <- boot.ci(boot_mass_max_both_thermo_adj_lmer,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 3) # CI for 2nd beta
print(bt_ci_mass_max_both_thermo_adj_lmer_2)

### Before and after thermo no outliers
## Unadjusted
mass_max_both_thermo_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp) + 
                                    scale(thermo_aft_max_temp) + (1|fsite) +
                                    (1|fnest_id), 
                                  data = subset(mass_outliers_removed,
                                                !is.na(x = thermo_bef_max_temp) &
                                                  !is.na(x = thermo_aft_max_temp)&
                                                  !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_both_thermo_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_max_both_thermo_noout_lmer))
  qqline(resid(mass_max_both_thermo_noout_lmer))}
# Histogram of residuals
hist(resid(mass_max_both_thermo_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_max_both_thermo_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_max_both_thermo_noout_lmer, vars=c("Studentized"))

summary(mass_max_both_thermo_noout_lmer)
confint(mass_max_both_thermo_noout_lmer)

# Calculate R squared 
r.squaredGLMM(mass_max_both_thermo_noout_lmer)

# Calcilate VIF
vif(mass_max_both_thermo_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_both_thermo_noout_lmer <- bootMer(x = mass_max_both_thermo_noout_lmer,
                                          FUN = fixef, nsim = 2000,
                                          seed = 632760,
                                          use.u = F, type = 'parametric')
tidy(boot_mass_max_both_thermo_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_both_thermo_noout_lmer <- boot.ci(boot_mass_max_both_thermo_noout_lmer,
                                           type = c('perc', 'norm', 'basic'),
                                           index = 2) # CI for 1st betas
print(bt_ci_mass_max_both_thermo_noout_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_max_both_thermo_noout_lmer_2 <- boot.ci(boot_mass_max_both_thermo_noout_lmer,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 3) # CI for 2nd beta
print(bt_ci_mass_max_both_thermo_noout_lmer_2)

## Adjusted
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

# Calculate R squared 
r.squaredGLMM(mass_max_both_thermo_adj_noout_lmer)

# Calcilate VIF
vif(mass_max_both_thermo_adj_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_both_thermo_adj_noout_lmer <- bootMer(x = mass_max_both_thermo_adj_noout_lmer,
                                                FUN = fixef, nsim = 2000,
                                                seed = 632760,
                                                use.u = F, type = 'parametric')
tidy(boot_mass_max_both_thermo_adj_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_both_thermo_adj_noout_lmer <- boot.ci(boot_mass_max_both_thermo_adj_noout_lmer,
                                                 type = c('perc', 'norm', 'basic'),
                                                 index = 2) # CI for 1st betas
print(bt_ci_mass_max_both_thermo_adj_noout_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_max_both_thermo_adj_noout_lmer_2 <- boot.ci(boot_mass_max_both_thermo_adj_noout_lmer,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 3) # CI for 2nd beta
print(bt_ci_mass_max_both_thermo_adj_noout_lmer_2)


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


### Outliers removed
## Unadjusted model
mass_iqr_before_thermo_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp) + 
                                      (1|fsite) +
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


## Adjusted model
mass_iqr_before_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp) + scale(nestling_number) +
                                          scale(days_summer) +  
                                          (1|fsite) +
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

### Outliers removed
## Unadjusted
mass_iqr_after_thermo_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_iqr_temp) + 
                                     (1|fsite) +
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


### After thermo adjusted 
mass_iqr_after_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_aft_iqr_temp) + scale(nestling_number) +
                                         scale(days_summer) + 
                                         (1|fsite) +
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

# Calculate R squared 
r.squaredGLMM(mass_iqr_both_thermo_lmer)

# Calculate VIF
vif(mass_iqr_both_thermo_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_both_thermo_lmer <- bootMer(x = mass_iqr_both_thermo_lmer,
                                                    FUN = fixef, nsim = 2000,
                                                    seed = 632760,
                                                    use.u = F, type = 'parametric')
tidy(boot_mass_iqr_both_thermo_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_both_thermo_lmer <- boot.ci(boot_mass_iqr_both_thermo_lmer,
                                                     type = c('perc', 'norm', 'basic'),
                                                     index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_both_thermo_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_iqr_both_thermo_lmer_2 <- boot.ci(boot_mass_iqr_both_thermo_lmer,
                                                       type = c('perc', 'norm', 'basic'),
                                                       index = 3) # CI for 2nd beta
print(bt_ci_mass_iqr_both_thermo_lmer_2)


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

# Calculate R squared 
r.squaredGLMM(mass_iqr_both_thermo_adj_lmer)

# Calculate VIF
vif(mass_iqr_both_thermo_adj_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_both_thermo_adj_lmer <- bootMer(x = mass_iqr_both_thermo_adj_lmer,
                                          FUN = fixef, nsim = 2000,
                                          seed = 632760,
                                          use.u = F, type = 'parametric')
tidy(boot_mass_iqr_both_thermo_adj_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_both_thermo_adj_lmer <- boot.ci(boot_mass_iqr_both_thermo_adj_lmer,
                                           type = c('perc', 'norm', 'basic'),
                                           index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_both_thermo_adj_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_iqr_both_thermo_adj_lmer_2 <- boot.ci(boot_mass_iqr_both_thermo_adj_lmer,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 3) # CI for 2nd beta
print(bt_ci_mass_iqr_both_thermo_adj_lmer_2)


### Outliers removed
## Unadjusted
mass_iqr_both_thermo_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp) + 
                                    scale(thermo_aft_iqr_temp) + (1|fsite) +
                                    (1|fnest_id), 
                                  data = subset(before_thermo_outliers_removed,
                                                !is.na(x = thermo_bef_iqr_temp) &
                                                  !is.na(x = thermo_aft_iqr_temp)&
                                                  !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_both_thermo_noout_lmer)
# Normal QQplot
{qqnorm(resid(mass_iqr_both_thermo_noout_lmer))
  qqline(resid(mass_iqr_both_thermo_noout_lmer))}
# Histogram of residuals
hist(resid(mass_iqr_both_thermo_noout_lmer))
# Checking for influential outliers
infIndexPlot(mass_iqr_both_thermo_noout_lmer, vars=c("Cook"))
infIndexPlot(mass_iqr_both_thermo_noout_lmer, vars=c("Studentized"))

summary(mass_iqr_both_thermo_noout_lmer)
confint(mass_iqr_both_thermo_noout_lmer)

# Calculate R squared 
r.squaredGLMM(mass_iqr_both_thermo_noout_lmer)

# Calculate VIF
vif(mass_iqr_both_thermo_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_both_thermo_noout_lmer <- bootMer(x = mass_iqr_both_thermo_noout_lmer,
                                          FUN = fixef, nsim = 2000,
                                          seed = 632760,
                                          use.u = F, type = 'parametric')
tidy(boot_mass_iqr_both_thermo_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_both_thermo_noout_lmer <- boot.ci(boot_mass_iqr_both_thermo_noout_lmer,
                                           type = c('perc', 'norm', 'basic'),
                                           index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_both_thermo_noout_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_iqr_both_thermo_noout_lmer_2 <- boot.ci(boot_mass_iqr_both_thermo_noout_lmer,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 3) # CI for 2nd beta
print(bt_ci_mass_iqr_both_thermo_noout_lmer_2)


## Before and after thermo adjusted
mass_iqr_both_thermo_adj_noout_lmer <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp) + 
                                        scale(thermo_aft_iqr_temp) + scale(nestling_number) +
                                        scale(days_summer) +  (1|fsite) +
                                        (1|fnest_id), 
                                      data = subset(mass_outliers_removed,
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

# Calculate R squared 
r.squaredGLMM(mass_iqr_both_thermo_adj_noout_lmer)

# Calculate VIF
vif(mass_iqr_both_thermo_adj_noout_lmer)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_both_thermo_adj_noout_lmer <- bootMer(x = mass_iqr_both_thermo_adj_noout_lmer,
                                                FUN = fixef, nsim = 2000,
                                                seed = 632760,
                                                use.u = F, type = 'parametric')
tidy(boot_mass_iqr_both_thermo_adj_noout_lmer) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_both_thermo_adj_noout_lmer <- boot.ci(boot_mass_iqr_both_thermo_adj_noout_lmer,
                                                 type = c('perc', 'norm', 'basic'),
                                                 index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_both_thermo_adj_noout_lmer)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_iqr_both_thermo_adj_noout_lmer_2 <- boot.ci(boot_mass_iqr_both_thermo_adj_noout_lmer,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 3) # CI for 2nd beta
print(bt_ci_mass_iqr_both_thermo_adj_noout_lmer_2)



############################ 5 days cut off ###################################
#### Minimum temp
### Before thermo
mass_min_before_thermo_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp_5) + (1|fsite) +
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


### Before thermo adjusted
mass_min_before_thermo_adj_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp_5) + scale(nestling_number) +
                                          scale(days_summer) + (1|fsite) +
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


### After thermo
mass_min_after_thermo_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_aft_min_temp_5) + (1|fsite) +
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


### After thermo adjusted
mass_min_after_thermo_adj_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_aft_min_temp_5) + scale(nestling_number) +
                                         scale(days_summer) + (1|fsite) +
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

### Before and after thermo
mass_min_both_thermo_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp_5) + 
                                    scale(thermo_aft_min_temp_5) + (1|fsite) +
                                    (1|fnest_id), 
                                  data = subset(late_nestling_parent_care,
                                                !is.na(x = thermo_bef_min_temp_5) &
                                                  !is.na(x = thermo_aft_min_temp_5)&
                                                  !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_both_thermo_lmer_5)
# Normal QQplot
{qqnorm(resid(mass_min_both_thermo_lmer_5))
  qqline(resid(mass_min_both_thermo_lmer_5))}
# Histogram of residuals
hist(resid(mass_min_both_thermo_lmer_5))
# Checking for influential outliers
infIndexPlot(mass_min_both_thermo_lmer_5, vars=c("Cook"))
infIndexPlot(mass_min_both_thermo_lmer_5, vars=c("Studentized"))

summary(mass_min_both_thermo_lmer_5)
confint(mass_min_both_thermo_lmer_5)

# Calculate R squared 
r.squaredGLMM(mass_min_both_thermo_lmer_5)

# Calculate VIF
vif(mass_min_both_thermo_lmer_5)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_both_thermo_lmer_5 <- bootMer(x = mass_min_both_thermo_lmer_5,
                                                    FUN = fixef, nsim = 2000,
                                                    seed = 632760,
                                                    use.u = F, type = 'parametric')
tidy(boot_mass_min_both_thermo_lmer_5) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_both_thermo_lmer_5 <- boot.ci(boot_mass_min_both_thermo_lmer_5,
                                                     type = c('perc', 'norm', 'basic'),
                                                     index = 2) # CI for 1st betas
print(bt_ci_mass_min_both_thermo_lmer_5)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_min_both_thermo_lmer_5_2 <- boot.ci(boot_mass_min_both_thermo_lmer_5,
                                                       type = c('perc', 'norm', 'basic'),
                                                       index = 3) # CI for 2nd beta
print(bt_ci_mass_min_both_thermo_lmer_5_2)


### Before and after thermo adjusted
mass_min_both_thermo_adj_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp_5) + 
                                        scale(thermo_aft_min_temp_5) + scale(nestling_number) +
                                        scale(days_summer) + (1|fsite) +
                                        (1|fnest_id), 
                                      data = subset(late_nestling_parent_care,
                                                    !is.na(x = thermo_bef_min_temp_5) &
                                                      !is.na(x = thermo_aft_min_temp_5)&
                                                      !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_both_thermo_adj_lmer_5)
# Normal QQplot
{qqnorm(resid(mass_min_both_thermo_adj_lmer_5))
  qqline(resid(mass_min_both_thermo_adj_lmer_5))}
# Histogram of residuals
hist(resid(mass_min_both_thermo_adj_lmer_5))
# Checking for influential outliers
infIndexPlot(mass_min_both_thermo_adj_lmer_5, vars=c("Cook"))
infIndexPlot(mass_min_both_thermo_adj_lmer_5, vars=c("Studentized"))

summary(mass_min_both_thermo_adj_lmer_5)
confint(mass_min_both_thermo_adj_lmer_5)

# Calculate R squared 
r.squaredGLMM(mass_min_both_thermo_adj_lmer_5)

# Calculate VIF
vif(mass_min_both_thermo_adj_lmer_5)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_both_thermo_adj_lmer_5 <- bootMer(x = mass_min_both_thermo_adj_lmer_5,
                                            FUN = fixef, nsim = 2000,
                                            seed = 632760,
                                            use.u = F, type = 'parametric')
tidy(boot_mass_min_both_thermo_adj_lmer_5) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_both_thermo_adj_lmer_5 <- boot.ci(boot_mass_min_both_thermo_adj_lmer_5,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 2) # CI for 1st betas
print(bt_ci_mass_min_both_thermo_adj_lmer_5)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_min_both_thermo_adj_lmer_5_2 <- boot.ci(boot_mass_min_both_thermo_adj_lmer_5,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 3) # CI for 2nd beta
print(bt_ci_mass_min_both_thermo_adj_lmer_5_2)


#### Maximum temp
### Before thermo
mass_max_before_thermo_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp_5) + (1|fsite) +
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

### Before thermo adjusted
mass_max_before_thermo_adj_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp_5) + scale(nestling_number) +
                                          scale(days_summer) + (1|fsite) +
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


### After thermo
mass_max_after_thermo_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_aft_max_temp_5) + (1|fsite) +
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


### After thermo adjusted
mass_max_after_thermo_adj_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_aft_max_temp_5) + scale(nestling_number) +
                                         scale(days_summer) + (1|fsite) +
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


### Before and after thermo
mass_max_both_thermo_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp_5) + 
                                    scale(thermo_aft_max_temp_5) + (1|fsite) +
                                    (1|fnest_id), 
                                  data = subset(late_nestling_parent_care,
                                                !is.na(x = thermo_bef_max_temp_5) &
                                                  !is.na(x = thermo_aft_max_temp_5)&
                                                  !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_both_thermo_lmer_5)
# Normal QQplot
{qqnorm(resid(mass_max_both_thermo_lmer_5))
  qqline(resid(mass_max_both_thermo_lmer_5))}
# Histogram of residuals
hist(resid(mass_max_both_thermo_lmer_5))
# Checking for influential outliers
infIndexPlot(mass_max_both_thermo_lmer_5, vars=c("Cook"))
infIndexPlot(mass_max_both_thermo_lmer_5, vars=c("Studentized"))

summary(mass_max_both_thermo_lmer_5)
confint(mass_max_both_thermo_lmer_5)

# Calculate R squared 
r.squaredGLMM(mass_max_both_thermo_lmer_5)

# Calculate VIF
vif(mass_max_both_thermo_lmer_5)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_both_thermo_lmer_5 <- bootMer(x = mass_max_both_thermo_lmer_5,
                                                FUN = fixef, nsim = 2000,
                                                seed = 632760,
                                                use.u = F, type = 'parametric')
tidy(boot_mass_max_both_thermo_lmer_5) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_both_thermo_lmer_5 <- boot.ci(boot_mass_max_both_thermo_lmer_5,
                                                 type = c('perc', 'norm', 'basic'),
                                                 index = 2) # CI for 1st betas
print(bt_ci_mass_max_both_thermo_lmer_5)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_max_both_thermo_lmer_5_2 <- boot.ci(boot_mass_max_both_thermo_lmer_5,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 3) # CI for 2nd beta
print(bt_ci_mass_max_both_thermo_lmer_5_2)


### Before and after thermo adjusted
mass_max_both_thermo_adj_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp_5) + 
                                        scale(thermo_aft_max_temp_5) + scale(nestling_number) +
                                        scale(days_summer) +  (1|fsite) +
                                        (1|fnest_id), 
                                      data = subset(late_nestling_parent_care,
                                                    !is.na(x = thermo_bef_max_temp_5) &
                                                      !is.na(x = thermo_aft_max_temp_5)&
                                                      !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_both_thermo_adj_lmer_5)
# Normal QQplot
{qqnorm(resid(mass_max_both_thermo_adj_lmer_5))
  qqline(resid(mass_max_both_thermo_adj_lmer_5))}
# Histogram of residuals
hist(resid(mass_max_both_thermo_adj_lmer_5))
# Checking for influential outliers
infIndexPlot(mass_max_both_thermo_adj_lmer_5, vars=c("Cook"))
infIndexPlot(mass_max_both_thermo_adj_lmer_5, vars=c("Studentized"))

summary(mass_max_both_thermo_adj_lmer_5)
confint(mass_max_both_thermo_adj_lmer_5)

# Calculate R squared 
r.squaredGLMM(mass_max_both_thermo_adj_lmer_5)

# Calculate VIF
vif(mass_max_both_thermo_adj_lmer_5)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_both_thermo_adj_lmer_5 <- bootMer(x = mass_max_both_thermo_adj_lmer_5,
                                            FUN = fixef, nsim = 2000,
                                            seed = 632760,
                                            use.u = F, type = 'parametric')
tidy(boot_mass_max_both_thermo_adj_lmer_5) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_both_thermo_adj_lmer_5 <- boot.ci(boot_mass_max_both_thermo_adj_lmer_5,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 2) # CI for 1st betas
print(bt_ci_mass_max_both_thermo_adj_lmer_5)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_max_both_thermo_adj_lmer_5_2 <- boot.ci(boot_mass_max_both_thermo_adj_lmer_5,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 3) # CI for 2nd beta
print(bt_ci_mass_max_both_thermo_adj_lmer_5_2)


#### IQR of temp
### Before thermo
mass_iqr_before_thermo_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp_5) + 
                                      (1|fsite) +
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


### Before thermo adjusted
mass_iqr_before_thermo_adj_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp_5) + scale(nestling_number) +
                                          scale(days_summer) +  
                                          (1|fsite) +
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



### After thermo
mass_iqr_after_thermo_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_aft_iqr_temp_5) + 
                                     (1|fsite) +
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


### After thermo adjusted
mass_iqr_after_thermo_adj_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_aft_iqr_temp_5) + scale(nestling_number) +
                                         scale(days_summer) + 
                                         (1|fsite) +
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


### Before and after thermo
mass_iqr_both_thermo_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp_5) + 
                                    scale(thermo_aft_iqr_temp_5) + (1|fsite) +
                                    (1|fnest_id), 
                                  data = subset(late_nestling_parent_care,
                                                !is.na(x = thermo_bef_iqr_temp_5) &
                                                  !is.na(x = thermo_aft_iqr_temp_5)&
                                                  !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_both_thermo_lmer_5)
# Normal QQplot
{qqnorm(resid(mass_iqr_both_thermo_lmer_5))
  qqline(resid(mass_iqr_both_thermo_lmer_5))}
# Histogram of residuals
hist(resid(mass_iqr_both_thermo_lmer_5))
# Checking for influential outliers
infIndexPlot(mass_iqr_both_thermo_lmer_5, vars=c("Cook"))
infIndexPlot(mass_iqr_both_thermo_lmer_5, vars=c("Studentized"))

summary(mass_iqr_both_thermo_lmer_5)
confint(mass_iqr_both_thermo_lmer_5)

# Calculate R squared 
r.squaredGLMM(mass_iqr_both_thermo_lmer_5)

# Calculate VIF
vif(mass_iqr_both_thermo_lmer_5)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_both_thermo_lmer_5 <- bootMer(x = mass_iqr_both_thermo_lmer_5,
                                            FUN = fixef, nsim = 2000,
                                            seed = 632760,
                                            use.u = F, type = 'parametric')
tidy(boot_mass_iqr_both_thermo_lmer_5) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_both_thermo_lmer_5 <- boot.ci(boot_mass_iqr_both_thermo_lmer_5,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_both_thermo_lmer_5)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_iqr_both_thermo_lmer_5_2 <- boot.ci(boot_mass_iqr_both_thermo_lmer_5,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 3) # CI for 2nd beta
print(bt_ci_mass_iqr_both_thermo_lmer_5_2)



### Before and after thermo adjusted
mass_iqr_both_thermo_adj_lmer_5 <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp_5) + 
                                        scale(thermo_aft_iqr_temp_5) + scale(nestling_number) +
                                        scale(days_summer) +  (1|fsite) +
                                        (1|fnest_id), 
                                      data = subset(late_nestling_parent_care,
                                                    !is.na(x = thermo_bef_iqr_temp_5) &
                                                      !is.na(x = thermo_aft_iqr_temp_5)&
                                                      !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_both_thermo_adj_lmer_5)
# Normal QQplot
{qqnorm(resid(mass_iqr_both_thermo_adj_lmer_5))
  qqline(resid(mass_iqr_both_thermo_adj_lmer_5))}
# Histogram of residuals
hist(resid(mass_iqr_both_thermo_adj_lmer_5))
# Checking for influential outliers
infIndexPlot(mass_iqr_both_thermo_adj_lmer_5, vars=c("Cook"))
infIndexPlot(mass_iqr_both_thermo_adj_lmer_5, vars=c("Studentized"))

summary(mass_iqr_both_thermo_adj_lmer_5)
confint(mass_iqr_both_thermo_adj_lmer_5)

# Calculate R squared 
r.squaredGLMM(mass_iqr_both_thermo_adj_lmer_5)

# Calculate VIF
vif(mass_iqr_both_thermo_adj_lmer_5)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_both_thermo_adj_lmer_5 <- bootMer(x = mass_iqr_both_thermo_adj_lmer_5,
                                            FUN = fixef, nsim = 2000,
                                            seed = 632760,
                                            use.u = F, type = 'parametric')
tidy(boot_mass_iqr_both_thermo_adj_lmer_5) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_both_thermo_adj_lmer_5 <- boot.ci(boot_mass_iqr_both_thermo_adj_lmer_5,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_both_thermo_adj_lmer_5)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_iqr_both_thermo_adj_lmer_5_2 <- boot.ci(boot_mass_iqr_both_thermo_adj_lmer_5,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 3) # CI for 2nd beta
print(bt_ci_mass_iqr_both_thermo_adj_lmer_5_2)


############################ 7 days cut off ###################################
#### Minimum temp
### Before thermo
mass_min_before_thermo_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp_7) + (1|fsite) +
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


### Before thermo adjusted
mass_min_before_thermo_adj_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp_7) + scale(nestling_number) +
                                            scale(days_summer) + (1|fsite) +
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


### After thermo
mass_min_after_thermo_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_aft_min_temp_7) + (1|fsite) +
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


### After thermo adjusted
mass_min_after_thermo_adj_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_aft_min_temp_7) + scale(nestling_number) +
                                           scale(days_summer) + (1|fsite) +
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

### Before and after thermo
mass_min_both_thermo_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp_7) + 
                                      scale(thermo_aft_min_temp_7) + (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  !is.na(x = thermo_bef_min_temp_7) &
                                                    !is.na(x = thermo_aft_min_temp_7)&
                                                    !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_both_thermo_lmer_7)
# Normal QQplot
{qqnorm(resid(mass_min_both_thermo_lmer_7))
  qqline(resid(mass_min_both_thermo_lmer_7))}
# Histogram of residuals
hist(resid(mass_min_both_thermo_lmer_7))
# Checking for influential outliers
infIndexPlot(mass_min_both_thermo_lmer_7, vars=c("Cook"))
infIndexPlot(mass_min_both_thermo_lmer_7, vars=c("Studentized"))

summary(mass_min_both_thermo_lmer_7)
confint(mass_min_both_thermo_lmer_7)

# Calculate R squared 
r.squaredGLMM(mass_min_both_thermo_lmer_7)

# Calculate VIF
vif(mass_min_both_thermo_lmer_7)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_both_thermo_lmer_7 <- bootMer(x = mass_min_both_thermo_lmer_7,
                                            FUN = fixef, nsim = 2000,
                                            seed = 632760,
                                            use.u = F, type = 'parametric')
tidy(boot_mass_min_both_thermo_lmer_7) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_both_thermo_lmer_7 <- boot.ci(boot_mass_min_both_thermo_lmer_7,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 2) # CI for 1st betas
print(bt_ci_mass_min_both_thermo_lmer_7)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_min_both_thermo_lmer_7_2 <- boot.ci(boot_mass_min_both_thermo_lmer_7,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 3) # CI for 2nd beta
print(bt_ci_mass_min_both_thermo_lmer_7_2)


### Before and after thermo adjusted
mass_min_both_thermo_adj_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_bef_min_temp_7) + 
                                          scale(thermo_aft_min_temp_7) + scale(nestling_number) +
                                          scale(days_summer) + (1|fsite) +
                                          (1|fnest_id), 
                                        data = subset(late_nestling_parent_care,
                                                      !is.na(x = thermo_bef_min_temp_7) &
                                                        !is.na(x = thermo_aft_min_temp_7)&
                                                        !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_min_both_thermo_adj_lmer_7)
# Normal QQplot
{qqnorm(resid(mass_min_both_thermo_adj_lmer_7))
  qqline(resid(mass_min_both_thermo_adj_lmer_7))}
# Histogram of residuals
hist(resid(mass_min_both_thermo_adj_lmer_7))
# Checking for influential outliers
infIndexPlot(mass_min_both_thermo_adj_lmer_7, vars=c("Cook"))
infIndexPlot(mass_min_both_thermo_adj_lmer_7, vars=c("Studentized"))

summary(mass_min_both_thermo_adj_lmer_7)
confint(mass_min_both_thermo_adj_lmer_7)

# Calculate R squared 
r.squaredGLMM(mass_min_both_thermo_adj_lmer_7)

# Calculate VIF
vif(mass_min_both_thermo_adj_lmer_7)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_min_both_thermo_adj_lmer_7 <- bootMer(x = mass_min_both_thermo_adj_lmer_7,
                                                FUN = fixef, nsim = 2000,
                                                seed = 632760,
                                                use.u = F, type = 'parametric')
tidy(boot_mass_min_both_thermo_adj_lmer_7) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_min_both_thermo_adj_lmer_7 <- boot.ci(boot_mass_min_both_thermo_adj_lmer_7,
                                                 type = c('perc', 'norm', 'basic'),
                                                 index = 2) # CI for 1st betas
print(bt_ci_mass_min_both_thermo_adj_lmer_7)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_min_both_thermo_adj_lmer_7_2 <- boot.ci(boot_mass_min_both_thermo_adj_lmer_7,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 3) # CI for 2nd beta
print(bt_ci_mass_min_both_thermo_adj_lmer_7_2)


#### Maximum temp
### Before thermo
mass_max_before_thermo_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp_7) + (1|fsite) +
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

### Before thermo adjusted
mass_max_before_thermo_adj_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp_7) + scale(nestling_number) +
                                            scale(days_summer) + (1|fsite) +
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


### After thermo
mass_max_after_thermo_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_aft_max_temp_7) + (1|fsite) +
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


### After thermo adjusted
mass_max_after_thermo_adj_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_aft_max_temp_7) + scale(nestling_number) +
                                           scale(days_summer) + (1|fsite) +
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


### Before and after thermo
mass_max_both_thermo_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp_7) + 
                                      scale(thermo_aft_max_temp_7) + (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  !is.na(x = thermo_bef_max_temp_7) &
                                                    !is.na(x = thermo_aft_max_temp_7)&
                                                    !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_both_thermo_lmer_7)
# Normal QQplot
{qqnorm(resid(mass_max_both_thermo_lmer_7))
  qqline(resid(mass_max_both_thermo_lmer_7))}
# Histogram of residuals
hist(resid(mass_max_both_thermo_lmer_7))
# Checking for influential outliers
infIndexPlot(mass_max_both_thermo_lmer_7, vars=c("Cook"))
infIndexPlot(mass_max_both_thermo_lmer_7, vars=c("Studentized"))

summary(mass_max_both_thermo_lmer_7)
confint(mass_max_both_thermo_lmer_7)

# Calculate R squared 
r.squaredGLMM(mass_max_both_thermo_lmer_7)

# Calculate VIF
vif(mass_max_both_thermo_lmer_7)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_both_thermo_lmer_7 <- bootMer(x = mass_max_both_thermo_lmer_7,
                                            FUN = fixef, nsim = 2000,
                                            seed = 632760,
                                            use.u = F, type = 'parametric')
tidy(boot_mass_max_both_thermo_lmer_7) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_both_thermo_lmer_7 <- boot.ci(boot_mass_max_both_thermo_lmer_7,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 2) # CI for 1st betas
print(bt_ci_mass_max_both_thermo_lmer_7)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_max_both_thermo_lmer_7_2 <- boot.ci(boot_mass_max_both_thermo_lmer_7,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 3) # CI for 2nd beta
print(bt_ci_mass_max_both_thermo_lmer_7_2)


### Before and after thermo adjusted
mass_max_both_thermo_adj_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_bef_max_temp_7) + 
                                          scale(thermo_aft_max_temp_7) + scale(nestling_number) +
                                          scale(days_summer) +  (1|fsite) +
                                          (1|fnest_id), 
                                        data = subset(late_nestling_parent_care,
                                                      !is.na(x = thermo_bef_max_temp_7) &
                                                        !is.na(x = thermo_aft_max_temp_7)&
                                                        !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_max_both_thermo_adj_lmer_7)
# Normal QQplot
{qqnorm(resid(mass_max_both_thermo_adj_lmer_7))
  qqline(resid(mass_max_both_thermo_adj_lmer_7))}
# Histogram of residuals
hist(resid(mass_max_both_thermo_adj_lmer_7))
# Checking for influential outliers
infIndexPlot(mass_max_both_thermo_adj_lmer_7, vars=c("Cook"))
infIndexPlot(mass_max_both_thermo_adj_lmer_7, vars=c("Studentized"))

summary(mass_max_both_thermo_adj_lmer_7)
confint(mass_max_both_thermo_adj_lmer_7)

# Calculate R squared 
r.squaredGLMM(mass_max_both_thermo_adj_lmer_7)

# Calculate VIF
vif(mass_max_both_thermo_adj_lmer_7)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_max_both_thermo_adj_lmer_7 <- bootMer(x = mass_max_both_thermo_adj_lmer_7,
                                                FUN = fixef, nsim = 2000,
                                                seed = 632760,
                                                use.u = F, type = 'parametric')
tidy(boot_mass_max_both_thermo_adj_lmer_7) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_max_both_thermo_adj_lmer_7 <- boot.ci(boot_mass_max_both_thermo_adj_lmer_7,
                                                 type = c('perc', 'norm', 'basic'),
                                                 index = 2) # CI for 1st betas
print(bt_ci_mass_max_both_thermo_adj_lmer_7)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_max_both_thermo_adj_lmer_7_2 <- boot.ci(boot_mass_max_both_thermo_adj_lmer_7,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 3) # CI for 2nd beta
print(bt_ci_mass_max_both_thermo_adj_lmer_7_2)


#### IQR of temp
### Before thermo
mass_iqr_before_thermo_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp_7) + 
                                        (1|fsite) +
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


### Before thermo adjusted
mass_iqr_before_thermo_adj_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp_7) + scale(nestling_number) +
                                            scale(days_summer) +  
                                            (1|fsite) +
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



### After thermo
mass_iqr_after_thermo_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_aft_iqr_temp_7) + 
                                       (1|fsite) +
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


### After thermo adjusted
mass_iqr_after_thermo_adj_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_aft_iqr_temp_7) + scale(nestling_number) +
                                           scale(days_summer) + 
                                           (1|fsite) +
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


### Before and after thermo
mass_iqr_both_thermo_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp_7) + 
                                      scale(thermo_aft_iqr_temp_7) + (1|fsite) +
                                      (1|fnest_id), 
                                    data = subset(late_nestling_parent_care,
                                                  !is.na(x = thermo_bef_iqr_temp_7) &
                                                    !is.na(x = thermo_aft_iqr_temp_7)&
                                                    !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_both_thermo_lmer_7)
# Normal QQplot
{qqnorm(resid(mass_iqr_both_thermo_lmer_7))
  qqline(resid(mass_iqr_both_thermo_lmer_7))}
# Histogram of residuals
hist(resid(mass_iqr_both_thermo_lmer_7))
# Checking for influential outliers
infIndexPlot(mass_iqr_both_thermo_lmer_7, vars=c("Cook"))
infIndexPlot(mass_iqr_both_thermo_lmer_7, vars=c("Studentized"))

summary(mass_iqr_both_thermo_lmer_7)
confint(mass_iqr_both_thermo_lmer_7)

# Calculate R squared 
r.squaredGLMM(mass_iqr_both_thermo_lmer_7)

# Calculate VIF
vif(mass_iqr_both_thermo_lmer_7)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_both_thermo_lmer_7 <- bootMer(x = mass_iqr_both_thermo_lmer_7,
                                            FUN = fixef, nsim = 2000,
                                            seed = 632760,
                                            use.u = F, type = 'parametric')
tidy(boot_mass_iqr_both_thermo_lmer_7) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_both_thermo_lmer_7 <- boot.ci(boot_mass_iqr_both_thermo_lmer_7,
                                             type = c('perc', 'norm', 'basic'),
                                             index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_both_thermo_lmer_7)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_iqr_both_thermo_lmer_7_2 <- boot.ci(boot_mass_iqr_both_thermo_lmer_7,
                                               type = c('perc', 'norm', 'basic'),
                                               index = 3) # CI for 2nd beta
print(bt_ci_mass_iqr_both_thermo_lmer_7_2)



### Before and after thermo adjusted
mass_iqr_both_thermo_adj_lmer_7 <- lmer(mass_pre_obs ~ scale(thermo_bef_iqr_temp_7) + 
                                          scale(thermo_aft_iqr_temp_7) + scale(nestling_number) +
                                          scale(days_summer) +  (1|fsite) +
                                          (1|fnest_id), 
                                        data = subset(late_nestling_parent_care,
                                                      !is.na(x = thermo_bef_iqr_temp_7) &
                                                        !is.na(x = thermo_aft_iqr_temp_7)&
                                                        !is.na(x = mass_pre_obs)))

## Check diagnostics for the full model
plot(mass_iqr_both_thermo_adj_lmer_7)
# Normal QQplot
{qqnorm(resid(mass_iqr_both_thermo_adj_lmer_7))
  qqline(resid(mass_iqr_both_thermo_adj_lmer_7))}
# Histogram of residuals
hist(resid(mass_iqr_both_thermo_adj_lmer_7))
# Checking for influential outliers
infIndexPlot(mass_iqr_both_thermo_adj_lmer_7, vars=c("Cook"))
infIndexPlot(mass_iqr_both_thermo_adj_lmer_7, vars=c("Studentized"))

summary(mass_iqr_both_thermo_adj_lmer_7)
confint(mass_iqr_both_thermo_adj_lmer_7)

# Calculate R squared 
r.squaredGLMM(mass_iqr_both_thermo_adj_lmer_7)

# Calculate VIF
vif(mass_iqr_both_thermo_adj_lmer_7)

## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_mass_iqr_both_thermo_adj_lmer_7 <- bootMer(x = mass_iqr_both_thermo_adj_lmer_7,
                                                FUN = fixef, nsim = 2000,
                                                seed = 632760,
                                                use.u = F, type = 'parametric')
tidy(boot_mass_iqr_both_thermo_adj_lmer_7) # beta estimates and SE
# use 'boot' package to generate 95% CI for 1st beta
bt_ci_mass_iqr_both_thermo_adj_lmer_7 <- boot.ci(boot_mass_iqr_both_thermo_adj_lmer_7,
                                                 type = c('perc', 'norm', 'basic'),
                                                 index = 2) # CI for 1st betas
print(bt_ci_mass_iqr_both_thermo_adj_lmer_7)

# use 'boot' package to generate 95% CI for 2nd beta
bt_ci_mass_iqr_both_thermo_adj_lmer_7_2 <- boot.ci(boot_mass_iqr_both_thermo_adj_lmer_7,
                                                   type = c('perc', 'norm', 'basic'),
                                                   index = 3) # CI for 2nd beta
print(bt_ci_mass_iqr_both_thermo_adj_lmer_7_2)



###############################################################################
##############                 Model visualizations              ##############
###############################################################################

################ Parental care models ######################
## Minimum temperature

pred_low <- ggpredict(mass_min_temp_blups_low_adj_lmer, c("nest_min_temp"))

plot(pred_low) 

pred_low$feeding <- "1"

pred_med <- ggpredict(mass_min_temp_blups_med_adj_lmer, c("nest_min_temp"))

plot(pred_med)

pred_med$feeding <- "2"

pred_high <- ggpredict(mass_min_temp_blups_high_adj_lmer, c("nest_min_temp"))

plot(pred_high)

pred_high$feeding <- "3"

pred_feeding <- rbind(pred_low, pred_med, pred_high)


sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_min_temp) &
                !is.na(x = feeding_expontd_blups))

cols <- c("#481567FF", "#20A387FF", "#95D840FF")
lines <- c(4, 5, 1)

temp_min_feeding_blups_predicted <- ggplot() +
  geom_line(data = pred_feeding, aes(x = x, y = predicted, col = feeding), size = 2) +
  # geom_ribbon(data = pred_feeding, aes(x = x, y = predicted, 
  # ymin = conf.low, ymax = conf.high, fill = feeding), 
  # alpha = .1, col = NA) +
  geom_point(data = sub, aes(x = nest_min_temp, y = mass_pre_obs,
                             col = as.character(feeding_expontd_blups_strat)), 
             show.legend = FALSE, size = 1.5, alpha = 0.5) +
  theme_classic() +
  labs(x = "Minimum temperature (C)", y = "Nestling mass (g)",
       colour = "Parent feeding level", fill = "Parent feeding level") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  scale_color_manual(values = cols, labels=c("Low", "Med", "High"),
                     aesthetics = c("colour", "fill"))

# Print
print(temp_min_feeding_blups_predicted)


## Maximum temperature
pred_low <- ggpredict(mass_max_temp_blups_low_adj_lmer, c("nest_max_temp"))

plot(pred_low) 

pred_low$feeding <- "1"

pred_med <- ggpredict(mass_max_temp_blups_med_adj_lmer, c("nest_max_temp"))

plot(pred_med)

pred_med$feeding <- "2"

pred_high <- ggpredict(mass_max_temp_blups_high_adj_lmer, c("nest_max_temp"))

plot(pred_high)

pred_high$feeding <- "3"

pred_feeding <- rbind(pred_low, pred_med, pred_high)


sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_max_temp) &
                !is.na(x = feeding_expontd_blups))

temp_max_feeding_blups_predicted <- ggplot() +
  geom_line(data = pred_feeding, aes(x = x, y = predicted, col = feeding), size = 2) +
  # geom_ribbon(data = pred_feeding, aes(x = x, y = predicted, 
  # ymin = conf.low, ymax = conf.high, fill = feeding), 
  # alpha = .1, col = NA) +
  geom_point(data = sub, aes(x = nest_max_temp, y = mass_pre_obs,
                             col = as.character(feeding_expontd_blups_strat)), 
             show.legend = FALSE, size = 1.5, alpha = 0.5) +
  theme_classic() +
  labs(x = "Maximum temperature (C)", y = "Nestling mass (g)",
       colour = "Parent feeding level", fill = "Parent feeding level") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  scale_color_manual(values = cols, labels=c("Low", "Med", "High"),
                     aesthetics = c("colour", "fill")) 


# Print
print(temp_max_feeding_blups_predicted)


## IQR of temperature
pred_low <- ggpredict(mass_iqr_temp_blups_low_adj_lmer, c("nest_iqr_temp"))

plot(pred_low) 

pred_low$feeding <- "1"

pred_med <- ggpredict(mass_iqr_temp_blups_med_adj_lmer, c("nest_iqr_temp"))

plot(pred_med)

pred_med$feeding <- "2"

pred_high <- ggpredict(mass_iqr_temp_blups_high_adj_lmer, c("nest_iqr_temp"))

plot(pred_high)

pred_high$feeding <- "3"

pred_feeding <- rbind(pred_low, pred_med, pred_high)


sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_iqr_temp) &
                !is.na(x = feeding_expontd_blups))


temp_iqr_feeding_blups_predicted <- ggplot() +
  geom_line(data = pred_feeding, aes(x = x, y = predicted, col = feeding), size = 2) +
  # geom_ribbon(data = pred_feeding, aes(x = x, y = predicted, 
  # ymin = conf.low, ymax = conf.high, fill = feeding), 
  # alpha = .1, col = NA) +
  geom_point(data = sub, aes(x = nest_iqr_temp, y = mass_pre_obs,
                             col = as.character(feeding_expontd_blups_strat)), 
             show.legend = FALSE, size = 1.5, alpha = 0.5) +
  theme_classic() +
  labs(x = "IQR of temperature (C)", y = "Nestling mass (g)",
       colour = "Parent feeding level", fill = "Parent feeding level") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  scale_color_manual(values = cols, labels=c("Low", "Med", "High"),
                     aesthetics = c("colour", "fill"))

# Print
print(temp_iqr_feeding_blups_predicted)


## Combined plot
combined_temp_feeding_predictions <- 
  ggarrange(temp_min_feeding_blups_predicted, temp_max_feeding_blups_predicted,
            temp_iqr_feeding_blups_predicted,
            ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom",
            labels = c("a", "b", "c"), font.label = list(size = 18))

# Print
print(combined_temp_feeding_predictions)

# Save
ggsave('combined_temp_feeding_predictions_solid.png', plot = combined_temp_feeding_predictions, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 6.5, 
       height = 12, 
       units = c('in'), dpi = 300, limitsize = TRUE)


################## Relative size models ##################

## Minimum temperature
pred_small <- ggpredict(mass_min_temp_mid_size_small_adj_lmer, c("nest_min_temp"))

plot(pred_small) 

pred_small$size <- "min"

pred_big <- ggpredict(mass_min_temp_mid_size_big_adj_lmer, c("nest_min_temp"))

plot(pred_big)

pred_big$size <- "other"

pred_size <- rbind(pred_small, pred_big)


sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_min_temp) &
                !is.na(x = mid_size_order))

cols <- c("#481567FF","#95D840FF")
lines <- c(4, 1)

temp_min_size_predicted <- ggplot() +
  geom_line(data = pred_size, aes(x = x, y = predicted, col = size), size = 2) +
  # geom_ribbon(data = pred_size, aes(x = x, y = predicted, 
  # ymin = conf.low, ymax = conf.high, fill = size), 
  # alpha = .1, col = NA) +
  geom_point(data = sub, aes(x = nest_min_temp, y = mass_pre_obs,
                             col = as.character(mid_size_order)), 
             show.legend = FALSE, size = 1.5, alpha = 0.5) +
  theme_classic() +
  labs(x = "Minimum temperature (C)", y = "Nestling mass (g)",
       colour = "Relative nestling size", fill = "Relative nestling size") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  scale_color_manual(values = cols, labels=c("Smallest", "Other"),
                     aesthetics = c("colour", "fill")) 


# Print
print(temp_min_size_predicted)


## Maximum temperature
pred_small <- ggpredict(mass_max_temp_mid_size_small_adj_lmer, c("nest_max_temp"))

plot(pred_small) 

pred_small$size <- "min"

pred_big <- ggpredict(mass_max_temp_mid_size_big_adj_lmer, c("nest_max_temp"))

plot(pred_big)

pred_big$size <- "other"

pred_size <- rbind(pred_small, pred_big)


sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_max_temp) &
                !is.na(x = mid_size_order))

cols <- c("#481567FF","#95D840FF")
lines <- c(4, 1)

temp_max_size_predicted <- ggplot() +
  geom_line(data = pred_size, aes(x = x, y = predicted, col = size), size = 2) +
  # geom_ribbon(data = pred_size, aes(x = x, y = predicted, 
  # ymin = conf.low, ymax = conf.high, fill = size), 
  # alpha = .1, col = NA) +
  geom_point(data = sub, aes(x = nest_max_temp, y = mass_pre_obs,
                             col = as.character(mid_size_order)), 
             show.legend = FALSE, size = 1.5, alpha = 0.5) +
  theme_classic() +
  labs(x = "Maximum temperature (C)", y = "Nestling mass (g)",
       colour = "Relative nestling size", fill = "Relative nestling size") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  scale_color_manual(values = cols, labels=c("Smallest", "Other"),
                     aesthetics = c("colour", "fill"))


# Print
print(temp_max_size_predicted)


## IQR of temperature
pred_small <- ggpredict(mass_iqr_temp_mid_size_small_adj_lmer, c("nest_iqr_temp"))

plot(pred_small) 

pred_small$size <- "min"

pred_big <- ggpredict(mass_iqr_temp_mid_size_big_adj_lmer, c("nest_iqr_temp"))

plot(pred_big)

pred_big$size <- "other"

pred_size <- rbind(pred_small, pred_big)


sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_iqr_temp) &
                !is.na(x = mid_size_order))


temp_iqr_size_predicted <- ggplot() +
  geom_line(data = pred_size, aes(x = x, y = predicted, col = size), size = 2) +
  # geom_ribbon(data = pred_size, aes(x = x, y = predicted, 
  # ymin = conf.low, ymax = conf.high, fill = size), 
  # alpha = .1, col = NA) +
  geom_point(data = sub, aes(x = nest_iqr_temp, y = mass_pre_obs,
                             col = as.character(mid_size_order)), 
             show.legend = FALSE, size = 1.5, alpha = 0.5) +
  theme_classic() +
  labs(x = "IQR of temperature (C)", y = "Nestling mass (g)",
       colour = "Relative nestling size", fill = "Relative nestling size") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  scale_color_manual(values = cols, labels=c("Smallest", "Other"),
                     aesthetics = c("colour", "fill")) 


# Print
print(temp_iqr_size_predicted)

## Combined plot
combined_temp_size_predictions <- 
  ggarrange(temp_min_size_predicted, temp_max_size_predicted,
            temp_iqr_size_predicted,
            ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom",
            labels = c("a", "b", "c"), font.label = list(size = 18))

# Print
print(combined_temp_size_predictions)

# Save
ggsave('combined_temp_size_predictions_solid.png', plot = combined_temp_size_predictions, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 6.5, 
       height = 12, 
       units = c('in'), dpi = 300, limitsize = TRUE)


################## Developmental stages model ################################

########### Dot and whisker plots
## Min temp before model: extract estimates and tidy the data frame
mass_min_before_thermo_adj_lmer_est <- tidy(mass_min_before_thermo_adj_lmer,
                                            conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(thermo_bef_min_temp)')

## Rename variables for estimates of interest
mass_min_before_thermo_adj_lmer_est['term'][mass_min_before_thermo_adj_lmer_est['term'] ==
                                              'scale(thermo_bef_min_temp)'] <- 'Before day 6'


## Min temp before model: Label the estimates in data frame
mass_min_before_thermo_adj_lmer_est$model <- c('Before model')

## Min temp after model: extract estimates and tidy the data frame
mass_min_after_thermo_adj_lmer_est <- tidy(mass_min_after_thermo_adj_lmer,
                                           conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(thermo_aft_min_temp)')


## Rename variables for estimates of interest
mass_min_after_thermo_adj_lmer_est['term'][mass_min_after_thermo_adj_lmer_est['term'] ==
                                             'scale(thermo_aft_min_temp)'] <- 'After day 6'


## Min temp before model: Label the estimates in data frame
mass_min_after_thermo_adj_lmer_est$model <- c('After model')


## Min temp both model: extract estimates and tidy the data frame
mass_min_both_thermo_adj_lmer_est <- tidy(mass_min_both_thermo_adj_lmer,
                                          conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(thermo_bef_min_temp)' | term == 'scale(thermo_aft_min_temp)')


## Rename variables for estimates of interest
mass_min_both_thermo_adj_lmer_est['term'][mass_min_both_thermo_adj_lmer_est['term'] ==
                                            'scale(thermo_bef_min_temp)'] <- 'Before day 6'


mass_min_both_thermo_adj_lmer_est['term'][mass_min_both_thermo_adj_lmer_est['term'] ==
                                            'scale(thermo_aft_min_temp)'] <- 'After day 6'


## Min temp before model: Label the estimates in data frame
mass_min_both_thermo_adj_lmer_est$model <- c('Both model')

## Combine regression estimates into a tidy table
mass_min_both_thermo_adj_lmer_tbl <- bind_rows(mass_min_before_thermo_adj_lmer_est,
                                               mass_min_after_thermo_adj_lmer_est,
                                               mass_min_both_thermo_adj_lmer_est)

## Re-code *nominal* factor (with ordered levels)  
# Set levels (odering) of 'model' variable
mass_min_both_thermo_adj_lmer_tbl <-
  transform(mass_min_both_thermo_adj_lmer_tbl,
            model = factor(model,
                           levels = c('Before model',
                                      'After model',
                                      'Both model')))

mass_min_both_thermo_adj_lmer_tbl <-
  transform(mass_min_both_thermo_adj_lmer_tbl,
            term = factor(term,
                          levels = c('Before day 6', 
                                     'After day 6')))

## Graph of estimates 
mass_min_both_thermo_adj_lmer_plot <-
  ggplot(mass_min_both_thermo_adj_lmer_tbl, aes(x = term, y = estimate,
                                                color = model)) +
  theme_classic() + 
  geom_hline(yintercept = 0, color = 'red',
             linetype = 2) + # line at null behind coefs
  geom_point(size = 6,
             position=position_dodge(width = 0.0)) +
  geom_errorbar(aes(ymin=(conf.low),
                    ymax=(conf.high)), width=.1,
                position=position_dodge(0.0), size = 1) +
  scale_color_manual(values=c('#481567FF', '#95D840FF','#20A387FF')) +
  facet_grid(cols = vars(model), scales = 'free', space = 'free') +
  theme(strip.background =element_rect(fill= 'white'))+
  theme(strip.text = element_text(colour = 'black')) +
  #coord_flip() + # flip x and y axes
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18)) +
  theme(text = element_text(size=18)) +
  # bold and size title and axes labels
  theme(legend.position = 'none') +
  theme(axis.ticks = element_blank()) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size = 16, angle=45,
                                   margin = margin(t = 50, r = 0,
                                                   b = 0, l = 0), face = 'bold'), 
        axis.text.y = element_text(size = 16, angle=0,
                                   margin = margin(t = 25, r = 0,
                                                   b = 0, l = 10)),
        axis.title = element_text(size = 18), legend.title=element_blank(),
        legend.text=element_text(size=18),
        legend.position = 'none', #c(0.91, 0.94),
        legend.key = element_blank()) +
  labs(colour = "Parent feeding level") +
  xlab(expression(italic('Standardized minimum temperature (SD)'))) +
  ylab(expression
       (atop(bold('Beta estimate and 95% CI'),
             paste(italic('Nestling mass (g)')))))
# remove axis ticks

print(mass_min_both_thermo_adj_lmer_plot)



## Max temp before model: extract estimates and tidy the data frame
mass_max_before_thermo_adj_lmer_est <- tidy(mass_max_before_thermo_adj_lmer,
                                            conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(thermo_bef_max_temp)')

## Rename variables for estimates of interest
mass_max_before_thermo_adj_lmer_est['term'][mass_max_before_thermo_adj_lmer_est['term'] ==
                                              'scale(thermo_bef_max_temp)'] <- 'Before day 6'


## Max temp before model: Label the estimates in data frame
mass_max_before_thermo_adj_lmer_est$model <- c('Before model')

## Max temp after model: extract estimates and tidy the data frame
mass_max_after_thermo_adj_lmer_est <- tidy(mass_max_after_thermo_adj_lmer,
                                           conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(thermo_aft_max_temp)')


## Rename variables for estimates of interest
mass_max_after_thermo_adj_lmer_est['term'][mass_max_after_thermo_adj_lmer_est['term'] ==
                                             'scale(thermo_aft_max_temp)'] <- 'After day 6'


## Max temp before model: Label the estimates in data frame
mass_max_after_thermo_adj_lmer_est$model <- c('After model')


## Max temp both model: extract estimates and tidy the data frame
mass_max_both_thermo_adj_lmer_est <- tidy(mass_max_both_thermo_adj_lmer,
                                          conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(thermo_bef_max_temp)' | term == 'scale(thermo_aft_max_temp)')


## Rename variables for estimates of interest
mass_max_both_thermo_adj_lmer_est['term'][mass_max_both_thermo_adj_lmer_est['term'] ==
                                            'scale(thermo_bef_max_temp)'] <- 'Before day 6'


mass_max_both_thermo_adj_lmer_est['term'][mass_max_both_thermo_adj_lmer_est['term'] ==
                                            'scale(thermo_aft_max_temp)'] <- 'After day 6'


## Max temp before model: Label the estimates in data frame
mass_max_both_thermo_adj_lmer_est$model <- c('Both model')

## Combine regression estimates into a tidy table
mass_max_both_thermo_adj_lmer_tbl <- bind_rows(mass_max_before_thermo_adj_lmer_est,
                                               mass_max_after_thermo_adj_lmer_est,
                                               mass_max_both_thermo_adj_lmer_est)

## Re-code *nominal* factor (with ordered levels)  
# Set levels (odering) of 'model' variable
mass_max_both_thermo_adj_lmer_tbl <-
  transform(mass_max_both_thermo_adj_lmer_tbl,
            model = factor(model,
                           levels = c('Before model',
                                      'After model',
                                      'Both model')))

mass_max_both_thermo_adj_lmer_tbl <-
  transform(mass_max_both_thermo_adj_lmer_tbl,
            term = factor(term,
                          levels = c('Before day 6', 
                                     'After day 6')))

## Graph of estimates
mass_max_both_thermo_adj_lmer_plot <-
  ggplot(mass_max_both_thermo_adj_lmer_tbl, aes(x = term, y = estimate,
                                                color = model)) +
  theme_classic() + 
  geom_hline(yintercept = 0, color = 'red',
             linetype = 2) + # line at null behind coefs
  geom_point(size = 6,
             position=position_dodge(width = 0.0)) +
  geom_errorbar(aes(ymin=(conf.low),
                    ymax=(conf.high)), width=.1,
                position=position_dodge(0.0), size = 1) +
  scale_color_manual(values=c('#481567FF', '#95D840FF','#20A387FF')) +
  facet_grid(cols = vars(model), scales = 'free', space = 'free') +
  theme(strip.background =element_rect(fill= 'white'))+
  theme(strip.text = element_text(colour = 'black')) +
  #coord_flip() + # flip x and y axes
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18)) +
  theme(text = element_text(size=18)) +
  # bold and size title and axes labels
  theme(legend.position = 'none') +
  theme(axis.ticks = element_blank()) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size = 16, angle=45,
                                   margin = margin(t = 50, r = 0,
                                                   b = 0, l = 0), face = 'bold'), 
        axis.text.y = element_text(size = 16, angle=0,
                                   margin = margin(t = 25, r = 0,
                                                   b = 0, l = 10)),
        axis.title = element_text(size = 18), legend.title=element_blank(),
        legend.text=element_text(size=18),
        legend.position = 'none', #c(0.91, 0.94),
        legend.key = element_blank()) +
  labs(colour = "Parent feeding level") +
  xlab(expression(italic('Standardized maximum temperature (SD)'))) +
  ylab(expression
       (atop(bold('Beta estimate and 95% CI'),
             paste(italic('Nestling mass (g)')))))


print(mass_max_both_thermo_adj_lmer_plot)



## IQR temp before model: extract estimates and tidy the data frame
mass_iqr_before_thermo_adj_lmer_est <- tidy(mass_iqr_before_thermo_adj_lmer,
                                            conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(thermo_bef_iqr_temp)')

## Rename variables for estimates of interest
mass_iqr_before_thermo_adj_lmer_est['term'][mass_iqr_before_thermo_adj_lmer_est['term'] ==
                                              'scale(thermo_bef_iqr_temp)'] <- 'Before day 6'


## Min temp before model: Label the estimates in data frame
mass_iqr_before_thermo_adj_lmer_est$model <- c('Before model')

## Min temp after model: extract estimates and tidy the data frame
mass_iqr_after_thermo_adj_lmer_est <- tidy(mass_iqr_after_thermo_adj_lmer,
                                           conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(thermo_aft_iqr_temp)')


## Rename variables for estimates of interest
mass_iqr_after_thermo_adj_lmer_est['term'][mass_iqr_after_thermo_adj_lmer_est['term'] ==
                                             'scale(thermo_aft_iqr_temp)'] <- 'After day 6'


## IQR temp before model: Label the estimates in data frame
mass_iqr_after_thermo_adj_lmer_est$model <- c('After model')


## IQR temp both model: extract estimates and tidy the data frame
mass_iqr_both_thermo_adj_lmer_est <- tidy(mass_iqr_both_thermo_adj_lmer,
                                          conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(thermo_bef_iqr_temp)' | term == 'scale(thermo_aft_iqr_temp)')


## Rename variables for estimates of interest
mass_iqr_both_thermo_adj_lmer_est['term'][mass_iqr_both_thermo_adj_lmer_est['term'] ==
                                            'scale(thermo_bef_iqr_temp)'] <- 'Before day 6'


mass_iqr_both_thermo_adj_lmer_est['term'][mass_iqr_both_thermo_adj_lmer_est['term'] ==
                                            'scale(thermo_aft_iqr_temp)'] <- 'After day 6'


## IQR temp before model: Label the estimates in data frame
mass_iqr_both_thermo_adj_lmer_est$model <- c('Both model')

## Combine regression estimates into a tidy table
mass_iqr_both_thermo_adj_lmer_tbl <- bind_rows(mass_iqr_before_thermo_adj_lmer_est,
                                               mass_iqr_after_thermo_adj_lmer_est,
                                               mass_iqr_both_thermo_adj_lmer_est)

## Re-code *nominal* factor (with ordered levels)  
# Set levels (odering) of 'model' variable
mass_iqr_both_thermo_adj_lmer_tbl <-
  transform(mass_iqr_both_thermo_adj_lmer_tbl,
            model = factor(model,
                           levels = c('Before model',
                                      'After model',
                                      'Both model')))

mass_iqr_both_thermo_adj_lmer_tbl <-
  transform(mass_iqr_both_thermo_adj_lmer_tbl,
            term = factor(term,
                          levels = c('Before day 6', 
                                     'After day 6')))

## Graph of estimates 
mass_iqr_both_thermo_adj_lmer_plot <-
  ggplot(mass_iqr_both_thermo_adj_lmer_tbl, aes(x = term, y = estimate,
                                                color = model)) +
  theme_classic() + 
  geom_hline(yintercept = 0, color = 'red',
             linetype = 2) + # line at null behind coefs
  geom_point(size = 6,
             position=position_dodge(width = 0.0)) +
  geom_errorbar(aes(ymin=(conf.low),
                    ymax=(conf.high)), width=.1,
                position=position_dodge(0.0), size = 1) +
  scale_color_manual(values=c('#440154', '#7ad151','#20A387FF')) +
  facet_grid(cols = vars(model), scales = 'free', space = 'free') +
  theme(strip.background =element_rect(fill= 'white'))+
  theme(strip.text = element_text(colour = 'black')) +
  #coord_flip() + # flip x and y axes
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18)) +
  theme(text = element_text(size=18)) +
  # bold and size title and axes labels
  theme(legend.position = 'none') +
  theme(axis.ticks = element_blank()) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size = 16, angle=45,
                                   margin = margin(t = 50, r = 0,
                                                   b = 0, l = 0), face = 'bold'), 
        axis.text.y = element_text(size = 16, angle=0,
                                   margin = margin(t = 25, r = 0,
                                                   b = 0, l = 10)),
        axis.title = element_text(size = 18), legend.title=element_blank(),
        legend.text=element_text(size=18),
        legend.position = 'none', #c(0.91, 0.94),
        legend.key = element_blank()) +
  labs(colour = "Parent feeding level") +
  xlab(expression(italic('Standardized IQR of temperature (SD)'))) +
  ylab(expression
       (atop(bold('Beta estimate and 95% CI'),
             paste(italic('Nestling mass (g)')))))
# remove axis ticks

print(mass_iqr_both_thermo_adj_lmer_plot)


##################### Predicted results graphs

## Make plots of temp effects before and after thermoreg
## Both thermo min temp
pred_bef <- ggpredict(mass_min_before_thermo_adj_lmer, c("thermo_bef_min_temp"))

plot(pred_bef) 

pred_aft <- ggpredict(mass_min_after_thermo_adj_lmer, c("thermo_aft_min_temp"))

plot(pred_aft) 


sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = thermo_bef_min_temp) &
                !is.na(x = thermo_aft_min_temp))

colors <- c("Before day 6" = "#481567FF", "After day 6" = "#95D840FF")
lines <- c("Before day 6" = 4, "After day 6" = 1)

temp_thermo_both_min_predicted <- ggplot() + 
  theme_classic() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  geom_line(data= pred_bef, size = 2, aes(x = x, y = predicted, col = "Before day 6")) +
  # geom_ribbon(data = pred, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high,
  #                              fill = "Before day 6"), alpha = .1, col = NA) +
  geom_line(data= pred_aft, size = 2, aes(x = x, y = predicted, col = "After day 6")) +
  # geom_ribbon(data = pred_aft, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, 
  #                                  fill = "After day 6"),  alpha = .1, col = NA) +
  geom_point(data = sub, aes(x = thermo_bef_min_temp, y = mass_pre_obs, col = "Before day 6"), 
             size = 1.5, alpha = 0.5, show.legend = FALSE) +
  geom_point(data = sub, aes(x = thermo_aft_min_temp, y = mass_pre_obs, col = "After day 6"), 
             size = 1.5, alpha = 0.5, show.legend = FALSE) +
  labs(x = "Minimum temperature (C)", y = "Nestling mass (g)", 
       color = "Before vs. after day 6", fill = "Before vs. after day 6") +
  scale_color_manual(values = colors, breaks = c("Before day 6", "After day 6"),
                     aesthetics = c("color", "fill"))
# Print
print(temp_thermo_both_min_predicted)


## Both thermo max temp
pred_bef <- ggpredict(mass_max_before_thermo_adj_lmer, c("thermo_bef_max_temp"))

plot(pred_bef) 

pred_aft <- ggpredict(mass_max_after_thermo_adj_lmer, c("thermo_aft_max_temp"))

plot(pred_aft) 


sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = thermo_bef_max_temp) &
                !is.na(x = thermo_aft_max_temp))

temp_thermo_both_max_predicted <- ggplot() + 
  theme_classic() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  geom_line(data= pred_bef, size = 2, aes(x = x, y = predicted, col = "Before day 6")) +
  # geom_ribbon(data = pred, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high,
  #                              fill = "Before day 6"), alpha = .1, col = NA) +
  geom_line(data= pred_aft, size = 2, aes(x = x, y = predicted, col = "After day 6")) +
  # geom_ribbon(data = pred_aft, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, 
  #                                  fill = "After day 6"),  alpha = .1, col = NA) +
  geom_point(data = sub, aes(x = thermo_bef_max_temp, y = mass_pre_obs, col = "Before day 6"), 
             size = 1.5, alpha = 0.5, show.legend = FALSE) +
  geom_point(data = sub, aes(x = thermo_aft_max_temp, y = mass_pre_obs, col = "After day 6"), 
             size = 1.5, alpha = 0.5, show.legend = FALSE) +
  labs(x = "Maximum temperature (C)", y = "Nestling mass (g)", 
       color = "Before vs. after day 6", fill = "Before vs. after day 6") +
  scale_color_manual(values = colors, breaks = c("Before day 6", "After day 6"),
                     aesthetics = c("color", "fill"))

# Print
print(temp_thermo_both_max_predicted)



## Both thermo IQR temp
pred_bef <- ggpredict(mass_iqr_before_thermo_adj_lmer, c("thermo_bef_iqr_temp"))

plot(pred_bef) 

pred_aft <- ggpredict(mass_iqr_after_thermo_adj_lmer, c("thermo_aft_iqr_temp"))

plot(pred_aft) 


sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = thermo_bef_iqr_temp) &
                !is.na(x = thermo_aft_iqr_temp))

temp_thermo_both_iqr_predicted <- ggplot() + 
  theme_classic() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  geom_line(data= pred_bef, size = 2, aes(x = x, y = predicted, col = "Before day 6")) +
  # geom_ribbon(data = pred, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high,
  #                              fill = "Before day 6"), alpha = .1, col = NA) +
  geom_line(data= pred_aft, size = 2, aes(x = x, y = predicted, col = "After day 6")) +
  # geom_ribbon(data = pred_aft, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, 
  #                                  fill = "After day 6"),  alpha = .1, col = NA) +
  geom_point(data = sub, aes(x = thermo_bef_iqr_temp, y = mass_pre_obs, col = "Before day 6"), 
             size = 1.5, alpha = 0.5, show.legend = FALSE) +
  geom_point(data = sub, aes(x = thermo_aft_iqr_temp, y = mass_pre_obs, col = "After day 6"), 
             size = 1.5, alpha = 0.5, show.legend = FALSE) +
  labs(x = "IQR of temperature (C)", y = "Nestling mass (g)", 
       color = "Before vs. after day 6", fill = "Before vs. after day 6") +
  scale_color_manual(values = colors, breaks = c("Before day 6", "After day 6"),
                     aesthetics = c("color", "fill"))

# Print
print(temp_thermo_both_iqr_predicted)


## Combined plot with points colored
combined_thermo_both_dw <- 
  ggarrange(mass_min_both_thermo_adj_lmer_plot, mass_max_both_thermo_adj_lmer_plot,
            mass_iqr_both_thermo_adj_lmer_plot,
            ncol = 3, nrow = 1, common.legend = TRUE,
            legend = "bottom", labels = c("a", "b", "c"), font.label = list(size = 18))

ggsave('combined_thermo_both_dw.png', plot = combined_thermo_both_dw, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 19, 
       height = 8, 
       units = c('in'), dpi = 300, limitsize = TRUE)

## Combined plot with points colored
combined_thermo_both_predictions <- 
  ggarrange(temp_thermo_both_min_predicted, 
            temp_thermo_both_max_predicted, temp_thermo_both_iqr_predicted,
            ncol = 3, nrow = 1, common.legend = TRUE,
            legend = "bottom", labels = c("d", "e", "f"), font.label = list(size = 18))

ggsave('combined_thermo_both_predictions.png', plot = combined_thermo_both_predictions, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 19, 
       height = 12, 
       units = c('in'), dpi = 300, limitsize = TRUE)

## Combined plot with both types of plot
combined_thermo_both_types <- 
  ggarrange(combined_thermo_both_dw, 
            combined_thermo_both_predictions,
            ncol = 1, nrow = 2, common.legend = FALSE,
            legend = "bottom")

# Print
print(combined_thermo_both_types)

# Save
ggsave('combined_thermo_both_types.png', plot = combined_thermo_both_types, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 19, 
       height = 12, 
       units = c('in'), dpi = 300, limitsize = TRUE)

