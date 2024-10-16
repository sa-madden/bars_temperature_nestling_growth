####### Purpose: calculate parental care index and format final
####### dataset for barn swallow nest microclimate and nestling 
####### growth dataset at nestling level
####### By: Sage Madden
####### Created: 12/16/2022
####### Last modified: 10/11/2024

# Code Blocks
# 1: Configure work space
# 2: Load data
# 3: Estimate BLUPs for parent feeding rate
# 4: Join BLUPs to full dataset
# 5: Save data

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

## Graph Plotting and Visualization Packages
library ('ggplot2')
library('hrbrthemes')
library('viridis')
library ('gridExtra')
library('broom')

## Modelling Packages
library('lme4')
library('DHARMa')

## boot packages
library ('boot')


### 1.3 Get Version and Session Info
R.Version()
sessionInfo()


###############################################################################
##############                        Load Data                  ##############
###############################################################################  

### Load Data
## Load RData tidy barn swallow data
load('Data/Tidy/tidy_parent_nestl_weather_data_10-4.RData')

###############################################################################
##############       Estimate BLUPs for parent feeding rate     ##############
############################################################################### 

### Feeding BLUPs
# NOTE: Use when there are repeated measuresments for a variable
# that is to be used as an explanatory variable in another analysis.
# Can control for other variables that bias estimates of explanatory
# variable

# NOTE: BLUPs are conditional modes from a generalized linear model
# (according to Doug Bates). 
# BLUP = fixef(intrcpt) + ranef
# Model  using glmmTMB (since these were zero inflated behavior data)
colnames(prim_merged)
feeding_blups_lmm <- lmer(total_feeding_visits ~ scale(nestling_age) + 
                             scale(nestling_number) + scale(obs_med_temp) +
                            scale(disturb_min) + 
                             offset(obs_duration/3600) +
                             (1|fnest_id) + (1|fsite),
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
                            (1|fnest_id) + (1|fsite),
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
                                    (1|fnest_id) + (1|fsite),
                                 data = subset(prim_merged,
                                               is.na(total_feeding_visits) == F &
                                                 is.na(nestling_age) == F &
                                                 is.na(nestling_number) == F &
                                                 is.na(obs_med_temp) == F &
                                                 is.na(disturb_min) == F))

simulation_output <- simulateResiduals(fittedModel = feeding_blups_glmm_nb, 
                                       plot = T)
# Looks pretty good




## b) Parameter estimates
summary(feeding_blups_glmm_nb) # model parameter estimates
# confint(feeding_blups_glmm_nb)  # 95% CIs
# confint gives an error


## Bootstrap parameter estimates
# bootstrapping number of resampling simulations
boot_feeding_blups_glmm_nb <- bootMer(x = feeding_blups_glmm_nb,
                                      FUN = fixef, nsim = 2000,
                                      seed = 632760,
                                      use.u = F, type = 'parametric')
tidy(boot_feeding_blups_glmm_nb) # beta estimates and SE
# use 'boot' package to generate 95% CI
bt_ci_feeding_blups_glmm_nb <- boot.ci(boot_feeding_blups_glmm_nb,
                            type = c('perc', 'norm', 'basic'),
                            index = 2) # CI for 1st betas
print(bt_ci_feeding_blups_glmm_nb)

# use 'boot' package to generate 95% CI
bt_ci_feeding_blups_glmm_nb_2 <- boot.ci(boot_feeding_blups_glmm_nb,
                                       type = c('perc', 'norm', 'basic'),
                                       index = 3) # CI for 2nd betas
print(bt_ci_feeding_blups_glmm_nb_2)

# use 'boot' package to generate 95% CI
bt_ci_feeding_blups_glmm_nb_3 <- boot.ci(boot_feeding_blups_glmm_nb,
                                         type = c('perc', 'norm', 'basic'),
                                         index = 4) # CI for 3rd betas
print(bt_ci_feeding_blups_glmm_nb_3)


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


###############################################################################
##############            Join BLUPs to full dataset             ##############
############################################################################### 

## Left join feeding BLUPs to late_nestling_parent_care
late_nestling_parent_care <- late_nestling_parent_care %>%
  left_join(feeding_blups, by = c('nest_id' = 'fnest_id'), 
            copy = F)


###############################################################################
##############                    Save data                      ##############
############################################################################### 

save(file = 'Data/Tidy/tidy_parent_nestl_weather_data_10-4_with_BLUPs.RData', 
     list = c('prim_merged', 'nestl_merged', 'govee_daily', 'late_nestling_parent_care'))



