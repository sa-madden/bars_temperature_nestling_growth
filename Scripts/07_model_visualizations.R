####### Purpose: visualize models for barn swallow
####### nest mircoclimate and nestling growth dataset
####### By: Sage Madden
####### Created: 12/20/2022
####### Last modified: 4/19/23

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
library('ggpubr')
library('ggeffects')
library('ggnewscale')

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
##############      Dot and whisker plots of coefs         ##############
###############################################################################
### Parental care * temperature
### FEEDING 

## Minimum temp 
# Remove terms I don't need 
mass_min_temp_blups_lmer_tidy <-
  broom::tidy(mass_min_temp_blups_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "With interaction")

mass_min_temp_blups_low_lmer_tidy <-
  broom::tidy(mass_min_temp_blups_low_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Stratified: LOW")

mass_min_temp_blups_med_lmer_tidy <-
  broom::tidy(mass_min_temp_blups_med_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Stratified: AVG")

mass_min_temp_blups_high_lmer_tidy <-
  broom::tidy(mass_min_temp_blups_high_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Stratified: HIGH")

four_models <- rbind(mass_min_temp_blups_lmer_tidy, mass_min_temp_blups_low_lmer_tidy,
                      mass_min_temp_blups_med_lmer_tidy, mass_min_temp_blups_high_lmer_tidy)

# Create plot 
mass_min_temp_dw <- dwplot(four_models)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ Min_temp * Feeding_blups")

## Minimum temp -- adjusted
# Remove terms I don't need 
mass_min_temp_blups_adj_lmer_tidy <-
  broom::tidy(mass_min_temp_blups_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "With interaction")

mass_min_temp_blups_low_adj_lmer_tidy <-
  broom::tidy(mass_min_temp_blups_low_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Stratified: LOW")

mass_min_temp_blups_med_adj_lmer_tidy <-
  broom::tidy(mass_min_temp_blups_med_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Stratified: AVG")

mass_min_temp_blups_high_adj_lmer_tidy <-
  broom::tidy(mass_min_temp_blups_high_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Stratified: HIGH")

four_models <- rbind(mass_min_temp_blups_adj_lmer_tidy, mass_min_temp_blups_low_adj_lmer_tidy,
                     mass_min_temp_blups_med_adj_lmer_tidy, mass_min_temp_blups_high_adj_lmer_tidy)

# Create plot 
mass_min_temp_adj_dw <- dwplot(four_models)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ Min_temp * Feeding_blups + Nestling_number + Hatch_date")


# Print plot 
print(mass_min_temp_adj_dw)

# Save plot
ggsave('mass_min_temp_adj_dw.png', plot = mass_min_temp_adj_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

## Maximum temp 
# Remove terms I don't need 
mass_max_temp_blups_lmer_tidy <-
  broom::tidy(mass_max_temp_blups_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "With interaction")

mass_max_temp_blups_noint_lmer_tidy <-
  broom::tidy(mass_max_temp_blups_noint_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Without interaction")

two_models <- rbind(mass_max_temp_blups_lmer_tidy, mass_max_temp_blups_noint_lmer_tidy)

# Create plot 
mass_max_temp_dw <- dwplot(two_models)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ Max_temp * Feeding_blups")


# Print plot 
print(mass_max_temp_dw)

# Save plot
ggsave('mass_max_temp_dw.png', plot = mass_max_temp_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 


## Maximum temp -- adjusted
# Remove terms I don't need 
mass_max_temp_blups_adj_lmer_tidy <-
  broom::tidy(mass_max_temp_blups_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "With interaction")

mass_max_temp_blups_noint_adj_lmer_tidy <-
  broom::tidy(mass_max_temp_blups_noint_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Without interaction")

two_models <- rbind(mass_max_temp_blups_adj_lmer_tidy, mass_max_temp_blups_noint_adj_lmer_tidy)

# Create plot 
mass_max_temp_adj_dw <- dwplot(two_models)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ Max_temp * Feeding_blups + Nestling_number + Hatch_date")


# Print plot 
print(mass_max_temp_adj_dw)

# Save plot
ggsave('mass_max_temp_adj_dw.png', plot = mass_max_temp_adj_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

## IQR of temp 
# Remove terms I don't need 
mass_iqr_temp_blups_lmer_tidy <-
  broom::tidy(mass_iqr_temp_blups_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "With interaction")

mass_iqr_temp_blups_noint_lmer_tidy <-
  broom::tidy(mass_iqr_temp_blups_noint_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Without interaction")

two_models <- rbind(mass_iqr_temp_blups_lmer_tidy, mass_iqr_temp_blups_noint_lmer_tidy)

# Create plot 
mass_iqr_temp_dw <- dwplot(two_models)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ IQR_temp * Feeding_blups")


# Print plot 
print(mass_iqr_temp_dw)

# Save plot
ggsave('mass_iqr_temp_dw.png', plot = mass_iqr_temp_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

## IQR of temp -- adjusted 
# Remove terms I don't need 
mass_iqr_temp_blups_adj_lmer_tidy <-
  broom::tidy(mass_iqr_temp_blups_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "With interaction")

mass_iqr_temp_blups_noint_adj_lmer_tidy <-
  broom::tidy(mass_iqr_temp_blups_noint_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Without interaction")

two_models <- rbind(mass_iqr_temp_blups_adj_lmer_tidy, mass_iqr_temp_blups_noint_adj_lmer_tidy)

# Create plot 
mass_iqr_temp_adj_dw <- dwplot(two_models)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ IQR_temp * Feeding_blups + Nestling_number + Hatch_date")


# Print plot 
print(mass_iqr_temp_adj_dw)

# Save plot
ggsave('mass_iqr_temp_adj_dw.png', plot = mass_iqr_temp_adj_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)


# Create combined plots
combined_temp_blups_mass <- 
  ggarrange(mass_min_temp_dw, mass_max_temp_dw, 
          mass_iqr_temp_dw, 
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)

print(combined_temp_blups_mass)

ggsave('combined_temp_blups_mass.png', plot = combined_temp_blups_mass, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 10, 
       height = 10, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# Create combined plots - adjusted 
combined_temp_blups_mass_adj <- 
  ggarrange(mass_min_temp_adj_dw, mass_max_temp_adj_dw, 
            mass_iqr_temp_adj_dw, 
            labels = c("A", "B", "C"),
            ncol = 1, nrow = 3)

print(combined_temp_blups_mass_adj)

ggsave('combined_temp_blups_mass_adj.png', plot = combined_temp_blups_mass_adj, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 10, 
       height = 10, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

## Relative nestling size 
# Min temp
# Remove terms I don't need 
mass_min_temp_size_lmer_tidy <-
  broom::tidy(mass_min_temp_size_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "With interaction")

# Create plot 
mass_min_temp_size_dw <- dwplot(mass_min_temp_size_lmer_tidy)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ Min_temp * Size_order")


# Print plot 
print(mass_min_temp_size_dw)

# Save plot
ggsave('mass_min_temp_size_dw.png', plot = mass_min_temp_size_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# Min temp - adjusted
# Remove terms I don't need 
mass_min_temp_size_adj_lmer_tidy <-
  broom::tidy(mass_min_temp_size_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "With interaction")

# Create plot 
mass_min_temp_size_adj_dw <- dwplot(mass_min_temp_size_adj_lmer_tidy)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ Min_temp * Size_order + Nestling_number + Hatch_date")


# Print plot 
print(mass_min_temp_size_adj_dw)

# Save plot
ggsave('mass_min_temp_size_adj_dw.png', plot = mass_min_temp_size_adj_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 


## Max temp
# Remove terms I don't need 
mass_max_temp_size_lmer_tidy <-
  broom::tidy(mass_max_temp_size_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "With interaction")

mass_max_temp_size_small_lmer_tidy <-
  broom::tidy(mass_max_temp_size_small_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Stratified: SMALL")

mass_max_temp_size_big_lmer_tidy <-
  broom::tidy(mass_max_temp_size_big_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Stratified: OTHER")

three_models <- rbind(mass_max_temp_size_lmer_tidy, mass_max_temp_size_small_lmer_tidy,
                      mass_max_temp_size_big_lmer_tidy)

# Create plot 
mass_max_temp_size_dw <- dwplot(three_models)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ Max_temp * Size_order")


# Print plot 
print(mass_max_temp_size_dw)

# Save plot
ggsave('mass_max_temp_size_dw.png', plot = mass_max_temp_size_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

## Max temp - adjusted
# Remove terms I don't need 
mass_max_temp_size_adj_lmer_tidy <-
  broom::tidy(mass_max_temp_size_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "With interaction")

mass_max_temp_size_small_adj_lmer_tidy <-
  broom::tidy(mass_max_temp_size_small_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Stratified: SMALL")

mass_max_temp_size_big_adj_lmer_tidy <-
  broom::tidy(mass_max_temp_size_big_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Stratified: OTHER")

three_models <- rbind(mass_max_temp_size_adj_lmer_tidy, mass_max_temp_size_small_adj_lmer_tidy,
                      mass_max_temp_size_big_adj_lmer_tidy)

# Create plot 
mass_max_temp_size_adj_dw <- dwplot(three_models)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ Max_temp * Size_order + Nestling_number + Hatch_date")


# Print plot 
print(mass_max_temp_size_adj_dw)

# Save plot
ggsave('mass_max_temp_size_adj_dw.png', plot = mass_max_temp_size_adj_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

## IQR 
# Remove terms I don't need 
mass_iqr_temp_size_lmer_tidy <-
  broom::tidy(mass_iqr_temp_size_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "With interaction")

mass_iqr_temp_size_small_lmer_tidy <-
  broom::tidy(mass_iqr_temp_size_small_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Stratified: SMALL")

mass_iqr_temp_size_big_lmer_tidy <-
  broom::tidy(mass_iqr_temp_size_big_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Stratified: OTHER")

three_models <- rbind(mass_iqr_temp_size_lmer_tidy, mass_iqr_temp_size_small_lmer_tidy,
                      mass_iqr_temp_size_big_lmer_tidy)

# Create plot 
mass_iqr_temp_size_dw <- dwplot(three_models)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ IQR_temp * Size_order")


# Print plot 
print(mass_iqr_temp_size_dw)

# Save plot
ggsave('mass_iqr_temp_size_dw.png', plot = mass_iqr_temp_size_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# IQR -adjusted
# Remove terms I don't need 
mass_iqr_temp_size_adj_lmer_tidy <-
  broom::tidy(mass_iqr_temp_size_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "With interaction")

mass_iqr_temp_size_small_adj_lmer_tidy <-
  broom::tidy(mass_iqr_temp_size_small_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Stratified: SMALL")

mass_iqr_temp_size_big_adj_lmer_tidy <-
  broom::tidy(mass_iqr_temp_size_big_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Stratified: OTHER")

three_models <- rbind(mass_iqr_temp_size_adj_lmer_tidy, mass_iqr_temp_size_small_adj_lmer_tidy,
                      mass_iqr_temp_size_big_adj_lmer_tidy)

# Create plot 
mass_iqr_temp_size_adj_dw <- dwplot(three_models)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ IQR_temp * Size_order + Nestling_number + Hatch_date")


# Print plot 
print(mass_iqr_temp_size_adj_dw)

# Save plot
ggsave('mass_iqr_temp_size_adj_dw.png', plot = mass_iqr_temp_size_adj_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# Combined plots
combined_temp_size_mass_adj <- 
  ggarrange(mass_min_temp_size_adj_dw, mass_max_temp_size_adj_dw, 
            mass_iqr_temp_size_adj_dw, 
            labels = c("A", "B", "C"),
            ncol = 1, nrow = 3)

print(combined_temp_size_mass_adj)

ggsave('combined_temp_size_mass_adj.png', plot = combined_temp_size_mass_adj, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 10, 
       height = 10, 
       units = c('in'), dpi = 300, limitsize = TRUE) 


## Before/after thermo
# Min temp
# Before thermo
# Remove terms I don't need 
mass_min_before_thermo_lmer_tidy <-
  broom::tidy(mass_min_before_thermo_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Before thermoreg. Indp.")

# Create plot 
mass_min_before_thermo_dw <- dwplot(mass_min_before_thermo_lmer_tidy)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ Min_temp_before_thermo")


# Print plot 
print(mass_min_before_thermo_dw)

# Save plot
ggsave('mass_min_before_thermo_dw.png', plot = mass_min_before_thermo_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# After thermo
# Remove terms I don't need 
mass_min_after_thermo_lmer_tidy <-
  broom::tidy(mass_min_after_thermo_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "After thermoreg. Indp.")

# Create plot 
mass_min_after_thermo_dw <- dwplot(mass_min_after_thermo_lmer_tidy)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ Min_temp_after_thermo")


# Print plot 
print(mass_min_after_thermo_dw)

# Save plot
ggsave('mass_min_after_thermo_dw.png', plot = mass_min_after_thermo_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# Both 
# Remove terms I don't need 
mass_min_both_thermo_lmer_tidy <-
  broom::tidy(mass_min_both_thermo_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Both stages")

# Create plot 
mass_min_both_thermo_dw <- dwplot(mass_min_both_thermo_lmer_tidy)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ Min_temp_before_thermo + Min_temp_after_thermo")


# Print plot 
print(mass_min_both_thermo_dw)

# Save plot
ggsave('mass_min_both_thermo_dw.png', plot = mass_min_both_thermo_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# Combined plot
combined_min_temp_bef_aft_thermo <- 
  ggarrange(mass_min_before_thermo_dw, mass_min_after_thermo_dw,
            mass_min_both_thermo_dw,
            labels = c("A", "B", "C"),
            ncol = 1, nrow = 3)

print(combined_min_temp_bef_aft_thermo)

ggsave('combined_min_temp_bef_aft_thermo.png', plot = combined_min_temp_bef_aft_thermo, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 10, 
       height = 10, 
       units = c('in'), dpi = 300, limitsize = TRUE) 


# Min temp adjusted
# Before thermo
# Remove terms I don't need 
mass_min_before_thermo_adj_lmer_tidy <-
  broom::tidy(mass_min_before_thermo_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Before thermoreg. Indp.")

# Create plot 
mass_min_before_thermo_dw <- dwplot(mass_min_before_thermo_adj_lmer_tidy)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ Min_temp_before_thermo")


# Print plot 
print(mass_min_before_thermo_dw)

# Save plot
ggsave('mass_min_before_thermo_adj_dw.png', plot = mass_min_before_thermo_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# After thermo
# Remove terms I don't need 
mass_min_after_thermo_adj_lmer_tidy <-
  broom::tidy(mass_min_after_thermo_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "After thermoreg. Indp.")

# Create plot 
mass_min_after_thermo_dw <- dwplot(mass_min_after_thermo_adj_lmer_tidy)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ Min_temp_after_thermo")


# Print plot 
print(mass_min_after_thermo_dw)

# Save plot
ggsave('mass_min_after_thermo_adj_dw.png', plot = mass_min_after_thermo_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# Both 
# Remove terms I don't need 
mass_min_both_thermo_adj_lmer_tidy <-
  broom::tidy(mass_min_both_thermo_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Both stages")

# Create plot 
mass_min_both_thermo_dw <- dwplot(mass_min_both_thermo_adj_lmer_tidy)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ Min_temp_before_thermo + Min_temp_after_thermo")


# Print plot 
print(mass_min_both_thermo_dw)

# Save plot
ggsave('mass_min_both_thermo_adj_dw.png', plot = mass_min_both_thermo_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# Combined plot
combined_min_temp_bef_aft_thermo <- 
  ggarrange(mass_min_before_thermo_dw, mass_min_after_thermo_dw,
            mass_min_both_thermo_dw,
            labels = c("A", "B", "C"),
            ncol = 1, nrow = 3)

print(combined_min_temp_bef_aft_thermo)

ggsave('combined_min_temp_bef_aft_thermo_adj.png', plot = combined_min_temp_bef_aft_thermo, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 10, 
       height = 10, 
       units = c('in'), dpi = 300, limitsize = TRUE) 



# Max temp
# Before thermo
# Remove terms I don't need 
mass_max_before_thermo_lmer_tidy <-
  broom::tidy(mass_max_before_thermo_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Before thermoreg. Indp.")

# Create plot 
mass_max_before_thermo_dw <- dwplot(mass_max_before_thermo_lmer_tidy)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ Max_temp_before_thermo")


# Print plot 
print(mass_max_before_thermo_dw)

# Save plot
ggsave('mass_max_before_thermo_dw.png', plot = mass_max_before_thermo_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# After thermo
# Remove terms I don't need 
mass_max_after_thermo_lmer_tidy <-
  broom::tidy(mass_max_after_thermo_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "After thermoreg. Indp.")

# Create plot 
mass_max_after_thermo_dw <- dwplot(mass_max_after_thermo_lmer_tidy)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ Max_temp_after_thermo")


# Print plot 
print(mass_max_after_thermo_dw)

# Save plot
ggsave('mass_max_after_thermo_dw.png', plot = mass_max_after_thermo_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# Both
# Remove terms I don't need 
mass_max_both_thermo_lmer_tidy <-
  broom::tidy(mass_max_both_thermo_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Both stages")

# Create plot 
mass_max_both_thermo_dw <- dwplot(mass_max_both_thermo_lmer_tidy)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ Max_temp_before_thermo + Max_temp_after_thermo")


# Print plot 
print(mass_max_both_thermo_dw)

# Save plot
ggsave('mass_max_both_thermo_dw.png', plot = mass_max_both_thermo_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# Combined plot
combined_max_temp_bef_aft_thermo <- 
  ggarrange(mass_max_before_thermo_dw, mass_max_after_thermo_dw,
            mass_max_both_thermo_dw,
            labels = c("A", "B", "C"),
            ncol = 1, nrow = 3)

print(combined_max_temp_bef_aft_thermo)

ggsave('combined_max_temp_bef_aft_thermo.png', plot = combined_max_temp_bef_aft_thermo, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 10, 
       height = 10, 
       units = c('in'), dpi = 300, limitsize = TRUE) 


# Max temp
# Before thermo
# Remove terms I don't need 
mass_max_before_thermo_adj_lmer_tidy <-
  broom::tidy(mass_max_before_thermo_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Before thermoreg. Indp.")

# Create plot 
mass_max_before_thermo_dw <- dwplot(mass_max_before_thermo_adj_lmer_tidy)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ Max_temp_before_thermo")


# Print plot 
print(mass_max_before_thermo_dw)

# Save plot
ggsave('mass_max_before_thermo_adj_dw.png', plot = mass_max_before_thermo_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# After thermo adjusted
# Remove terms I don't need 
mass_max_after_thermo_adj_lmer_tidy <-
  broom::tidy(mass_max_after_thermo_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "After thermoreg. Indp.")

# Create plot 
mass_max_after_thermo_dw <- dwplot(mass_max_after_thermo_adj_lmer_tidy)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ Max_temp_after_thermo")


# Print plot 
print(mass_max_after_thermo_dw)

# Save plot
ggsave('mass_max_after_thermo_adj_dw.png', plot = mass_max_after_thermo_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# Both
# Remove terms I don't need 
mass_max_both_thermo_adj_lmer_tidy <-
  broom::tidy(mass_max_both_thermo_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Both stages")

# Create plot 
mass_max_both_thermo_dw <- dwplot(mass_max_both_thermo_adj_lmer_tidy)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ Max_temp_before_thermo + Max_temp_after_thermo")


# Print plot 
print(mass_max_both_thermo_dw)

# Save plot
ggsave('mass_max_both_thermo_adj_dw.png', plot = mass_max_both_thermo_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# Combined plot
combined_max_temp_bef_aft_thermo <- 
  ggarrange(mass_max_before_thermo_dw, mass_max_after_thermo_dw,
            mass_max_both_thermo_dw,
            labels = c("A", "B", "C"),
            ncol = 1, nrow = 3)

print(combined_max_temp_bef_aft_thermo)

ggsave('combined_max_temp_bef_aft_thermo_adj.png', plot = combined_max_temp_bef_aft_thermo, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 10, 
       height = 10, 
       units = c('in'), dpi = 300, limitsize = TRUE) 


## IQR of temp
# Before thermo
# Remove terms I don't need 
mass_iqr_before_thermo_lmer_tidy <-
  broom::tidy(mass_iqr_before_thermo_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Before thermoreg. Indp.")

# Create plot 
mass_iqr_before_thermo_dw <- dwplot(mass_iqr_before_thermo_lmer_tidy)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ IQR_temp_before_thermo")


# Print plot 
print(mass_iqr_before_thermo_dw)

# Save plot
ggsave('mass_iqr_before_thermo_dw.png', plot = mass_iqr_before_thermo_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# After thermo
# Remove terms I don't need 
mass_iqr_after_thermo_lmer_tidy <-
  broom::tidy(mass_iqr_after_thermo_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "After thermoreg. Indp.")

# Create plot 
mass_iqr_after_thermo_dw <- dwplot(mass_iqr_after_thermo_lmer_tidy)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ IQR_temp_after_thermo")


# Print plot 
print(mass_iqr_after_thermo_dw)

# Save plot
ggsave('mass_iqr_after_thermo_dw.png', plot = mass_iqr_after_thermo_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# Both
# Remove terms I don't need 
mass_iqr_both_thermo_lmer_tidy <-
  broom::tidy(mass_iqr_both_thermo_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Both stages")

# Create plot 
mass_iqr_both_thermo_dw <- dwplot(mass_iqr_both_thermo_lmer_tidy)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ IQR_temp_before_thermo + IQR_temp_after_thermo")


# Print plot 
print(mass_iqr_both_thermo_dw)

# Save plot
ggsave('mass_iqr_both_thermo_dw.png', plot = mass_iqr_both_thermo_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# Combined plot
combined_iqr_temp_bef_aft_thermo <- 
  ggarrange(mass_iqr_before_thermo_dw, mass_iqr_after_thermo_dw,
            mass_iqr_both_thermo_dw,
            labels = c("A", "B", "C"),
            ncol = 1, nrow = 3)

print(combined_iqr_temp_bef_aft_thermo)

ggsave('combined_iqr_temp_bef_aft_thermo.png', plot = combined_iqr_temp_bef_aft_thermo, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 10, 
       height = 10, 
       units = c('in'), dpi = 300, limitsize = TRUE) 



## IQR of temp adjiusted
# Before thermo
# Remove terms I don't need 
mass_iqr_before_thermo_adj_lmer_tidy <-
  broom::tidy(mass_iqr_before_thermo_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Before thermoreg. Indp.")

# Create plot 
mass_iqr_before_thermo_dw <- dwplot(mass_iqr_before_thermo_adj_lmer_tidy)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ IQR_temp_before_thermo")


# Print plot 
print(mass_iqr_before_thermo_dw)

# Save plot
ggsave('mass_iqr_before_thermo_adj_dw.png', plot = mass_iqr_before_thermo_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# After thermo
# Remove terms I don't need 
mass_iqr_after_thermo_adj_lmer_tidy <-
  broom::tidy(mass_iqr_after_thermo_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "After thermoreg. Indp.")

# Create plot 
mass_iqr_after_thermo_dw <- dwplot(mass_iqr_after_thermo_adj_lmer_tidy)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ IQR_temp_after_thermo")


# Print plot 
print(mass_iqr_after_thermo_dw)

# Save plot
ggsave('mass_iqr_after_thermo_adj_dw.png', plot = mass_iqr_after_thermo_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# Both
# Remove terms I don't need 
mass_iqr_both_thermo_adj_lmer_tidy <-
  broom::tidy(mass_iqr_both_thermo_adj_lmer) %>% 
  filter(effect != "ran_pars") %>% mutate(model = "Both stages")

# Create plot 
mass_iqr_both_thermo_dw <- dwplot(mass_iqr_both_thermo_adj_lmer_tidy)+
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed")+
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(title = "Mass ~ IQR_temp_before_thermo + IQR_temp_after_thermo")


# Print plot 
print(mass_iqr_both_thermo_dw)

# Save plot
ggsave('mass_iqr_both_thermo_adj_dw.png', plot = mass_iqr_both_thermo_dw, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# Combined plot
combined_iqr_temp_bef_aft_thermo <- 
  ggarrange(mass_iqr_before_thermo_dw, mass_iqr_after_thermo_dw,
            mass_iqr_both_thermo_dw,
            labels = c("A", "B", "C"),
            ncol = 1, nrow = 3)

print(combined_iqr_temp_bef_aft_thermo)

ggsave('combined_iqr_temp_bef_aft_thermo_adj.png', plot = combined_iqr_temp_bef_aft_thermo, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 10, 
       height = 10, 
       units = c('in'), dpi = 300, limitsize = TRUE) 


################################################################


## Model predictions over raw data
## Minimum temp + parental care adjusted
pred <- ggpredict(mass_min_temp_blups_adj_lmer, c("nest_min_temp", "feeding_expontd_blups"))

plot(pred) 

sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_min_temp) &
                !is.na(x = feeding_expontd_blups))

temp_min_feeding_blups_predicted <- ggplot(pred, aes(x = x, y = predicted)) + 
  geom_line(size = 1.5, aes(col = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, col = group, fill = group), alpha = .1,
              col = NA) +
  geom_point(data = sub, aes(x = nest_min_temp, y = mass_pre_obs), size = 1.5) +
  theme_classic() + 
  scale_color_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                      labels = c("-1 SD (13.48 visits/hr)", "Mean (14.17 visits/hr)",
                                 "+1 SD (14.87 visits/hr)")) +
  scale_fill_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                     labels = c("-1 SD (13.48 visits/hr)", "Mean (14.17 visits/hr)",
                                "+1 SD (14.87 visits/hr)")) +
  labs(x = "Minimum temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       colour = "Parent feeding level (visits/hr)", fill = "Parent feeding level (visits/hr)") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18))

# Print
print(temp_min_feeding_blups_predicted)

# Save
ggsave('temp_min_feeding_blups_predicted.png', plot = temp_min_feeding_blups_predicted, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 8, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# With points colored 
temp_min_feeding_blups_col_predicted <- ggplot(pred, aes(x = x, y = predicted)) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  geom_line(size = 1.5, aes(col = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, col = group, fill = group), alpha = .1,
              col = NA) +
  scale_color_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                      labels = c("-1 SD (13.48 visits/hr)", "Mean (14.17 visits/hr)",
                                 "+1 SD (14.87 visits/hr)")) +
  scale_fill_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                     labels = c("-1 SD (13.48 visits/hr)", "Mean (14.17 visits/hr)",
                                "+1 SD (14.87 visits/hr)")) +
  labs(x = "Minimum temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       color = "Parent feeding level (visits/hr)", fill = "Parent feeding level (visits/hr)") +
  new_scale("color") + 
  geom_point(data = sub, aes(x = nest_min_temp, y = mass_pre_obs, col = feeding_expontd_blups), size = 1.5) +
  scale_color_viridis(option = "viridis", end = 0.8) +
  labs(color = "Parent feeding level (visits/hr)") 

## Maximum temp and parental care adjusted
pred <- ggpredict(mass_max_temp_blups_adj_lmer, c("nest_max_temp", "feeding_expontd_blups"))

plot(pred) 

sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_max_temp) &
                !is.na(x = feeding_expontd_blups))

temp_max_feeding_blups_predicted <- ggplot(pred, aes(x = x, y = predicted)) + 
  geom_line(size = 1.5, aes(col = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, col = group, fill = group), alpha = .1,
              col = NA) +
  geom_point(data = sub, aes(x = nest_max_temp, y = mass_pre_obs), size = 1.5) +
  theme_classic() + 
  scale_color_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                      labels = c("-1 SD (13.48 visits/hr)", "Mean (14.17 visits/hr)",
                                 "+1 SD (14.87 visits/hr)")) +
  scale_fill_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                     labels = c("-1 SD (13.48 visits/hr)", "Mean (14.17 visits/hr)",
                                "+1 SD (14.87 visits/hr)")) +
  labs(x = "Maximum temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       colour = "Parent feeding level (visits/hr)", fill = "Parent feeding level (visits/hr)") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18))

# Print
print(temp_max_feeding_blups_predicted)

# Save
ggsave('temp_max_feeding_blups_predicted.png', plot = temp_max_feeding_blups_predicted, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 8, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# With points colored 
temp_max_feeding_blups_col_predicted <- ggplot(pred, aes(x = x, y = predicted)) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  geom_line(size = 1.5, aes(col = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, col = group, fill = group), alpha = .1,
              col = NA) +
  scale_color_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                      labels = c("-1 SD (13.48 visits/hr)", "Mean (14.17 visits/hr)",
                                 "+1 SD (14.87 visits/hr)")) +
  scale_fill_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                     labels = c("-1 SD (13.48 visits/hr)", "Mean (14.17 visits/hr)",
                                "+1 SD (14.87 visits/hr)")) +
  labs(x = "Maximum temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       color = "Parent feeding level (visits/hr)", fill = "Parent feeding level (visits/hr)") +
  new_scale("color") + 
  geom_point(data = sub, aes(x = nest_max_temp, y = mass_pre_obs, col = feeding_expontd_blups), size = 1.5) +
  scale_color_viridis(option = "viridis", end = 0.8) +
  labs(color = "Parent feeding level (visits/hr)") 

## IQR temp and parental care adjusted
pred <- ggpredict(mass_iqr_temp_blups_adj_lmer, c("nest_iqr_temp", "feeding_expontd_blups"))

plot(pred) 

sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_iqr_temp) &
                !is.na(x = feeding_expontd_blups))

temp_iqr_feeding_blups_predicted <- ggplot(pred, aes(x = x, y = predicted)) + 
  geom_line(size = 1.5, aes(col = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, col = group, fill = group), alpha = .1,
              col = NA) +
  geom_point(data = sub, aes(x = nest_iqr_temp, y = mass_pre_obs), size = 1.5) +
  theme_classic() + 
  scale_color_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                      labels = c("-1 SD (13.48 visits/hr)", "Mean (14.17 visits/hr)",
                                 "+1 SD (14.87 visits/hr)")) +
  scale_fill_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                     labels = c("-1 SD (13.48 visits/hr)", "Mean (14.17 visits/hr)",
                                "+1 SD (14.87 visits/hr)")) +
  labs(x = "IQR of temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       colour = "Parent feeding level (visits/hr)", fill = "Parent feeding level (visits/hr)") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18))

# Print
print(temp_iqr_feeding_blups_predicted)

# Save
ggsave('temp_iqr_feeding_blups_predicted.png', plot = temp_iqr_feeding_blups_predicted, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 8, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# With points colored 
temp_iqr_feeding_blups_col_predicted <- ggplot(pred, aes(x = x, y = predicted)) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  geom_line(size = 1.5, aes(col = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, col = group, fill = group), alpha = .1,
              col = NA) +
  scale_color_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                      labels = c("-1 SD (13.48 visits/hr)", "Mean (14.17 visits/hr)",
                                 "+1 SD (14.87 visits/hr)")) +
  scale_fill_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                     labels = c("-1 SD (13.48 visits/hr)", "Mean (14.17 visits/hr)",
                                "+1 SD (14.87 visits/hr)")) +
  labs(x = "IQR of temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       color = "Parent feeding level (visits/hr)", fill = "Parent feeding level (visits/hr)") +
  new_scale("color") + 
  geom_point(data = sub, aes(x = nest_iqr_temp, y = mass_pre_obs, col = feeding_expontd_blups), size = 1.5) +
  scale_color_viridis(option = "viridis", end = 0.8) +
  labs(color = "Parent feeding level (visits/hr)") 


## Combined plot
combined_temp_parent_predictions <- 
  ggarrange(temp_min_feeding_blups_predicted, temp_max_feeding_blups_predicted,
            temp_iqr_feeding_blups_predicted,
            ncol = 3, nrow = 1, common.legend = TRUE)

# Print
print(combined_temp_parent_predictions)

# Save
ggsave('combined_temp_parent_predictions.png', plot = combined_temp_parent_predictions, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 18, 
       height = 8, 
       units = c('in'), dpi = 300, limitsize = TRUE)

## Combined plot with points colored
combined_temp_parent_col_predictions <- 
  ggarrange(temp_min_feeding_blups_col_predicted, temp_max_feeding_blups_col_predicted,
            temp_iqr_feeding_blups_col_predicted,
            ncol = 3, nrow = 1, common.legend = TRUE,
            legend = "right")

# Print
print(combined_temp_parent_col_predictions)

# Save
ggsave('combined_temp_parent_col_predictions.png', plot = combined_temp_parent_col_predictions, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 22, 
       height = 8, 
       units = c('in'), dpi = 300, limitsize = TRUE)



## Minimum temp + brooding adjusted
pred <- ggpredict(mass_min_temp_brooding_blups_adj_lmer, c("nest_min_temp", "brooding_blups"))

plot(pred) 

sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_min_temp) &
                !is.na(x = brooding_blups))

temp_min_brooding_blups_predicted <- ggplot(pred, aes(x = x, y = predicted)) + 
  geom_line(size = 1.5, aes(col = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, col = group, fill = group), alpha = .1,
              col = NA) +
  geom_point(data = sub, aes(x = nest_min_temp, y = mass_pre_obs), size = 1.5) +
  theme_classic() + 
  scale_color_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                      labels = c("-1 SD (22.32 minutes/hr)", "Mean (22.82 minutes/hr)",
                                 "+1 SD (23.32 minutes/hr)")) +
  scale_fill_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                     labels = c("-1 SD (22.32 minutes/hr)", "Mean (22.82 minutes/hr)",
                                "+1 SD (23.32 minutes/hr)")) +
  labs(x = "Minimum temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       colour = "Parent brooding level (minutes/hr)", fill = "Parent brooding level (minutes/hr)") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18))

# Print
print(temp_min_brooding_blups_predicted)

# Save
ggsave('temp_min_brooding_blups_predicted.png', plot = temp_min_brooding_blups_predicted, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 8, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# With points colored 
temp_min_brooding_blups_col_predicted <- ggplot(pred, aes(x = x, y = predicted)) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  geom_line(size = 1.5, aes(col = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, col = group, fill = group), alpha = .1,
              col = NA) +
  scale_color_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                      labels = c("-1 SD (22.32 minutes/hr)", "Mean (22.82 minutes/hr)",
                                 "+1 SD (23.32 minutes/hr)")) +
  scale_fill_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                     labels = c("-1 SD (22.32 minutes/hr)", "Mean (22.82 minutes/hr)",
                                "+1 SD (23.32 minutes/hr)")) +
  labs(x = "Minimum temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       color = "Parent brooding level (minutes/hr)", fill = "Parent brooding level (minutes/hr)") +
  new_scale("color") + 
  geom_point(data = sub, aes(x = nest_min_temp, y = mass_pre_obs, col = brooding_blups), size = 1.5) +
  scale_color_viridis(option = "viridis", end = 0.8) +
  labs(color = "Parent brooding level (minutes/hr)") 

## Maximum temp and parental care adjusted
pred <- ggpredict(mass_max_temp_brooding_blups_adj_lmer, c("nest_max_temp", "brooding_blups"))

plot(pred) 

sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_max_temp) &
                !is.na(x = brooding_blups))

temp_max_brooding_blups_predicted <- ggplot(pred, aes(x = x, y = predicted)) + 
  geom_line(size = 1.5, aes(col = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, col = group, fill = group), alpha = .1,
              col = NA) +
  geom_point(data = sub, aes(x = nest_max_temp, y = mass_pre_obs), size = 1.5) +
  theme_classic() + 
  scale_color_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                      labels = c("-1 SD (22.32 minutes/hr)", "Mean (22.82 minutes/hr)",
                                 "+1 SD (23.32 minutes/hr)")) +
  scale_fill_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                     labels = c("-1 SD (22.32 minutes/hr)", "Mean (22.82 minutes/hr)",
                                 "+1 SD (23.32 minutes/hr)")) +
  labs(x = "Maximum temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       colour = "Parent brooding level (minutes/hr)", fill = "Parent brooding level (minutes/hr)") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18))

# Print
print(temp_max_brooding_blups_predicted)

# Save
ggsave('temp_max_brooding_blups_predicted.png', plot = temp_max_brooding_blups_predicted, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 8, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# With points colored 
temp_max_brooding_blups_col_predicted <- ggplot(pred, aes(x = x, y = predicted)) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  geom_line(size = 1.5, aes(col = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, col = group, fill = group), alpha = .1,
              col = NA) +
  scale_color_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                      labels = c("-1 SD (22.32 minutes/hr)", "Mean (22.82 minutes/hr)",
                                 "+1 SD (23.32 minutes/hr)")) +
  scale_fill_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                     labels = c("-1 SD (22.32 minutes/hr)", "Mean (22.82 minutes/hr)",
                                "+1 SD (23.32 minutes/hr)")) +
  labs(x = "Maximum temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       color = "Parent brooding level (minutes/hr)", fill = "Parent brooding level (minutes/hr)") +
  new_scale("color") + 
  geom_point(data = sub, aes(x = nest_max_temp, y = mass_pre_obs, col = brooding_blups), size = 1.5) +
  scale_color_viridis(option = "viridis", end = 0.8) +
  labs(color = "Parent brooding level (minutes/hr)") 

## IQR temp and parental care adjusted
pred <- ggpredict(mass_iqr_temp_brooding_blups_adj_lmer, c("nest_iqr_temp", "brooding_blups"))

plot(pred) 

sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_iqr_temp) &
                !is.na(x = brooding_blups))

temp_iqr_brooding_blups_predicted <- ggplot(pred, aes(x = x, y = predicted)) + 
  geom_line(size = 1.5, aes(col = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, col = group, fill = group), alpha = .1,
              col = NA) +
  geom_point(data = sub, aes(x = nest_iqr_temp, y = mass_pre_obs), size = 1.5) +
  theme_classic() + 
  scale_color_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                      labels = c("-1 SD (22.32 minutes/hr)", "Mean (22.82 minutes/hr)",
                                 "+1 SD (23.32 minutes/hr)")) +
  scale_fill_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                     labels = c("-1 SD (22.32 minutes/hr)", "Mean (22.82 minutes/hr)",
                                "+1 SD (23.32 minutes/hr)")) +
  labs(x = "IQR of temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       colour = "Parent brooding level (minutes/hr)", fill = "Parent brooding level (minutes/hr)") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18))

# Print
print(temp_iqr_brooding_blups_predicted)

# Save
ggsave('temp_iqr_brooding_blups_predicted.png', plot = temp_iqr_brooding_blups_predicted, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 8, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# With points colored 
temp_iqr_brooding_blups_col_predicted <- ggplot(pred, aes(x = x, y = predicted)) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  geom_line(size = 1.5, aes(col = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, col = group, fill = group), alpha = .1,
              col = NA) +
  scale_color_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                      labels = c("-1 SD (22.32 minutes/hr)", "Mean (22.82 minutes/hr)",
                                 "+1 SD (23.32 minutes/hr)")) +
  scale_fill_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                     labels = c("-1 SD (22.32 minutes/hr)", "Mean (22.82 minutes/hr)",
                                "+1 SD (23.32 minutes/hr)")) +
  labs(x = "IQR of temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       color = "Parent brooding level (minutes/hr)", fill = "Parent brooding level (minutes/hr)") +
  new_scale("color") + 
  geom_point(data = sub, aes(x = nest_iqr_temp, y = mass_pre_obs, col = brooding_blups), size = 1.5) +
  scale_color_viridis(option = "viridis", end = 0.8) +
  labs(color = "Parent brooding level (minutes/hr)") 


## Combined plot
combined_temp_parent_brooding_predictions <- 
  ggarrange(temp_min_brooding_blups_predicted, temp_max_brooding_blups_predicted,
            temp_iqr_brooding_blups_predicted,
            ncol = 3, nrow = 1, common.legend = TRUE)

# Print
print(combined_temp_parent_brooding_predictions)

# Save
ggsave('combined_temp_parent_brooding_predictions.png', plot = combined_temp_parent_brooding_predictions, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 18, 
       height = 8, 
       units = c('in'), dpi = 300, limitsize = TRUE)

## Combined plot with points colored
combined_temp_parent_brooding_col_predictions <- 
  ggarrange(temp_min_brooding_blups_col_predicted, temp_max_brooding_blups_col_predicted,
            temp_iqr_brooding_blups_col_predicted,
            ncol = 3, nrow = 1, common.legend = TRUE,
            legend = "right")

# Print
print(combined_temp_parent_brooding_col_predictions)

# Save
ggsave('combined_temp_parent_brooding_col_predictions.png', plot = combined_temp_parent_brooding_col_predictions, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 22, 
       height = 8, 
       units = c('in'), dpi = 300, limitsize = TRUE)





# Minimum temp + nestling size adjusted
pred <- ggpredict(mass_min_temp_size_adj_lmer, c("nest_min_temp", "mid_size_order"))

plot(pred) 

sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_min_temp) &
                !is.na(x = mid_size_order))

temp_min_size_predicted <- ggplot(pred, aes(x = x, y = predicted)) + 
  geom_line(size = 1.5, aes(col = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, col = group, fill = group), alpha = .1,
              col = NA) +
  geom_point(data = sub, aes(x = nest_min_temp, y = mass_pre_obs), size = 1.5) +
  theme_classic() + 
  scale_color_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                      labels = c("Smallest", "Other")) +
  scale_fill_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                     labels = c("Smallest", "Other")) +
  labs(x = "Minimum temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       colour = "Nestling size", fill = "Nestling size")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18))

# Print
print(temp_min_size_predicted)

# Save
ggsave('temp_min_size_predicted.png', plot = temp_min_size_predicted, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 8, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 


## Maximum temp + nestling size adjusted
pred <- ggpredict(mass_max_temp_size_adj_lmer, c("nest_max_temp", "mid_size_order"))

plot(pred) 

sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_max_temp) &
                !is.na(x = mid_size_order))

temp_max_size_predicted <- ggplot(pred, aes(x = x, y = predicted)) + 
  geom_line(size = 1.5, aes(col = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, col = group, fill = group), alpha = .1,
              col = NA) +
  geom_point(data = sub, aes(x = nest_max_temp, y = mass_pre_obs), size = 1.5) +
  theme_classic() + 
  scale_color_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                      labels = c("Smallest", "Other")) +
  scale_fill_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                     labels = c("Smallest", "Other")) +
  labs(x = "Maximum temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       colour = "Nestling size", fill = "Nestling size") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18))

# Print
print(temp_max_size_predicted)

# Save
ggsave('temp_max_size_predicted.png', plot = temp_max_size_predicted, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 8, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

## IQR temp + nestling size adjusted
pred <- ggpredict(mass_iqr_temp_size_adj_lmer, c("nest_iqr_temp", "mid_size_order"))

plot(pred) 

sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_iqr_temp) &
                !is.na(x = mid_size_order))

temp_iqr_size_predicted <- ggplot(pred, aes(x = x, y = predicted)) + 
  geom_line(size = 1.5, aes(col = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, col = group, fill = group), alpha = .1,
              col = NA) +
  geom_point(data = sub, aes(x = nest_iqr_temp, y = mass_pre_obs), size = 1.5) +
  theme_classic() + 
  scale_color_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                      labels = c("Smallest", "Other")) +
  scale_fill_viridis(discrete = TRUE, option = "viridis", end = 0.8,
                     labels = c("Smallest", "Other")) +
  labs(x = "IQR of temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       colour = "Nestling size", fill = "Nestling size")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18))

# Print
print(temp_iqr_size_predicted)

# Save
ggsave('temp_iqr_size_predicted.png', plot = temp_iqr_size_predicted, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 8, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 


## Combined plot
combined_temp_size_predictions <- 
  ggarrange(temp_min_size_predicted, temp_max_size_predicted,
            temp_iqr_size_predicted,
            ncol = 1, nrow = 3)

# Print
print(combined_temp_size_predictions)

# Save
ggsave('combined_temp_size_predictions.png', plot = combined_temp_size_predictions, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 9.5, 
       height = 13.5, 
       units = c('in'), dpi = 300, limitsize = TRUE)






# Make model vis with stratified models - Minimum temp

pred_low <- ggpredict(mass_min_temp_blups_low_adj_lmer, c("nest_min_temp"))

plot(pred_low) 

pred_med <- ggpredict(mass_min_temp_blups_med_adj_lmer, c("nest_min_temp"))

plot(pred_med)

pred_high <- ggpredict(mass_min_temp_blups_high_adj_lmer, c("nest_min_temp"))

plot(pred_high)

sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_min_temp) &
                !is.na(x = feeding_expontd_blups))

cols <- c("#440154", "#2a788e", "#7ad151")

temp_min_feeding_blups_predicted <- ggplot() +
  geom_line(data = pred_low, aes(x = x, y = predicted), col = "#440154", size = 1.5) +
  geom_ribbon(data = pred_low, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = .1,
              fill = "#440154", col = NA) +
  geom_line(data = pred_med, aes(x = x, y = predicted), col = "#2a788e", size = 1.5) +
  geom_ribbon(data = pred_med, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = .1,
              fill = "#2a788e", col = NA) +
  geom_line(data = pred_high, aes(x = x, y = predicted), col = "#7ad151", size = 1.5) +
  geom_ribbon(data = pred_high, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = .1,
              fill = "#7ad151", col = NA) +
  geom_point(data = sub, aes(x = nest_min_temp, y = mass_pre_obs,
                             col = as.character(feeding_expontd_blups_strat)), size = 2) +
  theme_classic() + 
  labs(x = "Minimum temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       colour = "Parent feeding level", fill = "Parent feeding level") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  scale_color_manual(values = cols, labels=c("Low", "Avg", "High"),
                     aesthetics = c("colour", "fill"))


# Print
print(temp_min_feeding_blups_predicted)


# Save
ggsave('temp_min_feeding_blups_predicted.png', plot = temp_min_feeding_blups_predicted, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 8, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 




# Make model vis for no interaction model - max temp
pred <- ggpredict(mass_max_temp_blups_noint_adj_lmer, c("nest_max_temp"))

plot(pred) 


sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_max_temp) &
                !is.na(x = feeding_expontd_blups))


temp_max_feeding_blups_predicted <- ggplot() +
  geom_line(data = pred, aes(x = x, y = predicted), col = "black", size = 1.5) +
  geom_ribbon(data = pred, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = .1,
              fill = "black", col = NA) +
  geom_point(data = sub, aes(x = nest_max_temp, y = mass_pre_obs,
                             col = as.character(feeding_expontd_blups_strat)), size = 2) +
  theme_classic() + 
  labs(x = "Maximum temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       colour = "Parent feeding level") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  scale_color_manual(values = cols, labels=c("Low", "Avg", "High"),
                     aesthetics = c("colour"))


# Print
print(temp_max_feeding_blups_predicted)


# Save
ggsave('temp_max_feeding_blups_predicted.png', plot = temp_max_feeding_blups_predicted, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 8, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

# Make model vis for no interaction model - iqr temp
pred <- ggpredict(mass_iqr_temp_blups_noint_adj_lmer, c("nest_iqr_temp"))

plot(pred) 


sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_iqr_temp) &
                !is.na(x = feeding_expontd_blups))


temp_iqr_feeding_blups_predicted <- ggplot() +
  geom_line(data = pred, aes(x = x, y = predicted), col = "black", size = 1.5) +
    geom_ribbon(data = pred, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = .1,
                fill = "black", col = NA) +
  geom_point(data = sub, aes(x = nest_iqr_temp, y = mass_pre_obs,
                             col = as.character(feeding_expontd_blups_strat)), size = 2) +
  theme_classic() + 
  labs(x = "IQR of temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       colour = "Parent feeding level") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  scale_color_manual(values = cols, labels=c("Low", "Avg", "High"),
                     aesthetics = c("colour"))


# Print
print(temp_iqr_feeding_blups_predicted)


# Save
ggsave('temp_iqr_feeding_blups_predicted.png', plot = temp_iqr_feeding_blups_predicted, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 8, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 




# Make model vis for no interaction model - min temp and brooding
pred <- ggpredict(mass_min_temp_brooding_blups_noint_adj_lmer, c("nest_min_temp"))

plot(pred) 


sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_min_temp) &
                !is.na(x = brooding_blups))


temp_min_brooding_blups_predicted <- ggplot() +
  geom_line(data = pred, aes(x = x, y = predicted), col = "black", size = 1.5) +
    geom_ribbon(data = pred, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = .1,
                fill = "black", col = NA) +
  geom_point(data = sub, aes(x = nest_min_temp, y = mass_pre_obs,
                             col = as.character(brooding_blups_strat)), size = 2) +
  theme_classic() + 
  labs(x = "Minimum temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       colour = "Parent brooding level") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  scale_color_manual(values = cols, labels=c("Low", "Avg", "High"),
                     aesthetics = c("colour"))


# Print
print(temp_min_brooding_blups_predicted)


# Save
ggsave('temp_min_brooding_blups_predicted.png', plot = temp_min_brooding_blups_predicted, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 8, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 


# Make model vis for no interaction model - max temp and brooding
pred <- ggpredict(mass_max_temp_brooding_blups_noint_adj_lmer, c("nest_max_temp"))

plot(pred) 


sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_max_temp) &
                !is.na(x = brooding_blups))


temp_max_brooding_blups_predicted <- ggplot() +
  geom_line(data = pred, aes(x = x, y = predicted), col = "black", size = 1.5) +
    geom_ribbon(data = pred, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = .1,
                fill = "black", col = NA) +
  geom_point(data = sub, aes(x = nest_max_temp, y = mass_pre_obs,
                             col = as.character(brooding_blups_strat)), size = 2) +
  theme_classic() + 
  labs(x = "Maximum temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       colour = "Parent brooding level") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  scale_color_manual(values = cols, labels=c("Low", "Avg", "High"),
                     aesthetics = c("colour"))


# Print
print(temp_max_brooding_blups_predicted)


# Save
ggsave('temp_max_brooding_blups_predicted.png', plot = temp_max_brooding_blups_predicted, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 8, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 


# Make model vis for no interaction model - iqr temp and brooding
pred <- ggpredict(mass_iqr_temp_brooding_blups_noint_adj_lmer, c("nest_iqr_temp"))

plot(pred) 


sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = nest_iqr_temp) &
                !is.na(x = brooding_blups))


temp_iqr_brooding_blups_predicted <- ggplot() +
  geom_line(data = pred, aes(x = x, y = predicted), col = "black", size = 1.5) +
    geom_ribbon(data = pred, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = .1,
                fill = "black", col = NA) +
  geom_point(data = sub, aes(x = nest_iqr_temp, y = mass_pre_obs,
                             col = as.character(brooding_blups_strat)), size = 2) +
  theme_classic() + 
  labs(x = "IQR of temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       colour = "Parent brooding level") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  scale_color_manual(values = cols, labels=c("Low", "Avg", "High"),
                     aesthetics = c("colour"))


# Print
print(temp_iqr_brooding_blups_predicted)


# Save
ggsave('temp_iqr_brooding_blups_predicted.png', plot = temp_iqr_brooding_blups_predicted, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 8, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 



## Create dit and whisker coef plots for ABS
mass_min_temp_blups_lmer_tidy[1,4]


# Remove terms I don't need 
mass_min_temp_blups_lmer_tidy <-
  broom::tidy(mass_min_temp_blups_lmer) %>% 
  filter(effect != "ran_pars")

# Create plot 
mass_min_temp_blups_adj_dw_abs <- dwplot(mass_min_temp_blups_lmer_tidy) +
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  scale_color_manual(values = cols,
                     aesthetics = c("colour", "fill")) + 
  theme(legend.position = "none")


mass_min_temp_blups_adj_dw_abs 



## Make plots of temp effects before and after thermoreg
## Both thermo min temp
pred <- ggpredict(mass_min_both_thermo_adj_lmer, c("thermo_bef_min_temp"))

plot(pred) 

pred_aft <- ggpredict(mass_min_both_thermo_adj_lmer, c("thermo_aft_min_temp"))

plot(pred_aft) 

sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = thermo_bef_min_temp) &
                !is.na(x = thermo_aft_min_temp))

colors <- c("Before day 6" = "#440154", "After day 6" = "#7ad151")

temp_thermo_both_min_predicted <- ggplot() + 
  theme_classic() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  geom_line(data= pred, size = 1.5, aes(x = x, y = predicted, col = "Before day 6")) +
  geom_ribbon(data = pred, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high,
                               fill = "Before day 6"), alpha = .1, col = NA) +
  geom_line(data= pred_aft, size = 1.5, aes(x = x, y = predicted, col = "After day 6")) +
  geom_ribbon(data = pred_aft, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, 
                                   fill = "After day 6"),  alpha = .1, col = NA) +
  geom_point(data = sub, aes(x = thermo_bef_min_temp, y = mass_pre_obs, col = "Before day 6"), size = 1.5) +
  geom_point(data = sub, aes(x = thermo_aft_min_temp, y = mass_pre_obs, col = "After day 6"), size = 1.5) +
  labs(x = "Minimum temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       colour = "Before vs. after day 6", fill = "Before vs. after day 6") +
  scale_color_manual(values = colors, breaks = c("Before day 6", "After day 6")) +
  scale_fill_manual(values = colors, breaks = c("Before day 6", "After day 6"))

# Print
print(temp_thermo_both_min_predicted)

# Save
ggsave('temp_thermo_both_min_predicted.png', plot = temp_thermo_both_min_predicted, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 8, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 


## Both thermo max temp
pred <- ggpredict(mass_max_both_thermo_adj_lmer, c("thermo_bef_max_temp"))

plot(pred) 

pred_aft <- ggpredict(mass_max_both_thermo_adj_lmer, c("thermo_aft_max_temp"))

plot(pred_aft) 

sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = thermo_bef_max_temp) &
                !is.na(x = thermo_aft_max_temp))

colors <- c("Before day 6" = "#440154", "After day 6" = "#7ad151")

temp_thermo_both_max_predicted <- ggplot() + 
  theme_classic() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  geom_line(data= pred, size = 1.5, aes(x = x, y = predicted, col = "Before day 6")) +
  geom_ribbon(data = pred, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high,
                               fill = "Before day 6"), alpha = .1, col = NA) +
  geom_line(data= pred_aft, size = 1.5, aes(x = x, y = predicted, col = "After day 6")) +
  geom_ribbon(data = pred_aft, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, 
                                   fill = "After day 6"),  alpha = .1, col = NA) +
  geom_point(data = sub, aes(x = thermo_bef_max_temp, y = mass_pre_obs, col = "Before day 6"), size = 1.5) +
  geom_point(data = sub, aes(x = thermo_aft_max_temp, y = mass_pre_obs, col = "After day 6"), size = 1.5) +
  labs(x = "Maximum temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       colour = "Before vs. after day 6", fill = "Before vs. after day 6") +
  scale_color_manual(values = colors, breaks = c("Before day 6", "After day 6")) +
  scale_fill_manual(values = colors, breaks = c("Before day 6", "After day 6"))

# Print
print(temp_thermo_both_max_predicted)

# Save
ggsave('temp_thermo_both_max_predicted.png', plot = temp_thermo_both_max_predicted, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 8, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 



## Both thermo IQR temp
pred <- ggpredict(mass_iqr_both_thermo_adj_lmer, c("thermo_bef_iqr_temp"))

plot(pred) 

pred_aft <- ggpredict(mass_iqr_both_thermo_adj_lmer, c("thermo_aft_iqr_temp"))

plot(pred_aft) 

sub <- subset(late_nestling_parent_care,
              !is.na(x = mass_pre_obs) & 
                !is.na(x = thermo_bef_iqr_temp) &
                !is.na(x = thermo_aft_iqr_temp))

colors <- c("Before day 6" = "#440154", "After day 6" = "#7ad151")

temp_thermo_both_iqr_predicted <- ggplot() + 
  theme_classic() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  geom_line(data= pred, size = 1.5, aes(x = x, y = predicted, col = "Before day 6")) +
  geom_ribbon(data = pred, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high,
                               fill = "Before day 6"), alpha = .1, col = NA) +
  geom_line(data= pred_aft, size = 1.5, aes(x = x, y = predicted, col = "After day 6")) +
  geom_ribbon(data = pred_aft, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, 
                                   fill = "After day 6"),  alpha = .1, col = NA) +
  geom_point(data = sub, aes(x = thermo_bef_iqr_temp, y = mass_pre_obs, col = "Before day 6"), size = 1.5) +
  geom_point(data = sub, aes(x = thermo_aft_iqr_temp, y = mass_pre_obs, col = "After day 6"), size = 1.5) +
  labs(x = "IQR of temperature during nestling period (C)", y = "Mass on day 12 (g)", 
       colour = "Before vs. after day 6", fill = "Before vs. after day 6") +
  scale_color_manual(values = colors, breaks = c("Before day 6", "After day 6")) +
  scale_fill_manual(values = colors, breaks = c("Before day 6", "After day 6"))

# Print
print(temp_thermo_both_iqr_predicted)

# Save
ggsave('temp_thermo_both_iqr_predicted.png', plot = temp_thermo_both_iqr_predicted, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 8, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

## Combined plot with points colored
combined_thermo_both_predictions <- 
  ggarrange(temp_thermo_both_min_predicted, temp_thermo_both_max_predicted,
            temp_thermo_both_iqr_predicted,
            ncol = 3, nrow = 1, common.legend = TRUE,
            legend = "top")

# Print
print(combined_thermo_both_predictions)

# Save
ggsave('combined_thermo_both_predictions.png', plot = combined_thermo_both_predictions, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 18.5, 
       height = 8, 
       units = c('in'), dpi = 300, limitsize = TRUE)
