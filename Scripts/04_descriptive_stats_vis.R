####### Purpose: calculate descriptive stats and visualize barn swallow
####### nest mircoclimate and nestling growth dataset at nestling level
####### By: Sage Madden
####### Created: 12/19/2022
####### Last modified: 3/13/2025

##### Fill in code blocks + write down parental care stats


# Code Blocks
# 1: Configure work space
# 2: Load data
# 3: Univariate descriptive statistics
# 4: Univariate visualizations
# 5: Bivariate descriptive statistics
# 6: Bivariate visualizations


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
library('ggpubr')


### 1.3 Get Version and Session Info
R.Version()
sessionInfo()


###############################################################################
##############                        Load Data                  ##############
###############################################################################  

### Load RData
## Load RData tidy barn swallow data
load('Data/Tidy/tidy_parent_nestl_weather_data_10-4_with_BLUPs.RData')

# Create dataset at nest level (only one row per nest)
nest_dat <- late_nestling_parent_care
nest_dat$duplicate <- duplicated(nest_dat$nest_id)
nest_dat <- nest_dat %>% filter(duplicate == FALSE) %>%
  select(-duplicate)

# Remove nests excluded from analyses from Govee dataset
govee_daily <- filter(govee_daily, nest_id != "schaaps 110" & nest_id != "hayes 7")

###############################################################################
##############             Univariate descriptive stats          ##############
###############################################################################

### Univariate statistics
## Parental care
# Feeding BLUPs
univar_feeding_blups <- late_nestling_parent_care %>%
  summarise (n = sum(!is.na(feeding_expontd_blups)),
             avg = round (mean(feeding_expontd_blups, 
                               na.rm = T),2),
             stdev = round (sd(feeding_expontd_blups, 
                               na.rm = T), 2),
             med = round(median(feeding_expontd_blups,
                                na.rm = T), 2),
             min = round(min(feeding_expontd_blups,
                             na.rm = T), 2),
             max = round(max(feeding_expontd_blups,
                             na.rm = T), 2)
  )



# Feeding rate
univar_feeding_rate <- prim_merged %>%
  summarise (n = sum(!is.na(total_visits)),
             avg = round (mean(total_feeding_visits/(obs_duration/3600), 
                               na.rm = T),2),
             stdev = round (sd(total_feeding_visits/(obs_duration/3600), 
                               na.rm = T), 2),
             med = round(median(total_feeding_visits/(obs_duration/3600),
                                na.rm = T), 2),
             min = round(min(total_feeding_visits/(obs_duration/3600),
                             na.rm = T), 2),
             max = round(max(total_feeding_visits/(obs_duration/3600),
                             na.rm = T), 2)
  )


# Brooding duration
univar_brooding_dur <- prim_merged %>%
  summarise (n = sum(!is.na(total_brooding_duration)),
             avg = round (mean(total_brooding_duration/(obs_duration/3600), 
                               na.rm = T)/60,2),
             stdev = round (sd(total_brooding_duration/(obs_duration/3600), 
                               na.rm = T)/60, 2),
             med = round(median(total_brooding_duration/(obs_duration/3600),
                                na.rm = T)/60, 2),
             min = round(min(total_brooding_duration/(obs_duration/3600),
                             na.rm = T)/60, 2),
             max = round(max(total_brooding_duration/(obs_duration/3600),
                             na.rm = T)/60, 2)
  )

# Check correlation of feeding and brooding
prim_merged <- mutate(prim_merged, 
                      brooding_prop = total_brooding_duration/(obs_duration/3600),
                      feeding_rate = total_feeding_visits/(obs_duration/3600))


cor.test(prim_merged$brooding_prop, prim_merged$feeding_rate, method = "spearman")

# Bind the parental care stats together
univar_parental_care <- rbind(univar_feeding_blups, univar_feeding_rate, univar_brooding_dur)

univar_parental_care$variable_name <- c("feeding_blups", "feeding_rate", "brooding_duration")

## save the data frame of summary stats as a pdf into output file
pdf('Output/univar_feeding_blups.pdf', height = 3, width = 14)
grid.table(univar_feeding_blups)
dev.off()


## Brood size
univar_brood_size <- nest_dat %>%
  summarise (n = sum(!is.na(nestling_number)),
             avg = round (mean(nestling_number, 
                               na.rm = T),2),
             stdev = round (sd(nestling_number, 
                               na.rm = T), 2),
             med = round(median(nestling_number,
                                na.rm = T), 2),
             min = round(min(nestling_number,
                             na.rm = T), 2),
             max = round(max(nestling_number,
                             na.rm = T), 2)
  )

iqr_brood_size <- nest_dat %>%
  summarise (n = sum(!is.na(nestling_number)),
             iqr = round (IQR(nestling_number, 
                               na.rm = T),2)
  )

univar_brood_size$variable_name <- c("brood_size")

## save the data frame of summary stats as a pdf into output file
pdf('Output/univar_brood_size.pdf', height = 3, width = 14)
grid.table(univar_brood_size)
dev.off()

## Nestling age
univar_nestling_age <- nest_dat %>%
  summarise (n = sum(!is.na(nestling_age)),
             avg = round (mean(nestling_age, 
                               na.rm = T),2),
             stdev = round (sd(nestling_age, 
                               na.rm = T), 2),
             med = round(median(nestling_age,
                                na.rm = T), 2),
             min = round(min(nestling_age,
                             na.rm = T), 2),
             max = round(max(nestling_age,
                             na.rm = T), 2)
  )

univar_nestling_age$variable_name <- c("nestling_age")

## save the data frame of summary stats as a pdf into output file
pdf('Output/univar_nestling_age.pdf', height = 3, width = 14)
grid.table(univar_nestling_age)
dev.off()

## Colony size
nest_dat$num_pairs <- as.integer(nest_dat$num_pairs)
univar_colony_size <- nest_dat %>%
  summarise (n = sum(!is.na(num_pairs)),
             avg = round (mean(num_pairs, 
                               na.rm = T),2),
             stdev = round (sd(num_pairs, 
                               na.rm = T), 2),
             med = round(median(num_pairs,
                                na.rm = T), 2),
             min = round(min(num_pairs,
                             na.rm = T), 2),
             max = round(max(num_pairs,
                             na.rm = T), 2)
  )

univar_colony_size$variable_name <- c("colony_size")

## save the data frame of summary stats as a pdf into output file
pdf('Output/univar_colony_size.pdf', height = 3, width = 14)
grid.table(univar_colony_size)
dev.off()

## Hatch date
univar_hatch_date <- nest_dat %>%
  summarise (n = sum(!is.na(days_summer)),
             avg = round (mean(days_summer, 
                               na.rm = T),2),
             stdev = round (sd(days_summer, 
                               na.rm = T), 2),
             med = round(median(days_summer,
                                na.rm = T), 2),
             min = round(min(days_summer,
                             na.rm = T), 2),
             max = round(max(days_summer,
                             na.rm = T), 2)
  )

iqr_hatch_date <- nest_dat %>%
  summarise (n = sum(!is.na(hatch_date)),
             iqr = round (IQR(hatch_date, 
                              na.rm = T),2)
  )

univar_hatch_date$variable_name <- c("hatch_date")

## save the data frame of summary stats as a pdf into output file
pdf('Output/univar_hatch_date.pdf', height = 3, width = 14)
grid.table(univar_hatch_date)
dev.off()

## Nest temperature
# Median temp
univar_temp_med <- nest_dat %>%
  summarise (n = sum(!is.na(nest_med_temp)),
             avg = round (mean(nest_med_temp, 
                               na.rm = T),2),
             stdev = round (sd(nest_med_temp, 
                               na.rm = T), 2),
             med = round(median(nest_med_temp,
                                na.rm = T), 2),
             min = round(min(nest_med_temp,
                             na.rm = T), 2),
             max = round(max(nest_med_temp,
                             na.rm = T), 2)
  )

# Min temp
univar_temp_min <- nest_dat %>%
  summarise (n = sum(!is.na(nest_min_temp)),
             avg = round (mean(nest_min_temp, 
                               na.rm = T),2),
             stdev = round (sd(nest_min_temp, 
                               na.rm = T), 2),
             med = round(median(nest_min_temp,
                                na.rm = T), 2),
             min = round(min(nest_min_temp,
                             na.rm = T), 2),
             max = round(max(nest_min_temp,
                             na.rm = T), 2)
  )

# Max temp
univar_temp_max <- nest_dat %>%
  summarise (n = sum(!is.na(nest_max_temp)),
             avg = round (mean(nest_max_temp, 
                               na.rm = T),2),
             stdev = round (sd(nest_max_temp, 
                               na.rm = T), 2),
             med = round(median(nest_max_temp,
                                na.rm = T), 2),
             min = round(min(nest_max_temp,
                             na.rm = T), 2),
             max = round(max(nest_max_temp,
                             na.rm = T), 2)
  )

# IQR of temp
univar_temp_iqr <- nest_dat %>%
  summarise (n = sum(!is.na(nest_iqr_temp)),
             avg = round (mean(nest_iqr_temp, 
                               na.rm = T),2),
             stdev = round (sd(nest_iqr_temp, 
                               na.rm = T), 2),
             med = round(median(nest_iqr_temp,
                                na.rm = T), 2),
             min = round(min(nest_iqr_temp,
                             na.rm = T), 2),
             max = round(max(nest_iqr_temp,
                             na.rm = T), 2)
  )

# Before thermo median
univar_temp_bef_thermo <- nest_dat %>%
  summarise (n = sum(!is.na(thermo_bef_med_temp)),
             avg = round (mean(thermo_bef_med_temp, 
                               na.rm = T),2),
             stdev = round (sd(thermo_bef_med_temp, 
                               na.rm = T), 2),
             med = round(median(thermo_bef_med_temp,
                                na.rm = T), 2),
             min = round(min(thermo_bef_med_temp,
                             na.rm = T), 2),
             max = round(max(thermo_bef_med_temp,
                             na.rm = T), 2)
  )

# After thermo median
univar_temp_aft_thermo <- nest_dat %>%
  summarise (n = sum(!is.na(thermo_aft_med_temp)),
             avg = round (mean(thermo_aft_med_temp, 
                               na.rm = T),2),
             stdev = round (sd(thermo_aft_med_temp, 
                               na.rm = T), 2),
             med = round(median(thermo_aft_med_temp,
                                na.rm = T), 2),
             min = round(min(thermo_aft_med_temp,
                             na.rm = T), 2),
             max = round(max(thermo_aft_med_temp,
                             na.rm = T), 2)
  )

univar_temp <- rbind(univar_temp_med, univar_temp_min, univar_temp_max, 
                     univar_temp_iqr, univar_temp_bef_thermo, univar_temp_aft_thermo)

univar_temp$variable_name <- c("nest_med_temp", "nest_min_temp", 
                               "nest_max_temp", "nest_iqr_temp",
                               "thermo_bef_med_temp", "thermo_aft_med_temp")

## save the data frame of summary stats as a pdf into output file
pdf('Output/univar_temp.pdf', height = 3, width = 14)
grid.table(univar_temp)
dev.off()

## Nestling growth 
# Wing day 12
univar_growth_wing <- late_nestling_parent_care %>%
  summarise (n = sum(!is.na(rt_wing_length)),
             avg = round (mean(rt_wing_length, 
                               na.rm = T),2),
             stdev = round (sd(rt_wing_length, 
                               na.rm = T), 2),
             med = round(median(rt_wing_length,
                                na.rm = T), 2),
             min = round(min(rt_wing_length,
                             na.rm = T), 2),
             max = round(max(rt_wing_length,
                             na.rm = T), 2)
  )

# Mass day 12
univar_growth_mass <- late_nestling_parent_care %>%
  summarise (n = sum(!is.na(mass_pre_obs)),
             avg = round (mean(mass_pre_obs, 
                               na.rm = T),2),
             stdev = round (sd(mass_pre_obs, 
                               na.rm = T), 2),
             med = round(median(mass_pre_obs,
                                na.rm = T), 2),
             min = round(min(mass_pre_obs,
                             na.rm = T), 2),
             max = round(max(mass_pre_obs,
                             na.rm = T), 2)
  )

univar_growth <- rbind(univar_growth_wing, univar_growth_mass)

univar_growth$variable_name <- c("univar_growth_wing", "univar_growth_mass")

## save the data frame of summary stats as a pdf into output file
pdf('Output/univar_growth.pdf', height = 3, width = 14)
grid.table(univar_growth)
dev.off()


###############################################################################
##############             Univariate visualization               ##############
###############################################################################

## Parental care 
# Feeding visits across develpment
feeding_hist <- prim_merged %>%
  ggplot(aes(x = total_feeding_visits, fill = obs_state)) +
  facet_wrap(~obs_state, ncol = 1) +
  geom_histogram(alpha=0.6, col = "gray5",
                 position = 'identity', 
                 binwidth = 5) +
  labs(title = 'Histogram of total feeding visits across development',
       x ='Total feeding visits', 
       y ='Frequency') +
  theme_classic() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6)

# Print plot 
print(feeding_hist)

# Save plot
ggsave('feeding_hist.png', plot = feeding_hist, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 


# Brooding duration across develpment
brooding_hist <- prim_merged %>%
  ggplot(aes(x = (total_brooding_duration/60), fill = obs_state)) +
  facet_wrap(~obs_state, ncol = 1) +
  geom_histogram(alpha=0.6, col = "gray5",
                 position = 'identity', 
                 binwidth = 5) +
  labs(title = 'Histogram of total brooding duration across development',
       x ='Total brooding duration (minutes)', 
       y ='Frequency') +
  theme_classic() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6)

# Print plot 
print(brooding_hist)

# Save plot
ggsave('brooding_hist.png', plot = brooding_hist, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 


## Brood size
brood_hist <- nest_dat %>%
  ggplot(aes(x = nestling_number, fill = nest_id)) +
  geom_bar(alpha=0.6, col = "gray5") +
  labs(title = 'Bar plot of brood size colored by nest',
       x ='Brood size', 
       y ='Frequency') +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_viridis(discrete = TRUE, alpha=0.6)

# Print plot 
print(brood_hist)

# Save plot
ggsave('brood_hist.png', plot = brood_hist, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

## Nestling age
age_hist <- nest_dat %>%
  ggplot(aes(x = nestling_age, fill = nest_id)) +
  geom_bar(alpha=0.6, col = "gray5") +
  labs(title = 'Bar plot of nestling age at last measure colored by nest',
       x ='Nestling age at last measure', 
       y ='Frequency') +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_viridis(discrete = TRUE, alpha=0.6)

# Print plot 
print(age_hist)

# Save plot
ggsave('age_hist.png', plot = age_hist, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

## Hatch date
hatch_hist <- nest_dat %>%
  ggplot(aes(x = days_summer, fill = nest_id)) +
  geom_histogram(alpha=0.6, col = "gray5",
                 position = 'stack', 
                 binwidth = 4) +
  labs(title = 'Histogram of hatch dates colored by nest',
       x ='Days since June 1', 
       y ='Frequency') +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6)

# Print plot 
print(hatch_hist)

# Save plot
ggsave('hatch_hist.png', plot = hatch_hist, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

## Site
site_hist <- nest_dat %>%
  ggplot(aes(x = site, fill = nest_id)) +
  geom_bar(alpha=0.6, col = "gray5") +
  labs(title = 'Bar plot of breeding sites colored by nest',
       x ='Site', 
       y ='Frequency') +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_viridis(discrete = TRUE, alpha=0.6)

# Print plot 
print(site_hist)

# Save plot
ggsave('site_hist.png', plot = site_hist, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

## Median temperature at nest 
med_temp_hist <- nest_dat %>%
  ggplot(aes(x = nest_med_temp, fill = nest_id)) +
  geom_histogram(alpha=0.6, col = "gray5",
                 position = 'stack', 
                 binwidth = 1) +
  labs(title = 'Histogram of median nest temperature colored by nest',
       x ='Median nest temperature (deg C)', 
       y ='Frequency') +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6)

# Print plot 
print(med_temp_hist)

# Save plot
ggsave('med_temp_hist.png', plot = med_temp_hist, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

## Minimum temperature at nest 
min_temp_hist <- nest_dat %>%
  ggplot(aes(x = nest_min_temp, fill = nest_id)) +
  geom_histogram(alpha=0.6, col = "gray5",
                 position = 'stack', 
                 binwidth = 1) +
  labs(title = 'Histogram of minimum nest temperature colored by nest',
       x ='Minimum nest temperature (deg C)', 
       y ='Frequency') +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6)

# Print plot 
print(min_temp_hist)

# Save plot
ggsave('min_temp_hist.png', plot = min_temp_hist, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

## Maximum temperature at nest 
max_temp_hist <- nest_dat %>%
  ggplot(aes(x = nest_max_temp, fill = nest_id)) +
  geom_histogram(alpha=0.6, col = "gray5",
                 position = 'stack', 
                 binwidth = 1) +
  labs(title = 'Histogram of maximum nest temperature colored by nest',
       x ='Maximum nest temperature (deg C)', 
       y ='Frequency') +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6)

# Print plot 
print(max_temp_hist)

# Save plot
ggsave('max_temp_hist.png', plot = max_temp_hist, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

## IQR temperature at nest 
iqr_temp_hist <- nest_dat %>%
  ggplot(aes(x = nest_iqr_temp, fill = nest_id)) +
  geom_histogram(alpha=0.6, col = "gray5",
                 position = 'stack', 
                 binwidth = 1) +
  labs(title = 'Histogram of IQR of nest temperature colored by nest',
       x ='IQR of nest temperature (deg C)', 
       y ='Frequency') +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6)

# Print plot 
print(iqr_temp_hist)

# Save plot
ggsave('iqr_temp_hist.png', plot = iqr_temp_hist, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

## Before thermo temperature at nest 
bef_thermo_temp_hist <- nest_dat %>%
  ggplot(aes(x = thermo_bef_med_temp, fill = nest_id)) +
  geom_histogram(alpha=0.6, col = "gray5",
                 position = 'stack', 
                 binwidth = 1) +
  labs(title = 'Histogram of median nest temperature before thermoregulatory
       independence colored by nest',
       x ='Median temperature before thermoregulatory independence (deg C)', 
       y ='Frequency') +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6)

# Print plot 
print(bef_thermo_temp_hist)

# Save plot
ggsave('bef_thermo_temp_hist.png', plot = bef_thermo_temp_hist, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

## After thermo temperature at nest 
aft_thermo_temp_hist <- nest_dat %>%
  ggplot(aes(x = thermo_aft_med_temp, fill = nest_id)) +
  geom_histogram(alpha=0.6, col = "gray5",
                 position = 'stack', 
                 binwidth = 1) +
  labs(title = 'Histogram of median nest temperature after thermoregulatory
       independence colored by nest',
       x ='Median temperature after thermoregulatory independence (deg C)', 
       y ='Frequency') +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6)

# Print plot 
print(aft_thermo_temp_hist)

# Save plot
ggsave('aft_thermo_temp_hist.png', plot = aft_thermo_temp_hist, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)


## Wing chord day 12
wing_12_hist <- late_nestling_parent_care %>%
  ggplot(aes(x = rt_wing_length, fill = nest_id)) +
  geom_histogram(alpha=0.6, col = "gray5",
                 position = 'stack', 
                 binwidth = 4) +
  labs(title = 'Histogram of right wing length on day 12 colored by nest',
       x ='Right wing length (mm)', 
       y ='Frequency') +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6)

# Print plot 
print(wing_12_hist)

# Save plot
ggsave('wing_12_hist.png', plot = wing_12_hist, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

## Mass day 12
mass_12_hist <- late_nestling_parent_care %>%
  ggplot(aes(x = mass_pre_obs, fill = nest_id)) +
  geom_histogram(alpha=0.6, col = "gray5",
                 position = 'stack', 
                 binwidth = 2) +
  labs(title = 'Histogram of mass on day 12 colored by nest',
       x ='Mass (g)', 
       y ='Frequency') +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6)

# Print plot 
print(mass_12_hist)

# Save plot
ggsave('mass_12_hist.png', plot = mass_12_hist, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)



###############################################################################
##############              Bivariate descriptive stats          ##############
###############################################################################

### Bivariate descriptive stats of growth by various predictors
## Growth by brood size 
# Wing by brood size 
bivar_wing_brood_size <- late_nestling_parent_care %>%
  group_by(nestling_number) %>%
  filter(!is.na(nestling_number)) %>%
  summarise (n = sum(!is.na(rt_wing_length)),
             avg = round (mean(rt_wing_length, 
                               na.rm = T),2),
             stdev = round (sd(rt_wing_length, 
                               na.rm = T), 2),
             med = round(median(rt_wing_length,
                                na.rm = T), 2),
             min = round(min(rt_wing_length,
                             na.rm = T), 2),
             max = round(max(rt_wing_length,
                             na.rm = T), 2)
  )

# Mass by brood size 
bivar_mass_brood_size <- late_nestling_parent_care %>%
  group_by(nestling_number) %>%
  filter(!is.na(nestling_number)) %>%
  summarise (n = sum(!is.na(mass_pre_obs)),
             avg = round (mean(mass_pre_obs, 
                               na.rm = T),2),
             stdev = round (sd(mass_pre_obs, 
                               na.rm = T), 2),
             med = round(median(mass_pre_obs,
                                na.rm = T), 2),
             min = round(min(mass_pre_obs,
                             na.rm = T), 2),
             max = round(max(mass_pre_obs,
                             na.rm = T), 2)
  )

bivar_growth_brood_size <- rbind(bivar_wing_brood_size, 
                                 bivar_mass_brood_size)

bivar_growth_brood_size$variable_name <- c(rep("rt_wing_length", 5),
                                           rep("mass_pre_obs", 5))
## save the data frame of summary stats as a pdf into output file
pdf('Output/bivar_growth_brood_size.pdf', height = 3, width = 14)
grid.table(bivar_growth_brood_size)
dev.off()

## Wing by relative nestling size
bivar_wing_relative_size <- late_nestling_parent_care %>%
  group_by(size_order) %>%
  filter(!is.na(size_order)) %>%
  summarise (n = sum(!is.na(rt_wing_length)),
             avg = round (mean(rt_wing_length, 
                               na.rm = T),2),
             stdev = round (sd(rt_wing_length, 
                               na.rm = T), 2),
             med = round(median(rt_wing_length,
                                na.rm = T), 2),
             min = round(min(rt_wing_length,
                             na.rm = T), 2),
             max = round(max(rt_wing_length,
                             na.rm = T), 2)
  )

## Mass growth by relative nestling size
bivar_mass_relative_size <- late_nestling_parent_care %>%
  group_by(size_order) %>%
  filter(!is.na(size_order)) %>%
  summarise (n = sum(!is.na(mass_pre_obs)),
             avg = round (mean(mass_pre_obs, 
                               na.rm = T),2),
             stdev = round (sd(mass_pre_obs, 
                               na.rm = T), 2),
             med = round(median(mass_pre_obs,
                                na.rm = T), 2),
             min = round(min(mass_pre_obs,
                             na.rm = T), 2),
             max = round(max(mass_pre_obs,
                             na.rm = T), 2)
  )

bivar_growth_relative_size <- rbind(bivar_wing_relative_size, 
                                    bivar_mass_relative_size)

bivar_growth_relative_size$variable_name <- c(rep("rt_wing_length", 2),
                                              rep("mass_pre_obs", 2))

## save the data frame of summary stats as a pdf into output file
pdf('Output/bivar_growth_relative_size.pdf', height = 3, width = 14)
grid.table(bivar_growth_relative_size)
dev.off()


## Mid devel wing by relative nestling size 
bivar_wing_mid_relative_size <- nestl_merged %>%
  filter(sample_state == "mid") %>%
  group_by(size_order) %>%
  filter(!is.na(size_order)) %>%
  summarise (n = sum(!is.na(rt_wing_length)),
             avg = round (mean(rt_wing_length, 
                               na.rm = T),2),
             stdev = round (sd(rt_wing_length, 
                               na.rm = T), 2),
             med = round(median(rt_wing_length,
                                na.rm = T), 2),
             min = round(min(rt_wing_length,
                             na.rm = T), 2),
             max = round(max(rt_wing_length,
                             na.rm = T), 2)
  )

## Mass growth by relative nestling size
bivar_mass_mid_relative_size <- nestl_merged %>%
  filter(sample_state == "mid") %>%
  group_by(size_order) %>%
  filter(!is.na(size_order)) %>%
  summarise (n = sum(!is.na(mass_pre_obs)),
             avg = round (mean(mass_pre_obs, 
                               na.rm = T),2),
             stdev = round (sd(mass_pre_obs, 
                               na.rm = T), 2),
             med = round(median(mass_pre_obs,
                                na.rm = T), 2),
             min = round(min(mass_pre_obs,
                             na.rm = T), 2),
             max = round(max(mass_pre_obs,
                             na.rm = T), 2)
  )

bivar_growth_mid_relative_size <- rbind(bivar_wing_mid_relative_size, 
                                    bivar_mass_mid_relative_size)

bivar_growth_mid_relative_size$variable_name <- c(rep("rt_wing_length", 2),
                                              rep("mass_pre_obs", 2))

## save the data frame of summary stats as a pdf into output file
pdf('Output/bivar_growth_mid_relative_size.pdf', height = 3, width = 14)
grid.table(bivar_growth_mid_relative_size)
dev.off()


## Wing by age
bivar_wing_age <- late_nestling_parent_care %>%
  group_by(nestling_age) %>%
  filter(!is.na(nestling_age)) %>%
  summarise (n = sum(!is.na(rt_wing_length)),
             avg = round (mean(rt_wing_length, 
                               na.rm = T),2),
             stdev = round (sd(rt_wing_length, 
                               na.rm = T), 2),
             med = round(median(rt_wing_length,
                                na.rm = T), 2),
             min = round(min(rt_wing_length,
                             na.rm = T), 2),
             max = round(max(rt_wing_length,
                             na.rm = T), 2)
  )

## Mass growth by nestling age
bivar_mass_age <- late_nestling_parent_care %>%
  group_by(nestling_age) %>%
  filter(!is.na(nestling_age)) %>%
  summarise (n = sum(!is.na(mass_pre_obs)),
             avg = round (mean(mass_pre_obs, 
                               na.rm = T),2),
             stdev = round (sd(mass_pre_obs, 
                               na.rm = T), 2),
             med = round(median(mass_pre_obs,
                                na.rm = T), 2),
             min = round(min(mass_pre_obs,
                             na.rm = T), 2),
             max = round(max(mass_pre_obs,
                             na.rm = T), 2)
  )

bivar_growth_age <- rbind(bivar_wing_age, 
                          bivar_mass_age)

bivar_growth_age$variable_name <- c(rep("rt_wing_length", 3),
                                    rep("mass_pre_obs", 3))

## save the data frame of summary stats as a pdf into output file
pdf('Output/bivar_growth_age.pdf', height = 3, width = 14)
grid.table(bivar_growth_age)
dev.off()

## Wing by colony size
bivar_wing_colony <- late_nestling_parent_care %>%
  group_by(num_pairs) %>%
  filter(!is.na(num_pairs)) %>%
  summarise (n = sum(!is.na(rt_wing_length)),
             avg = round (mean(rt_wing_length, 
                               na.rm = T),2),
             stdev = round (sd(rt_wing_length, 
                               na.rm = T), 2),
             med = round(median(rt_wing_length,
                                na.rm = T), 2),
             min = round(min(rt_wing_length,
                             na.rm = T), 2),
             max = round(max(rt_wing_length,
                             na.rm = T), 2)
  )

bivar_mass_colony <- late_nestling_parent_care %>%
  group_by(num_pairs) %>%
  filter(!is.na(num_pairs)) %>%
  summarise (n = sum(!is.na(mass_pre_obs)),
             avg = round (mean(mass_pre_obs, 
                               na.rm = T),2),
             stdev = round (sd(mass_pre_obs, 
                               na.rm = T), 2),
             med = round(median(mass_pre_obs,
                                na.rm = T), 2),
             min = round(min(mass_pre_obs,
                             na.rm = T), 2),
             max = round(max(mass_pre_obs,
                             na.rm = T), 2)
  )

bivar_growth_colony <- rbind(bivar_wing_colony, 
                             bivar_mass_colony)

bivar_growth_colony$variable_name <- c(rep("rt_wing_length", 5),
                                    rep("mass_pre_obs", 5))

## save the data frame of summary stats as a pdf into output file
pdf('Output/bivar_growth_colony.pdf', height = 6, width = 14)
grid.table(bivar_growth_colony)
dev.off()

## Wing by developmental stage
bivar_wing_devel <- nestl_merged %>%
  group_by(sample_state) %>%
  filter(!is.na(sample_state)) %>%
  summarise (n = sum(!is.na(rt_wing_length)),
             avg = round (mean(rt_wing_length, 
                               na.rm = T),2),
             stdev = round (sd(rt_wing_length, 
                               na.rm = T), 2),
             med = round(median(rt_wing_length,
                                na.rm = T), 2),
             min = round(min(rt_wing_length,
                             na.rm = T), 2),
             max = round(max(rt_wing_length,
                             na.rm = T), 2)
  )

# Mass by developmental stage
bivar_mass_devel <- nestl_merged %>%
  group_by(sample_state) %>%
  filter(!is.na(sample_state)) %>%
  summarise (n = sum(!is.na(mass_pre_obs)),
             avg = round (mean(mass_pre_obs, 
                               na.rm = T),2),
             stdev = round (sd(mass_pre_obs, 
                               na.rm = T), 2),
             med = round(median(mass_pre_obs,
                                na.rm = T), 2),
             min = round(min(mass_pre_obs,
                             na.rm = T), 2),
             max = round(max(mass_pre_obs,
                             na.rm = T), 2)
  )

bivar_growth_devel <- rbind(bivar_wing_devel, 
                             bivar_mass_devel)

bivar_growth_devel$variable_name <- c(rep("rt_wing_length", 3),
                                       rep("mass_pre_obs", 3))

## save the data frame of summary stats as a pdf into output file
pdf('Output/bivar_growth_devel.pdf', height = 6, width = 14)
grid.table(bivar_growth_devel)
dev.off()


### Bivariate descriptive stats of parental care by developmental stage
# Feeding by developmental stage
bivar_feeding_devel <- prim_merged %>%
  group_by(obs_state) %>%
  filter(!is.na(obs_state)) %>%
  summarise (n = sum(!is.na(total_feeding_visits)),
             avg = round (mean(total_feeding_visits/(obs_duration/3600), 
                               na.rm = T),2),
             stdev = round (sd(total_feeding_visits/(obs_duration/3600), 
                               na.rm = T), 2),
             med = round(median(total_feeding_visits/(obs_duration/3600),
                                na.rm = T), 2),
             min = round(min(total_feeding_visits/(obs_duration/3600),
                             na.rm = T), 2),
             max = round(max(total_feeding_visits/(obs_duration/3600),
                             na.rm = T), 2)
  )

early <- filter(prim_merged, obs_state == "early")
mid <- filter(prim_merged, obs_state == "mid")
late <- filter(prim_merged, obs_state == "late")

unique(early$nest_id)
unique(late$nest_id)


test <- filter(prim_merged, is.na(total_feeding_visits) == T)
unique(test$nest_id)

bivar_feeding_devel$variable_name <- c("total_feeding_rate")

## save the data frame of summary stats as a pdf into output file
pdf('Output/bivar_feeding_devel.pdf', height = 6, width = 14)
grid.table(bivar_feeding_devel)
dev.off()

# Brooding by developmental stage
bivar_brooding_devel <- prim_merged %>%
  group_by(obs_state) %>%
  filter(!is.na(obs_state)) %>%
  summarise (n = sum(!is.na(total_brooding_duration)),
             avg = round (mean(total_brooding_duration/(obs_duration/3600), 
                               na.rm = T)/60,2),
             stdev = round (sd(total_brooding_duration/(obs_duration/3600), 
                               na.rm = T)/60, 2),
             med = round(median(total_brooding_duration/(obs_duration/3600),
                                na.rm = T)/60, 2),
             min = round(min(total_brooding_duration/(obs_duration/3600),
                             na.rm = T)/60, 2),
             max = round(max(total_brooding_duration/(obs_duration/3600),
                             na.rm = T)/60, 2)
  )



bivar_brooding_devel$variable_name <- c("total_brooding_dur")

## save the data frame of summary stats as a pdf into output file
pdf('Output/bivar_brooding_devel.pdf', height = 6, width = 14)
grid.table(bivar_brooding_devel)
dev.off()


### Create summary stats table

# Get summary stats for all different variables into the same format 

# Right wing length
bivar_wing_table <- bivar_wing_devel
bivar_wing_table$variable_name <- bivar_wing_table$sample_state
bivar_wing_table$sample_state <- NULL

# Mass
bivar_mass_table <- bivar_mass_devel
bivar_mass_table$variable_name <- bivar_mass_table$sample_state
bivar_mass_table$sample_state <- NULL

# Feeding rate
bivar_feeding_table <- bivar_feeding_devel
bivar_feeding_table$variable_name <- bivar_feeding_table$obs_state
bivar_feeding_table$obs_state <- NULL

# Temperature
univar_temp_min$variable_name <- "Minimum temperature (C)"

univar_temp_max$variable_name <- "Maximum temperature (C)"

univar_temp_iqr$variable_name <- "Temperature variability (C)"

univar_brood_size$variable_name <- "Brood size (number of nestlings)"

univar_hatch_date$variable_name <- "Hatch date (days since June 1)"

univar_nestling_age$variable_name <- "Nest age at last measure (days since hatch)"

# Bind all variables together
summary_stats_df <- rbind(univar_temp_min, univar_temp_max, univar_temp_iqr, univar_brood_size,
                            univar_hatch_date, univar_nestling_age, bivar_feeding_table,
                            bivar_wing_table, bivar_mass_table)

summary_stats_df$variable_name[summary_stats_df$variable_name == "early"] <- "Early"
summary_stats_df$variable_name[summary_stats_df$variable_name == "mid"] <- "Mid"
summary_stats_df$variable_name[summary_stats_df$variable_name == "late"] <- "Late"

summary_stats_df$med <- NULL

summary_stats_df <- rename(summary_stats_df,
                           n = n,
                           Mean = avg, 
                           SD = stdev,
                           Minimum = min, 
                           Maximum = max)

summary_stats_table <- gt(summary_stats_df, rowname_col = "variable_name") %>%
  tab_header(
    title = md("**Supplemental Table 2.** Background characteristics of wild nestling Barn Swallows in Boulder County, CO. Total number of feeding visits, nestling mass, and nestling wing length are separately estimated for three developmental stages (early, mid, and late). Temperature variability is defined as the interquartile range.")
  ) %>%
  tab_stubhead(
    label = md("Variable")
  ) %>%
  tab_row_group(
    label = md("**Nestling mass (g) across development**"), 
    rows = c(13:15)
  ) %>%
  tab_row_group(
    label = md("**Right wing length (mm) across development**"), 
    rows = c(10:12)
  ) %>%
  tab_row_group(
    label = md("**Total feeding rate (counts/hour) across development**"), 
    rows = c(7:9)
  ) %>%
  tab_row_group(
    label = md("**Other nest characteristics**"), 
    rows = c(4:6)
  ) %>%
  tab_row_group(
    label = md("**Near-nest temperature**"), 
    rows = c(1:3)
  ) %>%
  tab_footnote(
    footnote = " Parental feeding behaviors are measured at the level of the nest and include the totals for both social parents (maternal and paternal care).",
    locations = cells_row_groups(groups = "**Total feeding rate (counts/hour) across development**")
  ) %>%
  opt_table_font(font = "Arial", size = 12)  %>%
  tab_options(footnotes.font.size = 10)


summary_stats_table
gtsave(summary_stats_table, filename = "Output/summary_stats_table.docx")
gtsave(summary_stats_table, filename = "Output/summary_stats_table.html")


citation(package = "ggeffects")


###############################################################################
##############             Bivariate visualizations              ##############
###############################################################################

### Visualization of nestling growth by other predictors
## Parental care
# Wing growth by feeding BLUPs
wing_feeding_blups_scatter <- late_nestling_parent_care %>%
  ggplot(aes(x = feeding_expontd_blups, y = rt_wing_length)) + 
  geom_jitter(size=2, alpha=0.9, width = 0.2) +
  geom_smooth() + 
  labs(title = 'Scatter plot of right wing length on day 12 by feeding rate BLUPs',
       x ='Feeding rate BLUPs', 
       y ='Right wing length (mm)') +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8)

# Print plot 
print(wing_feeding_blups_scatter)

# Save plot
ggsave('wing_feeding_blups_scatter.png', plot = wing_feeding_blups_scatter, 
       device = NULL, 
       path = 'Output', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# Wing growth by feeding BLUPs
mass_feeding_blups_scatter <- late_nestling_parent_care %>%
  ggplot(aes(x = feeding_expontd_blups, y = mass_pre_obs)) + 
  geom_jitter(size=2, alpha=0.9, width = 0.2) +
  geom_smooth() + 
  labs(title = 'Scatter plot of mass on day 12 by feeding rate BLUPs',
       x ='Feeding rate BLUPs', 
       y ='Mass (g)') +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8)

# Print plot 
print(mass_feeding_blups_scatter)

# Save plot
ggsave('mass_feeding_blups_scatter.png', plot = mass_feeding_blups_scatter, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)


## Temperature 
# Wing growth by median temperature 
wing_med_temp_scatter <- late_nestling_parent_care %>%
  ggplot(aes(x = nest_med_temp, y = rt_wing_length)) + 
  geom_jitter(size=2, alpha=0.9, width = 0.2) +
  geom_smooth() + 
  labs(title = 'Scatter plot of right wing length on day 12 by median nest temperature',
       x ='Median nest temperature', 
       y ='Right wing length (mm)') +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8)

# Print plot 
print(wing_med_temp_scatter)

# Save plot
ggsave('wing_med_temp_scatter.png', plot = wing_med_temp_scatter, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# Mass growth by median temperature 
mass_med_temp_scatter <- late_nestling_parent_care %>%
  ggplot(aes(x = nest_med_temp, y = mass_pre_obs)) + 
  geom_jitter(size=2, alpha=0.9, width = 0.2) +
  geom_smooth() + 
  labs(title = 'Scatter plot of mass on day 12 by median nest temperature',
       x ='Median nest temperature', 
       y ='Mass (g)') +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8)

# Print plot 
print(mass_med_temp_scatter)

# Save plot
ggsave('mass_med_temp_scatter.png', plot = mass_med_temp_scatter, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# Wing growth by min temperature 
wing_min_temp_scatter <- late_nestling_parent_care %>%
  ggplot(aes(x = nest_min_temp, y = rt_wing_length)) + 
  geom_jitter(size=2, alpha=0.9, width = 0.2) +
  geom_smooth() + 
  labs(title = 'Scatter plot of right wing length on day 12 by minimum nest temperature',
       x ='Minimum nest temperature', 
       y ='Right wing length (mm)') +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8)

# Print plot 
print(wing_min_temp_scatter)

# Save plot
ggsave('wing_min_temp_scatter.png', plot = wing_min_temp_scatter, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# Mass growth by median temperature 
mass_min_temp_scatter <- late_nestling_parent_care %>%
  ggplot(aes(x = nest_min_temp, y = mass_pre_obs)) + 
  geom_jitter(size=2, alpha=0.9, width = 0.2) +
  geom_smooth() + 
  labs(title = 'Scatter plot of mass on day 12 by minimum nest temperature',
       x ='Minimum nest temperature', 
       y ='Mass (g)') +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8)

# Print plot 
print(mass_min_temp_scatter)

# Save plot
ggsave('mass_min_temp_scatter.png', plot = mass_min_temp_scatter, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# Wing growth by maximum temperature 
wing_max_temp_scatter <- late_nestling_parent_care %>%
  ggplot(aes(x = nest_max_temp, y = rt_wing_length)) + 
  geom_jitter(size=2, alpha=0.9, width = 0.2) +
  geom_smooth() + 
  labs(title = 'Scatter plot of right wing length on day 12 by maximum nest temperature',
       x ='Maximum nest temperature', 
       y ='Right wing length (mm)') +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8)

# Print plot 
print(wing_max_temp_scatter)

# Save plot
ggsave('wing_max_temp_scatter.png', plot = wing_max_temp_scatter, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# Mass growth by median temperature 
mass_max_temp_scatter <- late_nestling_parent_care %>%
  ggplot(aes(x = nest_max_temp, y = mass_pre_obs)) + 
  geom_jitter(size=2, alpha=0.9, width = 0.2) +
  geom_smooth() + 
  labs(title = 'Scatter plot of mass on day 12 by maximum nest temperature',
       x ='Maximum nest temperature', 
       y ='Mass (g)') +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8)

# Print plot 
print(mass_max_temp_scatter)

# Save plot
ggsave('mass_max_temp_scatter.png', plot = mass_max_temp_scatter, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)


# Wing growth by IQR of temperature 
wing_iqr_temp_scatter <- late_nestling_parent_care %>%
  ggplot(aes(x = nest_iqr_temp, y = rt_wing_length)) + 
  geom_jitter(size=2, alpha=0.9, width = 0.2) +
  geom_smooth() + 
  labs(title = 'Scatter plot of right wing length on day 12 by IQR of nest temperature',
       x ='IQR of nest temperature', 
       y ='Right wing length (mm)') +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8)

# Print plot 
print(wing_iqr_temp_scatter)

# Save plot
ggsave('wing_iqr_temp_scatter.png', plot = wing_iqr_temp_scatter, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# Mass growth by median temperature 
mass_iqr_temp_scatter <- late_nestling_parent_care %>%
  ggplot(aes(x = nest_iqr_temp, y = mass_pre_obs)) + 
  geom_jitter(size=2, alpha=0.9, width = 0.2) +
  geom_smooth() + 
  labs(title = 'Scatter plot of mass on day 12 by IQR of nest temperature',
       x ='IQR of nest temperature', 
       y ='Mass (g)') +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8)

# Print plot 
print(mass_iqr_temp_scatter)

# Save plot
ggsave('mass_iqr_temp_scatter.png', plot = mass_iqr_temp_scatter, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# Wing growth by temp before thermo independence
wing_bef_thermo_temp_scatter <- late_nestling_parent_care %>%
  ggplot(aes(x = thermo_bef_max_temp, y = rt_wing_length)) + 
  geom_jitter(size=2, alpha=0.9, width = 0.2) +
  geom_smooth() + 
  labs(title = 'Scatter plot of right wing length on day 12 by maximum
       nest temperature before thermoregulatory independence',
       x ='Maximum nest temperature', 
       y ='Right wing length (mm)') +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8)

# Print plot 
print(wing_bef_thermo_temp_scatter)

# Save plot
ggsave('wing_bef_thermo_temp_scatter.png', plot = wing_bef_thermo_temp_scatter, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# Mass growth by median temperature 
mass_before_thermo_temp_scatter <- late_nestling_parent_care %>%
  ggplot(aes(x = thermo_bef_max_temp, y = mass_pre_obs)) + 
  geom_jitter(size=2, alpha=0.9, width = 0.2) +
  geom_smooth() + 
  labs(title = 'Scatter plot of mass on day 12 by maximum
       nest temperature before thermoregulatory independence',
       x ='Maximum nest temperature', 
       y ='Mass (g)') +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8)

# Print plot 
print(mass_before_thermo_temp_scatter)

# Save plot
ggsave('mass_before_thermo_temp_scatter.png', plot = mass_before_thermo_temp_scatter, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)


# Wing growth by temp after thermo independence
wing_aft_thermo_temp_scatter <- late_nestling_parent_care %>%
  ggplot(aes(x = thermo_aft_max_temp, y = rt_wing_length)) + 
  geom_jitter(size=2, alpha=0.9, width = 0.2) +
  geom_smooth() + 
  labs(title = 'Scatter plot of right wing length on day 12 by maximum
       nest temperature after thermoregulatory independence',
       x ='Maximum nest temperature', 
       y ='Right wing length (mm)') +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8)

# Print plot 
print(wing_aft_thermo_temp_scatter)

# Save plot
ggsave('wing_aft_thermo_temp_scatter.png', plot = wing_aft_thermo_temp_scatter, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# Mass growth by median temperature 
mass_aft_thermo_temp_scatter <- late_nestling_parent_care %>%
  ggplot(aes(x = thermo_aft_max_temp, y = mass_pre_obs)) + 
  geom_jitter(size=2, alpha=0.9, width = 0.2) +
  geom_smooth() + 
  labs(title = 'Scatter plot of mass on day 12 by maximum
       nest temperature after thermoregulatory independence',
       x ='Maximum nest temperature', 
       y ='Mass (g)') +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8)

# Print plot 
print(mass_aft_thermo_temp_scatter)

# Save plot
ggsave('mass_after_thermo_temp_scatter.png', plot = mass_aft_thermo_temp_scatter, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)


## Brood size 
# Wing growth by brood size
wing_growth_brood_size_box <- late_nestling_parent_care %>%
  filter(!is.na(nestling_number)) %>%
  ggplot(aes(x = as.factor(nestling_number), y = rt_wing_length, 
             col = as.factor(nestling_number))) + 
  geom_boxplot() +
  geom_jitter(size=2, alpha=0.8, width = 0.2) +
  labs(title = 'Scatter plot of right wing length by brood size',
       x ='Brood size', 
       y ='Right wing length (mm)', col = "Brood size") +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8)

# Print plot 
print(wing_growth_brood_size_box)

# Save plot
ggsave('wing_growth_brood_size_box.png', plot = wing_growth_brood_size_box, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# Mass growth by brood size
mass_growth_brood_size_box <- late_nestling_parent_care %>%
  filter(!is.na(nestling_number)) %>%
  ggplot(aes(x = as.factor(nestling_number), y = mass_pre_obs, 
             col = as.factor(nestling_number))) + 
  geom_boxplot() +
  geom_jitter(size=2, alpha=0.8, width = 0.2) +
  labs(title = 'Scatter plot of mass by brood size',
       x ='Brood size', 
       y ='Mass (g)', col = "Brood size") +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8)

# Print plot 
print(mass_growth_brood_size_box)

# Save plot
ggsave('mass_growth_brood_size_box.png', plot = mass_growth_brood_size_box, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)


## Nestling age
# Wing growth by brood size
wing_growth_age_box <- late_nestling_parent_care %>%
  filter(!is.na(nestling_age)) %>%
  ggplot(aes(x = as.factor(nestling_age), y = rt_wing_length, 
             col = as.factor(nestling_age))) + 
  geom_boxplot() +
  geom_jitter(size=2, alpha=0.8, width = 0.2) +
  labs(title = 'Scatter plot of right wing length by age',
       x ='Nestling age', 
       y ='Right wing length (mm)', col = "Brood size") +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8)

# Print plot 
print(wing_growth_age_box)

# Save plot
ggsave('wing_growth_age_box.png', plot = wing_growth_age_box, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# Mass growth by age
mass_growth_brood_age <- late_nestling_parent_care %>%
  filter(!is.na(nestling_age)) %>%
  ggplot(aes(x = as.factor(nestling_age), y = mass_pre_obs, 
             col = as.factor(nestling_age))) + 
  geom_boxplot() +
  geom_jitter(size=2, alpha=0.8, width = 0.2) +
  labs(title = 'Scatter plot of mass by nestling age',
       x ='Nestling age', 
       y ='Mass (g)', col = "Brood size") +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8)

# Print plot 
print(mass_growth_brood_age)

# Save plot
ggsave('mass_growth_brood_age.png', plot = mass_growth_brood_age, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

## Relative nestling size
# Wing growth by relative nestling size
wing_growth_relative_size_box <- late_nestling_parent_care %>%
  filter(!is.na(size_order)) %>%
  ggplot(aes(x = as.factor(size_order), y = rt_wing_length, 
             col = as.factor(size_order))) + 
  geom_boxplot() +
  geom_jitter(size=2, alpha=0.8, width = 0.2) +
  labs(title = 'Scatter plot of right wing length by relative size',
       x ='Relative size', 
       y ='Right wing length (mm)', col = "Relative size") +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8)

# Print plot 
print(wing_growth_relative_size_box)

# Save plot
ggsave('wing_growth_relative_size_box.png', plot = wing_growth_relative_size_box, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# Mass growth by relative size
mass_growth_relative_size_box <- late_nestling_parent_care %>%
  filter(!is.na(size_order)) %>%
  ggplot(aes(x = as.factor(size_order), y = mass_pre_obs, 
             col = as.factor(size_order))) + 
  geom_boxplot() +
  geom_jitter(size=2, alpha=0.8, width = 0.2) +
  labs(title = 'Scatter plot of mass by relative size',
       x ='Relative size', 
       y ='Mass (g)', col = "Relative size") +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8)

# Print plot 
print(mass_growth_relative_size_box)

# Save plot
ggsave('mass_growth_relative_size_box.png', plot = mass_growth_relative_size_box, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# Wing growth by hatch date
wing_growth_hatch_box <- late_nestling_parent_care %>%
  filter(!is.na(days_summer)) %>%
  ggplot(aes(x = as.factor(days_summer), y = rt_wing_length, 
             col = as.factor(days_summer))) + 
  geom_boxplot() +
  geom_jitter(size=2, alpha=0.8, width = 0.2) +
  labs(title = 'Box plot of right wing length by hatch date',
       x ='Hatch date (days since 6/1)', 
       y ='Right wing length (mm)', col = "Hatch date (days since 6/1)") +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8,
                      option = "turbo")

# Print plot 
print(wing_growth_hatch_box)

# Save plot
ggsave('wing_growth_hatch_box.png', plot = wing_growth_hatch_box, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# Mass growth by hatch date
mass_growth_hatch_box <- late_nestling_parent_care %>%
  filter(!is.na(days_summer)) %>%
  ggplot(aes(x = as.factor(days_summer), y = mass_pre_obs, 
             col = as.factor(days_summer))) + 
  geom_boxplot() +
  geom_jitter(size=2, alpha=0.8, width = 0.2) +
  labs(title = 'Box plot of mass by hatch date',
       x ='Hatch date (days since 6/1)', 
       y ='Mass (g)', col = "Hatch date (days since 6/1)") +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8,
                      option = "turbo")

# Print plot 
print(mass_growth_hatch_box)

# Save plot
ggsave('mass_growth_hatch_box.png', plot = mass_growth_hatch_box, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)



### Visualization of temperature by site and nest
# Min temp by site
min_temp_site_box <- late_nestling_parent_care %>%
  filter(!is.na(nest_min_temp)) %>%
  ggplot(aes(x = site, y = nest_min_temp, 
             col = site)) + 
  geom_boxplot() +
  geom_jitter(size=2, alpha=0.8, width = 0.2) +
  labs(title = 'Box plot of minimum nest temperature by site',
       x ='Site', 
       y ='Minimum temperature (C)', col = "Site") +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8,
                      option = "turbo")

# Print plot 
print(min_temp_site_box)

# Save plot
ggsave('min_temp_site_box.png', plot = min_temp_site_box, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# Max temp by site
max_temp_site_box <- late_nestling_parent_care %>%
  filter(!is.na(nest_max_temp)) %>%
  ggplot(aes(x = site, y = nest_max_temp, 
             col = site)) + 
  geom_boxplot() +
  geom_jitter(size=2, alpha=0.8, width = 0.2) +
  labs(title = 'Box plot of maximum nest temperature by site',
       x ='Site', 
       y ='Maximum temperature (C)', col = "Site") +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8,
                      option = "turbo")

# Print plot 
print(max_temp_site_box)

# Save plot
ggsave('max_temp_site_box.png', plot = max_temp_site_box, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# IQR temp by site
iqr_temp_site_box <- late_nestling_parent_care %>%
  filter(!is.na(nest_iqr_temp)) %>%
  ggplot(aes(x = site, y = nest_iqr_temp, 
             col = site)) + 
  geom_boxplot() +
  geom_jitter(size=2, alpha=0.8, width = 0.2) +
  labs(title = 'Box plot of nest temperature IQR by site',
       x ='Site', 
       y ='IQR of temperature (C)', col = "Site") +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8,
                      option = "turbo")

# Print plot 
print(iqr_temp_site_box)

# Save plot
ggsave('iqr_temp_site_box.png', plot = iqr_temp_site_box, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# Combined plot 
combined_temp_by_site <- 
  ggarrange(min_temp_site_box, max_temp_site_box,
            iqr_temp_site_box,
            labels = c("A", "B", "C"),
            ncol = 1, nrow = 3)

print(combined_temp_by_site)

ggsave('combined_temp_by_site.png', plot = combined_temp_by_site, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 6, 
       height = 12, 
       units = c('in'), dpi = 300, limitsize = TRUE) 


## Temperature by nest
# Min temp by nest
min_temp_nest_box <- govee_daily %>%
  filter(!is.na(daily_min_temp)) %>%
  ggplot(aes(x = nest_id, y = daily_min_temp, 
             col = nest_id)) + 
  geom_boxplot() +
  geom_jitter(size=2, alpha=0.8, width = 0.2) +
  labs(title = 'Box plot of minimum daily temperature by nest',
       x ='Nest ID', 
       y ='Minimum daily temperature (C)', col = "Site") +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8,
                      option = "turbo")

# Print plot 
print(min_temp_nest_box)

# Save plot
ggsave('min_temp_nest_box.png', plot = min_temp_nest_box, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)


# Max temp by nest
max_temp_nest_box <- govee_daily %>%
  filter(!is.na(daily_max_temp)) %>%
  ggplot(aes(x = nest_id, y = daily_max_temp, 
             col = nest_id)) + 
  geom_boxplot() +
  geom_jitter(size=2, alpha=0.8, width = 0.2) +
  labs(title = 'Box plot of maximum daily temperature by nest',
       x ='Nest ID', 
       y ='Maximum daily temperature (C)', col = "Site") +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8,
                      option = "turbo")

# Print plot 
print(max_temp_nest_box)

# Save plot
ggsave('max_temp_nest_box.png', plot = max_temp_nest_box, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# IQR temp by nest
iqr_temp_nest_box <- govee_daily %>%
  filter(!is.na(daily_iqr_temp)) %>%
  ggplot(aes(x = nest_id, y = daily_iqr_temp, 
             col = nest_id)) + 
  geom_boxplot() +
  geom_jitter(size=2, alpha=0.8, width = 0.2) +
  labs(title = 'Box plot of IQR of daily temperature by nest',
       x ='Nest ID', 
       y ='Daily temperature IQR (C)', col = "Site") +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8,
                      option = "turbo")

# Print plot 
print(iqr_temp_nest_box)

# Save plot
ggsave('iqr_temp_nest_box.png', plot = iqr_temp_nest_box, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 12, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)

## Plot of daily max temps at each nest
govee_daily$nest_ids <- govee_daily$nest_id
govee_split <- separate(govee_daily,
                        nest_id, 
                        into = c("site", "nest"),
                        sep = " ")
unique(govee_split$nest)

for(i in 2:(length(govee_split$site)-1)){
  if(govee_split$nest[i] != govee_split$nest[i + 1]){
    govee_split$date[i] <- NA
  } else if(govee_split$nest[i] != govee_split$nest[i - 1]){
    govee_split$date[i] <- NA
  } else{
    govee_split$date[i] <- govee_split$date[i]
  }
}
govee_split <- filter(govee_split, is.na(date) == F)


daily_max_temps_site_scatter <- 
  ggplot() + 
  geom_point(aes(x = date, y = daily_max_temp, 
                 col = site), data = govee_split) +
  geom_line(aes(x = date, y = daily_max_temp, 
                col = site, group = nest_ids), data = govee_split,
            size = 0.8) +
  labs(x ='Date', 
       y ='Daily maximum temperature (C)',
       col = "Site") +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8,
                      option = "viridis") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18))
  

# Print plot 
print(daily_max_temps_site_scatter)

ggsave('daily_max_temps_site_scatter.png', plot = daily_max_temps_site_scatter, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 9, 
       height = 5.5, 
       units = c('in'), dpi = 300, limitsize = TRUE)


## Plot of daily min temps at each nest
daily_min_temps_site_scatter <- 
  ggplot() + 
  geom_point(aes(x = date, y = daily_min_temp, 
                 col = site), data = govee_split) +
  geom_line(aes(x = date, y = daily_min_temp, 
                col = site, group = nest_ids), data = govee_split,
            size = 0.8) +
  labs(x ='Date', 
       y ='Daily minimum temperature (C)',
       col = "Site") +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8,
                      option = "viridis") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18))


# Print plot 
print(daily_min_temps_site_scatter)

ggsave('daily_min_temps_site_scatter.png', plot = daily_min_temps_site_scatter, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 9, 
       height = 5.5, 
       units = c('in'), dpi = 300, limitsize = TRUE)

## Plot of daily IQR of temps at each nest
daily_iqr_temps_site_scatter <- 
  ggplot() + 
  geom_point(aes(x = date, y = daily_iqr_temp, 
                 col = site), data = govee_split) +
  geom_line(aes(x = date, y = daily_iqr_temp, 
                col = site, group = nest_ids), data = govee_split,
            size = 0.8) +
  labs(x ='Date', 
       y ='Daily temperature variability (C)',
       col = "Site") +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, alpha=0.8,
                      option = "viridis") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18))


# Print plot 
print(daily_iqr_temps_site_scatter)

ggsave('daily_iqr_temps_site_scatter.png', plot = daily_iqr_temps_site_scatter, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 9, 
       height = 5.5, 
       units = c('in'), dpi = 300, limitsize = TRUE)

# Combined plot 
combined_temp_date <- 
  ggarrange(daily_min_temps_site_scatter, daily_max_temps_site_scatter,
            daily_iqr_temps_site_scatter,
            labels = c("a", "b", "c"),
            ncol = 1, nrow = 3, common.legend = TRUE, 
            legend = "right")

print(combined_temp_date)

ggsave('combined_temp_date.png', plot = combined_temp_date, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 10, 
       height = 14, 
       units = c('in'), dpi = 300, limitsize = TRUE) 


### BLUPs by strata for statified analysis
## Create three categories for BLUPs to allow stratification 
late_nestling_parent_care <- late_nestling_parent_care %>%
  mutate(feeding_expontd_blups_strat =  as.integer(Hmisc::cut2(feeding_expontd_blups, g=3)))

late_nestling_parent_care$blups_factor <- as.factor(late_nestling_parent_care$feeding_expontd_blups_strat)

levels(late_nestling_parent_care$blups_factor) <- c("Low", "Med", "High")

cols <- c("#481567FF", "#20A387FF", "#95D840FF")



blups_stratified_boxplot <- ggplot(data = late_nestling_parent_care, 
                                   aes(x = blups_factor, y = feeding_expontd_blups, 
                                       fill = blups_factor), col = "black") +
  geom_jitter(aes(col = blups_factor), size = 1.5, alpha = 0.5, width = 0.2) +
  geom_boxplot(alpha = 0.5, size = 1) +
  theme_classic() +
  labs(x = "Parent feeding level", y = "Feeding BLUPs (visits/hr)", fill = "Parent feeding level") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18)) +
  theme(legend.position="none") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols)


# Print
print(blups_stratified_boxplot)

# Save
ggsave('blups_stratified_boxplot.png', plot = blups_stratified_boxplot, 
       device = NULL, 
       path = 'Output/', scale = 1, width = 6, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE)


