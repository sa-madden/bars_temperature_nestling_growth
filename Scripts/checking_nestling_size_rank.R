
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
library("dplyr")

## Graph Plotting and Visualization Packages
library ('ggplot2')

### 1.3 Get Version and Session Info
R.Version()
sessionInfo()

###############################################################################
##############                    2. Load RData                  ##############
###############################################################################  

### 2.1 Load RData
## Load RData tidy barn swallow data
load('Data/Tidy/tidy_parent_nestl_weather_data_10-4_with_BLUPs.RData')

# Make site a factor
late_nestling_parent_care$fsite <- as.factor(late_nestling_parent_care$site)
late_nestling_parent_care$fnest_id <- as.factor(late_nestling_parent_care$nest_id)

late_nestling_parent_care$num_pairs <- as.numeric(late_nestling_parent_care$num_pairs)

# Nestling size at day 8 and 12
nestl_ranked <- nestl_merged %>% 
  filter(sample_state != "early") %>%
  filter(nestling_band != "") %>%
  group_by(nest_id, sample_state) %>%
  mutate(size_rank = rank(rt_wing_length))

nestl_ranked[112, 3]

df <- nestl_ranked %>% select(nest_id, nestling_band, sample_state, size_rank) %>%
  spread(sample_state, size_rank) %>%
  filter(is.na(mid) == FALSE)


sum(df$mid == df$late)/length(df$mid)


nestl_ranked <- nestl_merged %>% 
  filter(sample_state != "early") %>%
  filter(nestling_band != "") %>%
  group_by(nest_id, sample_state) %>%
  mutate(size_rank = rank(mass_pre_obs))

nestl_ranked[112, 3]

df <- nestl_ranked %>% select(nest_id, nestling_band, sample_state, size_rank) %>%
  spread(sample_state, size_rank) %>%
  filter(is.na(mid) == FALSE)


sum(df$mid == df$late)/length(df$mid)


## Minimum size on day 8 vs. 12
late_nestling_parent_care_sub <- filter(late_nestling_parent_care, is.na(mid_size_order) == FALSE &
                                          is.na(size_order) == FALSE)

sum(late_nestling_parent_care_sub$size_order == late_nestling_parent_care_sub$mid_size_order)/length(late_nestling_parent_care_sub$female_band)


## Minimum size on day 8 vs. 12
late_nestling_parent_care_sub <- filter(late_nestling_parent_care, is.na(mid_size_order) == FALSE &
                                          is.na(size_order) == FALSE)

sum(late_nestling_parent_care_sub$size_order == late_nestling_parent_care_sub$mid_size_order)/length(late_nestling_parent_care_sub$female_band)



## relative size
mid_size_sub <- filter(late_nestling_parent_care, mid_size_order == "min")
late_size_sub <- filter(late_nestling_parent_care, size_order == "min")

relative_size_joined <- left_join(mid_size_sub, late_size_sub, by = "nest_id")

relative_size_joined_band <- left_join(mid_size_sub, late_size_sub, by = c("nest_id", "nestling_band"))
relative_size_joined_band <- filter(relative_size_joined_band, size_order.x == "min")
