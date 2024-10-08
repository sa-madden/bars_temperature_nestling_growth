####### Purpose: combine all govee datasheets into a single df and clean it
####### By: Sage Madden
####### Created: 3/1/2022
####### Last modified: 10/4/2024

# Code Blocks
# 1: Configure work space
# 2: Load data
# 3: Tidy data
# 4: Save data

###############################################################################
##############                Configure work space              ##############
###############################################################################

### Global options
# clear global environment
rm(list = ls())

# prevent R from automatically reading character strings as factors
options(stringsAsFactors = FALSE)


### Load relevant packages 
library(tidyverse)
library(lubridate)
library(data.table)


###############################################################################
##############                  Load data                        ##############
###############################################################################

### Read in govee .csv files and combine them into a single df

# Create a function to add the filename to the information read in for
# each .csv file
read_plus <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(filename = flnm)
}

# Creates a list, with each the name of each .csv file in the specified folder
# as an element, then create a dataframe by reading each csv file in the 
# list, row binding, and including the filename as part of each .csv
govee_tbl <-
  list.files(path = "./Data/used_nests_15min/", 
             pattern = "*.csv", 
             full.names = T) %>% 
  map_dfr(~read_plus(.))

str(govee_tbl)

# Changing to df because I'm more familiar with that format
govee_df <- data.frame(govee_tbl)
str(govee_df)



###############################################################################
##############                   Tidy data                      ##############
###############################################################################

### Clean govee data

# Extract nest and site info from each file name
govee_df_splt <- govee_df %>% separate(col = filename, 
                     into = c(NA, NA, NA, "file_name"), 
                     sep = "/")
unique(govee_df_splt$file_name)


govee_df_splt <- govee_df_splt %>%
  separate(col = file_name, 
           into = c("site", "nest"), 
           sep = "_") 
# Works okay despite warning. I did not fill in all the NAs because 
# there were differences in how files were named
unique(govee_df_splt$nest)
unique(govee_df_splt$site)
str(govee_df_splt)

# Rename columns
govee_df_splt <- rename(govee_df_splt, 
                        ymd_hms = Timestamp.for.sample.frequency.every.15.min.min, 
                        temp_c = Temperature_Celsius,
                        humid_perc = Relative_Humidity,
                        site = site, 
                        nest = nest)

# Make sure date time column is in correct lubridate format
govee_df_splt$ymd_hms <- ymd_hms(govee_df_splt$ymd_hms)
str(govee_df_splt)

# Split date and time
govee_df_splt$time <- format(as.POSIXct(govee_df_splt$ymd_hms, 
                                        format="%Y-%m-%d %H:%M:%S"),
                             "%H:%M:%S")

govee_df_splt$date <- format(as.POSIXct(govee_df_splt$ymd_hms, 
                                        format="%Y-%m-%d %H:%M:%S"),
                             "%Y-%m-%d")

#write.csv(govee_df_splt, file = "Output/govee_data_used_nests_unfiltered.csv")

## Filter to include dates and times from hatch day placement to 
## end time of last observation period 

# Read in parental care observation data
obs_dat <- read.csv("Data/2021_parental_care_data_entry.csv")
str(obs_dat)

# Format times and dates
obs_dat$obs_date <- mdy(obs_dat$obs_date)
obs_dat$obs_start_time <- hm(obs_dat$obs_start_time)
obs_dat$obs_end_time <- hm(obs_dat$obs_end_time)
str(obs_dat)

govee_df_splt$date <- ymd(govee_df_splt$date)
govee_df_splt$time <- hms(govee_df_splt$time)
str(govee_df_splt)

# Merge datasets 
govee_df_splt$nest <- as.character(govee_df_splt$nest)
obs_dat$nest <- as.character(obs_dat$nest)
obs_dat <- filter(obs_dat, site == "Hayes" | nestling_age > 9)
obs_dat <- select(obs_dat, site, nest, obs_date, obs_start_time, obs_end_time)
unique(obs_dat$site)
unique(govee_df_splt$site)

obs_dat$site <- gsub(" ", "", obs_dat$site)
obs_dat$site <- tolower(obs_dat$site)

# Now I've got the end obs date and time
# Join with govee data
joined_df <- left_join(obs_dat, govee_df_splt, by = c("nest", "site"))
unique(joined_df$site)
unique(joined_df$nest)

# Also need to get hatch date
# Read in nestling data, including hatch data
cal_dat <- read.csv("Data/2021_nestling_data_entry.csv")
str(cal_dat)

# Format date, site, and nest columns to match other datasheets
cal_dat$hatch_date <- mdy(cal_dat$hatch_date)

cal_dat$site <- gsub(" ", "", cal_dat$site)
cal_dat$site <- tolower(cal_dat$site)

cal_dat$nest <- as.character(cal_dat$nest)

cal_dat <- select(cal_dat, site, nest, hatch_date)

cal_dat$site[cal_dat$site == "urbanfarmgirlz"] <- "urbanfarm"

# Remove duplicates due to multiple nestlings per nest
cal_dat$duplicate <- duplicated(cal_dat) 
cal_dat <- cal_dat %>% filter(duplicate == FALSE) %>%
  select(-duplicate)

# Join hatch date with parental care and govee data
joined_df2 <- left_join(joined_df, cal_dat, 
                        by = c("site", "nest"))

str(joined_df2)

# Filter to only include measures on or after the hatch date
joined_df2 <- joined_df2 %>% filter(as.integer(date) >= 
                                      as.integer(hatch_date))
# Now starts at midnight on the hatch day

# Filter to end on the last obs day when the obs end
joined_df3 <- joined_df2 %>% filter(date <= obs_date)
joined_df3$filter_last_day <- NA
unique(joined_df3$site)

# Create an end time around the time of other observations for nests where we
# did not complete a final observation successfully
joined_df3$obs_end_time[is.na(joined_df3$obs_end_time) == TRUE] <- hms("7:0:0")

# Create a column indicating whether the logger time is before or after the obs time
# on the final day
for(i in 1:length(joined_df3$date)){
  if(joined_df3$date[i] == joined_df3$obs_date[i] & 
     joined_df3$time[i] > joined_df3$obs_end_time[i]){
    joined_df3$filter_last_day[i] <- "Y"
  } 
  else{
    joined_df3$filter_last_day[i] <- "N"
  }
}

# Remove logger data from after the observation on the last day
joined_df3 <- filter(joined_df3, filter_last_day == "N")

# ALL DONE YAY

# Now get rid of columns I don't need
joined_select <- joined_df3 %>% select(ymd_hms, temp_c, humid_perc, site, nest, 
                                       time, date)


###############################################################################
##############                       Save data                   ##############
###############################################################################
write.csv(joined_df3, file = "Data/Tidy/govee_data_for_filling_in_camera_obs.csv")
write.csv(joined_select, file = "Data/Tidy/govee_used_nests_filt.csv")


