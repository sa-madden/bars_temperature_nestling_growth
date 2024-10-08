####### Purpose: run functions to calculate parental care variables from
####### animal behavior pro logs and link up with parental care
####### metadata
####### By: Sage Madden
####### Created: 4/27/2021
####### Last modified: 10/4/2024

## Load relevant packages 
library(tidyverse)
library(lubridate)
library(data.table)

## Load parental care data
parent_data <- read.csv("./Data/parental_care_data.csv")

## Clean up parental care metadata to match AB pro logs
parent_data_lowercase <- parent_data
parent_data_lowercase$site <- tolower(parent_data_lowercase$site)
parent_data_lowercase$site <- gsub(" ", "", parent_data_lowercase$site)
parent_data_lowercase$nest <- as.character(parent_data_lowercase$nest)

## Create empty columns in parental care data sheet for each parental care variabke
## I want to calculate 
parent_data_lowercase$total_visits <- NA
parent_data_lowercase$female_visits <- NA
parent_data_lowercase$male_visits <- NA
parent_data_lowercase$total_feeding_visits <- NA
parent_data_lowercase$female_feeding_visits <- NA
parent_data_lowercase$male_feeding_visits <- NA
parent_data_lowercase$total_an_duration <- NA
parent_data_lowercase$female_an_duration <- NA
parent_data_lowercase$male_an_duration <- NA
parent_data_lowercase$total_brooding_duration <- NA
parent_data_lowercase$female_brooding_duration <- NA
parent_data_lowercase$male_brooding_duration <- NA
parent_data_lowercase$obs_duration <- NA
parent_data_lowercase$total_sanitizing_visits <- NA
parent_data_lowercase$female_sanitizing_visits <- NA
parent_data_lowercase$male_sanitizing_visits <- NA


## Read in AB pro logs, calculate parental care measures,
## and add to the relevant row of the parental care data
for(i in 1:length(parent_data_lowercase$nest)){
  file_name <-  paste("Data/", 
                      "revised_abpro_logs/", 
                      parent_data_lowercase$site[i], "_", 
                      as.character(parent_data_lowercase$nest[i]), "_", 
                      "d", as.character(parent_data_lowercase$nestling_age[i]), "_", 
                      "abpro", "_",
                      "rev",
                      ".csv", 
                      sep = "")
  if(file.exists(file_name)){
    abpro_df <- read.csv(file_name)
    abpro_df <- abpro_df %>% mutate(across(where(is.character), 
                                           str_trim, 
                                           side = c("both"))) 
    abpro_df[] <- lapply(abpro_df, gsub, pattern='\"', replacement='')
    abpro_df$Duration_s <- as.numeric(abpro_df$Duration_s)
    rbind()
    # Add all visits
    parent_data_lowercase$total_visits[i] <- calculate_visitation(abpro_df)[1,2]
    parent_data_lowercase$female_visits[i] <- calculate_visitation(abpro_df)[2,2]
    parent_data_lowercase$male_visits[i] <- calculate_visitation(abpro_df)[3,2]
    # Add feeding visits
    parent_data_lowercase$total_feeding_visits[i] <- calculate_feeding(abpro_df)[1,2]
    parent_data_lowercase$female_feeding_visits[i] <- calculate_feeding(abpro_df)[2,2]
    parent_data_lowercase$male_feeding_visits[i] <- calculate_feeding(abpro_df)[3,2]
    # Add duration at nest
    parent_data_lowercase$total_an_duration[i] <- calculate_an_duration(abpro_df)[1,2]
    parent_data_lowercase$female_an_duration[i] <- calculate_an_duration(abpro_df)[2,2]
    parent_data_lowercase$male_an_duration[i] <- calculate_an_duration(abpro_df)[3,2]
    # Add brooding duration
    parent_data_lowercase$total_brooding_duration[i] <- calculate_brooding_duration(abpro_df)[1,2]
    parent_data_lowercase$female_brooding_duration[i] <- calculate_brooding_duration(abpro_df)[2,2]
    parent_data_lowercase$male_brooding_duration[i] <- calculate_brooding_duration(abpro_df)[3,2]
    # Add sanitation visits 
    parent_data_lowercase$total_sanitizing_visits[i] <- calculate_sanitizing(abpro_df)[1,2]
    parent_data_lowercase$female_sanitizing_visits[i] <- calculate_sanitizing(abpro_df)[2,2]
    parent_data_lowercase$male_sanitizing_visits[i] <- calculate_sanitizing(abpro_df)[3,2]
    # Add total observation duration
    parent_data_lowercase$obs_duration[i] <- calculate_obs_duration(abpro_df)
  } 
}

# Remove rows with no total visits or duration, indicating no successful observation
filter(parent_data_lowercase, is.na(total_visits) == TRUE)
filter(parent_data_lowercase, is.na(total_an_duration) == TRUE)

colnames(parent_data_lowercase)

## SET durations to NA for MJA logs (these were determined to have been inaccurately recorded)
for(i in 1:length(parent_data_lowercase$female_band)){
  if(parent_data_lowercase$observer[i] == "MJA"){
    parent_data_lowercase$total_an_duration[i] <- NA
    parent_data_lowercase$female_an_duration[i] <- NA
    parent_data_lowercase$male_an_duration[i] <- NA
    parent_data_lowercase$total_brooding_duration[i] <- NA
    parent_data_lowercase$female_brooding_duration[i] <- NA
    parent_data_lowercase$male_brooding_duration[i] <- NA
  } 
}

# Select only the columns I need for the final dataset
parent_data_sel <- select(parent_data_lowercase, female_band, 
                                nest, site, obs_date, nestling_age, 
                           21:36)

# Write out the data
write.csv(parent_data_sel, file = "Data/Tidy/ab_pro_summary_data.csv")

