####### Purpose: write functions to calculate parental care variables from
####### animal behavior pro logs
####### By: XXX and XXX
####### Created: 7/15/2021
####### Last modified: 06/10/2025


# Code Blocks
# 1: Configure work space
# 2: Create functions


###############################################################################
##############                Configure work space              ##############
###############################################################################

### Load relevant packages 
library(tidyverse)


###############################################################################
##############                Create functions                  ##############
###############################################################################

### Create functions to quantify different types of parental care behavior

# Function to calculate total visitation rate from an animal behavior pro .csv
calculate_visitation <- function(abpro_df) {
  #Calculate female visitation
  #Filter to include only female visits
  female_df <- dplyr::filter(abpro_df, Actor == "Female")
  #Create a column encoding the first AN of each visit as 1 and everything
  #else as 0
  for(i in 1:length(female_df$Actor)){
    if(i == 1 & female_df$Behavior[i] == "AN"){
      female_df$first_AN[i] <- 1
    } else if(i == 1 & female_df$Behavior[i] == "GFN"){
      female_df$first_AN[i] <- 0
    } else if(female_df$Behavior[i] == "AN" & 
              female_df$Behavior[i-1] == "GFN"){
      female_df$first_AN[i] <- 1
    } else{
      female_df$first_AN[i] <- 0
    }
  }
  #Sum female visits
  female_visits <- sum(female_df$first_AN)
  #Repeat for males
  male_df <- dplyr::filter(abpro_df, Actor == "Male")
  for(i in 1:length(male_df$Actor)){
    if(i == 1 & male_df$Behavior[i] == "AN"){
      male_df$first_AN[i] <- 1
    } else if(i == 1 & male_df$Behavior[i] == "GFN"){
      male_df$first_AN[i] <- 0
    } else if(male_df$Behavior[i] == "AN" & 
              male_df$Behavior[i-1] == "GFN"){
      male_df$first_AN[i] <- 1
    } else{
      male_df$first_AN[i] <- 0
    }
  }
  male_visits <- sum(male_df$first_AN)
  #Repeat for unknown
  unknown_df <- dplyr::filter(abpro_df, Actor == "Unknown")
  for(i in 1:length(unknown_df$Actor)){
    if(is.na(unknown_df$Behavior[i]) == TRUE){
     #The break is in here to deal with cases where the observer forget to
      #hit GFN for unknown at the start of the observation
       break
    } else if(i == 1 & unknown_df$Behavior[i] == "AN"){
      unknown_df$first_AN[i] <- 1
    } else if(i == 1 & unknown_df$Behavior[i] == "GFN"){
      unknown_df$first_AN[i] <- 0
    } else if(unknown_df$Behavior[i] == "AN" & 
              unknown_df$Behavior[i-1] == "GFN"){
      unknown_df$first_AN[i] <- 1
    } else{
      unknown_df$first_AN[i] <- 0
    }
  }
  #If else statement sets unknown visits to 0 if unknown was never entered 
  #in the observation. Otherwise, it sums unknown visits the same way as
  #female and male
  if(is.na(unknown_df$Behavior[1]) == TRUE){
    unknown_visits <- 0
  } else{
    unknown_visits <- sum(unknown_df$first_AN)
  }
  #Sum female, male, unknown visits to get total visits count
  total_visits <- female_visits + male_visits + unknown_visits
  total_visits
  # Create data frame to output
  data.frame(name = c("total_vis", "fem_vis", "male_vis"), 
             calculation = c(total_visits, female_visits,
                             male_visits))
}

# Function to calculate feeding rate from an animal behavior pro .csv
calculate_feeding <- function(abpro_df) {
  #Calculate female visitation
  #Filter to include only female visits
  female_df <- dplyr::filter(abpro_df, Actor == "Female" & Event_Type == "State start")
  #Create a column encoding each feeding visit as 1 and everything
  #else as 0
  for(i in 1:length(female_df$Actor)){
    if(female_df$Behavior[i] == "AN" & female_df$Modifier_1[i] == "Feeding"){
      female_df$feeding_vis[i] <- 1
    } else{
      female_df$feeding_vis[i] <- 0
    }
  }
  #Sum female visits
  female_feeding_visits <- sum(female_df$feeding_vis)
  #Repeat for males
  male_df <- dplyr::filter(abpro_df, Actor == "Male" & Event_Type == "State start")
  for(i in 1:length(male_df$Actor)){
    if(male_df$Behavior[i] == "AN" & male_df$Modifier_1[i] == "Feeding"){
      male_df$feeding_vis[i] <- 1
    } else{
      male_df$feeding_vis[i] <- 0
    }
  }
  male_feeding_visits <- sum(male_df$feeding_vis)
  #Repeat for unknown
  unknown_df <- dplyr::filter(abpro_df, Actor == "Unknown" & Event_Type == "State start")
  for(i in 1:length(unknown_df$Actor)){
    if(is.na(unknown_df$Behavior[i]) == TRUE){
      #The break is in here to deal with cases where the observer forget to
      #hit GFN for unknown at the start of the observation
      break
    }
      else if(unknown_df$Behavior[i] == "AN" & unknown_df$Modifier_1[i] == "Feeding"){
        unknown_df$feeding_vis[i] <- 1
      } else{
        unknown_df$feeding_vis[i] <- 0
      }
    }
  #If else statement sets unknown visits to 0 if unknown was never entered 
  #in the observation. Otherwise, it sums unknown visits the same way as
  #female and male
  if(is.na(unknown_df$Behavior[1]) == TRUE){
    unknown_feeding_visits <- 0
  } else{
    unknown_feeding_visits <- sum(unknown_df$feeding_vis)
  }
  #Sum female, male, unknown visits to get total visits count
  total_feeding_visits <- female_feeding_visits + male_feeding_visits + unknown_feeding_visits
  total_feeding_visits
  # Create dataframe to output
  data.frame(name = c("total_feeding_vis", "fem_feeding_vis", "male_feeding_vis"), 
             calculation = c(total_feeding_visits, 
                             female_feeding_visits,
                             male_feeding_visits))
}


# Function to calculate sanitation rate from an animal behavior pro .csv
calculate_sanitizing <- function(abpro_df) {
  #Calculate female visitation
  #Filter to include only female visits
  female_df <- dplyr::filter(abpro_df, Actor == "Female" & Event_Type == "State start")
  #Create a column encoding each feeding visit as 1 and everything
  #else as 0
  for(i in 1:length(female_df$Actor)){
    if(female_df$Behavior[i] == "AN" & female_df$Modifier_1[i] == "Sanitizing"){
      female_df$sanitizing_vis[i] <- 1
    } else{
      female_df$sanitizing_vis[i] <- 0
    }
  }
  #Sum female visits
  female_sanitizing_visits <- sum(female_df$sanitizing_vis)
  #Repeat for males
  male_df <- dplyr::filter(abpro_df, Actor == "Male" & Event_Type == "State start")
  for(i in 1:length(male_df$Actor)){
    if(male_df$Behavior[i] == "AN" & male_df$Modifier_1[i] == "Sanitizing"){
      male_df$sanitizing_vis[i] <- 1
    } else{
      male_df$sanitizing_vis[i] <- 0
    }
  }
  male_sanitizing_visits <- sum(male_df$sanitizing_vis)
  #Repeat for unknown
  unknown_df <- dplyr::filter(abpro_df, Actor == "Unknown" & Event_Type == "State start")
  for(i in 1:length(unknown_df$Actor)){
    if(is.na(unknown_df$Behavior[i]) == TRUE){
      #The break is in here to deal with cases where the observer forget to
      #hit GFN for unknown at the start of the observation
      break
    }
    else if(unknown_df$Behavior[i] == "AN" & unknown_df$Modifier_1[i] == "Sanitizing"){
      unknown_df$sanitizing_vis[i] <- 1
    } else{
      unknown_df$sanitizing_vis[i] <- 0
    }
  }
  #If else statement sets unknown visits to 0 if unknown was never entered 
  #in the observation. Otherwise, it sums unknown visits the same way as
  #female and male
  if(is.na(unknown_df$Behavior[1]) == TRUE){
    unknown_sanitizing_visits <- 0
  } else{
    unknown_sanitizing_visits <- sum(unknown_df$sanitizing_vis)
  }
  #Sum female, male, unknown visits to get total visits count
  total_sanitizing_visits <- female_sanitizing_visits + male_sanitizing_visits + unknown_sanitizing_visits
  total_sanitizing_visits
  # Create dataframe to output
  data.frame(name = c("total_sanitizing_vis", "fem_sanitizing_vis", "male_sanitizing_vis"), 
             calculation = c(total_sanitizing_visits, 
                             female_sanitizing_visits,
                             male_sanitizing_visits))
}

# Function to calculate duration at the nest from an animal behavior pro .csv
calculate_an_duration <- function(abpro_df) {
  #Calculate female visitation
  #Filter to include only female visits
  female_df <- dplyr::filter(abpro_df, Actor == "Female" & Behavior == "AN")
  #Sum the duration column to get total duration AN
  female_an_duration <- sum(female_df$Duration_s)
  #Repeat for males
  male_df <- dplyr::filter(abpro_df, Actor == "Male" & Behavior == "AN")
  male_an_duration <- sum(male_df$Duration_s)
  #Repeat for unknown
  unknown_df <- dplyr::filter(abpro_df, Actor == "Unknown" & Behavior == "AN")
  unknown_an_duration <- sum(unknown_df$Duration_s)
  #Sum female, male, unknown durations to get total duration
  total_an_duration <- female_an_duration + male_an_duration + unknown_an_duration
  data.frame(name = c("total_an_duration", "fem_an_duration", "male_an_duration"), 
             calculation = c(total_an_duration, 
                             female_an_duration, 
                             male_an_duration))
}


# unction to calculate duration brooding at the nest from an animal behavior pro .csv
calculate_brooding_duration <- function(abpro_df) {
  #Calculate female visitation
  #Filter to include only female visits
  female_df <- dplyr::filter(abpro_df, Actor == "Female" & Behavior == "AN" & 
                               Modifier_1 == "Brooding")
  #Sum the duration column to get total duration AN
  female_an_duration <- sum(female_df$Duration_s)
  #Repeat for males
  male_df <- dplyr::filter(abpro_df, Actor == "Male" & Behavior == "AN" & 
                             Modifier_1 == "Brooding")
  male_an_duration <- sum(male_df$Duration_s)
  #Repeat for unknown
  unknown_df <- dplyr::filter(abpro_df, Actor == "Unknown" & Behavior == "AN" & 
                                Modifier_1 == "Brooding")
  unknown_an_duration <- sum(unknown_df$Duration_s)
  #Sum female, male, unknown durations to get total duration
  total_an_duration <- female_an_duration + male_an_duration + unknown_an_duration
  #Divide by observation duration (hours) to get visitation rate
  data.frame(name = c("total_brooding_duration", "fem_brooding_duration", 
                      "male_brooding_duration"), 
             calculation = c(total_an_duration, 
                             female_an_duration, 
                             male_an_duration))
}


# Function to calculate duration brooding at the nest from an animal behavior pro .csv
calculate_sanitizing_duration <- function(abpro_df) {
  #Calculate female visitation
  #Filter to include only female visits
  female_df <- dplyr::filter(abpro_df, Actor == "Female" & Behavior == "AN" & 
                               Modifier_1 == "Sanitizing")
  #Sum the duration column to get total duration AN
  female_an_duration <- sum(female_df$Duration_s)
  #Repeat for males
  male_df <- dplyr::filter(abpro_df, Actor == "Male" & Behavior == "AN" & 
                             Modifier_1 == "Sanitizing")
  male_an_duration <- sum(male_df$Duration_s)
  #Repeat for unknown
  unknown_df <- dplyr::filter(abpro_df, Actor == "Unknown" & Behavior == "AN" & 
                                Modifier_1 == "Sanitizing")
  unknown_an_duration <- sum(unknown_df$Duration_s)
  #Sum female, male, unknown durations to get total duration
  total_an_duration <- female_an_duration + male_an_duration + unknown_an_duration
  #Divide by observation duration (hours) to get visitation rate
  data.frame(name = c("total_sanitizing_duration", "fem_sanitizing_duration", 
                      "male_sanitizing_duration"), 
             calculation = c(total_an_duration, 
                             female_an_duration, 
                             male_an_duration))
}


# Function to calculate total length of obs in seconds
calculate_obs_duration <- function(abpro_df) {
  #Divide by observation duration (hours) to get visitation rate
  obs_dur <- (as.numeric(abpro_df$Time_Relative_s[length(abpro_df$Time_Relative_s)]) -
                     as.numeric(abpro_df$Time_Relative_s[1]))
  obs_dur
}


## Function to calculate avg feeding and sanitizing visit lengths to add to sheets
## with missing time stamps
calculate_vis_length <- function(abpro_df) {
  #Calculate female visitation
  start_only_df <- dplyr::filter(abpro_df, Event_Type == "State start")
  #Create a column encoding each feeding visit as 1 and everything
  #else as 0
  for(i in 1:length(start_only_df$Actor)){
    if(start_only_df$Behavior[i] == "AN" & 
       start_only_df$Modifier_1[i] == "Feeding"){
      start_only_df$feeding_vis[i] <- 1
    } else{
      start_only_df$feeding_vis[i] <- 0
    }
  }
  #Sum female visits
  feeding_visit_number <- sum(start_only_df$feeding_vis)
  #Filter to include only feeding visits
  length_df <- dplyr::filter(abpro_df, Behavior == "AN" & 
                               Modifier_1 == "Feeding")
  #Sum the duration column to get total feeding duration
  feeding_duration <- sum(length_df$Duration_s)
  #Create a column encoding each feeding visit as 1 and everything
  #else as 0
  for(i in 1:length(start_only_df$Actor)){
    if(start_only_df$Behavior[i] == "AN" & 
       start_only_df$Modifier_1[i] == "Sanitizing"){
      start_only_df$sanitizing_vis[i] <- 1
    } else{
      start_only_df$sanitizing_vis[i] <- 0
    }
  }
  #Sum female visits
  san_visit_number <- sum(start_only_df$sanitizing_vis)
  #Filter to include only feeding visits
  san_length_df <- dplyr::filter(abpro_df, Behavior == "AN" & 
                               Modifier_1 == "Sanitizing")
  #Sum the duration column to get total duration AN
  san_duration <- sum(san_length_df$Duration_s)

  # Create dataframe to output
  data.frame(name = c("avg_feeding_vis_length", "avg_san_visit_length"), 
             calculation = c(feeding_duration/feeding_visit_number,
                             san_duration/san_visit_number))
}


