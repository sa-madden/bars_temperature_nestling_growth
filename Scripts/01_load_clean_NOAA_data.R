####### Purpose: clean up NOAA weather dataset for Boulder, CO
####### By: Sage Madden
####### Created: 3/7/2022
####### Last modified: 3/7/2022

# Load packages
library(tidyverse)

# Manually deleted some text containing descriptive information from this 
# text file. 
# Original file left in the same folder
boco <- read.table("Data/boulderdaily.complete_ed.txt", 
                   header = FALSE, sep = "", dec = ".")
colnames(boco) <- c("year", "month", "day", "max_temp_f", "min_temp_f", 
                    "precip_in", "snowfall_in", "snowdepth_in")
# Important info: missing values indicated by -998.0
# Trace values indicated by -999.0

# Select only columns I want and filter to only year and months I want
str(boco)
boco <- boco %>% select(year, month, day, max_temp_f, 
                        min_temp_f, precip_in) %>% 
  filter(year == 2021 & month < 8 & month > 5)

# Change those -999.0s to NAs because I said so 
boco[boco == -999.0] <- NA

write.csv(boco, file = "Output/noaa_data_cleaned.csv")

# LETS LOOK AT THIS SHIT
boco$fmonth <- as.factor(boco$month)

ggplot(boco, aes(x = day, y = max_temp_f,
                 col = fmonth)) +
  facet_wrap(~ fmonth, 1) + 
  geom_point() + 
  theme_classic() +
  labs(x = "Date", y = "Max temp F", col = "Month")

ggplot(boco, aes(x = day, y = min_temp_f,
                 col = fmonth)) +
  facet_wrap(~ fmonth, 1) + 
  geom_point() + 
  theme_classic() +
  labs(x = "Date", y = "Min temp F", col = "Month")

ggplot(boco, aes(x = day, y = precip_in,
                 col = fmonth)) +
  facet_wrap(~ fmonth, 1) + 
  geom_point() + 
  theme_classic() +
  labs(x = "Date", y = "Precipitation (in)", col = "Month")






