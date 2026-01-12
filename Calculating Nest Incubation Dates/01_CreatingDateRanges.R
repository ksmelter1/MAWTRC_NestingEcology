#---
# title: Creating Date Ranges for ACC Downloads
# author: "V. Winter, K. Smelter
# date: "`r format(Sys.time(), '%d %B, %Y')`"
# output:
#   html_document: 
#     toc: true
#---
#  
# **Purpose**: This script creates specified date ranges to be used in 02DownloadingACCDataByBirdsandDates

################################################################################
## Load Packages and Data

# Initialize workspace
rm(list = ls())
gc()

library(dplyr)
library(lubridate)

# Load in data ---
nsts <- readRDS("Nest_data/Nesting_KJS/NJ/20260109_CleanedNests_2023_2024_2025_NJ.rds")
colnames(nsts)

################################################################################
## Data Prep

# Clean data -----
# We want birds with (nestfound=Y)
# Remove all nests with a fate of unconfirmed
# Obtain data from nests only, not available nests
# Year constraint to 2025 since data has been downloaded for prior years
nsts_clean <- nsts %>% 
  # Filter out where nests were found
  dplyr::filter(NestFound == "Y") %>%
  dplyr::filter(NestFate!= "Unconfirmed") %>%
  dplyr::filter(PlotType == "Nest") %>%
  dplyr::filter(CheckYr == "2025")

# Great, now we want to find a lists of bird IDs and start date
ids <- unique(nsts_clean$NestID)
ids
dates <- nsts_clean$CheckDate
dates
str(dates)
str(nsts_clean$CheckDate)

# Get daterange 
# Set range to 50 days prior to the CheckDate  and 3 days after
date_id <- nsts_clean %>%
  dplyr::select(NestID, BandID, CheckDate) %>%
  dplyr::mutate(
    CheckDate = as.Date(CheckDate, format = "%m/%d/%Y"),
    startsearchdate = CheckDate - days(50),
    endsearchdate = CheckDate + days(3),
    CheckDate = format(CheckDate, "%m/%d/%Y"),
    startsearchdate = format(startsearchdate, "%m/%d/%Y"),
    endsearchdate = format(endsearchdate, "%m/%d/%Y")
  )  

# Set date formatting to month day year
endsearchdate<-as.Date(date_id$endsearchdate, "%m/%d/%Y")
startsearchdate<-as.Date(date_id$startsearchdate, "%m/%d/%Y")

################################################################################
## Output Data

#write.csv(date_id, "Nest_data/Nesting_KJS/20250125_missingbirds_KS.csv")
saveRDS(date_id, "Nest_data/Nesting_KJS/MD/20260109_NestingRanges_2025_NJ.rds") # Switch for each state

################################################################################
###############################################################################X