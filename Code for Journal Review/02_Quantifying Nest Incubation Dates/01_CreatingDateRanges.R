#---
# title: Creating Date Ranges for ACC Downloads
# author: "V. Winter, K. Smelter
# date: "`r format(Sys.time(), '%d %B, %Y')`"
# output:
#   html_document: 
#     toc: true
#---
#  
# **Purpose**: This script creates specified date ranges to be used in 02_DownloadingACCDataByBirdsandDates

# Initialize workspace
rm(list = ls())
gc()

################################################################################
## Load Packages and Data

library(dplyr)
library(lubridate)

# Load in RDS file
nsts <- readRDS("Nest_data/Nesting_KJS/NJ/20260109_CleanedNests_2024_NJ.rds")
colnames(nsts)

################################################################################
## Data Prep

# We want birds with (nestfound=Y)
# Remove all nests with a fate of unconfirmed
# Obtain data from nests only, not available nests
# Do not include 2025 data
nsts_clean <- nsts %>% 
  # Ensure all nests were found (No Nests not included)
  dplyr::filter(NestFound == "Y") %>%
  dplyr::filter(NestFate!= "Unconfirmed") %>%
  dplyr::filter(PlotType == "Nest") %>%
  dplyr::filter(CheckYr != "2025")

# Find a lists of bird IDs and start date
ids <- unique(nsts_clean$NestID)
ids

# Get dates nests were checked
# Check structure
dates <- nsts_clean$CheckDate
dates
str(dates)
str(nsts_clean$CheckDate)

# Get date range 
# Set range to 50 days prior to the CheckDate  and 3 days after
# Format data mdy to match access databases for each state
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

# Output RDS files of birds and date ranges
saveRDS(date_id, "Nest_data/Nesting_KJS/NJ/20260109_NestingRanges_2025_NJ.rds") # Switch for each state

################################################################################
###############################################################################X