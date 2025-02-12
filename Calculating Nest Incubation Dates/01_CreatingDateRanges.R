
#'---
#' title: Calculating Start and End Incubation Dates
#' author: K. Smelter, F. Buderman, V. Winter, K. Lamp
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output: 
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script organizes the nesting dates for further processing 
#' **Last Updated**: 1/31/2025


################################################################################
## Download ACC Data from Specified Date Ranges 

rm(list = ls())
gc()

#' load in packages
library(dplyr)
library(lubridate)


#' Load in data ---
nsts <- read.csv("Nest_data/Nesting_KJS/20250121_CleanedNests_2022_2023_KJS.csv")
colnames(nsts)

#' Clean data -----
nsts_clean <- nsts %>% 
  dplyr::filter(NestBowlFound == "Y") %>%
  dplyr::filter(NestFate!= "Unknown") %>%
  dplyr::filter(NestFate!= "Unconfirmed") %>%
  dplyr::select(-X)


ids <- unique(nsts_clean$NestID)
ids
dates <- nsts_clean$CheckDate
dates
str(dates)
str(nsts_clean$CheckDate)

date_id <- nsts_clean %>%
  select(bandid, CheckDate) %>%
  mutate(
    CheckDate = as.Date(CheckDate),  
    startsearchdate = CheckDate - days(50),  
    endsearchdate = CheckDate + days(3),  
    CheckDate = format(CheckDate, "%m/%d/%Y"),  
    endsearchdate = format(endsearchdate, "%m/%d/%Y"),
    startsearchdate = format(startsearchdate, "%m/%d/%Y"))  

endsearchdate<-as.Date(date_id$endsearchdate, "%m/%d/%Y")
startsearchdate<-as.Date(date_id$startsearchdate, "%m/%d/%Y")

write.csv(date_id, "Nest_data/Nesting_KJS/20250125_missingbirds_KS.csv")
saveRDS(date_id, "Nest_data/20250103_NestingDates_KS.rds")