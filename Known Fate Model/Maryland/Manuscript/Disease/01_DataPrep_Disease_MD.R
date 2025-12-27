#'---
#' title: Daily Nest Survival Modeling of Female Wild Turkeys in Maryland
#' authors: K. Smelter
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script uses simulated data to fit a Bayesian known fate model for female wild turkeys in our study
#' **Key Changes**: This script uses the cloglog link instead of the logit link for modeling daily nest survival
#' **Last Updated**: 12/27/25

################################################################################
## Load Packages

library(tidyverse)
library(terra)
library(mapview)
library(sf)
library(stringr)
library(daymetr)

################################################################################
## Data Prep- Nest-Level Covs

# Nest start and end date csv
nests <- read_csv("Data Management/Csvs/Processed/Incubation Dates/Maryland/20250717_NestAttempts_allbirds_MD_Ready.csv")
nests 

# Nest veg csv
nests.veg <- read_csv("Data Management/Csvs/Processed/Nests/Nests/Maryland/20250219_CleanedNests_2022_2023_MD.csv") #%>%
  #dplyr::filter(CheckYr == "2023")

# Filter nests.veg to only include observations that have the same nestid as nests csv
nests.veg.filtered <- nests.veg %>%
  dplyr::filter(NestID%in% nests$NestID) 

# Merge the filtered nests.veg with nests by nestid
nests <- dplyr::right_join(nests, nests.veg.filtered, by = "NestID") %>%
  dplyr::select(-CheckDate.x) %>%
  dplyr::rename("NestFate" = NestFate.y)

# Rename and consolidate columns 
nests <- nests %>%
  dplyr::select(NestID, 
                BandID, 
                startI, 
                endI,
                NestFate, 
                Lat,
                Long) 

# Switch coding to UTF-8
nests <- nests %>% 
  dplyr::mutate(across(everything(), ~ iconv(., to = "UTF-8")))


################################################################################
## Data Prep- Individual Covariates

# Read in captures csv
captures <- read_csv("Data Management/Csvs/Raw/Captures/captures_md.csv")
captures
 
# Filter data to include only hens 
captures <- captures %>%
   dplyr::filter(sex == "F") %>%
   dplyr::select(bandid, sex, age, captyr) %>%
   dplyr::rename("BandID" = bandid)
  captures

# Merge columns 
 nests.scaled <- merge(captures, nests, by = "BandID") 

 # (?<=_): Only match is there is an underscore immediately before the number we are trying to extract
 # (\\d{4}): Match exactly 4 digits
 # (?=_) : Only match if the four digits are followed 
 # Create NestYr column 
nests.scaled <- nests.scaled %>%
  dplyr::mutate(NestYr = str_extract(NestID, "(?<=_)(\\d{4})(?=_)")) 

# Convert to numeric 
nests.scaled$NestYr <- as.numeric(nests.scaled$NestYr)
nests.scaled$captyr <- as.numeric(nests.scaled$captyr)

# Create a years since capture column
nests.scaled <- nests.scaled %>%
 dplyr::mutate(yrsincecap = NestYr-captyr)
 glimpse(nests.scaled)

# Assign Adult as the reference level
nests.scaled$age <- ifelse(nests.scaled$age == "J", 1, 
                             ifelse(nests.scaled$age == "A", 0, NA))

# Dealing with scaling age ad hoc
# If the bird is an adult and the years since capture is >1 assign it as an adult
# If not keep the age as juvenile because turkeys will nest the first year as a juvenile
nests.scaled$age <- ifelse(nests.scaled$age == 1 & nests.scaled$yrsincecap >= 1, 0, nests.scaled$age)
table(nests.scaled$age)
glimpse(nests.scaled)

################################################################################
 ## Data Prep- Disease Covariates

# Read in raw virus csv from Maryland
virus.2024 <- read_csv("Data Management/Csvs/Raw/Disease/LPDV_REV/Maryland/MD_virus.csv")
virus.2024

# Organizing data from SCWDS formatting 
virus.2024 <- virus.2024 %>%
  slice(-89:-94) %>%
  dplyr::rename("BandID" = `State ID`) %>%
  dplyr::select(BandID, LPDV, REV) %>%
  dplyr::mutate(LPDV = ifelse(LPDV == "+", 1, 0)) %>%
  dplyr::mutate(REV = ifelse(REV == "+", 1, 0))
  virus.2024

# Read in 2023 SCWDS data
virus.2023 <- read.csv("Data Management/Csvs/Raw/Disease/LPDV_REV/Maryland/MD_virus_2023.csv")
virus.2023

# Organizing data from SCWDS formatting 
virus.2023 <- virus.2023 %>%
  dplyr::select(-c(X.7:X.16)) %>%
  dplyr::select(-c(X:X.6)) %>%
  slice(-49:-1001) %>%
  dplyr::rename("BandID" = ID..) %>%
  dplyr::select(BandID, LPDV, REV) %>%
  dplyr::mutate(LPDV = ifelse(LPDV == "POS", 1, 0)) %>%
  dplyr::mutate(REV = ifelse(REV == "POS", 1, 0))

# Combine datasets
virus <- rbind(virus.2023, virus.2024)

# Filter data to contain only individuals within the nests.scaled df
virus.filter <- virus %>%
dplyr::filter(BandID %in% nests.scaled$BandID)

# Join virus data by all observations within nests.scaled
nests.scaled <- right_join(virus.filter, nests.scaled) %>%
  tidyr::drop_na()

# Assign 1 if the NestFate "Hatched", otherwise assign 0
nests.scaled$NestFate <- ifelse(nests.scaled$NestFate == "Hatched", 1, 0)

################################################################################
###############################################################################X
