#'---
#' title: Nest and Vegetation Survey Data Management
#' author: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script shows the amount of hens that nested across years and organizes nest data and shapefiles for analysis

#####################
## Load Packages 

require(tidyverse)
require(sf)

##############################
## Load in Data 

#' Load in raw nest csv file with data from 2022-2023
#' This is just nests, no veg sampling or no nest values
nests.raw <- read_csv("Data Management/Csvs/Raw/20250120_nests_raw.csv")
nests.raw

#' Load in raw nest vegetation csv file with data from 2022 surveys
#' This file also contains original nest data 
#' For nest fates we will have to merge the two 
veg.raw.2022 <- read_csv("Data Management/Csvs/Raw/VegPlotTable2022.csv") 
veg.raw.2022

#' Load in raw nest vegetation csv file with data from 2023 surveys
#' This file also contains original nest data 
#' For nest fates we will have to merge the two 
veg.raw.2023 <- read_csv("Data Management/Csvs/Raw/VegPlotTable2023.csv")
veg.raw.2023

#######################################################
## Check Nesting Attempts per Each Hen Across Years 

#' If all values within a row are missing remove the row
nests.raw <- nests.raw %>%
  filter(!apply(., 1, function(x) all(is.na(x))))

#' Create a dataframe and remove "No Nest" values as nest success models are reliant upon the nest being found
nest_df<- nests.raw %>%
  dplyr::filter(NestBowlFound == "Y")

#' Create unique identifier column 
nest_df$NestID <- with(nest_df, paste0(BandID,"_",NestNumber,"_",CheckYr))
glimpse(nest_df)

#' Subset nest dataframe by grouping it by nest check year, then keeping columns BandID and NestID
nest.year <- nest_df %>%
  dplyr::group_by(CheckYr)%>%
  dplyr::select(CheckYr, BandID,Lat, Long, WMU, CheckDay, NestID, NestFate)

#' Using dplyr to filter birds that nested in both 2022 and 2023
nests.both <- nest.year %>% dplyr::select(BandID, CheckYr)%>%
  dplyr::mutate(x=1) %>%
  distinct() %>%
  tidyr::pivot_wider( names_from = CheckYr, values_from = x, 
                      values_fill = 0 ) %>%
  dplyr::filter(`2022`==1 & `2023`==1) %>%
  pull(BandID) 

#' Create dataframe 
nests.year.both <- dplyr::filter(nest.year,bandid %in% nests.both)
unique(nests.year.both$bandid)

## 16 hens nested across 2022 and 2023
## We determined this was not as large of a sample as we needed to look into the dominant hen hypothesis
  
#####################################################
## Clean Nest Vegetation Data from 2022 and 2023 

#' Remove all NA values across the rows of 2 dataframes
veg.raw.2022 <- veg.raw.2022 %>%
  filter(!apply(., 1, function(x) all(is.na(x))))

veg.raw.2023 <- veg.raw.2023 %>%
  filter(!apply(., 1, function(x) all(is.na(x))))

#' Change all values of specified parameters to character strings so we can bind dfs
veg.raw.2022 <- veg.raw.2022 %>%
  dplyr::mutate(across(c(Lat, Long, PercWoody, PercFern, PercGrassForb, PercLitter, 
                  PercBare, PercBoulder, HtWoody, HtFern, HtGrassForb, 
                  HtLitter, HtBoulder, StemCount), as.character))

veg.raw.2023 <- veg.raw.2023 %>%
  dplyr::mutate(across(c(Lat, Long, PercWoody, PercFern, PercGrassForb, PercLitter, 
                  PercBare, PercBoulder, HtWoody, HtFern, HtGrassForb, 
                  HtLitter, HtBoulder, StemCount), as.character))

#' Now combine the data frames
veg_df <- bind_rows(veg.raw.2022, veg.raw.2023)


#' Fill empty cells with 99 as a proxy for NA
veg_df$VegComments[veg_df$VegComments== ""] <- as.character("99")
veg_df$WoodyType1[veg_df$WoodyType1== ""] <- as.character("99")
veg_df$WoodyType2[veg_df$WoodyType2== ""] <- as.character("99")
veg_df$GuardVO[veg_df$GuardVO== ""] <- as.character("99")

#' Remove no nest values from nestid column 
veg_clean <- veg_df %>% 
  dplyr::group_by(PlotType)
  
  veg_clean$Case <- ifelse(veg_clean$PlotType =="North"|
                             veg_clean$PlotType == "South"|
                             veg_clean$PlotType == "East"|
                             veg_clean$PlotType == "West", "0","1")

#' Create unique identifier column 
veg_clean$NestID <- with(veg_clean, paste0(BandID,"_",NestNumber,"_",SurveyYr))
glimpse(veg_clean)


####################################
## Plot Data for Further Cleaning
 
#' Fix coding inconsistency
 veg_clean$Lat <- iconv(veg_clean$Lat, from = "UTF-8", to = "ASCII", sub = "")
 veg_clean$Long <- iconv(veg_clean$Long, from = "UTF-8", to = "ASCII", sub = "")
 
 #' Change coords to numeric
 veg_clean$Lat <- as.numeric(veg_clean$Lat)
 veg_clean$Long <- as.numeric(veg_clean$Long)
 
 #' Filter out coordinate and NestID values with incorrect coords
 #' Create Sf object
 veg_clean.sf <- veg_clean %>%
   dplyr::mutate("Lat1" = Lat) %>%
   dplyr::mutate("Long1"= Long) %>%
   dplyr::filter(!is.na(Lat) & !is.na(Long)) %>%
   dplyr::filter(Lat >39) %>%
   dplyr::filter(Long <48) %>%
   dplyr::filter(NestID != "8173_1_2022") %>%
                sf::st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>%
                st_transform(5070) 
mapview(veg_clean.sf)
plot(st_geometry(veg_clean.sf))


#' Filter cleaned veg df to only include results from the veg sf object
veg_clean <- veg_clean.sf %>%
  dplyr::filter(NestID %in% veg_clean.sf$NestID) %>%
  dplyr::rename("Lat" = Lat1) %>%
  dplyr::rename("Long" = Long1) %>%
  dplyr::select(-geometry, -...36, -`Surveyed (Y/N)`,-`Surveyed (Y/N)`, -...37)

##########################
## Create Date Columns 

#' Build dataframe with information we need
veg_clean$checkdate <- paste0(veg_clean$SurveyYr,
                              veg_clean$SurveyMo,
                              veg_clean$SurveyDay) 
head(veg_clean$checkdate)

# Correct formatting for SurveyMo and SurveyDay (ensuring 2 digits)
veg_clean$CheckDate <- paste0(veg_clean$SurveyYr,
                              sprintf("%02d", veg_clean$SurveyMo),  # Format month as two digits
                              sprintf("%02d", veg_clean$SurveyDay)) # Format day as two digits
head(veg_clean$CheckDate)

#' Format timestamps correctly in year, month, day format and create start date column 
veg_clean$CheckDate<-as.Date(as.character(veg_clean$CheckDate),
                               format="%Y%m%d")
head(veg_clean$CheckDate)

#' Write cleaned vegetation csv
write.csv(veg_clean, "nests.veg_22_23.csv")

###################################
## Finalize Nesting Data

#' Filter veg sheet to only include used nests
nests_clean <- nest_df %>% 
dplyr::filter(NestID %in% veg_clean$NestID)

################################################################################