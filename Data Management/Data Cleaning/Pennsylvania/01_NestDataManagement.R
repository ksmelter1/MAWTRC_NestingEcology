#---
# title: Nest and Vegetation Survey Data Management
# author: "K. Smelter
# date: "`r format(Sys.time(), '%d %B, %Y')`"
# output:
#   html_document: 
#     toc: true
#---
#  
# **Purpose**: This script organizes nest and vegetation data for each of the three models

################################################################################
## Load Packages 

library(tidyverse)
library(sf)
library(mapview)
library(units)

################################################################################
## Load in Data 

# Load in raw nest csv file with data from 2022, 2023, and 2024
# This is just nests, no veg sampling or no nest values
nests.raw <- read_csv("Data Management/Csvs/Raw/Nests/nests.raw.pa.csv")
nests.raw

# Load in raw nest vegetation csv file with data from 2022, 2023, and 2024 surveys
# This file also contains original nest data 
# For nest fates we will have to merge the two 
veg.raw <- read_csv("Data Management/Csvs/Raw/VegetationSurveys_2024/nest_vegetation.csv")
veg.raw


################################################################################
## Clean Nest Vegetation Data from 2022 and 2023 

# Remove first row of data (NAs)
veg.raw <- veg.raw[-1,]
summary(veg.raw)

# Remove no nest values from nestid column 
veg_clean <- veg.raw %>% 
  dplyr::group_by(plottype)

# Define Case column for used and potential nests  
veg_clean$Case <- ifelse(veg_clean$plottype =="Nest", 1,0)


################################################################################
## Plot Data for Further Cleaning
 
 # Change coords to numeric
 veg_clean$lat_v <- as.numeric(veg_clean$lat_v)
 veg_clean$long_v <- as.numeric(veg_clean$long_v)
 
 # Filter out coordinate and NestID values with incorrect coords
 # Create Sf object
 veg_clean.sf <- veg_clean %>%
   dplyr::mutate("lat_v1" = lat_v) %>%
   dplyr::mutate("long_v1"= long_v) %>%
   dplyr::filter(!is.na(lat_v) & !is.na(long_v)) %>%
   dplyr::filter(lat_v >39) %>%
   dplyr::filter(lat_v <69) %>%
   dplyr::filter(long_v <48) %>%
   dplyr::filter(Case != "NA") %>%
                sf::st_as_sf(coords = c("long_v", "lat_v"), crs = 4326) %>%
                st_transform(5070) %>%
   dplyr::group_by(nestid_v)
 
 nestids_to_remove <- veg_clean.sf %>%
   group_by(nestid_v) %>%
   filter(n_distinct(Case) > 1) %>%  # Ensure both Case 0 and 1 are present
   summarise(
     geometry_0 = st_combine(geometry[Case == 0]),
     geometry_1 = st_combine(geometry[Case == 1]),
     .groups = "drop"
   ) %>%
   mutate(max_dist = st_distance(geometry_0, geometry_1, by_element = TRUE)) %>%
   filter(max_dist >= set_units(110, "m")) %>%
   pull(nestid_v)
 
 veg_clean.sf.filtered <- veg_clean.sf %>%
   filter(!nestid_v %in% nestids_to_remove)
 
 # Create the map by Case value
 mapview(veg_clean.sf.filtered, 
         zcol = "Case", 
         col.regions = c("blue", "red"), 
         legend = TRUE)
plot(st_geometry(veg_clean.sf.filtered))

# Remove nests that do not have any veg plots associated with them
# Also remove nests with outlier coordinates or with potential nests that didn't match spatially
veg_filtered <- veg_clean.sf.filtered %>%
  group_by(nestid_v) %>%
  filter(!(all(plottype == "Nest") & n_distinct(plottype) == 1)) %>%
  ungroup() %>%
   dplyr::filter(nestid_v != "8751_2023_1") %>%
   dplyr::filter(nestid_v != "8806_2023_1") %>%
   dplyr::filter(nestid_v != "8586_2023_1") %>%
   dplyr::filter(nestid_v != "8405_2022_1") %>%
   dplyr::filter(nestid_v != "8203_2022_1") %>%
   dplyr::filter(nestid_v != "8208_2022_1") %>%
   dplyr::filter(nestid_v != "8908_2023_1") %>%
   dplyr::filter(nestid_v != "9092_2023_2") %>%
   dplyr::filter(nestid_v != "9024_2023_3") %>%
   dplyr::filter(nestid_v != "9024_2024_1") %>%
   dplyr::filter(nestid_v != "9021_2023_1") 
   
# Create the map by Case value
mapview(veg_filtered, 
        zcol = "Case", 
        col.regions = c("blue", "red"), 
        legend = TRUE)
plot(st_geometry(veg_filtered))

# Filter cleaned veg df to only include results from the veg filtered object
veg_clean <- veg_clean.sf %>%
  arrange(nestid_v) %>%
  dplyr::filter(nestid_v %in% veg_filtered$nestid_v) 
head(veg_clean)

################################################################################
## Create Date Columns- Veg Data

# Build dataframe with information we need
veg_clean$CheckDate <- paste0(veg_clean$surveyyr,
                              veg_clean$surveymo,
                              veg_clean$surveyday) 
head(veg_clean$CheckDate)

# Correct formatting for SurveyMo and SurveyDay (ensuring 2 digits)
veg_clean$CheckDate <- paste0(veg_clean$surveyyr,
                              sprintf("%02d", veg_clean$surveymo),  
                              sprintf("%02d", veg_clean$surveyday)) 
head(veg_clean$CheckDate)

# Format timestamps correctly in year, month, day format and create start date column 
veg_clean$CheckDate<-as.Date(as.character(veg_clean$CheckDate),
                               format="%Y%m%d") 
head(veg_clean$CheckDate)
veg_clean

# Subset columns
veg_clean <- veg_clean %>%
  dplyr::select(-surveymo, -surveyday) %>%
  dplyr::rename("NestID" = nestid_v)
veg_clean

################################################################################
## Subset Nesting Data

# Create a dataframe and remove "No Nest" values, then replace 99s with NA
# 99s were a code for NAs in the database
nest_df <- nests.raw %>%
  dplyr::filter(nestfound == "Y") %>%
  dplyr::rename("NestID" = nestid) %>%
  dplyr::mutate(across(where(is.numeric), ~na_if(., 99))) %>%
  dplyr::mutate(across(where(is.character), ~na_if(., "99"))) %>%
  dplyr::filter(NestID %in% veg_clean$NestID)

# Merge cleaned Veg Data with Nest Data by keeping all columns in Veg Data
# Will remove unnecessary columns in the following pipeline
nests.veg <- right_join(veg_clean, nest_df)

# Arrange Dataframe by NestID and Case
# Remove columns not needed for analysis and BandIDs in which Tyler flagged
# Rename columns 
nests.veg <- nests.veg %>%
  dplyr::mutate(
    NestID = as.character(NestID),
    Case = as.character(Case)
  ) %>%
  arrange(NestID, Case) %>%
  dplyr::select(
    -checkyr, -checkmo, -checkday,
    -transmitterid_n, -landownership_n, -township_n, -county_v,
    -woodyht, -fernht, -grassforbht, -litterht, -boulderht,
    -lat_n, -long_n, -wmu_n, -county_n, -geometry,
    -unhatchids, -transmitterid_v, -landownership_v, -ageunhatched,
    -guardobject, -guardvo, -guardht, -secondaryguardobj, -averagemaxvo,
    -township_v, -nestnumber_v, -avgtotalcoverht, -avgtotalpercover
  ) %>%
  dplyr::rename(
    NestFate = nestfate,
    PlotType = plottype,
    WMU = wmu_v,
    BandID = bandid,
    SurveyYr = surveyyr,
    NestNumber = nestnumber,
    AvgVO = averagevo,
    PercWoody = percwoody,
    PercFern = percfern,
    PercGrassForb = percgrassforb,
    PercLitter = perclitter,
    PercBare = percbare,
    PercBoulder = percboulder,
    EggsUnhatched = eggsunhatched,
    EggsDestroyed = eggsdestroyed,
    EggsHatched = eggshatched,
    NestComments = nestcomments,
    VegComments = vegcomments,
    Lat = lat_v1,
    Long = long_v1,
    StemCount = stemcount,
    NestFound = nestfound
  ) 

# Remove hens that have NAs in EggsDestroyed, EggsHatched, and EggsUnhatched columns
# Create a column for nestfate where it is a binary outcome 0 or 1 (NestFate_Binary)
  nests.veg <- nests.veg %>%
    dplyr::filter(Case != "NA") %>%
  dplyr::filter(!BandID %in% c(
    "8619", "8811", "8969", "9069", "9074", "10012", "1056"
  )) %>%
  dplyr::mutate(NestFate_Binary = ifelse(NestFate == "Hatched", 1, 0)) %>%
  dplyr::select(NestID, Case, WMU, NestFate, Lat, Long, everything()) 
colnames(nests.veg)
summary(nests.veg)

# Create a clutch size column where for each row the eggsdestroyed, eggshatched, and eggsunhatched columns are summed
ClutchSize <- numeric(nrow(nests.veg))
for (i in 1:nrow(nests.veg)) {
  ClutchSize[i] <- sum(
    nests.veg$EggsHatched[i],
    nests.veg$EggsDestroyed[i],
    nests.veg$EggsUnhatched[i],
    na.rm = TRUE
  )
}
nests.veg$ClutchSize <- ClutchSize

# Filter out NestID values where there are potential nests but no used nests
# This would throw off our strata
nests.veg <- nests.veg %>%
  dplyr::group_by(NestID) %>%
  dplyr::filter(any(Case == "1")) %>%
  dplyr::ungroup()

# Create a nests object by filtering nests.veg to only include values where the case = 1
# Also include PlotType is equal to Nest to check for mismatch
nests <- nests.veg %>%
  dplyr::filter(Case == "1" & PlotType == "Nest")
table(nests$NestFate)
table(nests$NestFate_Binary)

# Make sure lengths match up
length(unique(nests$NestID))
length(unique(nests.veg$NestID))

# Final spatial checks
nests.veg.sf <- nests.veg %>%
  st_as_sf(coords = c("Lat", "Long"), crs = 4326) %>%
  st_transform(5070)
# Create the map by Case value
mapview(nests.veg.sf, 
        zcol = "Case", 
        col.regions = c("blue", "red"), 
        legend = TRUE)
plot(st_geometry(veg_filtered))


################################################################################
## Write Csvs

saveRDS(nests.veg, "nests.veg.RDS")
# Cleaned nest vegetation csv 
# Subset of nests that were marked as found
write_csv(nests.veg, "Data Management/Csvs/Processed/Nests/Vegetation Surveys/Pennsylvania/20250629_CleanedNestsVeg_2022_2023_2024.csv")

# Cleaned nests csv file
# Filtered out unknown fates
write.csv(nests, "Data Management/Csvs/Processed/Nests/Nests/Pennsylvania/20250629_CleanedNests_2022_2023_2024.csv")

################################################################################
################################################################################
