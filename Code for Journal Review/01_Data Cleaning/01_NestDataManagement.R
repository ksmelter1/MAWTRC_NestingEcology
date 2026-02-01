#---
#' title: Nest and Vegetation Survey Data Management for PA, MD, and NJ
#' author: "K. Smelter"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'  html_document: 
#'     toc: true
#---
#  
#' **Purpose**: This script organizes nest and vegetation data for Pennsylvania, Maryland, and New Jersey
#' **Last Updated**: 22 January 2026

################################################################################
## Load Packages 

library(tidyverse)
library(sf)
library(mapview)
library(units)


###################
## Pennsylvania ##
###################

################################################################################
## Load in Data 

# Load in raw nest csv file with data from 2022, 2023, 2024, and 2025
# This is just nests, no veg sampling or no nest values
nests.raw <- read_csv("Data Management/Csvs/Raw/Nests/nests.raw.pa.csv")
nests.raw

# 2025 data is not in this manuscript, may come in a revision
nests.raw <- nests.raw %>% dplyr::filter(checkyr != "2025")

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
 
 # Filter out coordinate and NestID values with incorrect coordinates
 # Create Sf object to visualize data
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
 
 # Pull observations of used and potential nests that aren't within 110 m of each other
 # Distance threshold was put in place to match spatial methods for Maryland and New Jersey (Potential nests 100m from used)
 # Ensure both used and available nests are present
 nestids_to_remove <- veg_clean.sf %>%
   dplyr::group_by(nestid_v) %>%
   summarise(
     geometry_0 = st_combine(geometry[Case == 0]),
     geometry_1 = st_combine(geometry[Case == 1]),
     .groups = "drop"
   ) %>%
  dplyr:: mutate(max_dist = st_distance(geometry_0, geometry_1, by_element = TRUE)) %>%
   dplyr::filter(max_dist >= set_units(110, "m")) %>%    # 110 m distance threshold 
   pull(nestid_v)
 
 # Subset sf object to only the removed nests
 removed_nests.sf <- veg_clean.sf %>%
   dplyr::filter(nestid_v %in% nestids_to_remove)
 
 # Map the flagged nests
 mapview(removed_nests.sf,
         zcol = "Case",
         col.regions = c("blue", "red"),
         legend = TRUE)
 
 # Remove nests with coordinate errors
 veg_clean.sf.filtered <- veg_clean.sf %>%
   dplyr::filter(!nestid_v %in% nestids_to_remove)
 
 # Change coding to UTF 8 to create maps
 veg_clean.sf.filtered <- veg_clean.sf.filtered %>%
   dplyr::mutate(across(where(is.character), ~ iconv(.x, from = "", to = "UTF-8")))
 mapview(veg_clean.sf.filtered)

# Remove nests that do not have any veg plots associated with them
# Also remove nests with outlier coordinates or with potential nests that didn't match spatially
# Some nests weren't flagged
veg_filtered <- veg_clean.sf.filtered %>%
  dplyr::group_by(nestid_v) %>%
 dplyr::filter(!(all(plottype == "Nest") & n_distinct(plottype) == 1)) %>%
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
# Check date formatting
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
  dplyr::filter(nestfound == "Y") %>%  # Remove No Nest Values (Nests that weren't found)
  dplyr::rename("NestID" = nestid) %>% # Rename NestID column
  dplyr::mutate(across(where(is.numeric), ~na_if(., 99))) %>% # Change 99s to NAs
  dplyr::mutate(across(where(is.character), ~na_if(., "99"))) %>% # Change to Character
  dplyr::filter(NestID %in% veg_clean$NestID) # Keep only NestID observations that are present within vegetation data

# Merge cleaned Veg Data with Nest Data by keeping all columns in Veg Data
# Will remove unnecessary columns in the following pipeline
nests.veg <- right_join(veg_clean, nest_df)

# Arrange Dataframe by NestID and Case
# Remove columns not needed for analysis 
nests.veg <- nests.veg %>%
  dplyr::mutate(
    NestID = as.character(NestID), # Convert columns to characters
    Case = as.character(Case)
  ) %>%
  arrange(NestID, Case) %>%
  dplyr::select(                   # Remove unnecessary columns (I have the nest check information elsewhere)
    -checkyr, -checkmo, -checkday,
    -transmitterid_n, -landownership_n, -township_n, -county_v,
    -woodyht, -fernht, -grassforbht, -litterht, -boulderht,
    -lat_n, -long_n, -wmu_n, -county_n, -geometry,
    -unhatchids, -transmitterid_v, -landownership_v, -ageunhatched,
    -guardobject, -guardvo, -guardht, -secondaryguardobj, -averagemaxvo,
    -township_v, -nestnumber_v, -avgtotalcoverht, -avgtotalpercover
  ) %>%                           
  dplyr::rename(                    # Rename columns        
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

# Create a column for nestfate where it is a binary outcome 0 or 1 (NestFate_Binary)
# Hens below were flagged by PSU so they were removed from the analysis
# Nests that were unknown failed from 2022-2024, will adjust for 2025
  nests.veg <- nests.veg %>%
  dplyr::filter(!BandID %in% c(
    "8619", "8811", "8969", "9069", "9074", "10012", "1056"    # Hens flagged by PSU
  )) %>%
  dplyr::mutate(NestFate_Binary = ifelse(NestFate == "Hatched", 1, 0)) %>%  # Binary Nest Fate column 
  dplyr::select(NestID, Case, WMU, NestFate, Lat, Long, everything())  # Select columns 

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

# Apply counts from clutch size column to dataframe
# Works because dataframes match up in length
nests.veg$ClutchSize <- ClutchSize

# Filter out NestID values where there are potential nests but no used nests
# This would throw off our strata
# We want the same nests in each analysis 
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

# Change coding to UTF 8 to create maps
nests.veg <- nests.veg%>%
  dplyr::mutate(across(where(is.character), ~ iconv(.x, from = "", to = "UTF-8")))

# Final spatial checks
nests.veg.sf <- nests.veg %>%
  st_as_sf(coords = c("Lat", "Long"), crs = 4326) %>%
  st_transform(5070)

# Create the map by Case value
mapview(nests.veg.sf, 
        zcol = "Case", 
        col.regions = c("blue", "red"), 
        legend = TRUE)


################################################################################
## Output Data

# Cleaned nest vegetation csv 
# Subset of nests that were marked as found
# write_csv(nests.veg, "Data Management/Csvs/Processed/Nests/Vegetation Surveys/Pennsylvania/20250629_CleanedNestsVeg_2022_2023_2024.csv")

# Cleaned nests csv file
# write.csv(nests, "Data Management/Csvs/Processed/Nests/Nests/Pennsylvania/20250629_CleanedNests_2022_2023_2024.csv")

# RDS of vegetation Data
# saveRDS(
#   nests.veg,
#   "Data Management/Data Management/Csvs/Processed/Nests/Vegetation Surveys/Pennsylvania/20250629_CleanedNestsVeg_2022_2023_2024.rds"
# )
################################################################################


##############
## Maryland ##
##############

################################################################################
## Read in Data

# Load raw nest csv file
nests.raw <- read.csv("Data Management/Csvs/Raw/Nests/nests.raw.md.csv")

################################################################################
## Data Prep

# Remove 2025 data
nests.raw <- nests.raw %>% dplyr::filter(checkyr!= "2025")

# Rename columns 
nests.raw.md <- nests.raw %>%
  dplyr::rename(
    BandID = bandid,
    CheckYr = checkyr,
    CheckDay = checkday,
    CheckMo = checkmo,
    NestID = nestid,
    NestFound = nestfound,
    NestNumber = nestnumber,
    WMU = wmu_n,
    County = county_n,
    Township = township_n,
    Landownership = landownership_n,
    EggsHatched = eggshatched,
    EggsUnhatched = eggsunhatched,
    EggsDestroyed = eggsdestroyed,
    NestFate = nestfate,
    NestSubFate = nestsubfate1,
    Lat = lat_n,
    Long = long_n
  )

# Remove "No Nests"
# These are nests that weren't found
nest_df <- nests.raw.md %>%
  dplyr::filter(NestFound == "Y") 

# Create a checkdate column
# Format a year, month, date
nest_df <- nest_df %>%
  dplyr::mutate(
    CheckDate = as.Date(
      paste0(
        CheckYr,
        sprintf("%02d", CheckMo),
        sprintf("%02d", CheckDay)
      ),
      format = "%Y%m%d"
    )
  )

# Remove unnecessary columns 
# Create an indicator for Nest under PlotType since all of these are used nests
# Will need the Nest indicator down the road when available nests are created
nests_clean <- nest_df %>%
  dplyr::select(-CheckMo, -CheckDay) %>%
  dplyr::mutate(PlotType = "Nest")

# Change coding to UTF-8
nests_clean <- nests_clean %>%
  dplyr::mutate(
    across(
      where(is.character),
      ~ iconv(.x, from = "latin1", to = "UTF-8", sub = "")
    )
  )

# Create secondary Long and Lat columns that won't be overwritten with geometry
nests_clean <- nests_clean %>%
  dplyr::mutate("Long1" = Long) %>%
  dplyr::mutate("Lat1" = Lat)

# Create spatial object
# Transform from WGS 84 to Albers
nests.sf <- st_as_sf(
  nests_clean,
  coords = c("Long1", "Lat1"),
  crs = 4326
) %>%
  st_transform(5070)

# Interactive map to inspect geometry locations
mapview(nests.sf)

# Remove nests with coordinate errors
# Removed via visual inspection of plotted data
nests.sf <- nests.sf %>%
  dplyr::filter(NestID != "1193_2023_1") %>%
  dplyr::filter(NestID != "1214_2024_2") %>%
  dplyr::filter(NestID != "1244_2025_1")

# Check again 
mapview(nests.sf)

################################################################################
## Generate Offset Points (100 m in each cardinal direction)

# Extract nest coordinates
coords <- st_coordinates(nests.sf)

# Function to generate available nests 100m away in each of the 4 cardinal directions
generate_points <- function(x, y, distance = 100) {
  tibble(
    Long = c(x, x, x + distance, x - distance),
    Lat  = c(y + distance, y - distance, y, y),
    PlotType = c("North", "South", "East", "West")
  )
}

valid_nest_ids <- nests.sf$NestID
new_points_df <- map_dfr(seq_along(valid_nest_ids), function(i) {
  generate_points(coords[i,1], coords[i,2]) %>%
    dplyr::mutate(NestID = valid_nest_ids[i])
})

# Create spatial object for coordinate extraction 
new_points_sf <- st_as_sf(
  new_points_df,
  coords = c("Long", "Lat"),
  crs = 5070
)

# Extract coordinates in meters
new.points <- new_points_sf %>%
  mutate(
    Long = st_coordinates(geometry)[,1],
    Lat  = st_coordinates(geometry)[,2]
  ) %>%
  st_drop_geometry()

# Create spatial object for final plot to make sure coords plot correct (Line 180)
new.points <- new.points %>%
  dplyr::mutate("Long1" = Long) %>%
  dplyr::mutate("Lat1" = Lat) %>%
  st_as_sf(coords = c("Long1", "Lat1"), crs = 5070)

# Check
mapview(new.points)

################################################################################
## Bring Data Together

# Combine used and available nests into a df
nests_clean <- bind_rows(
  nests.sf,
  new.points
) %>%
  dplyr::filter()

# If the PlotType is Nest it was a used nest
# If the PlotType isn't it is an available nest (100m away)
nests_clean <- nests_clean %>%
  dplyr::mutate(Case = if_else(PlotType == "Nest", 1, 0))

# Fix character encoding
nests_cleaned_md <- nests_clean %>%
  dplyr::mutate(
    across(
      where(is.character),
      ~ iconv(.x, from = "", to = "UTF-8", sub = "")
    )
  )

# Arrange data by NestID
nests_cleaned_md <- nests_cleaned_md %>%
  arrange(NestID)

# Check
mapview(nests_cleaned_md)

# Convert geometry to Lat and Long
# Put coordinates in the same crs for further analysis
nests_cleaned_md <- nests_cleaned_md %>%
  dplyr::mutate(
    Long = st_coordinates(geometry)[,1],  
    Lat  = st_coordinates(geometry)[,2]   
  )

################################################################################
## Finalize Data 

# Fill NAs in available nests with info from the used nest within each NestID
# This just completes the strata with all information 
nests_cleaned_md <- nests_cleaned_md %>%
  dplyr::group_by(NestID) %>%
  dplyr::mutate(across(
    .cols = c(BandID, CheckYr, NestFound, transmitterid_n, NestNumber, WMU, County, Township,
              Landownership, EggsHatched, EggsUnhatched, EggsDestroyed, NestFate,
              NestSubFate, nestsubfate2, nestcomments, CheckDate),
    .fns = ~ if_else(is.na(.x) & any(Case == 1), .x[Case == 1][1], .x)
  )) %>%
  ungroup() %>%
  st_drop_geometry()

################################################################################
## Output Data

# No vegetation data so just output nesting data
# Save RDS file
# saveRDS(
#   nests_cleaned_md,
#   "Data Management/Csvs/Processed/Nests/Nests/Maryland/20260109_CleanedNests_2023_2024_MD.rds"
# )

# # Write csv
# write.csv(nests_clean, "Data Management/Csvs/Processed/Nests/Nests/20260109_CleanedNests_2024_MD.csv")

################################################################################


#################
## New Jersey ##
#################


################################################################################
## Load Data

# Load raw nest csv file
nests.raw <- read_csv("Data Management/Csvs/Raw/Nests/nests.raw.nj.csv")
nests.raw 

################################################################################
## Data Prep

# Remove 2025 data
nests.raw <- nests.raw %>% dplyr::filter(checkyr!= "2025")

# Rename columns 
nests.raw.nj <- nests.raw %>%
  dplyr::rename(
    BandID = bandid,
    CheckYr = checkyr,
    CheckDay = checkday,
    CheckMo = checkmo,
    NestID = nestid,
    NestFound = nestfound,
    NestNumber = nestnumber,
    WMU = wmu_n,
    County = county_n,
    Township = township_n,
    Landownership = landownership_n,
    EggsHatched = eggshatched,
    EggsUnhatched = eggsunhatched,
    EggsDestroyed = eggsdestroyed,
    NestFate = nestfate,
    NestSubFate = nestsubfate1,
    Lat = lat_n,
    Long = long_n
  )


# Remove nests that weren't found ("No Nests")
nest_df <- nests.raw.nj %>%
  dplyr::filter(NestFound == "Y") 

# Create a checkdate column
# Format a year, month, date
nest_df <- nest_df %>%
  dplyr::mutate(
    CheckDate = as.Date(
      paste0(
        CheckYr,
        sprintf("%02d", CheckMo),
        sprintf("%02d", CheckDay)
      ),
      format = "%Y%m%d"
    )
  )

# Remove unnecessary columns 
# Create an indicator for Nest under PlotType since all of these are used nests
# Will need the Nest indicator down the road when available nests are created
nests_clean <- nest_df %>%
  dplyr::select(-CheckMo, -CheckDay) %>%
  dplyr::mutate(PlotType = "Nest")

# Change coding to UTF-8
nests_clean <- nests_clean %>%
  dplyr::mutate(
    across(
      where(is.character),
      ~ iconv(.x, from = "latin1", to = "UTF-8", sub = "")
    )
  )

# Create secondary Long and Lat columns that won't be overwritten with geometry
nests_clean <- nests_clean %>%
  dplyr::mutate("Long1" = Long) %>%
  dplyr::mutate("Lat1" = Lat)

# Create spatial object
# Transform from WGS 84 to Albers
nests.sf <- st_as_sf(
  nests_clean,
  coords = c("Long1", "Lat1"),
  crs = 4326
) %>%
  st_transform(5070)

# Plot geometry locations
mapview(nests.sf)

################################################################################
## Generate Offset Points (100 m in each cardinal direction)

# Extract nest coordinates
coords <- st_coordinates(nests.sf)

# Function to generate available nests 100m away in each of the 4 cardinal directions
generate_points <- function(x, y, distance = 100) {
  tibble(
    Long = c(x, x, x + distance, x - distance),
    Lat  = c(y + distance, y - distance, y, y),
    PlotType = c("North", "South", "East", "West")
  )
}

valid_nest_ids <- nests.sf$NestID
new_points_df <- map_dfr(seq_along(valid_nest_ids), function(i) {
  generate_points(coords[i,1], coords[i,2]) %>%
    dplyr::mutate(NestID = valid_nest_ids[i])
})

# Create spatial object for coordinate extraction 
new_points_sf <- st_as_sf(
  new_points_df,
  coords = c("Long", "Lat"),
  crs = 5070
)

# Extract coordinates in meters
new.points <- new_points_sf %>%
  dplyr::mutate(
    Long = st_coordinates(geometry)[,1],
    Lat  = st_coordinates(geometry)[,2]
  ) %>%
  st_drop_geometry()

# Create spatial object for final plot to make sure coords plot correct 
new.points <- new.points %>%
  dplyr::mutate("Long1" = Long) %>%
  dplyr::mutate("Lat1" = Lat) %>%
  st_as_sf(coords = c("Long1", "Lat1"), crs = 5070)

# Check
mapview(new.points)

################################################################################
## Bring Data Together

# Combine used and available nests into a df
nests_clean <- bind_rows(
  nests.sf,
  new.points
) %>%
  dplyr::filter()

# If the PlotType is Nest it was a used nest
# If the PlotType isn't it is an available nest (100m away)
nests_clean <- nests_clean %>%
  dplyr::mutate(Case = if_else(PlotType == "Nest", 1, 0))

# Fix character encoding
nests_cleaned_nj <- nests_clean %>%
  dplyr::mutate(
    across(
      where(is.character),
      ~ iconv(.x, from = "", to = "UTF-8", sub = "")
    )
  )

# Arrange data by NestID
nests_cleaned_nj <- nests_cleaned_nj %>%
  arrange(NestID)

# Check
mapview(nests_cleaned_nj)

# Convert geometry to Lat and Long
# Put coordinates in the same crs for further analysis
nests_cleaned_nj <- nests_cleaned_nj %>%
  dplyr::mutate(
    Long = st_coordinates(geometry)[,1],  
    Lat  = st_coordinates(geometry)[,2]   
  )

################################################################################
## Finalize Data 

# Fill NAs in available nests with info from the used nest within each NestID
# This just completes the strata with all information 
nests_cleaned_nj <- nests_cleaned_nj %>%
  dplyr::group_by(NestID) %>%
  dplyr::mutate(across(
    .cols = c(BandID, CheckYr, NestFound, transmitterid_n, NestNumber, WMU, County, Township,
              Landownership, EggsHatched, EggsUnhatched, EggsDestroyed, NestFate,
              NestSubFate, nestsubfate2, nestcomments, CheckDate),
    .fns = ~ if_else(is.na(.x) & any(Case == 1), .x[Case == 1][1], .x)
  )) %>%
  ungroup() %>%
  st_drop_geometry()

################################################################################
## Output Data

# Save RDS file
# saveRDS(
#   nests_cleaned_nj,
#   "Data Management/Csvs/Processed/Nests/Nests/New Jersey//20260109_CleanedNests_2024_NJ.rds"
# )

# Write csv
# write.csv(nests_clean, "Data Management/Csvs/Processed/Nests/Nests/20260109_CleanedNests_2024_NJ.csv")

################################################################################
###############################################################################X
