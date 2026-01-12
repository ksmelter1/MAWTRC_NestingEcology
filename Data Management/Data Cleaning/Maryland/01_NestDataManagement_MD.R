#'---
#' title: Nest and Vegetation Survey Data Management for Maryland 
#' author: "K. Smelter
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script shows the amount of hens that nested across years and organizes nest data and shapefiles for analysis


################################################################################
## Load Packages and Data

library(tidyverse)
library(sf)
library(mapview)

# Load raw nest csv file
nests.raw <- read.csv("Data Management/Csvs/Raw/Nests/nests.raw.md.csv")

################################################################################
## Data Prep

# Rename columns 
nests.raw.md <- nests.raw %>%
  rename(
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
# One no nest wasn't removed from earlier
nest_df <- nests.raw.md %>%
  dplyr::filter(NestFound == "Y") %>%
  dplyr::filter(NestID != "0320_NoNEST_06_09_2025") 

# Create a checkdate column
# Format a year, month, date
nest_df <- nest_df %>%
  mutate(
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
  mutate(
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
# Removed via visual inspection
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
    mutate(NestID = valid_nest_ids[i])
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
  mutate(Case = if_else(PlotType == "Nest", 1, 0))

# Fix character encoding
nests_cleaned_md <- nests_clean %>%
  mutate(
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
  mutate(
    Long = st_coordinates(geometry)[,1],  
    Lat  = st_coordinates(geometry)[,2]   
  )

################################################################################
## Finalize Data 

# Fill NAs in available nests with info from the used nest within each NestID
# This just completes the strata with all information 
nests_cleaned_md <- nests_cleaned_md %>%
  group_by(NestID) %>%
  mutate(across(
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
saveRDS(
  nests_cleaned_md,
  "Data Management/Csvs/Processed/Nests/Nests/Maryland/20260109_CleanedNests_2023_2024_2025_MD.rds"
)

################################################################################
###############################################################################X