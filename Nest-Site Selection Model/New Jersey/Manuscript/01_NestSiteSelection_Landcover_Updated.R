#'---
#' title: Nest-site selection of female wild turkeys in New Jersey
#' author: "K. Smelter
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#'---
#'  
#' **Purpose**: This script obtains NLCD values for used and available nests in each state
#' **Last Updated**: 12/27/2025

################################################################################
## Load Packages 

packages <- c("dplyr",
              "FedData",
              "mapview",
              "sf",
              "terra",
              "tidyr",
              "tigris",
              "tidyverse",
              "rasterVis")

load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

lapply(packages, load_packages)


################################################################################
## Load in NLCD

# Read in created NLCD file from above
pa.nlcd <- terra::rast("Data Management/Rasters/NLCD/Atlantic/atlantic.nlcd.tif")
nj.nlcd <- terra::rast("Data Management/Rasters/NLCD/Atlantic/atlantic.nlcd.tif")
md.nlcd <- terra::rast("Data Management/Rasters/NLCD/Atlantic/atlantic.nlcd.tif")


################################################################################
## Pennsylvania

# Load in nests csv
pa.nests <- read_csv("Data Management/Csvs/Processed/Nests/Vegetation Surveys/Pennsylvania/20250629_CleanedNestsVeg_2022_2023_2024.csv")
pa.nests

# Convert to sf object
pa.nests.sf <- pa.nests %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>%
  st_transform(5070)
mapview(pa.nests.sf)

# Read in csv with incubation start and end dates
nests.inc <- read_csv("Data Management/Csvs/Processed/Incubation Dates/Pennsylvania/20250709_NestAttempts_allbirds_PA_Ready.csv")  
nests.inc

# Merge the filtered nests.veg with nests by nestid
pa.nests<- inner_join(pa.nests.sf, nests.inc, by = "NestID") %>%
  dplyr::select(-NestFate.y, -CheckDate.y )

# Extract point value at each nest
landcov <-terra::extract(pa.nlcd, pa.nests)

# Bind columns together
pa.nests.landcov <- cbind(pa.nests, landcov) %>%
  dplyr::rename("landuse" = Class)


################################################################################
## New Jersey

# Read in nests csv
nj.nests <- read_csv("Data Management/Csvs/Processed/Nests/Nests/New Jersey/20250219_CleanedNests_2022_2023_NJ.csv")
nj.nests

# Convert to sf object
nj.nests.sf <- nj.nests %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>%
  st_transform(5070)
mapview(nj.nests.sf)

# Read in csv with incubation start and end dates
nests.inc.nj <- read_csv("Data Management/Csvs/Processed/Incubation Dates/New Jersey/20250717_NestAttempts_allbirds_NJ_Ready.csv")  
nests.inc.nj

# Merge the filtered nests.veg with nests by nestid
nj.nests<- inner_join(nj.nests.sf, nests.inc.nj, by = "NestID") %>%
  dplyr::select(-NestFate.y, -...1, -CheckDate.y )

# Extract point value at each nest
landcov <-terra::extract(nj.nlcd, nj.nests.sf)

# Bind columns together
nj.nests.landcov <- cbind(nj.nests.sf, landcov) %>%
  dplyr::rename("landuse" = Class)


################################################################################
## Maryland

# Read in nests csv
md.nests <- read_csv("Data Management/Csvs/Processed/Nests/Nests/Maryland/20250219_CleanedNests_2022_2023_MD.csv")
md.nests

# Convert to sf object
md.nests.sf <- md.nests %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>%
  st_transform(5070)
mapview(md.nests.sf)

# Read in csv with incubation start and end dates
nests.inc.md <- read_csv("Data Management/Csvs/Processed/Incubation Dates/Maryland/20250717_NestAttempts_allbirds_MD_Ready.csv")  
nests.inc.md

# Merge the filtered nests.veg with nests by nestid
md.nests<- inner_join(md.nests.sf, nests.inc.md, by = "NestID") %>%
  dplyr::select(-NestFate.y, -...1, -CheckDate.y )

# Extract point value at each nest
landcov <-terra::extract(md.nlcd, md.nests)

# Bind columns together
md.nests.landcov <- cbind(md.nests, landcov, by = "NestID") %>%
  dplyr::rename("landuse" = Class)


################################################################################
## Save RData file *InsertDate*_Landcover.RData

################################################################################
###############################################################################X