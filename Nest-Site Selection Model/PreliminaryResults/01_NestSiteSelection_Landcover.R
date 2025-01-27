
#'---
#' title: Nest-site selection of wild turkeys in Pennsylvania (an SSF analysis)
#' author: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output: 20250118_Landcover.RData
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script obtains Pennsylvania NLCD data and creates buffers for each used and available nest to extract the proportion of landcover within
#' **Last Updated**: 1/18/2025

#######################
## Load Packages 

#' Vector of package names
packages <- c("dplyr",
              "FedData",
              "mapview",
              "sf",
              "terra",
              "tidyr",
              "tigris",
              "tidyverse"
                      )

#' Function to load a package or install it if not already installed
load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

#' Apply the function to each package name
lapply(packages, load_packages)


##############
## NLCD Prep

#' GEOIDs above 60 are territories and islands so they are being removed for scaling
st <- tigris::states() 

#' Transform to albers equal area conic projection, epsg code is 5070
st <- st_transform(st, 5070)

#' Grab outline of PA
pa.outline <- subset(st, st$NAME=="Pennsylvania") %>%
  st_transform(5070)

#' Obtain Pennsylvania NLCD raster
pa.nlcd <- FedData::get_nlcd(template= pa.outline, year = 2019,
                             label = 'pa', 
                             force.redo = T)

#' Reclassify NLCD -- Skyrockets RAM
#' See covertyps for NLCD here: https://www.mrlc.gov/data/legends/national-land-cover-database-class-legend-and-description
terra::values(pa.nlcd) <- ifelse(terra::values(pa.nlcd) %in% c(21:24), yes = "Developed", ## Developed Open, low, medium, high intensity
                                 no = ifelse(terra::values(pa.nlcd) %in% c(41), yes = " Deciduous Forest", ## Deciduous
                                             no= ifelse(terra::values(pa.nlcd) %in% c(42), yes= "Evergreen Forest", ## Evergreen Forest
                                                        no= ifelse(terra::values(pa.nlcd) %in% c(43), yes= "Mixed Forest", ## Mixed Forest
                                                               no= ifelse(terra::values(pa.nlcd) %in% c(52,71), yes= "Grassland/Shrub", ## Grassland/Herbacaeous, Shrub/Scrub  
                                                                   no= ifelse(terra::values(pa.nlcd) %in% c(81:83), yes= "Agriculture", ## Pasture, Cultivated Crops
                                                                              no= "Other")))))) #Everything else
#' Check the reclassified values 
unique(terra::values(pa.nlcd))

######################
## Load in NLCD

#' Load in NLCD
pa.nlcd <- terra::rast("Data Management/Rasters/nlcd/paNLCD.tiff")

pa.nests <- read_csv("Data Management/Csvs/Processed/Nests/Vegetation Surveys/20250121_CleanedNestsVeg_2022_2023.csv")
pa.nests

pa.nests.sf <- pa.nests %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>%
  st_transform(5070)
mapview(pa.nests.sf)

#' Read in csv with incubation start and end dates
nests.inc <- read_csv("Data Management/Csvs/Processed/IncubationDates/Draft2/20250124_NestAttempts_allbirds.csv")  
nests.inc

#' Merge the filtered nests.veg with nests by nestid
pa.nests<- inner_join(pa.nests.sf, nests.inc, by = "NestID") %>%
  dplyr::select(-NestFate.y, -...1, -CheckDate.y )

#' Extract point value at each nest
landcov <-terra::extract(pa.nlcd, pa.nests)
  
#' Bind columns together
pa.nests.landcov <- cbind(pa.nests, landcov) %>%
  dplyr::rename("landuse" = Class)

############################################
## Save RData file 20250118_Landcover.RData

############################################################################################
############################################################################################
