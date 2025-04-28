#'---
#' title: Nest Success Modeling of Wild Turkeys in the Mid-Atlantic Region
#' authors: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script prepares data for the Habitat and Climate Known Fate Model for Maryland
#' **Note**: We are still waiting on 2024 weather data for Maryland and New Jersey.This will come in the update


################################################################################
## Load Packages

library(tidyverse)
library(terra)
library(mapview)
library(sf)
library(stringr)
library(daymetr)

################################################################################
## Data Prep- Nest-Scale Covs

#' Nest start and end date csv
nests <- read_csv("Data Management/Csvs/Processed/Incubation Dates/Maryland/20250131_NestAttempts_allbirds_MD.csv")
nests 

#' Nest veg csv
nests.veg <- read_csv("Data Management/Csvs/Processed/Nests/Nests/Maryland/20250219_CleanedNests_2022_2023_MD.csv") #%>%
  #dplyr::filter(CheckYr == "2023")

#' Filter nests.veg to only include observations that have the same nestid as nests csv
nests.veg.filtered <- nests.veg %>%
  dplyr::filter(NestID%in% nests$NestID) 

#' Merge the filtered nests.veg with nests by nestid
nests <- dplyr::right_join(nests, nests.veg.filtered, by = "NestID") %>%
  dplyr::select(-CheckDate.x) %>%
  dplyr::rename("NestFate" = NestFate.y)

#' Rename and consolidate columns 
nests <- nests %>%
  dplyr::select(NestID, 
                BandID, 
                startI, 
                endI,
                NestFate, 
                NestNumber,
                Lat,
                Long) 

#' Switch coding to UTF-8
nests <- nests %>% 
  dplyr::mutate(across(everything(), ~ iconv(., to = "UTF-8")))

################################################################################
## Data Prep- Individual Covariates

#' Read in captures csv
captures <- read_csv("Data Management/Csvs/Raw/Captures/captures_md.csv")
captures
 
#' Filter data to include only hens 
captures <- captures %>%
   dplyr::filter(sex == "F") %>%
   dplyr::select(bandid, sex, age, captyr) %>%
   dplyr::rename("BandID" = bandid)
  captures

#' Merge columns 
 nests.scaled <- merge(captures, nests, by = "BandID") 

 #' (?<=_): Only match is there is an underscore immediately before the number we are trying to extract
 #' (\\d{4}): Match exactly 4 digits
 #' (?=_) : Only match if the four digits are followed 
 #' Create NestYr column 
nests.scaled <- nests.scaled %>%
  dplyr::mutate(NestYr = str_extract(NestID, "(?<=_)(\\d{4})(?=_)")) 

#' Convert to numeric 
nests.scaled$NestYr <- as.numeric(nests.scaled$NestYr)
nests.scaled$captyr <- as.numeric(nests.scaled$captyr)

#' Create a years since capture column
nests.scaled <- nests.scaled %>%
 dplyr::mutate(yrsincecap = NestYr-captyr)
 glimpse(nests.scaled)

#' Assign Adult as the reference level
nests.scaled$age <- ifelse(nests.scaled$age == "J", 1, 
                             ifelse(nests.scaled$age == "A", 0, NA))

#' Dealing with scaling age ad hoc
#' If the bird is a juvenile and the years since capture is >1 assign it as an adult
#' If not keep the age as juvenile because turkeys will nest the first year as a juvenile
nests.scaled$age <- ifelse(nests.scaled$age == 1 & nests.scaled$yrsincecap >= 1, 0, nests.scaled$age)
table(nests.scaled$age)
glimpse(nests.scaled)

#' Nest Incubation Dates- Julian Date
nests.scaled <- nests.scaled %>%
  dplyr::mutate("Nest Incubation Date" = scale(lubridate::yday(startI))) 
nests.scaled$`Nest Incubation Date` <- as.numeric(nests.scaled$`Nest Incubation Date`)
nests.scaled$Renest <- ifelse(nests.scaled$NestNumber > 1, 1,0)


################################################################################
## Data Prep- Lanscape-Scale Covs

#' Create sf object and check projection using mapview
nests.sf <- st_as_sf(nests.scaled, coords = c("Long", "Lat"), crs = 4326) %>%
  st_transform(5070)
mapview(nests.sf)

#' Load in land cover information
md.nlcd <- terra::rast("Data Management/Rasters/NLCD/Maryland/md.nlcd.tif")
md.roads.prim <- terra::rast("Data Management/Rasters/Roads/Maryland/MdRoadRast.Prim.tiff")
md.roads.sec <- terra::rast("Data Management/Rasters/Roads/Maryland/MdRoadRast.Sec.tiff")

#' Extract landcover point value for each nest
nests.landcov <- terra::extract(md.nlcd, nests.sf) %>%
  dplyr::rename("landuse" = Class)

#' Create dummy variables for land cover classification
#' Zero nests were placed in water so I have removed it from the analysis
nests.scaled$Agriculture <- ifelse(nests.landcov$landuse == "Agriculture", 1, 0)
nests.scaled$Developed <- ifelse(nests.landcov$landuse == "Developed", 1, 0)
nests.scaled$Deciduous <- ifelse(nests.landcov$landuse == "Deciduous Forest", 1, 0)
nests.scaled$Evergreen <- ifelse(nests.landcov$landuse == "Evergreen Forest", 1, 0)
nests.scaled$Mixed <- ifelse(nests.landcov$landuse == "Mixed Forest", 1, 0)
nests.scaled$Grassland <- ifelse(nests.landcov$landuse == "Grassland/Shrub", 1, 0)
nests.scaled$Wetland <- ifelse(nests.landcov$landuse == "Wetland", 1, 0)

#' Extract distance from primary and secondary road structures
nests.prim.roads <- terra::extract(md.roads.prim, nests.sf) %>%
  dplyr::rename("Primary" = layer)
nests.sec.roads <- terra::extract(md.roads.sec, nests.sf) %>%
  dplyr::rename("Secondary" = layer)

#' Paste in distance from primary and secondary roads
nests.scaled$Primary <- nests.prim.roads$Primary
nests.scaled$Secondary <- nests.sec.roads$Secondary

#' Scale continous predictors
nests.scaled$Primary <- scale(nests.scaled$Primary) 
nests.scaled$Secondary <- scale(nests.scaled$Secondary)

#' Change predictors back to numeric 
nests.scaled$Primary <- as.numeric(nests.scaled$Primary)
nests.scaled$Secondary <- as.numeric(nests.scaled$Secondary)

#' Assign 1 if the NestFate "Hatched", otherwise assign 0
nests.scaled$NestFate <- ifelse(nests.scaled$NestFate == "Hatched", 1, 0)


################################################################################
## Data Prep- Behavior Covariates

#' Read in csv with behavior covariates
hens.behav.out <- readRDS("Data Management/Csvs/Processed/Covariates/Maryland/Behavior/hen.behav.covs.MD.RDS")
hens.behav.out

#' Only keep observations of nests.scaled that exist in nests.sample
#' Scale incubation constancy and Sum of step lengths
nests.scaled.ready <- right_join(hens.behav.out, nests.scaled) %>%
  dplyr::mutate(IncubationConstancy = scale(IncubationConstancy)) %>%
  dplyr::mutate(IncubationConstancy = as.numeric(IncubationConstancy)) %>%
  dplyr::mutate(sum_sl = scale(sum_sl)) %>%
  dplyr::mutate(sum_sl = as.numeric(sum_sl)) %>%
  dplyr::filter(TotalLocations != "NA") 
mapview(nests.scaled.ready)
nrow(nests.scaled.ready)

#' Our sample for the rest of our models
md.sample <- nests.scaled.ready %>%
  dplyr::select(NestID, BandID, NestYr)
#write_csv(md.sample,"Samples/Maryland/NestingSample_MD.csv")


################################################################################
## Data Prep- Weather Covariates - 2023

#' Create a csv file of needed information for Daymet download
#' rename and select needed columns 
nest.sites <- nests.scaled.ready %>%
  dplyr::rename("latitude" = Lat, "longitude" = Long, "site" = NestID) %>%
  dplyr::select(latitude, longitude, site) %>%
  st_drop_geometry()
write.csv(nest.sites,"Data Management/Csvs/Processed/Nests/Nests/nest.sites.weather_MD.csv")

#' Perform Daymet download  
w <- download_daymet_batch(
  file_location = 'Data Management/Csvs/Processed/Nests/Nests/nest.sites.weather_MD.csv',
  start = 2023,
  end = 2024,
  internal = TRUE)


w[[1]]$data[105:195,]

#' Check the total length of the data for nest 102 (Should be 730)
length(w[[102]]$data$tmin..deg.c)
length(w[[102]]$data$prcp..mm.day)

#' Number of weather covariates (minimum temperature and precipitation)
n.weather.cov <- 2

weather.array <- array(NA, dim = c(nrow(nests.scaled.ready),nrow(w[[1]]$data),n.weather.cov))

for (i in 1:102) {
  weather.array[i, 1:365, 1] <- w[[i]]$data$tmin..deg.c.
  weather.array[i, 1:365, 2] <- w[[i]]$data$prcp..mm.day.
}
# Assuming the 'weather.array' is pre-allocated for both 2023 and 2024
# Update the dimension to account for both years (730 days for each year)

weather.array[1,,1]
weather.array[1,,2]


################################################################################
## Weather Covariates - 2024

# Assuming that the weather.array is already pre-allocated with 730 days for both 2023 and 2024

# Assuming weather.array already has the 365 days data for 2023
# Create a new array to store the data for both years (2023 and 2024)
weather.array.copy <- array(NA, dim = c(nrow(nests.scaled.ready), 730, n.weather.cov))

# Copy the original weather.array (for 2023) to the first 365 rows of the new array
weather.array.copy[, 1:365, 1] <- weather.array[, 1:365, 1]  # tmin for 2023
weather.array.copy[, 1:365, 2] <- weather.array[, 1:365, 2]  # prcp for 2023

#' Precip 2024 extraction
precip.2024 <- terra::rast("Known Fate Model/Maryland/Weather Arrays/daymet_v4_daily_na_prcp_202400.nc")
precip.2024

#' Tmin 2024 extraction
tmin.2024 <- terra::rast("Known Fate Model/Maryland/Weather Arrays/daymet_v4_daily_na_tmin_202400.nc")
tmin.2024

# Match CRS directly from a raster layer
nests.scaled.ready <- nests.scaled.ready %>%
  st_drop_geometry() %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>%
  st_transform(crs = crs(precip.2024)) %>%
  dplyr::mutate("julianDayStart" = yday(startI),
                "julianDayEnd" = yday(endI))
mapview(nests.scaled.ready)

# Process 2024 Data (Precipitation and Tmin)
birds2024 <- which(nests.scaled.ready$NestYr == 2024)

# Loop through each bird for 2024
for(i in 1:length(birds2024)) {      # 72
  subdata <- nests.scaled.ready[birds2024[i],]  
  julians <- (subdata$julianDayStart-2):subdata$julianDayEnd
  
  # Initialize temp vectors for each nest
  precips <- c()
  tmins <- c()
  
  # Loop through Julian days for 2024 (make sure these align with 2024 weather data)
  for(j in 1:(length(julians))) {
    # Extract precipitation for 2024
    preciprast <- precip.2024[[julians[j]]]  # Precipitation raster for the corresponding Julian day
    precipvalue <- extract(preciprast, subdata)  # Extract the precipitation value
    precips <- c(precips, precipvalue[, 2])  # Append to precips vector
    
    # Extract Tmin for 2024
    tminrast <- tmin.2024[[julians[j]]]  # Tmin raster for the corresponding Julian day
    tminvalue <- extract(tminrast, subdata)  # Extract the Tmin value
    tmins <- c(tmins, tminvalue[, 2])  # Append to tmins vector
  }
  
  # Store the 2024 weather covariates in the weather.array (for the correct Julian days)
  # Offset by 365 to place data in 2024 section (366 to 730)
  weather.array.copy[birds2024[i], (julians + 365), 1] <- tmins
  weather.array.copy[birds2024[i], (julians + 365), 2] <- precips
}

weather.array.copy[1,62:90,1]

#Clean up temporary variables
#rm(precips)
#rm(tmins)

# Check the first few rows of the 2024 data (should be stored in 366-730)
print(weather.array.copy[4,,1])  # Tmin for both 2023 and 20245print(weather.array[102,,2])  # Precipitation for both 2023 and 2024

