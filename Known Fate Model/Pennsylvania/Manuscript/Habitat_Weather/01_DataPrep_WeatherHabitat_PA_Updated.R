#---
# title: Nest Success Modeling of Wild Turkeys in the Mid-Atlantic Region
# authors: "K. Smelter
# date: "`r format(Sys.time(), '%d %B, %Y')`"
# output:
#   html_document: 
#     toc: true
#---
#  
# **Purpose**: This script prepares Pennsylvania data for the nest survival analysis
# **Key Changes**: This script incorporates 2024 data

################################################################################
## Load Packages

library(tidyverse)
library(terra)
library(mapview)
library(sf)
library(stringr)
library(daymetr)
library(FedData)
library(GGally)

################################################################################
## Data Prep- Nest-Level Covs

# Nest start and end date csv
nests <- read_csv("Data Management/Csvs/Processed/Incubation Dates/Pennsylvania/20250709_NestAttempts_allbirds_PA_NoFilter.csv")
nests

# Nest veg csv
nests.veg <- read_csv("Data Management/Csvs/Processed/Nests/Vegetation Surveys/Pennsylvania/20250629_CleanedNestsVeg_2022_2023_2024.csv")
nests.veg

# Filter nests.veg to only include observations that have the same nestid as nests csv
nests.veg.filtered <- nests.veg %>%
  dplyr::filter(NestID%in% nests$NestID) %>%
  dplyr::filter(PlotType == "Nest")

# Merge the filtered nests.veg with nests by nestid
nests <- dplyr::right_join(nests, nests.veg.filtered, by = "NestID") %>%
  dplyr::select(-CheckDate.x) %>%
  dplyr::rename("NestFate" = NestFate.y)

# Create a NestID1 column that contains plottype to merge basal area data
nests$NestID1 <- paste(nests$NestID, nests$PlotType, sep = "_")

# Rename and consolidate columns 
nests <- nests %>%
  dplyr::select(NestID, 
                NestID1,
                BandID, 
                startI, 
                endI,
                NestFate, 
                PercWoody,
                PercGrassForb,
                AvgVO,
                PercFern,
                StemCount,
                PercLitter,
                PercBare,
                PercBoulder,
                PlotType,
                Lat,
                Long) 

# Read in basal area csv
basal <- read_csv("Data Management/Csvs/Processed/Covariates/Pennsylvania/BasalArea.csv")
basal

# Join nest data with basal area
nests <- right_join(basal, nests)

# Replace NA values of basal area with 0 because the NAs are zero values
nests$Basal <- tidyr::replace_na(nests$Basal, 0)

# Switch coding to UTF-8
nests <- nests %>% 
  dplyr::mutate(across(everything(), ~ iconv(., to = "UTF-8")))

# Change columns to numeric
nests$AvgVO <- as.numeric(nests$AvgVO)
nests$PercWoody <- as.numeric(nests$PercWoody)
nests$PercGrassForb <- as.numeric(nests$PercGrassForb)
nests$PercFern <- as.numeric(nests$PercFern)
nests$StemCount <- as.numeric(nests$StemCount)
nests$PercBoulder <- as.numeric(nests$PercBoulder)
nests$PercLitter <- as.numeric(nests$PercLitter)
nests$PercBare <- as.numeric(nests$PercBare)
nests$Basal <- as.numeric(nests$Basal)

# Scale continous predictors
nests.scaled <- nests %>% 
  dplyr::mutate(AvgVO = scale(AvgVO)) %>%
  dplyr::mutate(PercWoody = scale(PercWoody)) %>%
  dplyr::mutate(PercGrassForb = scale(PercGrassForb)) %>%
  dplyr::mutate(PercFern = scale(PercFern)) %>%
  dplyr::mutate(StemCount = scale(StemCount)) %>%
  dplyr::mutate(PercBoulder = scale(PercBoulder)) %>% 
  dplyr::mutate(PercLitter = scale(PercLitter)) %>%
  dplyr::mutate(PercBare = scale(PercBare)) %>%
  dplyr::mutate(Basal = scale(Basal)) 
nests.scaled

# Change columns back to numeric
nests.scaled$AvgVO <- as.numeric(nests.scaled$AvgVO)
nests.scaled$PercWoody <- as.numeric(nests.scaled$PercWoody)
nests.scaled$PercGrassForb <- as.numeric(nests.scaled$PercGrassForb)
nests.scaled$PercFern <- as.numeric(nests.scaled$PercFern)
nests.scaled$StemCount <- as.numeric(nests.scaled$StemCount)
nests.scaled$PercBoulder <- as.numeric(nests.scaled$PercBoulder)
nests.scaled$PercLitter <- as.numeric(nests.scaled$PercLitter)
nests.scaled$PercBare <- as.numeric(nests.scaled$PercBare)
nests.scaled$Basal <- as.numeric(nests.scaled$Basal)
nests.scaled

# Check Correlations
C <- cbind(
     nests.scaled$AvgVO,                       # Visual Obstruction
     nests.scaled$PercWoody,                   # Percent Woody Vegetation
     nests.scaled$PercGrassForb,               # Percent Grass/Forb
     nests.scaled$PercFern,                    # Percent Fern
     nests.scaled$PercLitter,                  # Percent Litter
     nests.scaled$PercBoulder,                 # Percent Boulder
     nests.scaled$PercBare,                    # Percent Bare Ground
     nests.scaled$StemCount,                   # Woody Stem Count
     nests.scaled$Basal                        # Basal Area
   ) 

# Use ggpairs to visualize correlations
ggpairs(as.data.frame(C), 
        upper = list(continuous = "cor"),
        diag = list(continuous = "barDiag"),
        lower = list(continuous = "points"))


################################################################################
## Data Prep- Individual Covariates

# Read in captures csv including 2024 data
captures <- read_csv("Data Management/Csvs/Raw/Captures/captures_pa.csv")
captures
 
# Filter data to include only hens 
captures <- captures %>%
   dplyr::filter(sex == "F") %>%
   dplyr::select(bandid, sex, age, captyr) %>%
   dplyr::rename("BandID" = bandid)
  captures

# Merge columns 
nests.scaled <- merge(captures, nests.scaled, by = "BandID") 

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
# If the bird is a juvenile and the years since capture is >1 assign it as an adult
# If not keep the age as juvenile because turkeys will nest the first year as a juvenile
nests.scaled$age <- ifelse(nests.scaled$age == 1 & nests.scaled$yrsincecap >= 1, 0, nests.scaled$age)
table(nests.scaled$age)
glimpse(nests.scaled)

# Nest Incubation Dates- Calendar Day
# Scale and keep as a continous variable
nests.scaled <- nests.scaled %>%
  dplyr::mutate("Nest Incubation Date" = scale(lubridate::yday(startI))) 
nests.scaled$`Nest Incubation Date` <- as.numeric(nests.scaled$`Nest Incubation Date`)


################################################################################
## Data Prep- Lanscape-Scale Covs

# Create sf object and check projection using mapview
nests.sf <- st_as_sf(nests.scaled, coords = c("Long", "Lat"), crs = 4326) %>%
  st_transform(5070)
mapview(nests.sf)

# Load in land cover information
pa.nlcd <- terra::rast("Data Management/Rasters/NLCD/Atlantic/atlantic.nlcd.tif")
pa.roads.prim <- terra::rast("Data Management/Rasters/Roads/Atlantic/atlantic.prim.tiff")
pa.roads.sec <- terra::rast("Data Management/Rasters/Roads/Atlantic/atlantic.sec.tiff")

# Extract landcover point value for each nest
nests.landcov <- terra::extract(pa.nlcd, nests.sf) %>%
  dplyr::rename("landuse" = Class)

# Create dummy variables for land cover classification
# Zero nests were placed in water so I have removed it from the analysis
nests.scaled$Pasture <- ifelse(nests.landcov$landuse == "Pasture", 1, 0)
nests.scaled$Crop <- ifelse(nests.landcov$landuse == "Crop", 1, 0)
nests.scaled$Developed <- ifelse(nests.landcov$landuse == "Developed", 1, 0)
nests.scaled$Deciduous <- ifelse(nests.landcov$landuse == "Deciduous Forest", 1, 0)
nests.scaled$Evergreen <- ifelse(nests.landcov$landuse == "Evergreen Forest", 1, 0)
nests.scaled$Mixed <- ifelse(nests.landcov$landuse == "Mixed Forest", 1, 0)
nests.scaled$Grassland <- ifelse(nests.landcov$landuse == "Grassland/Shrub", 1, 0)
nests.scaled$Wetland <- ifelse(nests.landcov$landuse == "Wetland", 1, 0)

# Extract distance from primary and secondary road structures
nests.prim.roads <- terra::extract(pa.roads.prim, nests.sf) %>%
  dplyr::rename("Primary" = layer)
nests.sec.roads <- terra::extract(pa.roads.sec, nests.sf) %>%
  dplyr::rename("Secondary" = layer)

# Paste in distance from primary and secondary roads
nests.scaled$Primary <- nests.prim.roads$Primary
nests.scaled$Secondary <- nests.sec.roads$Secondary

# Scale continous predictors
nests.scaled$Primary <- scale(nests.scaled$Primary) 
nests.scaled$Secondary <- scale(nests.scaled$Secondary)

# Change predictors back to numeric 
nests.scaled$Primary <- as.numeric(nests.scaled$Primary)
nests.scaled$Secondary <- as.numeric(nests.scaled$Secondary)

# Assign 1 if the NestFate "Hatched", otherwise assign 0
nests.scaled$NestFate <- ifelse(nests.scaled$NestFate == "Hatched", 1, 0)


################################################################################
## Data Prep- Behavior Covariates

# Read in RDS file from the behavior script
hens.behav.out <- readRDS("Data Management/Csvs/Processed/Covariates/Pennsylvania/Behavior/hens.behav.out.updated.RDS")

# Drop geometry
# Keep all observations of hens.behav.out that exist in nests.scaled 
nests.scaled <- nests.scaled %>%
  st_drop_geometry() %>%
  right_join(hens.behav.out, nests.scaled, by = "NestID")


# Calculate mean and standard deviation by age class
nests.scaled %>%
  group_by(age) %>%
  summarise(
    mean_constancy = mean(IncubationConstancy, na.rm = TRUE),
    sd_constancy = sd(IncubationConstancy, na.rm = TRUE)
  )

# Total number of nests by age class
nests.scaled %>%
  group_by(age) %>%
  summarise(count = n())

# Scale incubation constancy and change it back to numeric
# Scale cumulative step length and change it back to numeric 
# Remove all observations where the total locations of gps locations within a buffer is NA
nests.scaled <- nests.scaled %>%
  dplyr::mutate(IncubationConstancy = scale(IncubationConstancy)) %>%
  dplyr::mutate(IncubationConstancy = as.numeric(IncubationConstancy)) %>%
  dplyr::mutate(sum_sl = scale(sum_sl)) %>%
  dplyr::mutate(sum_sl = as.numeric(sum_sl)) %>%
  dplyr::filter(TotalLocations != "NA") 

################################################################################
## Data Prep- Weather Covariates

# Create a csv file of needed information for Daymet download
# rename and select needed columns 
nest.sites <- nests.scaled %>%
  dplyr::rename("latitude" = Lat, "longitude" = Long, "site" = NestID) %>%
  dplyr::select(latitude, longitude, site)
write.csv(nest.sites, "Data Management/Csvs/Processed/Covariates/Pennsylvania/Weather/nests.sites.updated.csv")

# Perform Daymet download  
# Download data from 2022 through the end of 2023
w <- download_daymet_batch(
  file_location = 'Data Management/Csvs/Processed/Covariates/Pennsylvania/Weather/nests.sites.updated.csv',
  start = 2022,
  end = 2024,
  internal = TRUE
)

# Check
w[[1]]$data[109:259,]

# Create a weather array with 2 covariates
# Precipitation and Minimum Temperature
n.weather.cov<-2
weather.array <- array(NA, dim = c(nrow(nests.scaled),nrow(w[[1]]$data),n.weather.cov))

# Fill the weather array with data from 2022-2023
for (i in 1:260) {
  weather.array[i, 1:730, 1] <- w[[i]]$data$tmin..deg.c.
  weather.array[i, 1:730, 2] <- w[[i]]$data$prcp..mm.day.
}


################################################################################
## Weather Covariates - 2024

# Assuming that the weather.array is already pre-allocated with 730 days for both 2023 and 2024

# Assuming weather.array already has the 365 days data for 2023
# Create a new array to store the data for both years (2023 and 2024)
weather.array.copy <- array(NA, dim = c(nrow(nests.scaled), 1095, n.weather.cov))

# Copy the original weather.array (for 2023) to the first 365 rows of the new array
weather.array.copy[, 1:730, 1] <- weather.array[, 1:730, 1]  # tmin for 2023
weather.array.copy[, 1:730, 2] <- weather.array[, 1:730, 2]  # prcp for 2023

# Precip 2024 extraction
# This is a raster of the whole Mid-Atlantic region
precip.2024 <- terra::rast("Known Fate Model/Maryland/Weather Arrays/daymet_v4_daily_na_prcp_202400.nc")
precip.2024

# Tmin 2024 extraction
# This is a raster of the whole Mid-Atlantic region
tmin.2024 <- terra::rast("Known Fate Model/Maryland/Weather Arrays/daymet_v4_daily_na_tmin_202400.nc")
tmin.2024

# Match CRS directly from a raster layer
nests.scaled.ready <- nests.scaled %>%
  st_drop_geometry() %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>%
  st_transform(crs = crs(precip.2024)) %>%
  dplyr::mutate("julianDayStart" = yday(startI),
                "julianDayEnd" = yday(endI))
mapview(nests.scaled.ready)

# Process 2024 Data (Precipitation and Tmin)
birds2024 <- which(nests.scaled.ready$NestYr == 2024)
birds2024

# Loop through each bird for 2024
for(i in 1:length(birds2024)) {      
  subdata <- nests.scaled.ready[birds2024[i],]  
  julians <- (subdata$julianDayStart-2):subdata$julianDayEnd
  
  # Initialize temp vectors for each nest
  precips <- c()
  tmins <- c()
  
  # Loop through Julian days for 2024 (make sure these align with 2024 weather data)
  for(j in 1:(length(julians))) {
    # Extract precipitation for 2024
    preciprast <- precip.2024[[julians[j]]]  # Precipitation raster for the corresponding Julian day
    precipvalue <- terra::extract(preciprast,subdata)  # Extract the precipitation value
    precips <- c(precips, precipvalue[, 2])  # Append to precips vector
    
    # Extract Tmin for 2024
    tminrast <- tmin.2024[[julians[j]]]  # Tmin raster for the corresponding Julian day
    tminvalue <- terra::extract(tminrast,subdata)  # Extract the Tmin value
    tmins <- c(tmins, tminvalue[, 2])  # Append to tmins vector
  }
  
  # Store the 2024 weather covariates in the weather.array (for the correct Julian days)
  # Offset by 730 to place data in 2024 section (731 to 1095)
  weather.array.copy[birds2024[i], (julians + 730), 1] <- tmins
  weather.array.copy[birds2024[i], (julians + 730), 2] <- precips
}

# Check to make sure that 2024 weather data is included in the array
weather.array.copy[13,731:1095,1]

# Check the first few rows of the 2024 data (should be stored in 731-1095)
print(weather.array.copy[13,731:1095,2])  


################################################################################
## Output sample for other models

# Output sample of used nests for other models
sample <- nests.scaled.ready %>%
  dplyr::select(NestID, BandID, startI, endI, NestFate)

# Write csv
#write_csv(sample, "Samples/Pennsylvania/PA_Sample_Updated.csv")

################################################################################
###############################################################################X