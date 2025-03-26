#'---
#' title: Nest Success Modeling of Wild Turkeys in the Mid-Atlantic Region
#' authors: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script uses simulated data to fit a Bayesian known fate model for female wild turkeys in our study
#' **Key Changes**: This script uses the cloglog link instead of the logit link for modeling daily nest survival

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
                Lat,
                Long) 

#' Switch coding to UTF-8
nests <- nests %>% 
  dplyr::mutate(across(everything(), ~ iconv(., to = "UTF-8")))

#' Change columns to numeric
# nests$AvgVO <- as.numeric(nests$AvgVO)
# nests$AvgMaxVO <- as.numeric(nests$AvgMaxVO)
# nests$PercWoody <- as.numeric(nests$PercWoody)
# nests$PercGrassForb <- as.numeric(nests$PercGrassForb)
# nests$PercFern <- as.numeric(nests$PercFern)

#' Scale continous predictors
# nests.scaled <- nests %>% 
#   dplyr::mutate(AvgVO = scale(AvgVO)) %>%
#   dplyr::mutate(PercWoody = scale(PercWoody)) %>%
#   dplyr::mutate(PercGrassForb = scale(PercGrassForb)) %>%
#   dplyr::mutate(AvgMaxVO = scale(AvgMaxVO)) %>%
#   dplyr::mutate(PercFern = scale(PercFern))
# nests.scaled

#' Change columns back to numeric
# nests.scaled$AvgVO <- as.numeric(nests.scaled$AvgVO)
# nests.scaled$PercWoody <- as.numeric(nests.scaled$PercWoody)
# nests.scaled$PercGrassForb <- as.numeric(nests.scaled$PercGrassForb)
# nests.scaled$PercFern <- as.numeric(nests.scaled$PercFern)
# nests.scaled$AvgMaxVO <- as.numeric(nests.scaled$AvgMaxVO)
# nests.scaled


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
#' If the bird is an adult and the years since capture is >1 assign it as an adult
#' If not keep the age as juvenile because turkeys will nest the first year as a juvenile
nests.scaled$age <- ifelse(nests.scaled$age == 1 & nests.scaled$yrsincecap >= 1, 0, nests.scaled$age)
table(nests.scaled$age)
glimpse(nests.scaled)

#' Nest Incubation Dates- Julian Date
nests.scaled <- nests.scaled %>%
  dplyr::mutate("Nest Incubation Date" = scale(lubridate::yday(startI))) 

nests.scaled$`Nest Incubation Date` <- as.numeric(nests.scaled$`Nest Incubation Date`)

#' ################################################################################
#' ## Data Prep- Disease Covariates
#' 
#' #' Read in raw virus csv
#' virus <- read_csv("Data Management/Csvs/Raw/Disease/LPDV_REV/virus_raw.csv")
#' virus
#' 
#' #' Rename columns 
#' virus <- virus %>%
#'   dplyr::rename("BandID" = bandid)
#' virus
#' 
#' #' Filter data to contain only individuals within the nests.scaled df
#' virus.filter <- virus %>%
#'   dplyr::filter(BandID %in% nests.scaled$BandID)
#' 
#' #' Read in raw parasite csv
#' parasite <- read_csv("Data Management/Csvs/Raw/Disease/Parasites/parasite_raw.csv")
#' parasite
#' 
#' #' Rename columns 
#' parasite <- parasite %>%
#'   dplyr::rename("BandID" = bandid)
#' parasite
#' 
#' #' Filter to only include individua;s within nests.scaled df
#' parasite.filter <- parasite %>%
#'   dplyr::filter(BandID %in% nests.scaled$BandID)
#' 
#' # Assuming 'parasite.filter' is your data frame
#' # Specify the columns you want to sum
#' columns_to_sum <- c("Capillaria sp", "Eimeria sp", "Ascarids", "Syngamus sp", "Tetrameres sp",
#'                     "Isospora sp", "Monocystis sp", "Raillietina sp", "Choanotaenia sp", "Strongylid")
#' 
#' # Create the 'ParasiteDiversity' column by summing the specified columns for each row
#' parasite.filter$ParasiteDiversity <- rowSums(parasite.filter[, columns_to_sum], na.rm = TRUE)
#' 
#' parasite.filter.1 <- left_join(parasite.filter, virus.filter)
#' 
#' #' Left join filtered viral data and nests together
#' nests.scaled <- left_join(parasite.filter.1, nests.scaled) %>%
#'   dplyr::select(-`Capillaria sp`, -`Eimeria sp`, -Ascarids, -`Syngamus sp`,
#'                 -`Tetrameres sp`, -`Isospora sp`, -`Monocystis sp`, -`Raillietina sp`, 
#'                 -`Choanotaenia sp`, -`Strongylid`)
#' 
#' nests.scaled$ParasiteDiversity <- scale(nests.scaled$ParasiteDiversity)
#' nests.scaled$ParasiteDiversity <- as.numeric(nests.scaled$ParasiteDiversity)


################################################################################
## Data Prep- Lanscape-Level Covs

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
## Data Prep- Weather Covariates

#' Create a csv file of needed information for Daymet download
#' rename and select needed columns 
nest.sites <- nests.scaled %>%
  dplyr::rename("latitude" = Lat, "longitude" = Long, "site" = NestID) %>%
  dplyr::select(latitude, longitude, site)
#write.csv(nest.sites,"Data Management/Csvs/Processed/Nests/Nests/nest.sites.weather_MD.csv")

#' Perform Daymet download  
w <- download_daymet_batch(
  file_location = 'Data Management/Csvs/Processed/Nests/Nests/nest.sites.weather_MD.csv',
  start = 2023,
  end = 2024,
  internal = TRUE)

w[[1]]$data[105:195,]

# Check the total length of the data for nest 104
length(w[[31]]$data$tmin..deg.c)
length(w[[31]]$data$prcp..mm.day)


# Number of weather covariates (minimum temperature and precipitation)
n.weather.cov <- 2

weather.array <- array(NA, dim = c(nrow(nests.scaled),nrow(w[[1]]$data),n.weather.cov))

# for (i in 1:31) {
#   weather.array[i, 1:365, 1] <- w[[i]]$data$tmin..deg.c.
#   weather.array[i, 1:365, 2] <- w[[i]]$data$prcp..mm.day.
# }

#' Array dimensions: We need to accommodate both 2023 (365 days) and 2024 (366 days)
days_in_2023 <- 365
days_in_2024 <- 366
total_days <- days_in_2023 + days_in_2024

#' Create a weather array to store data for both years
weather.array <- array(NA, dim = c(nrow(nests.scaled), total_days, n.weather.cov))

#' Loop through each nest data (assuming we have 104 nests)
 for (i in 1:104) {
 #' 2023: The first 365 days are for 2023
weather.array[i, 1:days_in_2023, 1] <- w[[i]]$data$tmin..deg.c[1:days_in_2023]
weather.array[i, 1:days_in_2023, 2] <- w[[i]]$data$prcp..mm.day[1:days_in_2023]

# 2024: The next 366 days are for 2024 (leap year)
weather.array[i, (days_in_2023+1):total_days, 1] <- w[[i]]$data$tmin..deg.c[(days_in_2023+1):total_days]
weather.array[i, (days_in_2023+1):total_days, 2] <- w[[i]]$data$prcp..mm.day[(days_in_2023+1):total_days]
}

#' Save weather.array object, weather.array2 will be used in the model (stored in model script)
#saveRDS(weather.array, "weather.RDS_MD.RDS")


################################################################################
## Data Prep- Behavior Covariates

#' Read in csv with behavior covariates
hens.behav.out <- readRDS("Data Management/Csvs/Processed/Covariates/Maryland/Behavior/hen.behav.covs.MD.RDS")
hens.behav.out

#' Only keep observations that exist in nessts.sample 
#nests.scaled <- right_join(hens.behav.out, nests.sample) 

#' Only keep observations of nests.scaled that exist in nests.sample
#nests.scaled.ready <- right_join(nests.scaled, nests.sample)


#' Only keep observations of nests.scaled that exist in nests.sample
#' Scale incubation constancy and Sum of step lengths
nests.scaled.ready <- right_join(hens.behav.out, nests.scaled) %>%
  dplyr::mutate(IncubationConstancy = scale(IncubationConstancy)) %>%
  dplyr::mutate(IncubationConstancy = as.numeric(IncubationConstancy)) %>%
  dplyr::mutate(sum_sl = scale(sum_sl)) %>%
  dplyr::mutate(sum_sl = as.numeric(sum_sl)) %>%
  dplyr::filter(TotalLocations != "NA") 
mapview(nests.scaled.ready)

md.sample <- nests.scaled.ready %>%
  dplyr::select(NestID, BandID, NestYr)
write_csv(md.sample,"Samples/Maryland/NestingSample_MD.csv")

################################################################################
################################################################################

