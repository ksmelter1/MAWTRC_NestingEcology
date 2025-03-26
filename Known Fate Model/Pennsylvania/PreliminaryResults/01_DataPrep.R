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
nests <- read_csv("Data Management/Csvs/Processed/Incubation Dates/Pennsylvania/20250131_NestAttempts_allbirds_PA.csv")
nests

#' Nest veg csv
nests.veg <- read_csv("Data Management/Csvs/Processed/Nests/Vegetation Surveys/Pennsylvania/20250121_CleanedNestsVeg_2022_2023.csv")
nests.veg

#' Filter nests.veg to only include observations that have the same nestid as nests csv
nests.veg.filtered <- nests.veg %>%
  dplyr::filter(NestID%in% nests$NestID) %>%
  dplyr::filter(PlotType == "Nest")

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
                PercWoody,
                PercGrassForb,
                AvgVO,
                AvgMaxVO, 
                PercFern,
                StemCount,
                Lat,
                Long) 

#' Switch coding to UTF-8
nests <- nests %>% 
  dplyr::mutate(across(everything(), ~ iconv(., to = "UTF-8")))

#' Change columns to numeric
nests$AvgVO <- as.numeric(nests$AvgVO)
nests$AvgMaxVO <- as.numeric(nests$AvgMaxVO)
nests$PercWoody <- as.numeric(nests$PercWoody)
nests$PercGrassForb <- as.numeric(nests$PercGrassForb)
nests$PercFern <- as.numeric(nests$PercFern)
nests$StemCount <- as.numeric(nests$StemCount)

#' Scale continous predictors
nests.scaled <- nests %>% 
  dplyr::mutate(AvgVO = scale(AvgVO)) %>%
  dplyr::mutate(PercWoody = scale(PercWoody)) %>%
  dplyr::mutate(PercGrassForb = scale(PercGrassForb)) %>%
  dplyr::mutate(AvgMaxVO = scale(AvgMaxVO)) %>%
  dplyr::mutate(PercFern = scale(PercFern)) %>%
  dplyr::mutate(StemCount = scale(StemCount)) 
nests.scaled

#' Change columns back to numeric
nests.scaled$AvgVO <- as.numeric(nests.scaled$AvgVO)
nests.scaled$PercWoody <- as.numeric(nests.scaled$PercWoody)
nests.scaled$PercGrassForb <- as.numeric(nests.scaled$PercGrassForb)
nests.scaled$PercFern <- as.numeric(nests.scaled$PercFern)
nests.scaled$AvgMaxVO <- as.numeric(nests.scaled$AvgMaxVO)
nests.scaled$StemCount <- as.numeric(nests.scaled$StemCount)
nests.scaled


################################################################################
## Data Prep- Individual Covariates

#' Read in captures csv
captures <- read_csv("Data Management/Csvs/Raw/Captures/captures.csv")
captures
 
#' Filter data to include only hens 
captures <- captures %>%
   dplyr::filter(sex == "F") %>%
   dplyr::select(bandid, sex, age, captyr) %>%
   dplyr::rename("BandID" = bandid)
  captures

#' Merge columns 
 nests.scaled <- merge(captures, nests.scaled, by = "BandID") 

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

#' Read in raw virus csv
virus <- read_csv("Data Management/Csvs/Raw/Disease/LPDV_REV/virus_raw.csv")
virus
 
#' Rename columns 
virus <- virus %>%
dplyr::rename("BandID" = bandid)
virus

#' Filter data to contain only individuals within the nests.scaled df
virus.filter <- virus %>%
dplyr::filter(BandID %in% nests.scaled$BandID)

#' Read in raw parasite csv
parasite <- read_csv("Data Management/Csvs/Raw/Disease/Parasites/parasite_raw.csv")
parasite

#' Rename columns 
parasite <- parasite %>%
dplyr::rename("BandID" = bandid)
parasite

parasite.filter <- right_join(parasite, virus.filter)
colnames(parasite.filter)

#' # Specify the summed columns
columns_to_sum <- c("Capillaria sp.", "Eimeria sp.", "Ascarids", "Syngamus sp.", "Tetrameres sp.",
                   "Isospora sp.", "Monocystis sp.", "Raillietina sp.", "Choanotaenia sp.", "Strongylid")
# Create the 'ParasiteDiversity' column by summing the specified columns for each row
parasite.filter$ParasiteDiversity <- rowSums(parasite.filter[, columns_to_sum], na.rm = TRUE)

 parasite.filter.1 <- left_join(parasite.filter, virus.filter)

#' Left join filtered viral data and nests together
nests.scaled <- left_join(parasite.filter.1, nests.scaled) %>%
dplyr::select(-`Capillaria sp.`, -`Eimeria sp.`, -Ascarids, -`Syngamus sp.`,
                -`Tetrameres sp.`, -`Isospora sp.`, -`Monocystis sp.`, -`Raillietina sp.`, 
                -`Choanotaenia sp.`, -`Strongylid`)

nests.scaled$ParasiteDiversity <- scale(nests.scaled$ParasiteDiversity)
nests.scaled$ParasiteDiversity <- as.numeric(nests.scaled$ParasiteDiversity)
nests.scaled$LPDV <- as.numeric(nests.scaled$LPDV)


################################################################################
## Data Prep- Lanscape-Level Covs

#' Create sf object and check projection using mapview
#' nests.sf <- st_as_sf(nests.scaled, coords = c("Long", "Lat"), crs = 4326) 
#' mapview(nests.sf)
#' 
#' #' Load in land cover information
#' pa.nlcd <- terra::rast("Data Management/Rasters/NLCD/Pennsylvania/pa.nlcd.tif")
#' pa.roads.prim <- terra::rast("Data Management/Rasters/Roads/Pennsylvania/PaRoadRast.Prim.tiff")
#' pa.roads.sec <- terra::rast("Data Management/Rasters/Roads/Pennsylvania/PaRoadRast.Sec.tiff")
#' 
#' #' Extract landcover point value for each nest
#' nests.landcov <- terra::extract(pa.nlcd, nests.sf) %>%
#'   dplyr::rename("landuse" = Class)
#' 
#' #' Create dummy variables for land cover classification
#' #' Zero nests were placed in water so I have removed it from the analysis
#' nests.scaled$Agriculture <- ifelse(nests.landcov$landuse == "Agriculture", 1, 0)
#' nests.scaled$Developed <- ifelse(nests.landcov$landuse == "Developed", 1, 0)
#' nests.scaled$Deciduous <- ifelse(nests.landcov$landuse == "Deciduous Forest", 1, 0)
#' nests.scaled$Evergreen <- ifelse(nests.landcov$landuse == "Evergreen Forest", 1, 0)
#' nests.scaled$Mixed <- ifelse(nests.landcov$landuse == "Mixed Forest", 1, 0)
#' nests.scaled$Grassland <- ifelse(nests.landcov$landuse == "Grassland/Shrub", 1, 0)
#' 
#' 
#' #' Extract distance from primary and secondary road structures
#' nests.prim.roads <- terra::extract(pa.roads.prim, nests.sf) %>%
#'   dplyr::rename("Primary" = layer)
#' 
#' nests.sec.roads <- terra::extract(pa.roads.sec, nests.sf) %>%
#'   dplyr::rename("Secondary" = layer)

#' Paste in distance from primary and secondary roads
#' nests.scaled$Primary <- nests.prim.roads$Primary
#' nests.scaled$Secondary <- nests.sec.roads$Secondary
#' 
#' #' Scale continous predictors
#' nests.scaled$Primary <- scale(nests.scaled$Primary) 
#' nests.scaled$Secondary <- scale(nests.scaled$Secondary)
#' 
#' #' Change predictors back to numeric 
#' nests.scaled$Primary <- as.numeric(nests.scaled$Primary)
#' nests.scaled$Secondary <- as.numeric(nests.scaled$Secondary)

#' Assign 1 if the NestFate "Hatched", otherwise assign 0
nests.scaled$NestFate <- ifelse(nests.scaled$NestFate == "Hatched", 1, 0)



################################################################################
## Data Prep- Weather Covariates

#' Create a csv file of needed information for Daymet download
#' rename and select needed columns 
#' nest.sites <- nests.scaled %>%
#'   dplyr::rename("latitude" = Lat, "longitude" = Long, "site" = NestID) %>%
#'   dplyr::select(latitude, longitude, site)
#' #write.csv(nests.sites, "nests.sites.csv")
#' 
#' #' Perform Daymet download  
#' w <- download_daymet_batch(
#'   file_location = 'Data Management/Csvs/Processed/Covariates/Pennsylvania/Weather/nest.sites.weather.csv',
#'   start = 2022,
#'   end = 2023,
#'   internal = TRUE
#' )
#' 
#' w[[1]]$data[109:259,]
#' 
#' n.weather.cov<-2
#' weather.array <- array(NA, dim = c(nrow(nests.scaled),nrow(w[[1]]$data),n.weather.cov))
#' 
#' for (i in 1:169) {
#'   weather.array[i, 1:730, 1] <- w[[i]]$data$tmin..deg.c.
#'   weather.array[i, 1:730, 2] <- w[[i]]$data$prcp..mm.day.
#}
#saveRDS(nests.scaled, "nests.scaled.RDS")
#saveRDS(weather.array, "weather.RDS")

################################################################################
## Data Prep- Behavior Covariates

#nests.scaled <- readRDS("nests.scaled.RDS")
# 
# hens.behav.out <- left_join(hens.all.sf.buffer_summary, dist.traveled.cov)
# 
# nests.scaled <- nests.scaled %>%
#   st_drop_geometry() %>%
#   left_join(hens.behav.out, nests.scaled, by = "NestID") %>%
#   dplyr::mutate(IncubationConstancy = scale(IncubationConstancy)) %>%
#   dplyr::mutate(IncubationConstancy = as.numeric(IncubationConstancy)) %>%
#   dplyr::mutate(sum_sl = scale(sum_sl)) %>%
#   dplyr::mutate(sum_sl = as.numeric(sum_sl)) %>%
#   dplyr::filter(TotalLocations != "NA") %>%
#   dplyr::filter(NestID != "4255_2022_1")

