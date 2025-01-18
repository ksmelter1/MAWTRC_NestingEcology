
#'---
#' title: Habitat selection of female wild turkeys during incubation (an SSF analysis)
#' author: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'
#+ include = FALSE
knitr::opts_chunk$set(warning = FALSE, message = FALSE, cache = TRUE)
#'  
#' **Purpose**: This script creates tracks and extracts covariates using the amt package
#'

#####################
## Load Packages ##
#####################

require(sf)
require(amt)
require(tigris)
require(FedData)
require(dplyr)
require(terra)

###############################
## Prepare NLCD Covariates ##
###############################

#' GEOIDs above 60 are territories and islands so they are being removed for scaling
st <- tigris::states() %>%
  dplyr::filter(GEOID < "60") %>% 
  tigris::shift_geometry()

#' Transform to albers equal area conic projection, epsg code is 5070
st <- st_transform(st, 5070)

#' Grab outline of PA
pa.outline <- subset(st, st$NAME=="Pennsylvania")

#' Obtain Pennsylvania NLCD raster
pa.nlcd <- FedData::get_nlcd(template= pa.outline, year = 2021, 
                             label = 'pa', 
                            force.redo = T)

#' Reclassify NLCD -- Skyrockets RAM
#' See covertyps for NLCD here: https://www.mrlc.gov/data/legends/national-land-cover-database-class-legend-and-description
terra::values(pa.nlcd) <- ifelse(terra::values(pa.nlcd) %in% c(21:24), yes = "Developed", ## Developed Open, low, medium, high intensity
                                no = ifelse(terra::values(pa.nlcd) %in% c(41:43), yes = "Forest", ## Deciduous, Evergreen, Mixed Forest -- Everything else is Non
                                            no = ifelse(terra::values(pa.nlcd) %in% c(51:74), yes = "Shrub", #Shrub scrub and herbaceous
                                                        no= ifelse(terra::values(pa.nlcd) %in% c(81:83), yes= "Agriculture", ## Pasture, Cultivated Crops
                                                                  no= "NonForest")))) #Everything else
#' Crop PA NLCD raster to the Pennsylvania outline
pa.nlcd <- crop(pa.nlcd, vect(pa.outline))
#' Mask the PA NLCD raster according to the Pennsylvania outline 
pa.nlcd <- mask(pa.nlcd, vect(pa.outline))

#' Save raster
#writeRaster(pa.nlcd, "paNLCD.tiff", overwrite=T)

pa.nlcd <- terra::rast("Data Management/Rasters/nlcd/paNLCD.tiff")

#################
## Elevation ##
################

#' Download elevation data following the pa.outline 
DEM <- FedData::get_ned(template = pa.outline, 
                        res = 1, label = 'dem', force.redo = T )
plot(DEM) #Check

#' Project DEM with pa.nlcd, this allows for a consistent crs and layer 1 is just the elevation
DEM.proj <- project(DEM, pa.nlcd, method= 'bilinear')

#' Save raster
#writeRaster(DEM,"paDEM.tiff", overwrite=T)

DEM.proj <- terra::rast("Data Management/Rasters/nlcd/paDEM.tiff")

##############################
## Create and Save Tracks ##
##############################

#' Open storage folder
turkey.data.files<-list.files(paste(getwd(),
                                    "Data Management/RData/HenMovementData_Nesting",sep="/"))
random_steps.sf<-vector(mode="list", length(turkey.data.files))

#' For loop to extract turkey nest movement data
for(i in 1:10){
#for (i in 1:length(turkey.data.files)){
  #' Read in RDS file and create an sf object
  turkey.ind<- readRDS(paste(getwd(),"Data Management/RData/HenMovementData_Nesting",
                             turkey.data.files[i],sep="/")) #%>%
    #st_as_sf(coords= c("long", "lat"),crs=5070)
  
  #' Create a track
  #' New amt track function is just amt::track
  gobble.ind <-amt::make_track(turkey.ind, .x= long, .y=lat, .t=timestamp, BirdID, 
                               crs= 4326) %>%
                                amt::transform_coords(5070)
  #' Check sampling rate
  amt::summarize_sampling_rate(gobble.ind)
  
  #' Check cumulative distance
  amt::cum_dist(gobble.ind)
  
  #' Check straightness (Directional persistence)
  amt::straightness(gobble.ind)
  
  #' Total distance
  amt::tot_dist(gobble.ind)
  
  #' Set criteria for steps 
  #' Rate: sampling rate
  #' Tolerance: the tolerance of deviations of the sampling rate
  #' Steps by burst: Returns NA for zero step lengths and calculates step lengths for points on a path
  #' Time of day: Was the fix taken during the day or night? we could filter out night locations
  stps <- amt::track_resample(gobble.ind,
                              rate = minutes(60),
                               tolerance = seconds(60))%>%
                              amt::steps_by_burst(keep_cols="end") %>%
                              amt::time_of_day(include.crepuscule=T) 
  
  #' Create random steps
  #' Exponential step length used due to issues with formatting gamma
  #' Extract covariates at the end of each used and available step
  #' Random step lengths drawn from gamma distribution
  #' Random turning angles drawn from a vonmises distribution
  #' Include_observed: Include all used steps in the analysis
  random_steps <-amt::random_steps(
    stps,
    n_control = 10,
    sl_distr = fit_distr(stps$sl_, "gamma"),
    ta_distr = fit_distr(stps$ta_, "vonmises"),
    include_observed = T) %>%
    extract_covariates(pa.nlcd, where="end")%>%
  extract_covariates(DEM.proj, where="end")
  
  random_steps.sf[[i]] <- st_as_sf(random_steps, coords=c("x2_", "y2_"),crs=5070) 
}

#' Create object with all steps for analysis
dat_2 <- do.call(rbind,(random_steps.sf)) %>%
  dplyr::rename("landuse"= Class) %>%
  dplyr::rename("elev"= Layer_1)

#' Convert land cover classifications to a categorical variable and create separate columns
dat_2$Agriculture <- ifelse(dat_2$landuse == "Agriculture", 1, 0)
dat_2$Shrub <- ifelse(dat_2$landuse == "Shrub", 1, 0)
dat_2$Developed <- ifelse(dat_2$landuse == "Developed", 1, 0)
dat_2$Forest <- ifelse(dat_2$landuse == "Forest", 1, 0)
dat_2$NonForest <- ifelse(dat_2$landuse == "NonForest", 1, 0)


#' Change case to numeric
dat_2$case_ <- as.numeric(dat_2$case_)


#' Save RDS file
#saveRDS(dat_2, "Working_SSF_Landcover.RDS")


dat_2 <-readRDS("Data Management/RData/Individual-Specific Movement Process/Working_SSF_Landcover.RDS")
which(dat_2$step_id_==5)

