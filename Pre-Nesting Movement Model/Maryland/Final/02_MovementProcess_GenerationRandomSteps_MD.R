
#'---
#' title: Habitat selection of female wild turkeys during pre-nesting (an SSF analysis)
#' author: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output: RandomSteps_Prep(R Workspace)
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script creates tracks and extracts covariates using the amt package
#' **Last Updated**: 1/25/25

######################################
## Load Packages 


#' Vector of package names
packages <- c("sf",
              "amt",
              "tigris",
              "FedData",
              "dplyr",
              "terra",
              "ggplot2",
              "mapview")

#' Function to load a package or install it if not already installed
load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

#' Apply the function to each package name
lapply(packages, load_packages)


#load("Data Management/RData/Pre-Nesting Movement Model/RData Files/Draft3/20250131_Movebank.RData")

#' Load NLCD Raster
md.nlcd <- terra::rast("Data Management/Rasters/nlcd/md.nlcd.tif")
plot(md.nlcd)


#####################
## Road Rasters 

#' Distance to primary roads raster
primary <- terra::rast("Data Management/Rasters/PA Roads/MdRoadRast.Prim.tiff")

#' Distance to secondary roads Raster
secondary <- terra::rast("Data Management/Rasters/PA Roads/MdRoadRast.Sec.tiff")

##############################
## Create and Save Tracks 

#' I follow the AMT Vignette with multiple individuals
#' https://cran.r-project.org/web/packages/amt/vignettes/p1_getting_started.html
#' 1. Load in data from movebank prep script
#' 2. Make a track using dataset with multiple individuals
#' 3. Nest the track by id
#' 4. Take one individual's data, view its sampling rate and adjust using steps by burst
#' 5. Use the map function to apply the same steps by burst parameters across the marked population

#' Temporary Data
dat <- hens.all %>%
  dplyr::select(BirdID, timestamp,long, lat) 

#' Create a track
#' New amt track function is just amt::track
trk <- amt::make_track(tbl=dat, .x= long, .y=lat, .t=timestamp, id=BirdID,
                             crs= 4326) %>% 
                                    amt::transform_coords(5070)
#' Check
class(trk)

#' Group the track by id and nest the track
trk1 <- trk %>% nest(data = -"id")
trk1

#' get the data for the first animal
x <- trk1$data[[1]]

#' Set criteria for steps 
#' Rate: sampling rate
#' Tolerance: the tolerance of deviations of the sampling rate
#' Steps by burst: Returns NA for zero step lengths and calculates step lengths for points on a path
#' Time of day: Was the fix taken during the day or night? we could filter out night locations
x %>% 
  track_resample(rate = minutes(60), tolerance = minutes(5)) %>%
  amt::steps_by_burst()

#' Check
class(x)

#' Summarize sampling rate
amt::summarize_sampling_rate(x)

#' Apply the same track resampling format to each hen within the dataset
#' This done by using the map function
#' Steps dataframe is the newly created column
trk2 <- trk1 %>%
  dplyr::mutate(steps = map(data, function(x) 
    x %>% amt::track_resample(rate = minutes(60), tolerance = minutes(5))
    %>% amt::steps_by_burst())) 

#' Check
class(trk2)
glimpse(trk2)

#' Create object with all used steps for analysis
stps<- trk2 %>% 
  dplyr::select(id, steps) %>% 
  unnest(cols = steps)

#' Check
glimpse(stps)

#' Cumulative step length during pre-nesting period for each bird
sl.sum <- stps %>%
  group_by(id) %>%
  summarise(sum_sl = sum(sl_, na.rm = TRUE)) %>%
  ungroup()


#' Create random steps
#' Exponential step length used due to issues with formatting gamma
#' Extract covariates at the end of each used and available step
#' Random step lengths drawn from gamma distribution
#' Random turning angles drawn from a vonmises distribution
#' Include_observed: Include all used steps in the analysis
#' Rename and organize columns
names(secondary) <- "secondary"
names(primary) <- "primary"
random_steps<- amt::random_steps(
  stps,
  n_control = 10,
  sl_distr = amt::fit_distr(stps$sl_, "gamma"),
  ta_distr = amt::fit_distr(stps$ta_, "vonmises"),
  include_observed = T) %>%
  amt::extract_covariates(md.nlcd, where ="end") %>%
  amt::extract_covariates(secondary, where ="end") %>%
  amt::extract_covariates(primary, where ="end") %>%
  dplyr::rename("landuse"= Class)
glimpse(random_steps)
table(random_steps$landuse, random_steps$case_)

#' Check
which(random_steps$step_id_==5)
table(random_steps$case_)
table(random_steps$landuse)

################################################################################
##  Consolidate Open Water Category

#' Identify step_ids to remove (where case is 1 and landuse is 'Open Water')
step_ids_to_remove <- unique(random_steps$step_id_[random_steps$case_ == 1 & random_steps$landuse == 'Open Water'])

#' Filter the data frame by removing rows with the identified step_ids
filtered_random_steps <- random_steps[!(random_steps$step_id_ %in% step_ids_to_remove), ]

#' Remove rows where 'landuse' is 'Open Water'
final_filtered_random_steps <- filtered_random_steps[filtered_random_steps$landuse != 'Open Water', ]
table(final_filtered_random_steps$landuse)


################################################################################
## Create dummy landcover variables

#' Create dummy variables for land cover classification in final_filtered_random_steps
final_filtered_random_steps$Agriculture <- ifelse(final_filtered_random_steps$landuse == "Agriculture", 1, 0)
final_filtered_random_steps$Developed <- ifelse(final_filtered_random_steps$landuse == "Developed", 1, 0)
final_filtered_random_steps$Deciduous <- ifelse(final_filtered_random_steps$landuse == "Deciduous Forest", 1, 0)
final_filtered_random_steps$Evergreen <- ifelse(final_filtered_random_steps$landuse == "Evergreen Forest", 1, 0)
final_filtered_random_steps$Mixed <- ifelse(final_filtered_random_steps$landuse == "Mixed Forest", 1, 0)
final_filtered_random_steps$Grassland <- ifelse(final_filtered_random_steps$landuse == "Grassland/Shrub", 1, 0)
final_filtered_random_steps$Wetland <- ifelse(final_filtered_random_steps$landuse == "Wetland", 1, 0)

#' Change case to numeric
final_filtered_random_steps$case_ <- as.numeric(final_filtered_random_steps$case_)


################################################################################
################################################################################
