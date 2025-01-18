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

#######################
## Load Packages

library(tidyverse)
library(terra)
library(mapview)


##############################
## Data Prep- Nest-Level Covs

#' Nest start and end date csv
nests <- read_csv("Data Management/Csvs/processed data/IncubationDates/NestAttempts_allbirds.csv")
nests

#' Nest veg csv
nests.veg <- read_csv("Data Management/Csvs/raw data/nest_vegetation_raw.csv")
nests.veg

#' Filter nests.veg to only include observations that have the same nestid as nests csv
nests.veg.filtered <- nests.veg %>%
  dplyr::rename("nestid"= nestid_v) %>%
  dplyr::filter(nestid%in% nests$nestid) %>%
  dplyr::filter(plottype == "Nest")

#' Merge the filtered nests.veg with nests by nestid
nests <- inner_join(nests, nests.veg.filtered, by = "nestid") 

#' Rename and consolidate columns 
nests <- nests %>%
  dplyr::select(nestid, bandid.x, startI, endI, nestfate, averagemaxvo, percwoody, percgrassforb, stemcount, lat_v, long_v) %>%
  dplyr::rename("Visual Obstruction" = averagemaxvo) %>%
  dplyr::rename("Percent Woody Vegetation" = percwoody) %>%
  dplyr::rename("Percent Grass/Forb" = percgrassforb) %>%
  dplyr::rename("Basal Area" = stemcount) %>%
  dplyr::rename("BandID" = bandid.x)

#' Scale continous predictors
nests.scaled <- nests %>% 
  dplyr::mutate(`Visual Obstruction` = scale(`Visual Obstruction`)) %>%
  dplyr::mutate(`Percent Woody Vegetation` = scale(`Percent Woody Vegetation`)) %>%
  dplyr::mutate(`Percent Grass/Forb` = scale(`Percent Grass/Forb`)) %>%
  dplyr::mutate(`Basal Area` = scale(`Basal Area`))
  
#' Convert scaled columns to numeric
nests.scaled <- nests.scaled %>% 
  dplyr::mutate(`Visual Obstruction` = as.numeric(`Visual Obstruction`)) %>%
  dplyr::mutate(`Percent Woody Vegetation` = as.numeric(`Percent Woody Vegetation`)) %>%
  dplyr::mutate(`Percent Grass/Forb` = as.numeric(`Percent Grass/Forb`)) %>%
  dplyr::mutate(`Basal Area` = as.numeric(`Basal Area`))


###################################
## Data Prep- Lanscape-Level Covs

#' Create sf object and check projection using mapview
nests.sf <- st_as_sf(nests.scaled, coords = c("long_v", "lat_v"), crs = 4326) %>%
  st_transform(5070)
  mapview(nests.sf)

#' Load in land cover information
pa.nlcd <- terra::rast("Data Management/Rasters/nlcd/paNLCD.tiff")
pa.roads.prim <- terra::rast("Data Management/Rasters/PA Roads/paroadrast.prim.tiff")
pa.roads.sec <- terra::rast("Data Management/Rasters/PA Roads/paroadrast.sec.tiff")

#' Extract landcover point value for each nest
nests.landcov <- terra::extract(pa.nlcd, nests.sf) %>%
  dplyr::rename("landuse" = Class)

#' Create dummy variables for land cover classification
nests.scaled$Agriculture <- ifelse(nests.landcov$landuse == "Agriculture", 1, 0)
nests.scaled$Developed <- ifelse(nests.landcov$landuse == "Developed", 1, 0)
nests.scaled$Deciduous <- ifelse(nests.landcov$landuse == "Deciduous Forest", 1, 0)
nests.scaled$Evergreen <- ifelse(nests.landcov$landuse == "Evergreen Forest", 1, 0)
nests.scaled$Mixed <- ifelse(nests.landcov$landuse == "Mixed Forest", 1, 0)
nests.scaled$Grassland <- ifelse(nests.landcov$landuse == "Grassland/Shrub", 1, 0)
nests.scaled$Other <- ifelse(nests.landcov$landuse == "Other", 1, 0)

#' Extract distance from primary and secondary road structures
nests.prim.roads <- terra::extract(pa.roads.prim, nests.sf) %>%
  dplyr::rename("Primary" = layer)

nests.sec.roads <- terra::extract(pa.roads.sec, nests.sf) %>%
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


####################################
## Data Prep- Weather Covariates




####################################
## Data Prep- Individual Covariates

#' Nest Incubation Dates- Julian Date
nests.scaled <- nests.scaled %>%
  dplyr::mutate("Nest Incubation Date" = lubridate::yday(startI))


####################################
## Data Prep- Encounter Histories

#' Create a sequence of all unique dates between the min start and max end date
inc.dates <- sort(unique(c(nests$startI, nests$endI)))
inc.dates <- seq(inc.dates[1], inc.dates[length(inc.dates)], by = 1)

#' Calculate the number of encounter days
occs <- length(inc.dates)

#' Initialize variables
first <- last <- array(NA, dim = nrow(nests))
surv.caps <- matrix(data = NA, nrow = nrow(nests), ncol = occs)

#' Loop through each row of the dataset
for(i in 1:nrow(nests)){ 
  first[i] <- which(inc.dates == nests$startI[i]) 
  for (j in first[i]:occs){  
    if (surv.caps[i, j]==0) {
      break
    }
  } 
  
  last[i]<-ifelse(surv.caps[i,occs] %in% c(1),occs,which(surv.caps[i, ] == 0))
}


############################################
## Compile Parameters into matrix format

#' Function to fill NA values with 0 throughout dataframe
fill_NA_with_value <- function(df, value = 0) {
  #' Ensure that all columns are of numeric type to handle the replacement
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      #' Replace NA values with the specified value
      replace(col, is.na(col), value)
    } else {
      #' Return the column unchanged if it's not numeric
      col
    }
  })
  return(df)
}

#' Apply the function to fill NA values with 0
nests.ready <- fill_NA_with_value(nests.scaled)
summary(nests.ready)

#' Model parameters
X <- cbind(
  rep(1, nrow(nests.ready)),               # Intercept (1)
  nests.ready$`Visual Obstruction`,        # Visual Obstruction
  nests.ready$`Percent Woody Vegetation`,  # Percent Woody Vegetation
  nests.ready$`Percent Grass/Forb`,        # Percent Grass Forb
  nests.ready$`Basal Area`,                # Basal Area
  nests.ready$Primary,                     # Distance to Primary Road
  nests.ready$Secondary,                   # Distance to Secondary Road
  nests.ready$Mixed,                       # Mixed Forest
  nests.ready$Evergreen,                   # Evergreen Forest
  nests.ready$Deciduous,                   # Deciduous Forest
  nests.ready$Agriculture,                 # Agriculture
  nests.ready$Grassland,                   # Grassland  
  nests.ready$`Nest Incubation Date`       # Nest Incubation Date
  )       

################################################################################
################################################################################

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



#####################
## Load Packages

library(nimble)
library(MCMCvis)


################################
## Nimble Model

nestsuccess <- nimbleCode({
  
  for (i in 1:n.ind) {
   for(t in first[i]:last[i]) {
       y[i, t] ~ dbern(mu[i, t])  
       cloglog(mu[i, t]) <- inprod(beta[1:J], X[i, 1:J])  
    }
  }
  
  for (j in 1:J){
    beta[j] ~ dnorm(0, sd = sqrt(1 / 0.0001))  
  }
  
  phi ~ dbeta(1, 1)
  
})

#' Constants
Consts <- list(
  n.ind = nrow(nests.ready),   # Number of rows of nests
  X = X,                       # Covariate matrix (temperature)
  J = ncol(X),                 # Columns of X
  first = first,               # First incubation day for each nest
  last = last                  # Last incubation day for each nest
)

#' Data (survival data: 1 = survived, 0 = failed)
Data <- list(
  y = surv.caps               # Survival data
)

#' Initial values for parameters
Inits <- list(
  beta = rep(0, Consts$J),    # Initial values for beta (length should match the number of covariates in X)
  phi = 0.1                   # Initial value for daily survival probability
)

start <- Sys.time()

#' Run MCMC sampling with Nimble
nimbleMCMC_samples <- nimbleMCMC(
  code = nestsuccess,          # Model code
  constants = Consts,          # Constants defined earlier
  inits = Inits,               # Initial values for the parameters
  data = Data,                 # Data (survival data)
  nburnin = 10000,             # Burn-in period
  niter = 40000,               # Number of iterations
  thin = 1                     # Thinning factor
)

#' View traceplots
MCMCtrace(nimbleMCMC_samples, pdf = FALSE)
