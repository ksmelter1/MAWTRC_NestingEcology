#'---
#' title: 2D Study Area Analysis
#' author: "K. Smelter
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#'---
#'  
#' **Purpose**: This script analyzes Pennsylvania 2D Data 
#' **Last Updated**: 6/12/2025
#' 
#' 

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
              "rasterVis",
              "move2",
              "amt",
              "purrr",
              "lubridate",
              "nimble",
              "ggplot2",
              "GGally",
              "daymetr",
              "matrixStats",
              "MCMCvis")

load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

lapply(packages, load_packages)


################################################################################
## Pre-Nesting Movement Model

# This section fits a pre-nesting movement model for WMU 2D

################################################################################
## Data Management

# Read in nests csv
pa.nests <- read_csv("Data Management/Csvs/Processed/Nests/Nests/Pennsylvania/2025_CleanedNests_2022_2023_PA.csv")
pa.nests

# Create clutch size column and remove unnecessary column
# Clutch size is a minimum count 
pa.nests <- pa.nests %>%
  dplyr::mutate(clutchsize = rowSums(select(., EggsHatch, EggsDestroyed, EggsUnhatch, EggsDepred), na.rm = TRUE)) 
glimpse(pa.nests)

#' Subset nesting data for 2D in year 2022
pa.nests.2D <- dplyr::filter(pa.nests, WMU =="2D")%>%
  dplyr::select(BandID, CheckDate, NestID, WMU, clutchsize)

#' Csv from incubation start and end script
nests.inc <- read_csv("Data Management/Csvs/Processed/Incubation Dates/Pennsylvania/20250131_NestAttempts_allbirds_PA.csv")
nests.inc

#' Sample from known fate model to ensure consistency
pa.sample <- read_csv("Samples/Pennsylvania/NestingSample_PA.csv")
pa.sample
nests.inc <- right_join(nests.inc, pa.sample)

#' Merge pa.nests.2D and nests.inc, only keep nests that exist in both pa.nests.2D and nests.inc
pa.nests.2D1 <- dplyr::inner_join(pa.nests.2D, nests.inc, by = "NestID") %>%
  dplyr::select(-CheckDate.y) %>%
  dplyr::rename("CheckDate" = CheckDate.x)
glimpse(pa.nests.2D1)

#' Run this twice down to as.Date
#' Iterate over the rows and subtract clutchsize days
for (i in 1:nrow(pa.nests.2D1)) {
  clutchsize <- pa.nests.2D1$clutchsize[i]  
  startI <- pa.nests.2D1$startI[i]  
  
  #' Create the enddate by subtracting clutchsize (in days) from startI
  #' 5 day buffer
  pa.nests.2D1$enddate[i] <- startI - clutchsize - days(5)
  
  #' Create the startdate by subtracting 14 days from the enddate
  pa.nests.2D1$startdate[i] <- pa.nests.2D1$enddate[i]- days(14)
  
  #' Check the calculated dates
  print(paste("Row", i, ": Startdate =", pa.nests.2D1$startdate[i], ", Enddate =", pa.nests.2D1$enddate[i]))
}

#' Convert 'startdate' and 'enddate' columns to Date format
pa.nests.2D1$startdate <- as.Date(pa.nests.2D1$startdate)
pa.nests.2D1$enddate <- as.Date(pa.nests.2D1$enddate)

#' View the data frame structure
glimpse(pa.nests.2D1)


################################################################################
## Connect to Movebank

login <- movebank_store_credentials(username = "Kyle.Smelter",
                                    password="Rayshawks5!",
                                    key="Kyle",
                                    force= T)


################################################################################
## WMU 2D

unique.ID.2D<-unique(pa.nests.2D1$NestID)

for (j in 1:length(unique.ID.2D)){
  tmp.subset.2D<-pa.nests.2D1[which(pa.nests.2D1$NestID==unique.ID.2D[j]),]
  tmp.subset.2D$TrackID<-paste(unique.ID.2D[j],seq(1,nrow(tmp.subset.2D),1),sep="_")
  
  for(i in 1:nrow(tmp.subset.2D)){
    BirdID<- as.character(tmp.subset.2D[i,1])
    EndDate <- gsub("\\D","", tmp.subset.2D$enddate[i]) 
    
    StartDate <- gsub("\\D","", tmp.subset.2D$startdate[i]) 
    
    Year <- lubridate::year(tmp.subset.2D$startI[i])
    
    
    dat.2D<- movebank_download_study(study ="Wild Turkey Pennsylvania WMU 2D", 
                                     login = login,
                                     individual_local_identifier= BirdID,
                                     timestamp_start= StartDate,
                                     timestamp_end= EndDate,
                                     removeDuplicatedTimestamps=T)
    mt_track_id(dat.2D)<-rep(tmp.subset.2D$TrackID[i],nrow(dat.2D))
    
    if(exists("full_all_2D")){ 
      full_all_2D <- rbind(full_all_2D, dat.2D)
      
      
    }else{
      full_all_2D <- dat.2D
    }
  }
}


################################################################################
## Process GPS Data for Pre-Nesting Movement Model

# Function to process GPS data
# Convert from move2 object to dataframe
# Rename the individual_local_identifier
# Extract lattitude and longitude
process_gps_data <- function(gps_data) {
  # Ensure the input is a dataframe
  gps_data <- as.data.frame(gps_data)
  
  # Rename the individual identifier to BirdID
  gps_data <- gps_data %>%
    dplyr::rename("BirdID" = individual_local_identifier)
  
  # Extract longitude and latitude from geometry column
  gps_data <- gps_data %>%
    mutate(long = unlist(map(geometry, 1)),
           lat = unlist(map(geometry, 2))) %>%
    dplyr::select(BirdID, timestamp, long, lat)
  
  return(gps_data)
}

hens_2D <- process_gps_data(full_all_2D)


################################################################################
## Creating Tracks for Pre-Nesting Movement Model

#' Load NLCD Raster
pa.nlcd <- terra::rast("Data Management/Rasters/NLCD/Pennsylvania/pa.nlcd.tif")
plot(pa.nlcd)

#' Distance to primary roads raster
primary <- terra::rast("Data Management/Rasters/Roads/Pennsylvania/PaRoadRast.Prim.tiff")

#' Distance to secondary roads Raster
secondary <- terra::rast("Data Management/Rasters/Roads/Pennsylvania/PaRoadRast.Sec.tiff")

################################################################################
## Create and Save Tracks 

#' I follow the AMT Vignette with multiple individuals
#' https://cran.r-project.org/web/packages/amt/vignettes/p1_getting_started.html
#' 1. Load in data from movebank prep script
#' 2. Make a track using dataset with multiple individuals
#' 3. Nest the track by id
#' 4. Take one individual's data, view its sampling rate and adjust using steps by burst
#' 5. Use the map function to apply the same steps by burst parameters across the marked population

#' Temporary Data
dat <- hens_2D %>%
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
  amt::extract_covariates(pa.nlcd, where ="end") %>%
  amt::extract_covariates(secondary, where ="end") %>%
  amt::extract_covariates(primary, where ="end") %>%
  dplyr::rename("landuse"= Class)
glimpse(random_steps)

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
final_filtered_random_steps$Pasture <- ifelse(final_filtered_random_steps$landuse == "Pasture", 1, 0)
final_filtered_random_steps$Crop <- ifelse(final_filtered_random_steps$landuse == "Crop", 1, 0)
final_filtered_random_steps$Developed <- ifelse(final_filtered_random_steps$landuse == "Developed", 1, 0)
final_filtered_random_steps$Deciduous <- ifelse(final_filtered_random_steps$landuse == "Deciduous Forest", 1, 0)
final_filtered_random_steps$Evergreen <- ifelse(final_filtered_random_steps$landuse == "Evergreen Forest", 1, 0)
final_filtered_random_steps$Mixed <- ifelse(final_filtered_random_steps$landuse == "Mixed Forest", 1, 0)
final_filtered_random_steps$Grassland <- ifelse(final_filtered_random_steps$landuse == "Grassland/Shrub", 1, 0)
final_filtered_random_steps$Wetland <- ifelse(final_filtered_random_steps$landuse == "Wetland", 1, 0)

#' Change case to numeric
final_filtered_random_steps$case_ <- as.numeric(final_filtered_random_steps$case_)


################################################################################
## Final Data Prep for Pre-Nesting Movement Model

#' Function to prepare step-selection function data 
#' Scales data
#' Converts back to numeric
#' Creates animal and stratum-id columns
#' Order by NAID
prepare_step_selection_data <- function(data) {
  
  # Scale continuous predictors and convert to numeric
  data <- data %>%
    mutate(
      secondary = scale(secondary),
      primary = scale(primary),
      secondary = as.numeric(secondary),
      primary = as.numeric(primary)
    )
  
  # Clean and convert ID to numeric
  data$id <- gsub("_", "", data$id)
  data$id <- as.numeric(data$id)
  
  # Create unique NA_ID using id and step_id_
  data$NA_ID <- paste(data$id, data$step_id_, sep = "_")
  data$NA_ID <- gsub("_", "", data$NA_ID)
  data$NA_ID <- as.numeric(data$NA_ID)
  
  # Add numerical ANIMAL_ID
  data$ANIMAL_ID <- as.numeric(as.factor(data$id))
  
  # Create stratum ID map
  d.map <- data.frame(
    NA_ID = unique(data$NA_ID),
    str_ID = seq_along(unique(data$NA_ID))
  )
  
  # Match stratum IDs back to main data
  data$str_ID <- d.map[match(data$NA_ID, d.map$NA_ID), "str_ID"]
  
  # Order by stratum ID
  data <- data[order(data$str_ID), ]
  
  return(data)
}

dat_2.ready <- prepare_step_selection_data(final_filtered_random_steps)


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

#' Apply the function to fill NA values with 0 (Mean for scaled continous variables and wasn't used for categorical)
dat_2.ready <- fill_NA_with_value(dat_2.ready)

#' Check
which(dat_2.ready$step_id_==3)
table(dat_2.ready$case_)
class(dat_2.ready)
summary(dat_2.ready)

#' Model parameters
X <- cbind(
  rep(1, nrow(dat_2.ready)),   # Intercept (1)
  dat_2.ready$primary,         # Distance to Primary Road
  dat_2.ready$secondary,       # Distance to Secondary Road
  dat_2.ready$Mixed,           # Mixed Forest
  dat_2.ready$Evergreen,       # Evergreen Forest
  dat_2.ready$Developed,       # Developed
  dat_2.ready$Pasture,         # Pasture
  dat_2.ready$Crop,            # Crop
  dat_2.ready$Grassland,       # Grassland    
  dat_2.ready$Wetland)         # Wetland

#' Keep only design matrix (X) and data columns
#' Cuts down on processing
rm(list = setdiff(ls(), c("dat_2.ready", "X", "pa.nests.2D1")))


################################################################################
## Nimble  Pre-Nesting Movement Model

#' Step-selection function in a Bayesian Framework
#' Conditional logistic regression
#' Set prior for stratum-specific intercept to 10^6
#' Set prior for betas to 0,1 which is a standard uninformative prior for a logistic regression

movemodel <- nimbleCode({
  for (i in 1:I){
    use[i] ~ dpois(lambda[i])  
    
    log(lambda[i]) <- inprod(beta[1:J], X[i, 1:J]) + alpha[str_ID[i]]  
  }
  
  #' Priors for beta 
  for (j in 1:J){
    beta[j] ~ dnorm(0, sd = 1)  
  }
  
  #' Priors for alpha 
  for (k in 1:K){
    alpha[k] ~ dnorm(0, sd = sqrt(1 / 0.000001))  
  }
})

#' Constants 
Consts <- list(
  str_ID = as.numeric(dat_2.ready$str_ID), 
  J = ncol(X),                             
  I = nrow(dat_2.ready),                   
  K = length(unique(dat_2.ready$str_ID))   
)

#' Data
Data <- list(
  X = X,                                   
  use = dat_2.ready$case_                  
)

#' Initial values 
Inits <- list(
  beta = rep(0, ncol(X)),                  
  alpha = rep(0, length(unique(dat_2.ready$str_ID)))  
)

start <- Sys.time()

#' Run MCMC sampling with Nimble
nimbleMCMC_samples <- nimbleMCMC(
  code = movemodel, 
  constants = Consts, 
  inits = Inits,
  data = Data,
  nburnin = 10000, 
  niter = 40000,
  thin = 1
)

end <- Sys.time()

#' Print the means and standard deviations of the posterior samples for the beta coefficients
colMeans(nimbleMCMC_samples[, 12959:12968])
colSds(nimbleMCMC_samples[, 12959:12968])

#' View traceplots
MCMCtrace(nimbleMCMC_samples[, 12959:12968], pdf = F)

#' Extract the posterior samples for the 'beta' parameters (columns 219 to 230)
beta_samples <- nimbleMCMC_samples[, 12959:12968]

# Convert mcmc.list to matrix
samples_matrix <- as.matrix(beta_samples)

# Convert the matrix to a data frame
samples_df <- as.data.frame(samples_matrix) 

#' Create a vector of new names
new_names <- c("Intercept", 
               "Distance to Primary Road",
               "Distance to Secondary Road", 
               "Mixed Forest",
               "Evergreen Forest",
               "Developed",
               "Pasture",
               "Crop",
               "Grassland/Shrub",
               "Wetland")

#' Assign the variable names to the columns of the beta_samples
colnames(samples_df) <- new_names

#' View the renamed data frame
head(samples_df)


################################################################################
## Visualize Pre-Nesting Outputs

################################################################################
## Visualize Pre-Nesting Outputs

# Reshape to long format
samples_long <- samples_df %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "estimate")

# Compute mean and 90% credible intervals
mean_estimates <- samples_long %>%
  group_by(parameter) %>%
  summarise(
    mean_estimate = mean(estimate),
    lower = quantile(estimate, 0.05),
    upper = quantile(estimate, 0.95),
    .groups = 'drop'
  ) %>%
  filter(!parameter %in% c("Intercept", "Elevation")) 

# Assign Scales
mean_estimates1 <- mean_estimates %>%
  mutate(Scale = case_when(
    parameter %in% c("Percent Grass/Forb", 
                     "Percent Woody Vegetation", 
                     "Visual Obstruction",
                     "Woody Stem Count",
                     "Basal Area",
                     "Percent Fern") ~ "Nest",
    parameter %in% c( "Grassland/Shrub",
                      "Mixed Forest",
                      "Evergreen Forest",
                      "Developed",
                      "Distance to Primary Road",
                      "Distance to Secondary Road",
                      "Pasture",
                      "Crop",
                      "Wetland"
    ) ~ "Landscape",
    parameter %in% c("Minimum Temperature",
                     "Precipitation",
                     "Precipitation * Minimum Temperature"
    ) ~ "Weather",
    TRUE ~ "Individual"
  ))

# Define factor levels for display order
mean_estimates1$parameter <- factor(mean_estimates1$parameter,
                                    levels = c("Wetland",
                                               "Grassland/Shrub",
                                               "Mixed Forest",
                                               "Evergreen Forest",
                                               "Developed",
                                               "Crop",
                                               "Pasture",
                                               "Distance to Primary Road",
                                               "Distance to Secondary Road"))

# Create ggplot
p1 <- ggplot(mean_estimates1, aes(x = parameter, y = mean_estimate, color = Scale, shape = Scale)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Parameter", y = "Selection Relative to Deciduous Forest") +
  theme_minimal() +
  coord_flip() +
  scale_color_manual(values = c("Landscape" = "#D65F5F")) +
  scale_shape_manual(values = c("Landscape" = 16)) +
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 10), hjust = 0.45),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  )

p1

mean_estimates1

save(mean_estimates1,
     samples_df,
     file = "MAWTRC Nesting Ecology Manuscript/Figures/RData/Pre-Nesting Movement/Study Area Analysis/pre.2D.RData")


################################################################################
###############################################################################X



################################################################################
## Nest-Site Selection

# This section creates a nest-site selection model for WMU 2D

################################################################################
## Organize 2D Nesting Data

#' Read in Pennsylvania NLCD
pa.nlcd <- terra::rast("Data Management/Rasters/NLCD/Pennsylvania/pa.nlcd.tif")

# Read in nests csv
pa.nests <- read_csv("Data Management/Csvs/Processed/Nests/Vegetation Surveys/Pennsylvania/20250121_CleanedNestsVeg_2022_2023.csv")
pa.nests

#' Subset nesting data for 2D in year 2022
pa.nests.2D <- dplyr::filter(pa.nests, WMU =="2D")

#' Csv from incubation start and end script
nests.inc <- read_csv("Data Management/Csvs/Processed/Incubation Dates/Pennsylvania/20250131_NestAttempts_allbirds_PA.csv")
nests.inc

#' Sample from known fate model to ensure consistency
pa.sample <- read_csv("Samples/Pennsylvania/NestingSample_PA.csv")
pa.sample
nests.inc <- right_join(nests.inc, pa.sample)

#' Merge the filtered nests.veg with nests by nestid
pa.nests<- inner_join(pa.nests.2D, nests.inc, by = "NestID") 

pa.nests.sf <- pa.nests %>%
  st_as_sf(., coords = c("Long", "Lat"), crs = 4326) %>%
  st_transform(5070)
mapview(pa.nests.sf)

#' Extract point value at each nest
landcov <-terra::extract(pa.nlcd, pa.nests.sf)

#' Bind columns together
pa.nests.landcov <- cbind(pa.nests.sf, landcov) %>%
  dplyr::rename("landuse" = Class)


################################################################################
## Extract Distances to Roads- Pennsylvania

#' PASDA Raster
pasda.dist.prim <- terra::rast("Data Management/Rasters/Roads/Pennsylvania/PaRoadRast.Prim.tiff")
pasda.dist.sec <- terra::rast("Data Management/Rasters/Roads/Pennsylvania/PaRoadRast.Sec.tiff")

#' Read in nest data with landcover
pa.nests.landcov.sf <- pa.nests.landcov

#' Extract point value at each nest
dist.prim.out <-terra::extract(pasda.dist.prim, pa.nests.landcov.sf) %>%
  dplyr::select(-ID)%>%
  dplyr::rename("primary" = layer)
dist.sec.out <-terra::extract(pasda.dist.sec, pa.nests.landcov.sf) %>%
  dplyr::select(-ID) %>%
  dplyr::rename("secondary" = layer)


pa.nests.landcov.sf.roads <- cbind(pa.nests.landcov.sf, dist.prim.out)
pa.nests.landcov.sf.roads <- cbind(pa.nests.landcov.sf.roads, dist.sec.out) %>%
  st_drop_geometry()


pa.nests.covs <- left_join(pa.nests.landcov.sf, pa.nests.landcov.sf.roads) %>%
  dplyr::rename("BandID" = BandID.x) %>%
  dplyr::select(landuse, 
                Case,
                NestID, 
                BandID, 
                PercGrassForb,
                PercWoody, 
                PercLitter,
                PercBare,
                PercBoulder,
                AvgMaxVO, 
                primary,
                secondary,
                StemCount,
                WoodyType1,
                AvgVO,
                PercFern, 
                GuardHt,
                HtWoody,
                HtFern,
                HtGrassForb,
                HtLitter, 
                HtBoulder,
                PlotType) 


#' Convert land cover classifications to a categorical variable and create separate columns
pa.nests.covs$Pasture <- ifelse(pa.nests.covs$landuse == "Pasture", 1, 0)
pa.nests.covs$Crop <- ifelse(pa.nests.covs$landuse == "Crop", 1, 0)
pa.nests.covs$Developed <- ifelse(pa.nests.covs$landuse == "Developed", 1, 0)
pa.nests.covs$Deciduous <- ifelse(pa.nests.covs$landuse == "Deciduous Forest", 1, 0)
pa.nests.covs$Evergreen <- ifelse(pa.nests.covs$landuse == "Evergreen Forest", 1, 0)
pa.nests.covs$Mixed <- ifelse(pa.nests.covs$landuse == "Mixed Forest", 1, 0)
pa.nests.covs$Grassland <- ifelse(pa.nests.covs$landuse == "Grassland/Shrub", 1, 0)
pa.nests.covs$Wetland <- ifelse(pa.nests.covs$landuse == "Wetland", 1, 0)
pa.nests.covs$Water <- ifelse(pa.nests.covs$landuse == "Open Water", 1, 0)


################################################################################
## Fit Nest-Site Selection Model

#' This sample matches the Pre-Nesting and known fate model
Sample_PA <- read_csv("Samples/Pennsylvania/NestingSample_PA.csv")
Sample_PA

#' Drop geometry column and create unique identifier for NestID
pa.nests.covs <- sf::st_drop_geometry(pa.nests.covs)
pa.nests.covs$NestID1 <- paste(pa.nests.covs$NestID, pa.nests.covs$PlotType, sep = "_")

#' Read in basal area data
basal_summary <- read_csv("Data Management/Csvs/Processed/Covariates/Pennsylvania/BasalArea.csv")

#' Keep only samples that exist in Sample_PA in pa.nests.covs
pa.nests.covs.1 <- right_join(pa.nests.covs, Sample_PA)
pa.nests.covs.2 <- right_join(basal_summary, pa.nests.covs.1) 

pa.nests.covs.2 <- pa.nests.covs.2[-c(219:332), ]

#' Zero values represent sites with no trees for basal area
pa.nests.covs.2 <- pa.nests.covs.2 %>%
  replace_na(list(Basal = 0))

#' Rename object
nest.data <- pa.nests.covs.2 
str(nest.data)
summary(nest.data)

#' Subset dataframe and rename columns 
nest.data <- nest.data %>%
  dplyr::select(NestID, BandID, , PercGrassForb, PercWoody, AvgMaxVO,
                Case, Developed, Deciduous, Mixed, Evergreen, Pasture, Crop,
                primary, secondary, StemCount, Grassland, AvgVO, Water, Wetland, 
                PercFern, Basal) %>%
  dplyr::rename("Primary" = primary) %>%
  dplyr::rename("Secondary" = secondary)
nest.data

#' Switch coding to UTF-8
nest.data <- nest.data %>% 
  dplyr::mutate(across(everything(), ~ iconv(., to = "UTF-8")))

#' Change variables to be scaled to numeric
nest.data <- nest.data %>%
  dplyr::mutate(PercWoody = as.numeric(PercWoody)) %>%
  dplyr::mutate(PercGrassForb = as.numeric(PercGrassForb)) %>%
  dplyr::mutate(AvgMaxVO = as.numeric(AvgMaxVO)) %>%
  dplyr::mutate(Primary = as.numeric(Primary)) %>%
  dplyr::mutate(Secondary = as.numeric(Secondary)) %>%
  dplyr::mutate(PercFern = as.numeric(PercFern)) %>%
  dplyr::mutate(AvgVO = as.numeric(AvgVO)) %>%
  dplyr::mutate(StemCount = as.numeric(StemCount)) %>%
  dplyr::mutate(Basal = as.numeric(Basal)) 
glimpse(nest.data)
str(nest.data)

#' Scale continous predictor variables
nest.data <- nest.data %>%
  dplyr::mutate(StemCount = scale(StemCount)) %>%
  dplyr::mutate(PercWoody = scale( PercWoody)) %>%
  dplyr::mutate(PercGrassForb = scale(PercGrassForb)) %>%
  dplyr::mutate(AvgMaxVO = scale(AvgMaxVO)) %>%
  dplyr::mutate(Primary = scale(Primary)) %>%
  dplyr::mutate(Secondary = scale(Secondary)) %>%
  dplyr::mutate(PercFern = scale(PercFern)) %>%
  dplyr::mutate(AvgVO = scale(AvgVO)) %>%
  dplyr::mutate(Basal = scale(Basal)) 
str(nest.data)
glimpse(nest.data)

#' Convert all variables to numeric including all scaled predictors
nest.data <- nest.data %>%
  dplyr::mutate(PercWoody = as.numeric(PercWoody)) %>%
  dplyr::mutate(PercGrassForb = as.numeric(PercGrassForb)) %>%
  dplyr::mutate(AvgMaxVO = as.numeric(AvgMaxVO)) %>%
  dplyr::mutate(Primary = as.numeric(Primary)) %>%
  dplyr::mutate(Secondary = as.numeric(Secondary)) %>%
  dplyr::mutate(PercFern = as.numeric(PercFern)) %>%
  dplyr::mutate(AvgVO = as.numeric(AvgVO)) %>%
  dplyr::mutate(Developed = as.numeric(Developed)) %>%
  dplyr::mutate(Pasture = as.numeric(Pasture)) %>%
  dplyr::mutate(Crop = as.numeric(Crop)) %>%
  dplyr::mutate(Deciduous = as.numeric(Deciduous)) %>%
  dplyr::mutate(Mixed = as.numeric(Mixed)) %>%
  dplyr::mutate(Evergreen = as.numeric(Evergreen)) %>%
  dplyr::mutate(Grassland = as.numeric(Grassland)) %>%
  dplyr::mutate(Wetland = as.numeric(Wetland)) %>%
  dplyr::mutate(Water = as.numeric(Water)) %>%
  dplyr::mutate(StemCount = as.numeric(StemCount)) %>%
  dplyr::mutate(Basal = as.numeric(Basal))
glimpse(nest.data)
str(nest.data)
summary(nest.data)

#' Drop geometry column and change case to numeric
#' 1 the nest was used, 0 it was available to the hen 
nest.data <- nest.data %>%
  st_drop_geometry() %>%
  dplyr::mutate(Case = as.numeric(Case))

#' Order dataframe by NestID
nest.data <- nest.data[order(nest.data$NestID),]

#' Change Nest_ID_V to numeric 
#' First must remove underlines
nest.data$NestID <- gsub("_", "", nest.data$NestID)
nest.data$NestID<- as.numeric(nest.data$NestID)

#' Create str_id column 
#' This allows the loop to iterate through the steps associated with each bird 
#' cur_group_id() gives a unique numeric identifier for the current group.
nest.data <- nest.data %>%
  group_by(NestID) %>%
  mutate(str_ID=cur_group_id())

#' Check data
str(nest.data)
glimpse(nest.data)
summary(nest.data)

#' Function to replace all NAs with 0 in a dataframe
#' Since we scaled the data, the mean is zero
fill_na_with_zero <- function(df) {
  df[is.na(df)] <- 0
  return(df)
}

#' Apply the function
nest.data.ready <- fill_na_with_zero(nest.data)

#' Check
summary(nest.data.ready)
glimpse(nest.data.ready)
cor(nest.data.ready$Primary, nest.data.ready$Secondary)


################################################################################
## Nimble Nest-Site Selection Model

#' Conditional logistic regression
#' Prior for stratum-specific intercept variance set to 10^6
#' Prior for betas set to 0,1

nestmodel<-nimbleCode({
  for (i in 1:I){
    use[i]~dpois(lambda[i])
    
    log(lambda[i])<-inprod(beta[1:J],X[i,1:J]) + alpha[str_ID[i]] 
  }
  
  ### Priors for betas ###
  for(j in 1:J){
    beta[j]~dnorm(0, sd= 1)
  }
  
  ### Priors for random effect ###
  for(k in 1:K){
    alpha[k]~dnorm(0,sd=sqrt(1/0.000001))
  }
}
)

#' Model parameters
X<- cbind(
  rep(1, nrow(nest.data.ready)),   # Intercept (1)
  nest.data.ready$Primary,         # Distance to Primary Road
  nest.data.ready$Secondary,       # Distance to Secondary Road
  nest.data.ready$Mixed,           # Mixed Forest
  nest.data.ready$Evergreen,       # Evergreen Forest
  nest.data.ready$Developed,       # Developed
  nest.data.ready$Pasture,         # Pasture
  nest.data.ready$Crop,            # Crop
  nest.data.ready$Wetland,         # Wetland
  nest.data.ready$PercGrassForb,   # Percent Grass/Forb
  nest.data.ready$PercWoody,       # Percent Woody
  nest.data.ready$AvgVO,           # Visual Obstruction
  nest.data.ready$PercFern,        # Percent Fern
  nest.data.ready$Grassland,       # Grassland/Shrub
  nest.data.ready$Basal,           # Basal Area
  nest.data.ready$StemCount        # Woody Stem Count
)

ggpairs(
  as.data.frame(X), 
  upper = list(continuous = wrap("cor", size = 5)),
  diag = list(continuous = "barDiag"),
  lower = list(continuous = "points")
)

Consts <- list(str_ID = as.numeric(nest.data.ready$str_ID),
               J=ncol(X),
               I = nrow(nest.data.ready),
               K =length(unique(nest.data.ready$str_ID)))

#' Data list for nimble
Data<-list(X=X,use = nest.data.ready$Case)

Inits<-list(beta=rep(0,ncol(X)),alpha=rep(0,length(unique(nest.data.ready$str_ID))))

nimbleMCMC_samples <- nimbleMCMC(code = nestmodel, 
                                 constants = Consts, 
                                 inits= Inits,
                                 data = Data,
                                 nburnin = 10000,
                                 niter = 20000,
                                 thin = 3,
                                 nchains = 1)

#' View traceplots
MCMCtrace(nimbleMCMC_samples, pdf = F)
colMeans(nimbleMCMC_samples[,45:60])
colSds(nimbleMCMC_samples[,45:60])

#' Extract the posterior samples for the 'beta' parameters (columns 219 to 230)
beta_samples <- nimbleMCMC_samples[, 45:60]

# Convert mcmc.list to matrix
samples_matrix <- as.matrix(beta_samples)

# Convert the matrix to a data frame
samples_df <- as.data.frame(samples_matrix) 

#' Create a vector of new names
new_names <- c("Intercept", 
               "Distance to Primary Road",
               "Distance to Secondary Road",
               "Mixed Forest",
               "Evergreen Forest", 
               "Developed",
               "Pasture", 
               "Crop",
               "Wetland",
               "Percent Grass/Forb",
               "Percent Woody Vegetation", 
               "Visual Obstruction",
               "Percent Fern", 
               "Grassland/Shrub",
               "Basal Area",
               "Woody Stem Count")

#' Assign the variable names to the columns of the beta_samples
colnames(samples_df) <- new_names

#' View the renamed data frame
head(samples_df)


################################################################################
## Visualize Nest-Site Selection Results for Pennsylvania WMU 4D

# Reshape to long format
samples_long <- samples_df %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "estimate")

# Compute mean and 90% credible intervals
mean_estimates <- samples_long %>%
  group_by(parameter) %>%
  summarise(
    mean_estimate = mean(estimate),
    lower = quantile(estimate, 0.05),
    upper = quantile(estimate, 0.95),
    .groups = 'drop'
  ) %>%
  filter(!parameter %in% c("Intercept", "Elevation")) 

# Assign nest and landscape scales
mean_estimates2 <- mean_estimates %>%
  mutate(Scale = case_when(
    parameter %in% c("Percent Grass/Forb", 
                     "Percent Woody Vegetation", 
                     "Visual Obstruction",
                     "Woody Stem Count",
                     "Basal Area",
                     "Percent Fern") ~ "Nest",
    parameter %in% c( "Grassland/Shrub",
                      "Mixed Forest",
                      "Evergreen Forest",
                      "Developed",
                      "Distance to Primary Road",
                      "Distance to Secondary Road",
                      "Pasture",
                      "Crop",
                      "Wetland"
    ) ~ "Landscape",
    parameter %in% c("Minimum Temperature",
                     "Precipitation",
                     "Precipitation * Minimum Temperature"
    ) ~ "Weather",
    TRUE ~ "Individual"
  ))

# Define factor levels for display order
mean_estimates2$parameter <- factor(mean_estimates2$parameter,
                                    levels = c("Wetland",
                                               "Grassland/Shrub",
                                               "Mixed Forest",
                                               "Evergreen Forest",
                                               "Developed",
                                               "Pasture",
                                               "Crop",
                                               "Distance to Secondary Road",
                                               "Distance to Primary Road",
                                               "Visual Obstruction",
                                               "Basal Area",
                                               "Percent Woody Vegetation", 
                                               "Percent Grass/Forb",
                                               "Percent Fern",
                                               "Woody Stem Count"))

# Create ggplot
p <- ggplot(mean_estimates2, aes(x = parameter, y = mean_estimate, color = Scale, shape = Scale)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Parameter", y = "Selection Relative to Deciduous Forest") +
  theme_minimal() +
  coord_flip() +
  scale_color_manual(values = c("Individual" = "#DA6509",
                                "Nest" = "#A44200",
                                "Landscape" = "#D65F5F",
                                "Weather" = "#197278")) +  
  scale_shape_manual(values = c("Individual" = 15, 
                                "Nest" = 17,
                                "Landscape" = 16,
                                "Weather" = 18)) +
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 10), hjust = 0.45),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  )

p

mean_estimates2

save(mean_estimates2,
     samples_df,
     file = "MAWTRC Nesting Ecology Manuscript/Figures/RData/Pre-Nesting Movement/Study Area Analysis/NSS.2D.RData")


################################################################################
###############################################################################X


################################################################################
## Known Fate Model

# This section fits a known fate model for Pennsylvania WMU 2D

################################################################################
## Data Prep- Nest-Level Covs

# Read in nests csv
pa.nests <- read_csv("Data Management/Csvs/Processed/Nests/Vegetation Surveys/Pennsylvania/20250121_CleanedNestsVeg_2022_2023.csv")
pa.nests

#' Subset nesting data for 2D in year 2022
pa.nests.2D <- dplyr::filter(pa.nests, WMU =="2D")

#' Csv from incubation start and end script
nests.inc <- read_csv("Data Management/Csvs/Processed/Incubation Dates/Pennsylvania/20250131_NestAttempts_allbirds_PA.csv")
nests.inc

#' Sample from known fate model to ensure consistency
pa.sample <- read_csv("Samples/Pennsylvania/NestingSample_PA.csv")
pa.sample
nests.inc <- right_join(nests.inc, pa.sample)

#' Merge the filtered nests.veg with nests by nestid
pa.nests<- inner_join(pa.nests.2D, nests.inc, by = "NestID") 

pa.nests.sf <- pa.nests %>%
  dplyr::mutate("Long1" = Long) %>%
  dplyr::mutate("Lat1" = Lat) %>%
  st_as_sf(., coords = c("Long1", "Lat1"), crs = 4326) %>%
  st_transform(5070)
mapview(pa.nests.sf)

#' Extract point value at each nest
landcov <-terra::extract(pa.nlcd, pa.nests.sf)

#' Bind columns together
pa.nests.landcov <- cbind(pa.nests.sf, landcov) %>%
  dplyr::rename("landuse" = Class)


################################################################################
## Extract Distances to Roads- Pennsylvania

#' PASDA Raster
pasda.dist.prim <- terra::rast("Data Management/Rasters/Roads/Pennsylvania/PaRoadRast.Prim.tiff")
pasda.dist.sec <- terra::rast("Data Management/Rasters/Roads/Pennsylvania/PaRoadRast.Sec.tiff")

#' Read in nest data with landcover
pa.nests.landcov.sf <- pa.nests.landcov

#' Extract point value at each nest
dist.prim.out <-terra::extract(pasda.dist.prim, pa.nests.landcov.sf) %>%
  dplyr::select(-ID)%>%
  dplyr::rename("primary" = layer)
dist.sec.out <-terra::extract(pasda.dist.sec, pa.nests.landcov.sf) %>%
  dplyr::select(-ID) %>%
  dplyr::rename("secondary" = layer)


pa.nests.landcov.sf.roads <- cbind(pa.nests.landcov.sf, dist.prim.out)
pa.nests.landcov.sf.roads <- cbind(pa.nests.landcov.sf.roads, dist.sec.out) %>%
  st_drop_geometry()

pa.nests.covs <- left_join(pa.nests.landcov.sf, pa.nests.landcov.sf.roads) %>%
  dplyr::rename("BandID" = BandID.x) %>%
  dplyr::rename("NestFate" = NestFate.y) %>%
  dplyr::select(landuse, 
                Case,
                NestID, 
                BandID, 
                PercGrassForb,
                PercWoody, 
                PercLitter,
                PercBare,
                PercBoulder,
                AvgMaxVO, 
                primary,
                secondary,
                StemCount,
                WoodyType1,
                AvgVO,
                PercFern, 
                GuardHt,
                HtWoody,
                HtFern,
                HtGrassForb,
                HtLitter, 
                HtBoulder,
                PlotType,
                startI,
                endI,
                NestFate,
                Long,
                Lat) 
#' Drop geometry column and create unique identifier for NestID
pa.nests.covs <- sf::st_drop_geometry(pa.nests.covs)
pa.nests.covs$NestID1 <- paste(pa.nests.covs$NestID, pa.nests.covs$PlotType, sep = "_")

#' Convert land cover classifications to a categorical variable and create separate columns
pa.nests.covs$Pasture <- ifelse(pa.nests.covs$landuse == "Pasture", 1, 0)
pa.nests.covs$Crop <- ifelse(pa.nests.covs$landuse == "Crop", 1, 0)
pa.nests.covs$Developed <- ifelse(pa.nests.covs$landuse == "Developed", 1, 0)
pa.nests.covs$Deciduous <- ifelse(pa.nests.covs$landuse == "Deciduous Forest", 1, 0)
pa.nests.covs$Evergreen <- ifelse(pa.nests.covs$landuse == "Evergreen Forest", 1, 0)
pa.nests.covs$Mixed <- ifelse(pa.nests.covs$landuse == "Mixed Forest", 1, 0)
pa.nests.covs$Grassland <- ifelse(pa.nests.covs$landuse == "Grassland/Shrub", 1, 0)
pa.nests.covs$Wetland <- ifelse(pa.nests.covs$landuse == "Wetland", 1, 0)
pa.nests.covs$Water <- ifelse(pa.nests.covs$landuse == "Open Water", 1, 0)

#' Rename and consolidate columns 
nests <- pa.nests.covs %>%
  dplyr::select(NestID, 
                NestID1,
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
                PercLitter,
                PercBare,
                PercBoulder,
                Deciduous,
                Pasture,
                Crop,
                Developed,
                Evergreen,
                Mixed,
                Grassland,
                Wetland,
                primary,
                secondary,
                PlotType,
                Case,
                Long,
                Lat) %>%
  dplyr::filter(Case == "1")

#' Read in basal area data
basal <- read_csv("Data Management/Csvs/Processed/Covariates/Pennsylvania/BasalArea.csv")

nests <- right_join(basal, nests, by = "NestID1")

#' Zero values represent sites with no trees for basal area
nests <- nests %>%
  replace_na(list(basal_summary = 0)) 

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
nests$PercBoulder <- as.numeric(nests$PercBoulder)
nests$PercLitter <- as.numeric(nests$PercLitter)
nests$PercBare <- as.numeric(nests$PercBare)
nests$Basal <- as.numeric(nests$Basal)
nests$primary <- as.numeric(nests$primary)
nests$secondary <- as.numeric(nests$secondary)

#' Scale continous predictors
nests.scaled <- nests %>% 
  dplyr::mutate(AvgVO = scale(AvgVO)) %>%
  dplyr::mutate(PercWoody = scale(PercWoody)) %>%
  dplyr::mutate(PercGrassForb = scale(PercGrassForb)) %>%
  dplyr::mutate(AvgMaxVO = scale(AvgMaxVO)) %>%
  dplyr::mutate(PercFern = scale(PercFern)) %>%
  dplyr::mutate(StemCount = scale(StemCount)) %>%
  dplyr::mutate(PercBoulder = scale(PercBoulder)) %>% 
  dplyr::mutate(PercLitter = scale(PercLitter)) %>%
  dplyr::mutate(PercBare = scale(PercBare)) %>%
  dplyr::mutate(Basal = scale(Basal)) %>%
  dplyr::mutate(secondary = scale(secondary)) %>%
  dplyr::mutate(primary = scale(primary))
nests.scaled

#' Change columns back to numeric
nests.scaled$AvgVO <- as.numeric(nests.scaled$AvgVO)
nests.scaled$PercWoody <- as.numeric(nests.scaled$PercWoody)
nests.scaled$PercGrassForb <- as.numeric(nests.scaled$PercGrassForb)
nests.scaled$PercFern <- as.numeric(nests.scaled$PercFern)
nests.scaled$AvgMaxVO <- as.numeric(nests.scaled$AvgMaxVO)
nests.scaled$StemCount <- as.numeric(nests.scaled$StemCount)
nests.scaled$PercBoulder <- as.numeric(nests.scaled$PercBoulder)
nests.scaled$PercLitter <- as.numeric(nests.scaled$PercLitter)
nests.scaled$PercBare <- as.numeric(nests.scaled$PercBare)
nests.scaled$Basal <- as.numeric(nests.scaled$Basal)
nests.scaled$primary <- as.numeric(nests.scaled$primary)
nests.scaled$secondary <- as.numeric(nests.scaled$secondary)
nests.scaled

#' Check Correlations
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

#' Use ggpairs to visualize correlations
ggpairs(as.data.frame(C), 
        upper = list(continuous = "cor"),
        diag = list(continuous = "barDiag"),
        lower = list(continuous = "points"))


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
nests.scaled <- nests.scaled %>% dplyr::rename("BandID" = BandID.y) %>% dplyr::select(-BandID.x)
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

#' Nest Incubation Dates- Calendar Day
nests.scaled <- nests.scaled %>%
  dplyr::mutate("Nest Incubation Date" = scale(lubridate::yday(startI))) 
nests.scaled$`Nest Incubation Date` <- as.numeric(nests.scaled$`Nest Incubation Date`)


################################################################################
## Data Prep- Behavior Covariates

#' Read in RDS file from the behavior script
hens.behav.out <- readRDS("Data Management/Csvs/Processed/Covariates/Pennsylvania/Behavior/hens.behav.out.RDS")

#' Drop geometry
#' Keep all observations of hens.behav.out that exist in nests.scaled 
nests.scaled <- nests.scaled %>%
  st_drop_geometry() %>%
  right_join(hens.behav.out, nests.scaled, by = "NestID")

nests.scaled <- nests.scaled[-c(45:159), ]


#' Calculate mean and standard deviation by age class
nests.scaled %>%
  group_by(age) %>%
  summarise(
    mean_constancy = mean(IncubationConstancy, na.rm = TRUE),
    sd_constancy = sd(IncubationConstancy, na.rm = TRUE)
  )

#' Total number of nests by age class
nests.scaled %>%
  group_by(age) %>%
  summarise(count = n())

#' Scale incubation constancy and change it back to numeric
#' Scale cumulative step length and change it back to numeric 
#' Remove all observations where the total locations of gps locations within a buffer is NA
nests.scaled <- nests.scaled %>%
  dplyr::mutate(IncubationConstancy = scale(IncubationConstancy)) %>%
  dplyr::mutate(IncubationConstancy = as.numeric(IncubationConstancy)) %>%
  dplyr::mutate(sum_sl = scale(sum_sl)) %>%
  dplyr::mutate(sum_sl = as.numeric(sum_sl)) %>%
  dplyr::filter(TotalLocations != "NA") %>%
  dplyr::filter(NestID != "4255_2022_1") 



################################################################################
## Data Prep- Weather Covariates

#' Create a csv file of needed information for Daymet download
#' Rename and select needed columns 
nest.sites <- nests.scaled %>%
  dplyr::rename("latitude" = Lat, "longitude" = Long, "site" = NestID) %>%
  dplyr::select(latitude, longitude, site)
write.csv(nest.sites, "Data Management/Csvs/Processed/Covariates/Pennsylvania/Weather/nests.sites.2D.csv")

#' Perform Daymet download  
w <- download_daymet_batch(
  file_location = 'Data Management/Csvs/Processed/Covariates/Pennsylvania/Weather/nests.sites.2D.csv',
  start = 2022,
  end = 2023,
  internal = TRUE
)

w[[1]]$data[109:259,]

n.weather.cov<-2
weather.array <- array(NA, dim = c(nrow(nests.scaled),nrow(w[[1]]$data),n.weather.cov))

for (i in 1:44) {
  weather.array[i, 1:730, 1] <- w[[i]]$data$tmin..deg.c.
  weather.array[i, 1:730, 2] <- w[[i]]$data$prcp..mm.day.
}
#saveRDS(nests.scaled, "nests.scaled.RDS")
#saveRDS(weather.array, "weather.RDS")


################################################################################
## Fit Weather and Habitat Known Fate Model

################################################################################
## Data Prep- Encounter Histories

#' Convert start and end of incubation days to date objects
nests.scaled$startI <- as.Date(nests.scaled$startI)
nests.scaled$endI <- as.Date(nests.scaled$endI)

#' Create a sequence of all unique dates between the min start and max end date
nests.scaled %>% 
  sf::st_drop_geometry() %>%
  group_by(NestYr) %>% summarise(mindate = min(startI),
                                 maxdate = max(endI))
#' Get minimum and maximum days within a year
min(nests.scaled$startI)
max(nests.scaled$endI)
inc.dates <- sort(unique(c(nests.scaled$startI, nests.scaled$endI)))
inc.dates <- seq(as.Date("2022-04-19"), as.Date("2022-09-09"), by = 1)
inc.dates <- format(inc.dates, "%m-%d")
inc.dates

#' Calculate the number of encounter days
occs <- length(inc.dates)

#' Initialize variables
#' Create matrix for encounter histories
first <- last <- c()
surv.caps <- matrix(data = NA, nrow = nrow(nests.scaled), ncol = occs)

#' Create encounter histories within surv.caps
#' If a nest fails it gets a zero at the end 
#' If it reaches hatch it ends with a 1
for(i in 1:nrow(nests.scaled)){ 
  first[i] <- which(inc.dates == format(nests.scaled$startI[i], "%m-%d")) 
  last[i] <- which(inc.dates == format(nests.scaled$endI[i], "%m-%d")) 
  surv.caps[i,first[i]:last[i]] <- 1 
  if(nests.scaled$NestFate[i] == 0)
  {surv.caps[i,last[i]] <- 0} 
}

#' Check work on first individual
surv.caps[11,]
first;
last;

#' Change surv.caps to a df
surv.caps = as.data.frame(surv.caps)
colnames(surv.caps) <- inc.dates
surv.caps$Year = nests.scaled$NestYr-2021


################################################################################
## Bring in Weather Data

#' Functions from Kery and Schaub to get the first and last encounter occasions
#' Make sure that f is always greater than k (First encounter always comes before last encounter)
getFirst <- function(x) {min(which(!is.na(x)))}
getLast <- function(x) {max(which(!is.na(x)))}
f <- apply(surv.caps[,1:144],1,getFirst)
k <- apply(surv.caps[,1:144],1,getLast)
f;k
f<k

#' Store weather data within arrays for 2022 and 2023
#' Rows = 158 (# of nests)
#' Columns = 365 (# of days within a year)
#' Dimensions = 1, or 2 (Temp = 1, Precip = 2)
temperatureC2022 <- weather.array[1:44,1:365,1]
temperatureC2023 <- weather.array[1:44,366:730,1]
precip2022 <- weather.array[1:44,1:365,2]
precip2023 <- weather.array[1:44,366:730,2]

################################################################################
## Get 3 day moving average

#' Convert from matrix to dataframe
#' For each row in surv.caps if the year is the first year of the study, collect the rolling mean for 2022 variables
#' For each row in surv.caps if the year is the second year of the study, collect the rolling mean for 2023 variables
testTempC3Roll = data.frame(matrix(NA, nrow = 44, ncol = 144))   
testPrecip3Roll = data.frame(matrix(NA, nrow = 44, ncol = 144))
for(i in 1:nrow(surv.caps)){
  if(surv.caps$Year[i] == 1){
    for(j in (f[i]+109):(k[i]+109)){
      testTempC3Roll[i,j-109] <- mean(temperatureC2022[i,((j-2):j)])
      testPrecip3Roll[i,j-109] <- mean(precip2022[i,((j-2):j)])
    }
  }
  else if(surv.caps$Year[i] == 2){
    for(j in (f[i]+109):(k[i]+109)){
      testTempC3Roll[i,j-109] <- mean(temperatureC2023[i,((j-2):j)])
      testPrecip3Roll[i,j-109] <- mean(precip2023[i,((j-2):j)])
    }
  }
  else{
    testTempC3Roll[i,j] <- "You messed up"
    testPrecip3Roll[i,] <- "You messed up"
  }
}

#' Check work
range(testPrecip3Roll, na.rm = T) 
range(testTempC3Roll, na.rm = T)

#' Scale weather data
#' Scale function had problems with the NAs within the matrices
weather.array2 <- array(data = NA, dim = c(nrow(testTempC3Roll), ncol(testTempC3Roll),2))
testTempC.mat <- as.matrix(testTempC3Roll)
testTempC.scale <- (testTempC.mat - mean(testTempC.mat,na.rm=T))/sd(testTempC.mat,na.rm=T)
testTempC.scale
testPrecip.mat <- as.matrix(testPrecip3Roll)
testPrecip.scale <- (testPrecip.mat - mean(testPrecip.mat,na.rm=T))/sd(testPrecip.mat,na.rm=T)
testPrecip.scale
weather.array2[,,1] <- testTempC.scale
weather.array2[,,2] <- testPrecip.scale

#' Check work
weather.array2[1,1:144,1]
weather.array2[1,1:144,2]

#' Save output RDS files
#saveRDS(weather.array2, "Known Fate Model/Pennsylvania/Weather Arrays/weather.pa_2022_2023_2Daymean.RDS")

#' Check correlations
temp_flat <- as.vector(weather.array2[,,1])
precip_flat <- as.vector(weather.array2[,,2])

#' Compute correlation, excluding NA pairs
cor_result <- cor(temp_flat, precip_flat, use = "complete.obs")

#' Print correlation
cor_result


################################################################################
## Compile Parameters into matrix format

#' Function to fill NA values with 0 throughout dataframe
#' These parameters are scaled so we are replacing NA values with the mean
fill_NA_with_value <- function(df, value = 0) {
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      replace(col, is.na(col), value)
    } else {
      col
    }
  })
  return(df)
}
nests.ready <- fill_NA_with_value(nests.scaled)
summary(nests.ready)


################################################################################
## Nimble Known Fate Model

# Estimates daily nest survival probability as a function of relevant variables

knownfate.habitat <- nimbleCode({
  
  #### Loop over individuals 
  for (i in 1:n.ind) {
    
    #### Get the first and last days a nest was active
    for(t in (first[i]+1):last[i]) {
      
      #### mu and z change if an animal is alive or dead 
      mu[i,t] <- phi[i,t]*z[i,t-1]
      
      #### Estimating daily nest survival in matrix format
      cloglog(phi[i,t]) <- inprod(beta[1:J], X[i, 1:J]) +
        weather.beta[1]*weather.array[i,t,1] + weather.beta[2]*weather.array[i,t,2] + 
        weather.beta[3]*weather.array[i,t,1]*weather.array[i,t,2]
      
      #### Likelihood of Nest Survival (Bernoulli)
      z[i,t]~dbern(mu[i, t])
      
    }
  }
  
  #### Priors for beta coefficients
  for (j in 1:J){
    beta[j] ~ dnorm(0, sd = 1)  
  }
  
  #### Priors for weather coefficients
  for (k in 1:n.weather.cov) {
    weather.beta[k] ~ dnorm(0, sd = 1)
  }
  
})

################################################################################

#' Model parameters
# Model parameters
X <- cbind(
  rep(1, nrow(nests.ready)),                                # Intercept (1)
  as.numeric(nests.ready$AvgVO),                            # Visual Obstruction
  as.numeric(nests.ready$StemCount),                        # Woody Stem Count
  as.numeric(nests.ready$Basal),                            # Basal Area
  as.numeric(nests.ready$PercWoody),                        # Percent Woody Vegetation
  as.numeric(nests.ready$PercGrassForb),                    # Percent Grass Forb
  as.numeric(nests.ready$primary),                          # Distance to Primary Road
  as.numeric(nests.ready$secondary),                        # Distance to Secondary Road
  as.numeric(nests.ready$Mixed),                            # Mixed Forest
  as.numeric(nests.ready$Evergreen),                        # Evergreen Forest
  as.numeric(nests.ready$Developed),                        # Developed
  as.numeric(nests.ready$Pasture),                          # Pasture
  as.numeric(nests.ready$Crop),                             # Crop
  as.numeric(nests.ready$Grassland),                        # Grassland
  as.numeric(nests.ready$Wetland),                          # Wetland
  as.numeric(nests.ready$PercFern),                         # Percent Fern
  as.numeric(nests.ready$`Nest Incubation Date`),           # Nest Incubation Date
  as.numeric(nests.ready$IncubationConstancy),              # Incubation Constancy
  as.numeric(nests.ready$age)                               # Age Class
)

#' Use ggpairs to visualize correlations
ggpairs(as.data.frame(X), 
        upper = list(continuous = "cor"),
        diag = list(continuous = "barDiag"),
        lower = list(continuous = "points"))

#' Constants
Consts <- list(
  n.ind = nrow(nests.ready),   
  X = X,                       
  J = ncol(X),                 
  first = first,               
  last = last,  
  n.weather.cov = dim(weather.array2)[3]+1,
  weather.array = weather.array2
)
surv.caps.matrix <- as.matrix(surv.caps)
colnames(surv.caps.matrix) <- NULL

#' Data (survival data: 1 = survived, 0 = failed)
Data <- list(
  z = surv.caps.matrix               
)

#' Initial values for parameters
Inits <- list(
  beta = rep(0, Consts$J),
  weather.beta = rep(0,Consts$n.weather.cov)
)

start <- Sys.time()

#' Run MCMC sampling with Nimble
nimbleMCMC_samples <- nimbleMCMC(
  code = knownfate.habitat,          
  constants = Consts,          
  inits = Inits,              
  data = Data,                 
  nburnin = 3000,             
  niter = 20000,               
  thin = 1                     
)

end <- Sys.time()

summary(nimbleMCMC_samples)
colMeans(nimbleMCMC_samples[,1:22])
colSds(nimbleMCMC_samples[,1:22])
MCMCtrace(nimbleMCMC_samples, pdf = FALSE)

#' Extract the posterior samples for the 'beta' parameters (columns 219 to 230)
beta_samples <- nimbleMCMC_samples[, 1:22]

# Convert mcmc.list to matrix
samples_matrix <- as.matrix(beta_samples)

# Convert the matrix to a data frame
samples_df <- as.data.frame(samples_matrix)

#' Create a vector of new names
new_names <- c("Intercept", 
               "Visual Obstruction",
               "Woody Stem Count",
               "Basal Area",
               "Percent Woody Vegetation",
               "Percent Grass/Forb",
               "Distance to Primary Road",
               "Distance to Secondary Road",
               "Mixed Forest",
               "Evergreen Forest",
               "Developed",
               "Pasture",
               "Crop",
               "Grassland/Shrub",
               "Wetland",
               "Percent Fern",
               "Nest Incubation Date",
               "Incubation Constancy",
               "Age Class",
               "Minimum Temperature",
               "Precipitation",
               "Precipitation * Minimum Temperature")

#' Assign the variable names to the columns of the beta_samples
colnames(samples_df) <- new_names

#' View the renamed data frame
head(samples_df)


################################################################################
## Data Prep for Beta Plot

#' Reshape the data into long format for ggplot
samples_long <- samples_df %>%
  pivot_longer(cols = everything(), 
               names_to = "parameter", 
               values_to = "estimate") 

#' View the reshaped data
head(samples_long)

credible_intervals <- samples_long %>%
  group_by(parameter) %>%
  summarise(
    lower = quantile(estimate, 0.05),   
    upper = quantile(estimate, 0.95),   
    .groups = 'drop'
  )

#' Calculate posterior mean for each parameter
#' Calculate 95% credible intervals using quantiles
mean_estimates <- samples_long %>%
  dplyr::group_by(parameter) %>%
  summarise(
    mean_estimate = mean(estimate),  
    lower = quantile(estimate, 0.05),  
    upper = quantile(estimate, 0.95),  
    .groups = 'drop'
  ) %>%
  dplyr::filter(parameter != "Intercept")
mean_estimates

mean_estimates <- mean_estimates %>%
  mutate(Scale = case_when(
    parameter %in% c("Percent Grass/Forb", 
                     "Percent Woody Vegetation", 
                     "Visual Obstruction",
                     "Woody Stem Count",
                     "Basal Area",
                     "Percent Fern") ~ "Nest",
    parameter %in% c( "Grassland/Shrub",
                      "Mixed Forest",
                      "Evergreen Forest",
                      "Developed",
                      "Distance to Primary Road",
                      "Distance to Secondary Road",
                      "Pasture",
                      "Crop",
                      "Wetland"
    ) ~ "Landscape",
    parameter %in% c("Minimum Temperature",
                     "Precipitation",
                     "Precipitation * Minimum Temperature"
    ) ~ "Weather",
    TRUE ~ "Individual"
  ))

mean_estimates

#' Organize variables into levels to be displayed
#' Filter out the intercept
mean_estimates3 <- mean_estimates %>%
  dplyr::mutate(parameter = factor(parameter, 
                                   levels = c( 
                                     "Precipitation * Minimum Temperature",
                                     "Minimum Temperature",
                                     "Precipitation",
                                     "Wetland",
                                     "Grassland/Shrub",
                                     "Mixed Forest",
                                     "Evergreen Forest",
                                     "Developed",
                                     "Pasture",
                                     "Crop",
                                     "Distance to Primary Road",
                                     "Distance to Secondary Road",
                                     "Percent Fern",
                                     "Percent Grass/Forb",
                                     "Percent Woody Vegetation",
                                     "Visual Obstruction",
                                     "Basal Area",
                                     "Woody Stem Count",
                                     "Incubation Constancy",
                                     "Nest Incubation Date",
                                     "Age Class"
                                   )))
mean_estimates3

mean_estimates3 <- mean_estimates3 %>%
  mutate(Scale = factor(Scale, levels = c("Individual", "Nest", "Landscape", "Weather")))


################################################################################
## Beta Plot 

#' Beta estimates and associated 90% credible intervals 
p1.betas <- ggplot(mean_estimates3, 
                   aes(x = parameter, 
                       y = mean_estimate, 
                       color = Scale, 
                       shape = Scale)) +
  geom_point(size = 3.5) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1.1) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  labs(x = "Parameter", y = "Posterior Mean") +
  theme_minimal() + 
  coord_flip()+
  scale_color_manual(values = c("Individual" = "#DA6509",
                                "Nest" = "#A44200",
                                "Landscape" = "#D65F5F",
                                "Weather" = "#197278")) +  
  scale_shape_manual(values = c("Individual" = 15, 
                                "Nest" = 17,
                                "Landscape" = 16,
                                "Weather" = 18)) +
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 10), hjust = 0.45),  
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  )
p1.betas

save(mean_estimates3,
     samples_df,
     file = "MAWTRC Nesting Ecology Manuscript/Figures/RData/Pre-Nesting Movement/Study Area Analysis/kf.2D.RData")


################################################################################
###############################################################################X
