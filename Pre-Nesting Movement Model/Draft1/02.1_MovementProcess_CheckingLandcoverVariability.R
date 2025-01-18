
library(terra)
library(sf)
library(dplyr)
library(mapview)
library(amt)
library(move2)

#' Read in raster of Pennsylvania NLCD
pa.nlcd <- terra::rast("Data Management/Rasters/nlcd/paNLCD.tiff")

#' Read in nests csv
pa.nests <- read_csv("Data Management/Csvs/processed data/Nests/nests_22_23_clean.csv")
pa.nests

######################
## Data Management

#' Change 99s into NA Values
pa.nests$eggshatched[pa.nests$eggshatched == 99] <- NA
pa.nests$eggsunhatched[pa.nests$eggsunhatched == 99] <- NA
pa.nests$eggsdestroyed[pa.nests$eggsdestroyed == 99] <- NA

#' Create clutch size column and remove unnecessary column
#' Clutch size is a minimum count 
pa.nests <- pa.nests %>%
  dplyr::mutate(clutchsize = rowSums(select(., eggshatched, eggsdestroyed, eggsunhatched), na.rm = TRUE)) %>%
  dplyr::select(-lat_n, -long_n)
glimpse(pa.nests)

#' Csv from incubation start and end script
nests.inc <- read_csv("Data Management/Csvs/processed data/IncubationDates/NestAttempts_25birdsfiltered.csv")
nests.inc

############################
## Prepare 2D Nest Data 

#' Subset nesting data for 4D in year 2022
pa.nests.2D <- dplyr::filter(pa.nests, wmu_n =="2D")%>%
  dplyr::select(bandid, checkyr, checkmo, checkday, nestid, wmu_n, clutchsize)


#' Clean up issue with zeros in dates column
pa.nests.2D$checkday<-ifelse(nchar(pa.nests.2D$checkday)==1,paste(0,pa.nests.2D$checkday,sep=""),pa.nests.2D$checkday)
pa.nests.2D$checkmo<-ifelse(nchar(pa.nests.2D$checkmo)==1,paste(0,pa.nests.2D$checkmo,sep=""),pa.nests.2D$checkmo)

#' Build df with information we need
pa.nests.2D$checkdate <- paste0(pa.nests.2D$checkyr, 
                                pa.nests.2D$checkmo, 
                                pa.nests.2D$checkday) 

pa.nests.2D$checkdate <- as.Date(pa.nests.2D$checkdate, format = "%Y%m%d")

#############################
## Incubation Data 2D

#' Merge pa.nests.2D and nests.inc, only keep nests that exist in both pa.nests.2D and nests.inc
pa.nests.2D1 <- dplyr::inner_join(pa.nests.2D, nests.inc, by = "nestid") %>%
  dplyr::rename("checkdate" = checkdate.x) %>%
  dplyr::rename("bandid" = bandid.x) %>%
  dplyr::rename("wmu" = wmu_n) %>%
  dplyr::select(-bandid.y, -checkdate.y)
glimpse(pa.nests.2D1)

#' Ensure startI is in Date format
pa.nests.2D1$startI <- as.Date(pa.nests.2D1$startI) 

#' Now, iterate over the rows and subtract clutchsize days
for (i in 1:nrow(pa.nests.2D1)) {
  clutchsize <- pa.nests.2D1$clutchsize[i]  
  startI <- pa.nests.2D1$startI[i]  
  
  #' Subtracting the clutch size (in days) from the startI and creating the new 'startdate' column
  pa.nests.2D1$startdate[i] <- startI - clutchsize - days(14)
}

glimpse(pa.nests.2D1)


pa.nests.2D1$startdate <- as.Date(pa.nests.2D1$startdate, format = "%Y%m%d")
glimpse(pa.nests.2D1)

############################################################
## Pull Movebank Data from Movebank for Specified Dates 

#' Login to movebank
login <- movebank_store_credentials(username = "Kyle.Smelter",
                                    password="Rayshawks5!",
                                    key="Kyle",
                                    force= T)

#################
##  WMU 2D 

#' List of unique identifier
unique.ID.2d<-unique(pa.nests.2D1$nestid)

for (j in 1:length(unique.ID.2d)){
  tmp.subset.2d<-pa.nests.2D1[which(pa.nests.2D1$nestid==unique.ID.2d[j]),]
  tmp.subset.2d$TrackID<-paste(unique.ID.2d[j],seq(1,nrow(tmp.subset.2d),1),sep="_")
  
  for(i in 1:nrow(tmp.subset.2d)){
    BirdID<- as.character(tmp.subset.2d[i,1])
    EndDate <- gsub("\\D","", tmp.subset.2d$startI[i]) 
    #' (Format for time is YYYYMMDDHHSSMM000)
    StartDate <- gsub("\\D","", tmp.subset.2d$startdate[i]) 
    #' 30 days earlier than check date 
    #' This will be the birds exact incubation period 
    Year <-tmp.subset.2d$checkyr[i]
    #' track_id <- as.character(tmp.subset$band_nestid[i])
    
    dat.2d<- movebank_download_study(study ="Wild Turkey Pennsylvania WMU 2D", 
                                     login = login,
                                     individual_local_identifier= BirdID,
                                     timestamp_start= StartDate,
                                     timestamp_end= EndDate,
                                     removeDuplicatedTimestamps=T)
    mt_track_id(dat.2d)<-rep(tmp.subset.2d$TrackID[i],nrow(dat.2d))
    
    if(exists("full_all_2d")){ #' rbind ind bird data to create one large df
      full_all_2d <- rbind(full_all_2d, dat.2d)
      
      
    }else{
      full_all_2d <- dat.2d
    }
  }
}

#######################################
## Generate Tracks and Random Steps

hens.all <- full_all_2d %>%
  as.data.frame()

#' Complete hen movement dataset from movebank
dat <- hens.all %>%
  dplyr::mutate("BirdID"=individual_local_identifier) %>%
  mutate(long = unlist(map(hens.all$geometry,1)),
         lat = unlist(map(hens.all$geometry,2))) %>%
  dplyr::select(BirdID, timestamp,long, lat) 

#' Create a track
#' New amt track function is just amt::track
trk <- amt::make_track(tbl=dat, .x= long, .y=lat, .t=timestamp, id=BirdID,
                       crs= 4326) %>% amt::transform_coords(5070)
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
x %>% track_resample(rate = minutes(30), tolerance = minutes(5)) %>%
  amt::steps_by_burst()

#' Check
class(x)

#' Summarize sampling rate
amt::summarize_sampling_rate(x)

#' Apply the same track resampling format to each hen within the dataset
#' This done by using the map function
#' Steps dataframe is the newly created column
trk2 <- trk1 %>%
  mutate(steps = map(data, function(x) 
    x %>% amt::track_resample(rate = minutes(30), tolerance = minutes(5))
    %>% amt::steps_by_burst())) 

#' Check
class(trk2)
glimpse(trk2)

#' Visualize step length distribution following vignette
# trk2 %>% dplyr::select(id, steps) %>% unnest(cols = steps) %>% dplyr::filter(id=="8202_2022_1_1") %>%
#   ggplot(aes(sl_, fill = factor(id))) + geom_density(alpha = 0.4)

#' Create object with all used steps for analysis
stps<- trk2 %>% dplyr::select(id, steps) %>% unnest(cols = steps)

#' Check
glimpse(stps)

#' Create random steps
#' Exponential step length used due to issues with formatting gamma
#' Extract covariates at the end of each used and available step
#' Random step lengths drawn from gamma distribution
#' Random turning angles drawn from a vonmises distribution
#' Include_observed: Include all used steps in the analysis
#' Rename and organize columns
random_steps<- amt::random_steps(
  stps,
  n_control = 10,
  sl_distr = amt::fit_distr(stps$sl_, "gamma"),
  ta_distr = amt::fit_distr(stps$ta_, "vonmises"),
  include_observed = T) %>%
  amt::extract_covariates(pa.nlcd, where ="end") 

glimpse(random_steps)

##################################
## Check Tracks on top of landcover

onestep <- random_steps %>%
  filter(id == "8176_2022_1_1", step_id_ >= 768, step_id_ <= 3000) %>%
  select(x1_, y1_, x2_, y2_, case_)  # Keep the 'case_' column for later

# Handle numeric columns separately
hold <- apply(onestep, 1, function(x) {
  # Extract only numeric columns (x1_, y1_, x2_, y2_)
  v <- as.numeric(x[1:4])  # Assume first four columns are numeric
  m <- matrix(v, byrow = TRUE, nrow = 2)
  
  # Return the line geometry with additional case_ information
  line <- st_as_sf(st_sfc(st_linestring(m)))
  line$case_ <- x["case_"]  # Add the 'case_' value to the line
  return(line)
})

lines.hold <- do.call("rbind", hold)

# Ensure lines.hold is an sf object
lines.hold.sf <- st_as_sf(lines.hold)

# Set the CRS for the lines (using CRS of the raster)
crs_pa <- crs(pa.nlcd)
lines.hold.sf <- st_set_crs(lines.hold.sf, crs_pa)

#' Create a 1000m buffer around the point
hen_buffer <- st_buffer(lines.hold.sf, dist = 1000) 

#' Crop the raster to the extent of the 1000m buffer
pa_cropped <- crop(pa.nlcd, hen_buffer)

#' Ensure the CRS of the buffer matches the raster's CRS
hen_buffer <- st_transform(hen_buffer, crs = crs(pa.nlcd))

# Fix invalid geometries
lines.hold.sf <- st_make_valid(lines.hold.sf)

# Add the 'case_' column to the lines.hold.sf object
lines.hold.sf$case_ <- onestep$case_

#' Check 
unique(lines.hold.sf$case_)  # Check unique values of case_

# Ensure case_ is a factor (this explicitly sets TRUE/FALSE levels)
lines.hold.sf$case_ <- factor(lines.hold.sf$case_, levels = c(TRUE, FALSE))

# Check unique values in 'case_' to ensure the factor levels are correctly set
unique_case <- unique(lines.hold.sf$case_)
print(unique_case)

# Create a vector of colors based on the 'case_' column
color_vector <- ifelse(lines.hold.sf$case_ == "TRUE", "red", "blue")

# Plot the raster alone first
mapview(pa_cropped) +
mapview(lines.hold.sf, color = color_vector, legend = TRUE)

################################################################################
onestep <- random_steps %>%
  filter(id == "8168_2022_1_1", step_id_ >= 768, step_id_ <= 3000) %>%
  select(x1_, y1_, x2_, y2_, case_)  # Keep the 'case_' column for later

# Handle numeric columns separately
hold <- apply(onestep, 1, function(x) {
  # Extract only numeric columns (x1_, y1_, x2_, y2_)
  v <- as.numeric(x[1:4])  # Assume first four columns are numeric
  m <- matrix(v, byrow = TRUE, nrow = 2)
  
  # Return the line geometry with additional case_ information
  line <- st_as_sf(st_sfc(st_linestring(m)))
  line$case_ <- x["case_"]  # Add the 'case_' value to the line
  return(line)
})

lines.hold <- do.call("rbind", hold)

# Ensure lines.hold is an sf object
lines.hold.sf <- st_as_sf(lines.hold)

# Set the CRS for the lines (using CRS of the raster)
crs_pa <- crs(pa.nlcd)
lines.hold.sf <- st_set_crs(lines.hold.sf, crs_pa)

#' Create a 1000m buffer around the point
hen_buffer <- st_buffer(lines.hold.sf, dist = 1000) 

#' Crop the raster to the extent of the 1000m buffer
pa_cropped <- crop(pa.nlcd, hen_buffer)

#' Ensure the CRS of the buffer matches the raster's CRS
hen_buffer <- st_transform(hen_buffer, crs = crs(pa.nlcd))

# Fix invalid geometries
lines.hold.sf <- st_make_valid(lines.hold.sf)

# Add the 'case_' column to the lines.hold.sf object
lines.hold.sf$case_ <- onestep$case_

#' Check 
unique(lines.hold.sf$case_)  # Check unique values of case_

# Ensure case_ is a factor (this explicitly sets TRUE/FALSE levels)
lines.hold.sf$case_ <- factor(lines.hold.sf$case_, levels = c(TRUE, FALSE))

# Check unique values in 'case_' to ensure the factor levels are correctly set
unique_case <- unique(lines.hold.sf$case_)
print(unique_case)

# Create a vector of colors based on the 'case_' column
color_vector <- ifelse(lines.hold.sf$case_ == "TRUE", "red", "blue")

# Plot the raster alone first
mapview(pa_cropped) +
mapview(lines.hold.sf, color = color_vector, legend = TRUE)

################################################################################
onestep <- random_steps %>%
  filter(id == "8173_2022_1_1", step_id_ >= 768, step_id_ <= 1500) %>%
  select(x1_, y1_, x2_, y2_, case_)  # Keep the 'case_' column for later

# Handle numeric columns separately
hold <- apply(onestep, 1, function(x) {
  # Extract only numeric columns (x1_, y1_, x2_, y2_)
  v <- as.numeric(x[1:4])  # Assume first four columns are numeric
  m <- matrix(v, byrow = TRUE, nrow = 2)
  
  # Return the line geometry with additional case_ information
  line <- st_as_sf(st_sfc(st_linestring(m)))
  line$case_ <- x["case_"]  # Add the 'case_' value to the line
  return(line)
})

lines.hold <- do.call("rbind", hold)

# Ensure lines.hold is an sf object
lines.hold.sf <- st_as_sf(lines.hold)

# Set the CRS for the lines (using CRS of the raster)
crs_pa <- crs(pa.nlcd)
lines.hold.sf <- st_set_crs(lines.hold.sf, crs_pa)

#' Create a 1000m buffer around the point
hen_buffer <- st_buffer(lines.hold.sf, dist = 1000) 

#' Crop the raster to the extent of the 1000m buffer
pa_cropped <- crop(pa.nlcd, hen_buffer)

#' Ensure the CRS of the buffer matches the raster's CRS
hen_buffer <- st_transform(hen_buffer, crs = crs(pa.nlcd))

# Fix invalid geometries
lines.hold.sf <- st_make_valid(lines.hold.sf)

# Add the 'case_' column to the lines.hold.sf object
lines.hold.sf$case_ <- onestep$case_

#' Check 
unique(lines.hold.sf$case_)  # Check unique values of case_

# Ensure case_ is a factor (this explicitly sets TRUE/FALSE levels)
lines.hold.sf$case_ <- factor(lines.hold.sf$case_, levels = c(TRUE, FALSE))

# Check unique values in 'case_' to ensure the factor levels are correctly set
unique_case <- unique(lines.hold.sf$case_)
print(unique_case)

# Create a vector of colors based on the 'case_' column
color_vector <- ifelse(lines.hold.sf$case_ == "TRUE", "red", "blue")

# Plot the raster alone first
mapview(pa_cropped) +
  mapview(lines.hold.sf, color = color_vector, legend = TRUE)
################################################################################
onestep <- random_steps %>%
  filter(id == "8169_2022_1_1") %>%
  select(x1_, y1_, x2_, y2_, case_)  # Keep the 'case_' column for later

# Handle numeric columns separately
hold <- apply(onestep, 1, function(x) {
  # Extract only numeric columns (x1_, y1_, x2_, y2_)
  v <- as.numeric(x[1:4])  # Assume first four columns are numeric
  m <- matrix(v, byrow = TRUE, nrow = 2)
  
  # Return the line geometry with additional case_ information
  line <- st_as_sf(st_sfc(st_linestring(m)))
  line$case_ <- x["case_"]  # Add the 'case_' value to the line
  return(line)
})

lines.hold <- do.call("rbind", hold)

# Ensure lines.hold is an sf object
lines.hold.sf <- st_as_sf(lines.hold)

# Set the CRS for the lines (using CRS of the raster)
crs_pa <- crs(pa.nlcd)
lines.hold.sf <- st_set_crs(lines.hold.sf, crs_pa)

#' Create a 1000m buffer around the point
hen_buffer <- st_buffer(lines.hold.sf, dist = 1000) 

#' Crop the raster to the extent of the 1000m buffer
pa_cropped <- crop(pa.nlcd, hen_buffer)

#' Ensure the CRS of the buffer matches the raster's CRS
hen_buffer <- st_transform(hen_buffer, crs = crs(pa.nlcd))

# Fix invalid geometries
lines.hold.sf <- st_make_valid(lines.hold.sf)

# Add the 'case_' column to the lines.hold.sf object
lines.hold.sf$case_ <- onestep$case_

#' Check 
unique(lines.hold.sf$case_)  # Check unique values of case_

# Ensure case_ is a factor (this explicitly sets TRUE/FALSE levels)
lines.hold.sf$case_ <- factor(lines.hold.sf$case_, levels = c(TRUE, FALSE))

# Check unique values in 'case_' to ensure the factor levels are correctly set
unique_case <- unique(lines.hold.sf$case_)
print(unique_case)

# Create a vector of colors based on the 'case_' column
color_vector <- ifelse(lines.hold.sf$case_ == "TRUE", "red", "blue")

# Plot the raster alone first
mapview(pa_cropped) +
  mapview(lines.hold.sf, color = color_vector, legend = TRUE)
################################################################################
onestep <- random_steps %>%
  filter(id == "8175_2023_1_1", step_id_ >= 5158, step_id_ <= 5400) %>%
  select(x1_, y1_, x2_, y2_, case_, step_id_)  # Keep the 'case_' column for later

# Handle numeric columns separately
hold <- apply(onestep, 1, function(x) {
  # Extract only numeric columns (x1_, y1_, x2_, y2_)
  v <- as.numeric(x[1:4])  # Assume first four columns are numeric
  m <- matrix(v, byrow = TRUE, nrow = 2)
  
  # Return the line geometry with additional case_ information
  line <- st_as_sf(st_sfc(st_linestring(m)))
  line$case_ <- x["case_"]  # Add the 'case_' value to the line
  return(line)
})

lines.hold <- do.call("rbind", hold)

# Ensure lines.hold is an sf object
lines.hold.sf <- st_as_sf(lines.hold)

# Set the CRS for the lines (using CRS of the raster)
crs_pa <- crs(pa.nlcd)
lines.hold.sf <- st_set_crs(lines.hold.sf, crs_pa)

#' Create a 1000m buffer around the point
hen_buffer <- st_buffer(lines.hold.sf, dist = 1000) 

#' Crop the raster to the extent of the 1000m buffer
pa_cropped <- crop(pa.nlcd, hen_buffer)

#' Ensure the CRS of the buffer matches the raster's CRS
hen_buffer <- st_transform(hen_buffer, crs = crs(pa.nlcd))

# Fix invalid geometries
lines.hold.sf <- st_make_valid(lines.hold.sf)

# Add the 'case_' column to the lines.hold.sf object
lines.hold.sf$case_ <- onestep$case_

#' Check 
unique(lines.hold.sf$case_)  # Check unique values of case_

# Ensure case_ is a factor (this explicitly sets TRUE/FALSE levels)
lines.hold.sf$case_ <- factor(lines.hold.sf$case_, levels = c(TRUE, FALSE))

# Check unique values in 'case_' to ensure the factor levels are correctly set
unique_case <- unique(lines.hold.sf$case_)
print(unique_case)

# Create a vector of colors based on the 'case_' column
color_vector <- ifelse(lines.hold.sf$case_ == "TRUE", "red", "blue")

# Plot the raster alone first
mapview(pa_cropped) +
  mapview(lines.hold.sf, color = color_vector, legend = TRUE)

