
#'---
#' title: Habitat selection of female wild turkeys during pre-nesting (an SSF analysis)
#' author: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output: MovementProcess_Prep.RData (R workspace)
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script downloads movement data associated with nesting hens and extracts cumulative step length during incubation
#' **Last Updated**: 1/25/25


#############################################################################
## Load Packages 

#' Vector of package names
packages <- c("purrr",
              "lubridate",
              "dplyr",
              "move2",
              "tidyverse",
              "amt")

#' Function to load a package or install it if not already installed
load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

#' Apply the function to each package name
lapply(packages, load_packages)


###############################################################################
## Data Management

#' Csv from incubation start and end script
nests.inc <- read_csv("Data Management/Csvs/Processed/IncubationDates/Draft2/20250124_NestAttempts_allbirds.csv") 
nests.inc

#' Read in nests csv
pa.nests <- read_csv("Data Management/Csvs/Processed/Nests/Nests/2025_CleanedNests_2022_2023.csv")
pa.nests

#' Subset nesting data for 4D in year 2022
pa.nests.4D <- dplyr::filter(pa.nests, WMU =="4D")%>%
  dplyr::select(BandID, CheckDate, NestID, WMU) %>%
  dplyr::filter(NestID != "8204_2023_1")

#' Merge pa.nests.4D and nests.inc, only keep nests that exist in both pa.nests.4D and nests.inc
pa.nests.4D1 <- dplyr::inner_join(pa.nests.4D, nests.inc, by = "NestID") %>%
  dplyr::select(-CheckDate.y) %>%
  dplyr::rename("CheckDate" = CheckDate.x)
glimpse(pa.nests.4D1)

#' Subset nesting data for 3D in year 2022
pa.nests.3D <- dplyr::filter(pa.nests, WMU =="3D")%>%
  dplyr::select(BandID, CheckDate, NestID, WMU)

#' Merge pa.nests.3D and nests.inc, only keep nests that exist in both pa.nests.3D and nests.inc
pa.nests.3D1 <- dplyr::inner_join(pa.nests.3D, nests.inc, by = "NestID") %>%
  dplyr::select(-CheckDate.y) %>%
  dplyr::rename("CheckDate" = CheckDate.x)
glimpse(pa.nests.3D1)

#' Subset nesting data for 2D in year 2022
pa.nests.2D <- dplyr::filter(pa.nests, WMU =="2D")%>%
  dplyr::select(BandID, CheckDate, NestID, WMU)

#' Merge pa.nests.2D and nests.inc, only keep nests that exist in both pa.nests.2D and nests.inc
pa.nests.2D1 <- dplyr::inner_join(pa.nests.2D, nests.inc, by = "NestID") %>%
  dplyr::select(-CheckDate.y) %>%
  dplyr::rename("CheckDate" = CheckDate.x)
glimpse(pa.nests.2D1)

#' Subset nesting data for 5C in year 2022
pa.nests.5C <- dplyr::filter(pa.nests, WMU =="5C")%>%
  dplyr::select(BandID, CheckDate, NestID, WMU)

#' Merge pa.nests.5C and nests.inc, only keep nests that exist in both pa.nests.5C and nests.inc
pa.nests.5C1 <- dplyr::inner_join(pa.nests.5C, nests.inc, by = "NestID") %>%
  dplyr::select(-CheckDate.y) %>%
  dplyr::rename("CheckDate" = CheckDate.x)
glimpse(pa.nests.5C1)


#############################################################################
## Login to Movebank

login <- movebank_store_credentials(username = "Kyle.Smelter",
                                    password="Rayshawks5!",
                                    key="Kyle",
                                    force= T)

#############################################################################
## WMU 4D 

unique.ID.4d<-unique(pa.nests.4D1$NestID)

for (j in 1:length(unique.ID.4d)){
  tmp.subset.4d<-pa.nests.4D1[which(pa.nests.4D1$NestID==unique.ID.4d[j]),]
  tmp.subset.4d$TrackID<-paste(unique.ID.4d[j],seq(1,nrow(tmp.subset.4d),1),sep="_")
  
  for(i in 1:nrow(tmp.subset.4d)){
    BirdID<- as.character(tmp.subset.4d[i,1])
    EndDate <- gsub("\\D","", tmp.subset.4d$endI[i]) 
    
    StartDate <- gsub("\\D","", tmp.subset.4d$startI[i]) 
    
    Year <- lubridate::year(tmp.subset.4d$startI[i])
    
    
    dat.4d<- movebank_download_study(study ="Wild Turkey Pennsylvania WMU 4D", 
                                     login = login,
                                     individual_local_identifier= BirdID,
                                     timestamp_start= StartDate,
                                     timestamp_end= EndDate,
                                     removeDuplicatedTimestamps=T)
    mt_track_id(dat.4d)<-rep(tmp.subset.4d$TrackID[i],nrow(dat.4d))
    
    if(exists("full_all_4d")){ 
      full_all_4d <- rbind(full_all_4d, dat.4d)
      
      
    }else{
      full_all_4d <- dat.4d
    }
  }
}


##########################################################################
## WMU 3D 


unique.ID.3d<-unique(pa.nests.3D1$NestID)

for (j in 1:length(unique.ID.3d)){
  tmp.subset.3d<-pa.nests.3D1[which(pa.nests.3D1$NestID==unique.ID.3d[j]),]
  tmp.subset.3d$TrackID<-paste(unique.ID.3d[j],seq(1,nrow(tmp.subset.3d),1),sep="_")
  
  for(i in 1:nrow(tmp.subset.3d)){
    BirdID<- as.character(tmp.subset.3d[i,1])
    EndDate <- gsub("\\D","", tmp.subset.3d$endI[i]) 
    
    StartDate <- gsub("\\D","", tmp.subset.3d$startI[i]) 
    
    Year <- lubridate::year(tmp.subset.3d$startI[i])
    
    
    dat.3d<- movebank_download_study(study ="Wild Turkey Pennsylvania WMU 3D", 
                                     login = login,
                                     individual_local_identifier= BirdID,
                                     timestamp_start= StartDate,
                                     timestamp_end= EndDate,
                                     removeDuplicatedTimestamps=T)
    mt_track_id(dat.3d)<-rep(tmp.subset.3d$TrackID[i],nrow(dat.3d))
    
    if(exists("full_all_3d")){ 
      full_all_3d <- rbind(full_all_3d, dat.3d)
      
      
    }else{
      full_all_3d <- dat.3d
    }
  }
}


############################################################################
##  WMU 2D 


unique.ID.2d<-unique(pa.nests.2D1$NestID)

for (j in 1:length(unique.ID.2d)){
  tmp.subset.2d<-pa.nests.2D1[which(pa.nests.2D1$NestID==unique.ID.2d[j]),]
  tmp.subset.2d$TrackID<-paste(unique.ID.2d[j],seq(1,nrow(tmp.subset.2d),1),sep="_")
  
  for(i in 1:nrow(tmp.subset.2d)){
    BirdID<- as.character(tmp.subset.2d[i,1])
    EndDate <- gsub("\\D","", tmp.subset.2d$endI[i]) 
    
    StartDate <- gsub("\\D","", tmp.subset.2d$startI[i]) 
    
    Year <- lubridate::year(tmp.subset.2d$startI[i])
    
    dat.2d<- movebank_download_study(study ="Wild Turkey Pennsylvania WMU 2D", 
                                     login = login,
                                     individual_local_identifier= BirdID,
                                     timestamp_start= StartDate,
                                     timestamp_end= EndDate,
                                     removeDuplicatedTimestamps=T)
    mt_track_id(dat.2d)<-rep(tmp.subset.2d$TrackID[i],nrow(dat.2d))
    
    if(exists("full_all_2d")){ 
      full_all_2d <- rbind(full_all_2d, dat.2d)
      
      
    }else{
      full_all_2d <- dat.2d
    }
  }
}


###########################################################################
##  WMU 5C 


unique.ID.5c<-unique(pa.nests.5C1$NestID)

for (j in 1:length(unique.ID.5c)){
  tmp.subset.5c<-pa.nests.5C1[which(pa.nests.5C$NestID==unique.ID.5c[j]),]
  tmp.subset.5c$TrackID<-paste(unique.ID.5c[j],seq(1,nrow(tmp.subset.5c),1),sep="_")
  
  for(i in 1:nrow(tmp.subset.5c)){
    BirdID<- as.character(tmp.subset.5c[i,1])
    EndDate <- gsub("\\D","", tmp.subset.5c$endI[i]) 
    
    StartDate <- gsub("\\D","", tmp.subset.5c$startI[i]) 
    
    Year <-lubridate::year(tmp.subset.5c$startI[i])
    
    
    dat.5c<- movebank_download_study(study ="Wild Turkey Pennsylvania WMU 5C", 
                                     login = login,
                                     individual_local_identifier= BirdID,
                                     timestamp_start= StartDate,
                                     timestamp_end= EndDate,
                                     removeDuplicatedTimestamps=T)
    mt_track_id(dat.5c)<-rep(tmp.subset.5c$TrackID[i],nrow(dat.5c))
    
    if(exists("full_all_5c")){ 
      full_all_5c <- rbind(full_all_5c, dat.5c)
      
      
    }else{
      full_all_5c <- dat.5c
    }
  }
}


############################################################################
## Extract cumulative distance traveled 

#' Convert move objects to dataframes
full_all_3d <- as.data.frame(full_all_3d)
full_all_4d <- as.data.frame(full_all_4d)
full_all_2d <- as.data.frame(full_all_2d)
full_all_5c <- as.data.frame(full_all_5c)

#' Create df with all study areas
df.all <- rbind(full_all_5c, 
                full_all_3d, 
                full_all_2d, 
                full_all_4d) %>%
  dplyr::rename("BirdID"= individual_local_identifier) 

#' Separate geometry lat and longs into separate columns and create new dataframe
#' Organize timestamp to be formatted in year, month, day, hour, minutes, seconds
#' Map function applies a function to each element of a vector
hens.all <- df.all%>%
  mutate(long = unlist(map(df.all$geometry,1)),
         lat = unlist(map(df.all$geometry,2))) %>%
  dplyr::select(BirdID, timestamp,long, lat) 

#' Temporary Data
dat <- hens.all %>%
  dplyr::select(BirdID, timestamp,long, lat) 

#' Create a track
trk <- amt::make_track(tbl=dat, .x= long, .y=lat, .t=timestamp, id=BirdID,
                       crs= 4326) %>% 
  amt::transform_coords(5070)
class(trk)
glimpse(trk)

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
glimpse(stps)

#' Cumulative step length during pre-nesting period for each bird
dist.traveled.cov <- stps %>%
  group_by(id) %>%
  summarise(sum_sl = sum(sl_, na.rm = TRUE)) %>%
  ungroup()
dist.traveled.cov
summary(dist.traveled.cov)

#' Write csv
write.csv(dist.traveled.cov, "Data Management/Csvs/Processed/Covariates/Behavior/dist.traveled.csv")

################################################################################
################################################################################

