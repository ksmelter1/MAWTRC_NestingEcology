
#'---
#' title: Habitat selection of female wild turkeys during pre-nesting (an SSF analysis)
#' author: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output: MovementProcess_Prep.RData (R workspace)
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script downloads movement data associated with each hens nesting attempt from movebank and exports hen movement data as RDS files.
#' **Last Updated**: 12/9/24


#####################
## Load Packages 

#' Vector of package names
packages <- c("purrr",
              "lubridate",
              "dplyr",
              "move2",
              "tidyverse")

#' Function to load a package or install it if not already installed
load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

#' Apply the function to each package name
lapply(packages, load_packages)


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


############################
## Prepare 4D Nest Data 

#' Subset nesting data for 4D in year 2022
pa.nests.4D <- dplyr::filter(pa.nests, wmu_n =="4D")%>%
  dplyr::select(bandid, checkyr, checkmo, checkday, nestid, wmu_n, clutchsize)


#' Clean up issue with zeros in dates column
pa.nests.4D$checkday<-ifelse(nchar(pa.nests.4D$checkday)==1,paste(0,pa.nests.4D$checkday,sep=""),pa.nests.4D$checkday)
pa.nests.4D$checkmo<-ifelse(nchar(pa.nests.4D$checkmo)==1,paste(0,pa.nests.4D$checkmo,sep=""),pa.nests.4D$checkmo)

#' Build df with information we need
pa.nests.4D$checkdate <- paste0(pa.nests.4D$checkyr, 
                                pa.nests.4D$checkmo, 
                                pa.nests.4D$checkday) 
# Convert checkdate to Date object
pa.nests.4D$checkdate <- as.Date(pa.nests.4D$checkdate, format = "%Y%m%d")
glimpse(pa.nests.4D)

##############################
## Incubation Data 4D

#' Csv from incubation start and end script
nests.inc <- read_csv("Data Management/Csvs/processed data/IncubationDates/NestAttempts_allbirds.csv")
nests.inc

#' Merge pa.nests.4D and nests.inc, only keep nests that exist in both pa.nests.4D and nests.inc
pa.nests.4D1 <- dplyr::inner_join(pa.nests.4D, nests.inc, by = "nestid") %>%
  dplyr::rename("checkdate" = checkdate.x) %>%
  dplyr::rename("bandid" = bandid.x) %>%
  dplyr::rename("wmu" = wmu_n) %>%
  dplyr::select(-bandid.y, -checkdate.y)
glimpse(pa.nests.4D1)

#' Ensure startI is in Date format
pa.nests.4D1$startI <- as.Date(pa.nests.4D1$startI) 

#' Now, iterate over the rows and subtract clutchsize days
for (i in 1:nrow(pa.nests.4D1)) {
  clutchsize <- pa.nests.4D1$clutchsize[i]  
  startI <- pa.nests.4D1$startI[i]  
  
  #' Subtracting the clutch size (in days) from the startI and creating the new 'startdate' column
  pa.nests.4D1$startdate[i] <- startI - clutchsize - days(14)
}

pa.nests.4D1$startdate <- as.Date(pa.nests.4D1$startdate, format = "%Y%m%d")
glimpse(pa.nests.4D1)

############################
## Prepare 3D Nest Data 

#' Subset nesting data for 4D in year 2022
pa.nests.3D <- dplyr::filter(pa.nests, wmu_n =="3D")%>%
  dplyr::select(bandid, checkyr, checkmo, checkday, nestid, wmu_n, clutchsize)


#' Clean up issue with zeros in dates column
pa.nests.3D$checkday<-ifelse(nchar(pa.nests.3D$checkday)==1,paste(0,pa.nests.3D$checkday,sep=""),pa.nests.3D$checkday)
pa.nests.3D$checkmo<-ifelse(nchar(pa.nests.3D$checkmo)==1,paste(0,pa.nests.3D$checkmo,sep=""),pa.nests.3D$checkmo)

#' Build df with information we need
pa.nests.3D$checkdate <- paste0(pa.nests.3D$checkyr, 
                                pa.nests.3D$checkmo, 
                                pa.nests.3D$checkday) 

pa.nests.3D$checkdate <- as.Date(pa.nests.3D$checkdate, format = "%Y%m%d")
glimpse(pa.nests.3D)

##############################
## Incubation Data 3D

#' Merge pa.nests.3D and nests.inc, only keep nests that exist in both pa.nests.4D and nests.inc
pa.nests.3D1 <- dplyr::inner_join(pa.nests.3D, nests.inc, by = "nestid") %>%
  dplyr::rename("checkdate" = checkdate.x) %>%
  dplyr::rename("bandid" = bandid.x) %>%
  dplyr::rename("wmu" = wmu_n) %>%
  dplyr::select(-bandid.y, -checkdate.y)
glimpse(pa.nests.3D1)


#' Ensure startI is in Date format
pa.nests.3D1$startI <- as.Date(pa.nests.3D1$startI) 

#' Now, iterate over the rows and subtract clutchsize days
for (i in 1:nrow(pa.nests.3D1)) {
  clutchsize <- pa.nests.3D1$clutchsize[i]  
  startI <- pa.nests.3D1$startI[i]  
  
  #' Subtracting the clutch size (in days) from the startI and creating the new 'startdate' column
  pa.nests.3D1$startdate[i] <- startI - clutchsize - days(14)
}

glimpse(pa.nests.3D1)


pa.nests.3D1$startdate <- as.Date(pa.nests.3D1$startdate, format = "%Y%m%d")
glimpse(pa.nests.3D1)

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

############################
## Prepare 5C Nest Data 

#' Subset nesting data for 4D in year 2022
pa.nests.5C <- dplyr::filter(pa.nests, wmu_n =="5C")%>%
  dplyr::select(bandid, checkyr, checkmo, checkday, nestid, wmu_n, clutchsize)

#' Clean up issue with zeros in dates column
pa.nests.5C$checkday<-ifelse(nchar(pa.nests.5C$checkday)==1,paste(0,pa.nests.5C$checkday,sep=""),pa.nests.5C$checkday)
pa.nests.5C$checkmo<-ifelse(nchar(pa.nests.5C$checkmo)==1,paste(0,pa.nests.5C$checkmo,sep=""),pa.nests.5C$checkmo)

#' Build df with information we need
pa.nests.5C$checkdate <- paste0(pa.nests.5C$checkyr, 
                                pa.nests.5C$checkmo, 
                                pa.nests.5C$checkday) 

pa.nests.5C$checkdate <- as.Date(pa.nests.5C$checkdate, format = "%Y%m%d")
glimpse(pa.nests.5C)

#############################
## Incubation Data 5C

#' Merge pa.nests.5C and nests.inc, only keep nests that exist in both pa.nests.5C and nests.inc
pa.nests.5C1 <- dplyr::inner_join(pa.nests.5C, nests.inc, by = "nestid") %>%
  dplyr::rename("checkdate" = checkdate.x) %>%
  dplyr::rename("bandid" = bandid.x) %>%
  dplyr::rename("wmu" = wmu_n) %>%
  dplyr::select(-bandid.y, -checkdate.y)
glimpse(pa.nests.5C1)

#' Ensure startI is in Date format
pa.nests.5C1$startI <- as.Date(pa.nests.5C1$startI) 

#' Now, iterate over the rows and subtract clutchsize days
for (i in 1:nrow(pa.nests.5C1)) {
  clutchsize <- pa.nests.5C1$clutchsize[i]  
  startI <- pa.nests.5C1$startI[i]  
  
  #' Subtracting the clutch size (in days) from the startI and creating the new 'startdate' column
  pa.nests.5C1$startdate[i] <- startI - clutchsize - days(14)
}

glimpse(pa.nests.5C1)

#' Change to date object
pa.nests.5C1$startdate <- as.Date(pa.nests.5C1$startdate, format = "%Y%m%d")
glimpse(pa.nests.5C1)

############################################################
## Pull Movebank Data from Movebank for Specified Dates 

#' Login to movebank
login <- movebank_store_credentials(username = "Kyle.Smelter",
                                    password="Rayshawks5!",
                                    key="Kyle",
                                    force= T)

####################
## WMU 4D 

#' List of unique identifier
unique.ID.4d<-unique(pa.nests.4D1$nestid)

for (j in 1:length(unique.ID.4d)){
  tmp.subset.4d<-pa.nests.4D1[which(pa.nests.4D1$nestid==unique.ID.4d[j]),]
  tmp.subset.4d$TrackID<-paste(unique.ID.4d[j],seq(1,nrow(tmp.subset.4d),1),sep="_")
  
for(i in 1:nrow(tmp.subset.4d)){
  BirdID<- as.character(tmp.subset.4d[i,1])
  EndDate <- gsub("\\D","", tmp.subset.4d$startI[i]) 
  #' (Format for time is YYYYMMDDHHSSMM000)
  StartDate <- gsub("\\D","", tmp.subset.4d$startdate[i]) 
  #' This will be the birds exact incubation period 
  Year <-tmp.subset.4d$checkyr[i]
 #' track_id <- as.character(tmp.subset$band_nestid[i])
  
  dat.4d<- movebank_download_study(study ="Wild Turkey Pennsylvania WMU 4D", 
                                   login = login,
                                   individual_local_identifier= BirdID,
                                   timestamp_start= StartDate,
                                   timestamp_end= EndDate,
                                   removeDuplicatedTimestamps=T)
  mt_track_id(dat.4d)<-rep(tmp.subset.4d$TrackID[i],nrow(dat.4d))
  
  if(exists("full_all_4d")){ #' rbind ind bird data to create one large df
    full_all_4d <- rbind(full_all_4d, dat.4d)
    
    
  }else{
    full_all_4d <- dat.4d
  }
}
}


##################
## WMU 3D 

#' List of unique identifier
unique.ID.3d<-unique(pa.nests.3D1$nestid)

for (j in 1:length(unique.ID.3d)){
  tmp.subset.3d<-pa.nests.3D1[which(pa.nests.3D1$nestid==unique.ID.3d[j]),]
  tmp.subset.3d$TrackID<-paste(unique.ID.3d[j],seq(1,nrow(tmp.subset.3d),1),sep="_")
  
  for(i in 1:nrow(tmp.subset.3d)){
    BirdID<- as.character(tmp.subset.3d[i,1])
    EndDate <- gsub("\\D","", tmp.subset.3d$checkdate[i]) 
    #' (Format for time is YYYYMMDDHHSSMM000)
    StartDate <- gsub("\\D","", tmp.subset.3d$startdate[i]) 
    #' 30 days earlier than check date 
    #' This will be the birds exact incubation period 
    #' This will be the birds exact incubation period 
    Year <-tmp.subset.3d$checkyr[i]
    #' track_id <- as.character(tmp.subset$band_nestid[i])
    
    dat.3d<- movebank_download_study(study ="Wild Turkey Pennsylvania WMU 3D", 
                                     login = login,
                                     individual_local_identifier= BirdID,
                                     timestamp_start= StartDate,
                                     timestamp_end= EndDate,
                                     removeDuplicatedTimestamps=T)
    mt_track_id(dat.3d)<-rep(tmp.subset.3d$TrackID[i],nrow(dat.3d))
    
    if(exists("full_all_3d")){ #' rbind ind bird data to create one large df
      full_all_3d <- rbind(full_all_3d, dat.3d)
      
      
    }else{
      full_all_3d <- dat.3d
    }
  }
}

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


###################
##  WMU 5C 

#' List of unique identifier
unique.ID.5c<-unique(pa.nests.5C1$nestid)

for (j in 1:length(unique.ID.5c)){
  tmp.subset.5c<-pa.nests.5C1[which(pa.nests.5C$nestid==unique.ID.5c[j]),]
  tmp.subset.5c$TrackID<-paste(unique.ID.5c[j],seq(1,nrow(tmp.subset.5c),1),sep="_")
  
  for(i in 1:nrow(tmp.subset.5c)){
    BirdID<- as.character(tmp.subset.5c[i,1])
    EndDate <- gsub("\\D","", tmp.subset.5c$checkdate[i]) 
    #' (Format for time is YYYYMMDDHHSSMM000)
    StartDate <- gsub("\\D","", tmp.subset.5c$startdate[i]) 
    #' 30 days earlier than check date 
    #' This will be the birds exact incubation period 
    Year <-tmp.subset.5c$checkyr[i]
    #' track_id <- as.character(tmp.subset$band_nestid[i])
    
    dat.5c<- movebank_download_study(study ="Wild Turkey Pennsylvania WMU 5C", 
                                     login = login,
                                     individual_local_identifier= BirdID,
                                     timestamp_start= StartDate,
                                     timestamp_end= EndDate,
                                     removeDuplicatedTimestamps=T)
    mt_track_id(dat.5c)<-rep(tmp.subset.5c$TrackID[i],nrow(dat.5c))
    
    if(exists("full_all_5c")){ #' rbind ind bird data to create one large df
      full_all_5c <- rbind(full_all_5c, dat.5c)
      
      
    }else{
      full_all_5c <- dat.5c
    }
  }
}


#################################
## Prep Data for Next Script

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

#' Save RData file for hens.all
#' Output all hens 
#save(hens.all, "Data Management/RData/Individual-Specific Movement Process/RData Files/01_MovementProcess_Prep.RData")
