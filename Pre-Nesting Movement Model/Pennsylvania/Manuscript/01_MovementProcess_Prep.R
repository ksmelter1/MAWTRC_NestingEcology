
#'---
#' title: Habitat selection of female wild turkeys during pre-nesting (an SSF analysis)
#' author: "K. Smelter
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#'---
#'  
#' **Purpose**: This script downloads movement data associated with each hens nesting attempt from movebank 
#' **Last Updated**: 1/25/25


################################################################################
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
pa.nests <- read_csv("Data Management/Csvs/Processed/Nests/Nests/Pennsylvania/2025_CleanedNests_2022_2023_PA.csv")
pa.nests

################################################################################
## Data Management

#' Create clutch size column and remove unnecessary column
#' Clutch size is a minimum count 
pa.nests <- pa.nests %>%
  dplyr::mutate(clutchsize = rowSums(select(., EggsHatch, EggsDestroyed, EggsUnhatch, EggsDepred), na.rm = TRUE)) 
glimpse(pa.nests)


################################################################################
## Prepare 4D Nest Data 

#' Subset nesting data for 4D in year 2022
pa.nests.4D <- dplyr::filter(pa.nests, WMU =="4D")%>%
  dplyr::select(BandID, CheckDate, NestID, WMU, clutchsize)

################################################################################
## Incubation Data 4D

#' Csv from incubation start and end script
nests.inc <- read_csv("Data Management/Csvs/Processed/Incubation Dates/Pennsylvania/20250131_NestAttempts_allbirds_PA.csv")
nests.inc

#' Sample from known fate model to ensure consistency
pa.sample <- read_csv("Samples/Pennsylvania/NestingSample_PA.csv")
pa.sample
nests.inc <- right_join(nests.inc, pa.sample)

#' Merge pa.nests.4D and nests.inc, only keep nests that exist in both pa.nests.4D and nests.inc
pa.nests.4D1 <- dplyr::inner_join(pa.nests.4D, nests.inc, by = "NestID") %>%
  dplyr::select(-CheckDate.y) %>%
  dplyr::rename("CheckDate" = CheckDate.x)
glimpse(pa.nests.4D1)

#' Iterate over the rows and subtract clutchsize days
for (i in 1:nrow(pa.nests.4D1)) {
  clutchsize <- pa.nests.4D1$clutchsize[i]  
  startI <- pa.nests.4D1$startI[i]  
  
  #' Create the enddate by subtracting clutchsize (in days) from startI
  #' 5 day buffer
  pa.nests.4D1$enddate[i] <- startI - clutchsize - days(5)
  
  #' Create the startdate by subtracting 14 days from the enddate
  pa.nests.4D1$startdate[i] <- pa.nests.4D1$enddate[i]- days(14)
  
  #' Check the calculated dates
  print(paste("Row", i, ": Startdate =", pa.nests.4D1$startdate[i], ", Enddate =", pa.nests.4D1$enddate[i]))
}

#' Convert 'startdate' and 'enddate' columns to Date format
pa.nests.4D1$startdate <- as.Date(pa.nests.4D1$startdate)
pa.nests.4D1$enddate <- as.Date(pa.nests.4D1$enddate)

#' View the data frame structure
glimpse(pa.nests.4D1)


################################################################################
## Prepare 3D Nest Data 

#' Subset nesting data for 4D in year 2022
pa.nests.3D <- dplyr::filter(pa.nests, WMU =="3D")%>%
  dplyr::select(BandID, CheckDate, NestID, WMU, clutchsize)


################################################################################
## Incubation Data 3D

#' Merge pa.nests.4D and nests.inc, only keep nests that exist in both pa.nests.4D and nests.inc
pa.nests.3D1 <- dplyr::inner_join(pa.nests.3D, nests.inc, by = "NestID") %>%
  dplyr::select(-CheckDate.y) %>%
  dplyr::rename("CheckDate" = CheckDate.x)
glimpse(pa.nests.3D1)

#' Iterate over the rows and subtract clutchsize days
for (i in 1:nrow(pa.nests.3D1)) {
  clutchsize <- pa.nests.3D1$clutchsize[i]  
  startI <- pa.nests.3D1$startI[i]  
  
  #' Create the enddate by subtracting clutchsize (in days) from startI
  pa.nests.3D1$enddate[i] <- startI - clutchsize - days(5)
  
  #' Create the startdate by subtracting 14 days from the enddate
  pa.nests.3D1$startdate[i] <- pa.nests.3D1$enddate[i] - days(14)
  
  #' Check the calculated dates
  print(paste("Row", i, ": Startdate =", pa.nests.3D1$startdate[i], ", Enddate =", pa.nests.3D1$enddate[i]))
}

#' Convert 'startdate' and 'enddate' columns to Date format
pa.nests.3D1$startdate <- as.Date(pa.nests.3D1$startdate)
pa.nests.3D1$enddate <- as.Date(pa.nests.3D1$enddate)

#' View the data frame structure
glimpse(pa.nests.3D1)

################################################################################
## Prepare 2D Nest Data 

#' Subset nesting data for 4D in year 2022
pa.nests.2D <- dplyr::filter(pa.nests, WMU =="2D")%>%
  dplyr::select(BandID, CheckDate, NestID, WMU, clutchsize)


################################################################################
## Incubation Data 2D

#' Merge pa.nests.4D and nests.inc, only keep nests that exist in both pa.nests.4D and nests.inc
pa.nests.2D1 <- dplyr::inner_join(pa.nests.2D, nests.inc, by = "NestID") %>%
  dplyr::select(-CheckDate.y) %>%
  dplyr::rename("CheckDate" = CheckDate.x)
glimpse(pa.nests.2D1)

#' Iterate over the rows and subtract clutchsize days
for (i in 1:nrow(pa.nests.2D1)) {
  clutchsize <- pa.nests.2D1$clutchsize[i]  
  startI <- pa.nests.2D1$startI[i]  
  
  #' Create the enddate by subtracting clutchsize (in days) from startI
  pa.nests.2D1$enddate[i] <- startI - clutchsize - days(5)
  
  #' Create the startdate by subtracting 14 days from the enddate
  pa.nests.2D1$startdate[i] <- pa.nests.2D1$enddate[i] - days(14)
  
  #' Check the calculated dates
  print(paste("Row", i, ": Startdate =", pa.nests.2D1$startdate[i], ", Enddate =", pa.nests.2D1$enddate[i]))
}

#' Convert 'startdate' and 'enddate' columns to Date format
pa.nests.2D1$startdate <- as.Date(pa.nests.2D1$startdate)
pa.nests.2D1$enddate <- as.Date(pa.nests.2D1$enddate)

#' View the data frame structure
glimpse(pa.nests.2D1)


################################################################################
## Prepare 5C Nest Data 

#' Subset nesting data for 4D in year 2022
pa.nests.5C <- dplyr::filter(pa.nests, WMU =="5C")%>%
  dplyr::select(BandID, CheckDate, NestID, WMU, clutchsize)


################################################################################
## Incubation Data 5C

#' Merge pa.nests.4D and nests.inc, only keep nests that exist in both pa.nests.4D and nests.inc
pa.nests.5C1 <- dplyr::inner_join(pa.nests.5C, nests.inc, by = "NestID") %>%
  dplyr::select(-CheckDate.y) %>%
  dplyr::rename("CheckDate" = CheckDate.x)
glimpse(pa.nests.5C1)

#' Iterate over the rows and subtract clutchsize days
for (i in 1:nrow(pa.nests.5C1)) {
  clutchsize <- pa.nests.5C1$clutchsize[i]  
  startI <- pa.nests.5C1$startI[i]  
  
  #' Create the enddate by subtracting clutchsize (in days) from startI
  pa.nests.5C1$enddate[i] <- startI - clutchsize - days(5)
  
  #' Create the startdate by subtracting 14 days from the enddate
  pa.nests.5C1$startdate[i] <- pa.nests.5C1$enddate[i] - days(14)
  
  #' Check the calculated dates
  print(paste("Row", i, ": Startdate =", pa.nests.5C1$startdate[i], ", Enddate =", pa.nests.5C1$enddate[i]))
}

#' Convert 'startdate' and 'enddate' columns to Date format
pa.nests.5C1$startdate <- as.Date(pa.nests.5C1$startdate)
pa.nests.5C1$enddate <- as.Date(pa.nests.5C1$enddate)

#' View the data frame structure
glimpse(pa.nests.5C1)


################################################################################
## Connect to Movebank

login <- movebank_store_credentials(username = "Kyle.Smelter",
                                    password="Rayshawks5!",
                                    key="Kyle",
                                    force= T)


################################################################################
## Loops to Download Pre-Nesting GPS Data in Pennsylvania


###########
## WMU 4D

unique.ID.4d<-unique(pa.nests.4D1$NestID)

for (j in 1:length(unique.ID.4d)){
  tmp.subset.4d<-pa.nests.4D1[which(pa.nests.4D1$NestID==unique.ID.4d[j]),]
  tmp.subset.4d$TrackID<-paste(unique.ID.4d[j],seq(1,nrow(tmp.subset.4d),1),sep="_")
  
for(i in 1:nrow(tmp.subset.4d)){
  BirdID<- as.character(tmp.subset.4d[i,1])
  EndDate <- gsub("\\D","", tmp.subset.4d$enddate[i]) 

  StartDate <- gsub("\\D","", tmp.subset.4d$startdate[i]) 
  
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


###########
## WMU 3D 


unique.ID.3d<-unique(pa.nests.3D1$NestID)

for (j in 1:length(unique.ID.3d)){
  tmp.subset.3d<-pa.nests.3D1[which(pa.nests.3D1$NestID==unique.ID.3d[j]),]
  tmp.subset.3d$TrackID<-paste(unique.ID.3d[j],seq(1,nrow(tmp.subset.3d),1),sep="_")
  
  for(i in 1:nrow(tmp.subset.3d)){
    BirdID<- as.character(tmp.subset.3d[i,1])
    EndDate <- gsub("\\D","", tmp.subset.3d$enddate[i]) 
 
    StartDate <- gsub("\\D","", tmp.subset.3d$startdate[i]) 

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


###########
##  WMU 2D 


unique.ID.2d<-unique(pa.nests.2D1$NestID)

for (j in 1:length(unique.ID.2d)){
  tmp.subset.2d<-pa.nests.2D1[which(pa.nests.2D1$NestID==unique.ID.2d[j]),]
  tmp.subset.2d$TrackID<-paste(unique.ID.2d[j],seq(1,nrow(tmp.subset.2d),1),sep="_")
  
  for(i in 1:nrow(tmp.subset.2d)){
    BirdID<- as.character(tmp.subset.2d[i,1])
    EndDate <- gsub("\\D","", tmp.subset.2d$enddate[i]) 

    StartDate <- gsub("\\D","", tmp.subset.2d$startdate[i]) 
   
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


###########
##  WMU 5C 


unique.ID.5c<-unique(pa.nests.5C1$NestID) 

for (j in 1:length(unique.ID.5c)){
  tmp.subset.5c<-pa.nests.5C1[which(pa.nests.5C1$NestID==unique.ID.5c[j]),]
  tmp.subset.5c$TrackID<-paste(unique.ID.5c[j],seq(1,nrow(tmp.subset.5c),1),sep="_")
  
  for(i in 1:nrow(tmp.subset.5c)){
    BirdID<- as.character(tmp.subset.5c[i,1])
    EndDate <- gsub("\\D","", tmp.subset.5c$enddate[i]) 

    StartDate <- gsub("\\D","", tmp.subset.5c$startdate[i]) 

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


################################################################################
## Organize GPS Data

#' Convert move objects to dataframes
full_all_3d <- as.data.frame(full_all_3d)
full_all_4d <- as.data.frame(full_all_4d)
full_all_2d <- as.data.frame(full_all_2d)
full_all_5c <- as.data.frame(full_all_5c)

unique(full_all_5c$individual_local_identifier)

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

################################################################################
###############################################################################X
