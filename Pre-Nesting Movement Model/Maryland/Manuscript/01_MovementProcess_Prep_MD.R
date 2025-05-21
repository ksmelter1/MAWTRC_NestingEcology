
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
#' **Last Updated**: 2/25/25


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
md.nests <- read_csv("Data Management/Csvs/Processed/Nests/Nests/Maryland/20250219_CleanedNests_2022_2023_MD.csv")
md.nests

md.nests <- md.nests %>%
  dplyr::filter(Case == "1")
md.nests

captures <- read_csv("Data Management/Csvs/Raw/Captures/captures_md.csv")
captures

################################################################################
## Data Management

#' Change 99s into NA Values
md.nests$EggsHatched[md.nests$EggsHatched == 99] <- NA
md.nests$EggsUnhatched[md.nests$EggsUnhatched == 99] <- NA
md.nests$EggsDestroyed[md.nests$EggsDestroyed == 99] <- NA

#' Create clutch size column and remove unnecessary column
#' Clutch size is a minimum count 
md.nests <- md.nests %>%
  dplyr::mutate(clutchsize = rowSums(select(., EggsHatched, EggsDestroyed, EggsUnhatched), na.rm = TRUE)) 
glimpse(md.nests)

################################################################################
## Incubation Data Maryland

#' Csv from incubation start and end script
nests.inc <- read_csv("Data Management/Csvs/Processed/Incubation Dates/Maryland/20250131_NestAttempts_allbirds_MD.csv")

sample <- read_csv("Samples/Maryland/NestingSample_MD.csv")
nests.inc <- right_join(nests.inc, sample)

md.nests.test <- dplyr::inner_join(nests.inc, captures, by = "bandid")

#' Merge pa.nests.4D and nests.inc, only keep nests that exist in both pa.nests.4D and nests.inc
md.nests1 <- dplyr::inner_join(md.nests,md.nests.test, by = "NestID") %>%
  dplyr::select(-CheckDate.y, -...1) %>%
  dplyr::rename("CheckDate" = CheckDate.x)
glimpse(md.nests1)

# Iterate over the rows and subtract clutchsize days
for (i in 1:nrow(md.nests1)) {
  clutchsize <- md.nests1$clutchsize[i]  
  startI <- md.nests1$startI[i]  
  
  # Create the enddate by subtracting clutchsize (in days) from startI
  md.nests1$enddate[i] <- startI - clutchsize - days(5)
  
  # Create the startdate by subtracting 14 days from the enddate
  md.nests1$startdate[i] <- md.nests1$enddate[i] - days(14)
  
  # Check the calculated dates
  print(paste("Row", i, ": Startdate =", md.nests1$startdate[i], ", Enddate =", md.nests1$enddate[i]))
}

# Convert 'startdate' and 'enddate' columns to Date format
md.nests1$startdate <- as.Date(md.nests1$startdate)
md.nests1$enddate <- as.Date(md.nests1$enddate)
glimpse(md.nests1)

md.nests1.EM <- md.nests1 %>%
  dplyr::filter(studyarea == "EM")

md.nests1.WM <- md.nests1 %>%
  dplyr::filter(studyarea == "WM")


################################################################################
## Connect to Movebank

login <- movebank_store_credentials(username = "Kyle.Smelter",
                                    password="Rayshawks5!",
                                    key="Kyle",
                                    force= T)


################################################################################
## Maryland GPS EM

unique.ID.EM<-unique(md.nests1.EM$NestID)

for (j in 1:length(unique.ID.EM)){
  tmp.subset.EM<-md.nests1.EM[which(md.nests1.EM$NestID==unique.ID.EM[j]),]
  tmp.subset.EM$TrackID<-paste(unique.ID.EM[j],seq(1,nrow(tmp.subset.EM),1),sep="_")
  
for(i in 1:nrow(tmp.subset.EM)){
  BirdID<- as.character(tmp.subset.EM[i,1])
  EndDate <- gsub("\\D","", tmp.subset.EM$enddate[i]) 

  StartDate <- gsub("\\D","", tmp.subset.EM$startdate[i]) 
  
  Year <- lubridate::year(tmp.subset.EM$startI[i])

  
  dat.EM<- movebank_download_study(study ="Wild Turkey Maryland East", 
                                   login = login,
                                   individual_local_identifier= BirdID,
                                   timestamp_start= StartDate,
                                   timestamp_end= EndDate,
                                   removeDuplicatedTimestamps=T)
  mt_track_id(dat.EM)<-rep(tmp.subset.EM$TrackID[i],nrow(dat.EM))
  
  if(exists("full_all_EM")){ 
    full_all_EM <- rbind(full_all_EM, dat.EM)
    
    
  }else{
    full_all_EM <- dat.EM
  }
}
}


################################################################################
## Maryland GPS WM

unique.ID.WM<-unique(md.nests1.WM$NestID)

for (j in 1:length(unique.ID.WM)){
  tmp.subset.WM<-md.nests1.WM[which(md.nests1.WM$NestID==unique.ID.WM[j]),]
  tmp.subset.WM$TrackID<-paste(unique.ID.WM[j],seq(1,nrow(tmp.subset.WM),1),sep="_")
  
  for(i in 1:nrow(tmp.subset.WM)){
    BirdID<- as.character(tmp.subset.WM[i,1])
    EndDate <- gsub("\\D","", tmp.subset.WM$enddate[i]) 
    
    StartDate <- gsub("\\D","", tmp.subset.WM$startdate[i]) 
    
    Year <- lubridate::year(tmp.subset.WM$startI[i])
    
    
    dat.WM<- movebank_download_study(study ="Wild Turkey Maryland West", 
                                     login = login,
                                     individual_local_identifier= BirdID,
                                     timestamp_start= StartDate,
                                     timestamp_end= EndDate,
                                     removeDuplicatedTimestamps=T)
    mt_track_id(dat.WM)<-rep(tmp.subset.WM$TrackID[i],nrow(dat.WM))
    
    if(exists("full_all_WM")){ 
      full_all_WM <- rbind(full_all_WM, dat.WM)
      
      
    }else{
      full_all_WM <- dat.WM
    }
  }
}


################################################################################
## Organize GPS Data

#' Convert move objects to dataframes
full_all_EM <- as.data.frame(full_all_EM)
full_all_WM <- as.data.frame(full_all_WM)

  #' Create df with all study areas
  df.all <- rbind(full_all_EM, 
                  full_all_WM) %>%
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
