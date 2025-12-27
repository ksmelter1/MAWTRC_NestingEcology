#'---
#' title: Habitat selection of female wild turkeys during pre-nesting (an SSF analysis)
#' author: K. Smelter
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#'---
#'  
#' **Purpose**: This script downloads movement data associated with each hens nesting attempt from movebank 
#' **Last Updated**: 12/27/25


################################################################################
## Load Packages 

# Vector of package names
packages <- c("purrr",
              "lubridate",
              "dplyr",
              "move2",
              "tidyverse")

# Function to load a package or install it if not already installed
load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

# Apply the function to each package name
lapply(packages, load_packages)


# Read in nests csv
nj.nests <- read_csv("Data Management/Csvs/Processed/Nests/Nests/New Jersey/20250219_CleanedNests_2022_2023_NJ.csv")
nj.nests

# Constrain to nests that were used (Case = 1) and not randomly generated nests available to hens
nj.nests <- nj.nests %>%
  dplyr::filter(Case == "1")
nj.nests

# Read in captures csv
captures <- read_csv("Data Management/Csvs/Raw/Captures/captures_nj.csv")
captures

################################################################################
## Data Management

# Change 99s into NA Values (There were 99s in this data unlike PA and MD)
# Originally the database was coded to fill NAs with 99s
nj.nests$EggsHatched[nj.nests$EggsHatched == 99] <- NA
nj.nests$EggsUnhatched[nj.nests$EggsUnhatched == 99] <- NA
nj.nests$EggsDestroyed[nj.nests$EggsDestroyed == 99] <- NA

# Create clutch size column and remove unnecessary column
# Clutch size is a minimum count 
nj.nests <- nj.nests %>%
  dplyr::mutate(clutchsize = rowSums(select(., EggsHatched, EggsDestroyed, EggsUnhatched), na.rm = TRUE)) 
glimpse(nj.nests)


################################################################################
## Incubation Data NJ

# Csv from incubation start and end script
nests.inc <- read_csv("Data Management/Csvs/Processed/Incubation Dates/New Jersey/20250717_NestAttempts_allbirds_NJ_Ready.csv")
nests.inc

# Read in sample from known fate model
sample <- read_csv("Samples/New Jersey/NestingSample_NJ.updated.csv")
sample  

# Bind columns together
# Create BirdID column
nests.inc <- right_join(nests.inc, sample) %>%
  dplyr::mutate(BirdID = str_extract(NestID, "\\d{4}")) 

# Rename columns
captures <- captures %>%
  dplyr::rename("BandID" = bandid) %>%
  dplyr::mutate("BandID" = as.character(BandID))

# Join columns together 
nj.nests.test <- dplyr::inner_join(nests.inc, captures, by = "BandID")

# Merge pa.nests.4D and nests.inc, only keep nests that exist in both pa.nests.4D and nests.inc
nj.nests1 <- dplyr::inner_join(nj.nests,nj.nests.test, by = "NestID") %>%
  dplyr::select(-CheckDate.y, -...1) %>%
  dplyr::rename("CheckDate" = CheckDate.x)
glimpse(nj.nests1)

# Iterate over the rows and subtract clutchsize days
for (i in 1:nrow(nj.nests1)) {
  clutchsize <- nj.nests1$clutchsize[i]  
  startI <- nj.nests1$startI[i]  
  
  # Create the enddate by subtracting clutchsize (in days) from startI
  nj.nests1$enddate[i] <- startI - clutchsize - days(5)
  
  # Create the startdate by subtracting 14 days from the enddate
  nj.nests1$startdate[i] <- nj.nests1$enddate[i] - days(14)
  
  # Check the calculated dates
  print(paste("Row", i, ": Startdate =", nj.nests1$startdate[i], ", Enddate =", nj.nests1$enddate[i]))
}

# Convert 'startdate' and 'enddate' columns to Date format
nj.nests1$startdate <- as.Date(nj.nests1$startdate)
nj.nests1$enddate <- as.Date(nj.nests1$enddate)

# Rename columns
nj.nests1 <- nj.nests1 %>%
  dplyr::rename("BandID" = BandID.y)

# Group nests by study area for loops below
nj.nests1.SJ <- nj.nests1 %>%
  dplyr::filter(studyarea == "SJ") %>%
  dplyr::select(-BirdID) %>%
  dplyr::mutate(BirdID = BandID) %>%
  dplyr::select(BirdID,BandID, startI, endI, clutchsize, NestID, startdate, enddate)
nj.nests1.NJ <- nj.nests1 %>%
  dplyr::filter(studyarea == "NJ") %>%
  dplyr::select(BirdID, BandID, startI, endI, clutchsize, NestID, startdate, enddate)


################################################################################
## Connect to Movebank

login <- movebank_store_credentials(username = "Kyle.Smelter",
                                    password="Rayshawks5!",
                                    key="Kyle",
                                    force= T)


################################################################################
## Loops to Download Pre-Nesting GPS Data in New Jersey


################################################################################
## NJ South

unique.ID.SJ<-unique(nj.nests1.SJ$NestID)

for (j in 1:length(unique.ID.SJ)){
  tmp.subset.SJ<-nj.nests1.SJ[which(nj.nests1.SJ$NestID==unique.ID.SJ[j]),]
  tmp.subset.SJ$TrackID<-paste(unique.ID.SJ[j],seq(1,nrow(tmp.subset.SJ),1),sep="_")
  
for(i in 1:nrow(tmp.subset.SJ)){
  BirdID<- as.character(tmp.subset.SJ[i,1])
  EndDate <- gsub("\\D","", tmp.subset.SJ$enddate[i]) 

  StartDate <- gsub("\\D","", tmp.subset.SJ$startdate[i]) 
  
  Year <- lubridate::year(tmp.subset.SJ$startI[i])

  
  dat.SJ<- movebank_download_study(study ="Wild Turkey New Jersey South", 
                                   login = login,
                                   individual_local_identifier= BirdID,
                                   timestamp_start= StartDate,
                                   timestamp_end= EndDate,
                                   removeDuplicatedTimestamps=T)
  mt_track_id(dat.SJ)<-rep(tmp.subset.SJ$TrackID[i],nrow(dat.SJ))
  
  if(exists("full_all_SJ")){ 
    full_all_SJ <- rbind(full_all_SJ, dat.SJ)
    
    
  }else{
    full_all_SJ <- dat.SJ
  }
}
}


################################################################################
## NJ North

unique.ID.NJ<-unique(nj.nests1.NJ$NestID)

for (j in 1:length(unique.ID.NJ)){
  tmp.subset.NJ<-md.nests1.WM[which(nj.nests1.NJ$NestID==unique.ID.NJ[j]),]
  tmp.subset.NJ$TrackID<-paste(unique.ID.NJ[j],seq(1,nrow(tmp.subset.NJ),1),sep="_")
  
  for(i in 1:nrow(tmp.subset.NJ)){
    BirdID<- as.character(tmp.subset.NJ[i,1])
    EndDate <- gsub("\\D","", tmp.subset.NJ$enddate[i]) 
    
    StartDate <- gsub("\\D","", tmp.subset.NJ$startdate[i]) 
    
    Year <- lubridate::year(tmp.subset.NJ$startI[i])
    
    
    dat.NJ<- movebank_download_study(study ="Wild Turkey Maryland West", 
                                     login = login,
                                     individual_local_identifier= BirdID,
                                     timestamp_start= StartDate,
                                     timestamp_end= EndDate,
                                     removeDuplicatedTimestamps=T)
    mt_track_id(dat.NJ)<-rep(tmp.subset.NJ$TrackID[i],nrow(dat.NJ))
    
    if(exists("full_all_NJ")){ 
      full_all_NJ <- rbind(full_all_NJ, dat.NJ)
      
      
    }else{
      full_all_NJ <- dat.NJ
    }
  }
}


################################################################################
## Organize GPS Data

# Convert move objects to dataframes
# No data for NJ North yet
full_all_SJ <- as.data.frame(full_all_SJ)
full_all_NJ <- as.data.frame(full_all_NJ)

# No NJ North data so proceed with this
df.all <- full_all_SJ
  
  
  # Create df with all study areas
  df.all <- rbind(full_all_SJ, 
                  full_all_NJ) %>%
    dplyr::rename("BirdID"= individual_local_identifier) 
  
# Separate geometry lat and longs into separate columns and create new dataframe
# Organize timestamp to be formatted in year, month, day, hour, minutes, seconds
# Map function applies a function to each element of a vector
hens.all <- df.all%>%
  mutate(long = unlist(map(df.all$geometry,1)),
         lat = unlist(map(df.all$geometry,2))) %>%
  dplyr::rename("BirdID" = individual_local_identifier) %>%
  dplyr::select(BirdID, timestamp,long, lat) %>%
  st_drop_geometry()

################################################################################

# Save GPS data for step generation 
# Saved here: Data Management/RData/Pre-Nesting Movement Model/New Jersey/Covariates/Movebank_NJ_buffer.RData

################################################################################
###############################################################################X
