#'---
#' title: Habitat selection of female wild turkeys during pre-nesting (an SSF analysis)
#' author: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output: 
#'   html_document: 
#'     toc: true
#'---
#' **Purpose**: This script identifies incubation behavior using daily average standard deviation calculations 
#' **Last Updated**: 1/31/2025

###########################################
## Load Packages 

library(data.table)
library(matrixStats)
library(lubridate)
library(tidyr)
library(tidyverse)

###############################################
## Data Prep- Process Csv Files for Each Hen

#' This loop writes a csv file for each bird filtered by bandid that has daily z-axis sd calculated during daylight hours
#' The issue is it doesn't filter properly by bandid and when it iterates, it'll copy over the same data and name the file differently (Checkpoint is files are all the same size)
#' I took the timezone portion from the second part of the script and added it to the loop, this saves time with having to read in files again
#' Other than the issue with subsetting properly this loop works
#' I will run this over the weekend and get separate csvs for each unique band id, I will probably need help with the identifying nest attempts (bottom of script) next week

working.dir <- getwd()


file_paths <- list.files(path = "E:/ACC/Raw", pattern = "\\.csv$", full.names = TRUE)


for (file in file_paths) {
  temp_data <- fread(file)
  unique_bandids <- unique(temp_data$bandid)
  
  for (i in 1:length(unique_bandids)) {
    bandid.tmp <- unique_bandids[i]
    bandid_data <- subset(temp_data, bandid == bandid.tmp)
    bandid_data$transmitter <- as.factor(substr(bandid_data$transmitter, 4, 8))
    bandid_data$timeofday <- ymd_hms(paste(bandid_data$date, bandid_data$time), tz = "UTC")
    bandid_data$timeofday <- with_tz(bandid_data$timeofday, "America/New_York")
    bandid_data <- bandid_data[hour(bandid_data$timeofday) >= 6 & hour(bandid_data$timeofday) <= 19, ]
    bandid_data[, z.sd := rowSds(as.matrix(.SD)), .SDcols = seq(11, 128, 3)]
    bandid_data <- data.table(transmitter = bandid_data$transmitter,
                              z.sd = bandid_data$z.sd,
                              date = bandid_data$date,
                              bandid = bandid_data$bandid)
    output_file <- paste0("E:/ACC/Draft2/", bandid.tmp, "_data.csv")
    counter <- 1
    original_output_file <- output_file
    while (file.exists(output_file)) {
      output_file <- paste0(tools::file_path_sans_ext(original_output_file), "_", 
                            counter,
                            ".",
                            tools::file_ext(original_output_file))
      counter <- counter + 1
    }
    
    fwrite(bandid_data, output_file)
    message("Finished processing bandid: ", bandid.tmp, " in file: ", file)
    
    rm(bandid.tmp, bandid_data)
    gc()  
    
  }
  
  rm(temp_data)
  gc()  
}


################################################################################
## Data Prep-Consolidate Nesting Data 

#' Read in nest csv 
nests <- read.csv("Data Management/Csvs/Processed/Nests/Nests/2025_CleanedNests_2022_2023.csv") 
unique(nests$NestID)


################################################################################
## Check to see if there are missing bandids that didn't download

#' ids <- (unique(nests$bandid)) %>% as.data.frame()
#' files <- (list.files("E:/ACC/Draft2/", pattern = "\\.csv$", full.names = TRUE))
#' bandids <- str_extract(files, "(?<=bandid_)\\d+") %>% as.data.frame()
#' missing_bandids <- setdiff(ids$., bandids$.)
#' nests.missing <- nests %>%
#'   dplyr::filter(bandid %in% missing_bandids) 
#' 
#' #' Use this file to download the rest of the bandids
#' write.csv(nests.missing, "Calculating Nest Incubation Dates/20250124_Nests.Missing.csv")

################################################################################
## Process Data-Get Incubation Start and End Dates

#' Get incubation start and end dates using daily z-axis standard deviation calculations
#' Make sure that only nest attempts where the nest was found are included
#' Create a begin column which truncates the ACC data to 30 days prior to the checkdate + 3 days
nests.1 <- nests %>%
  dplyr::filter(NestBowlFound == "Y") %>%
  dplyr::mutate(CheckDate = as.Date(CheckDate)) %>%
  dplyr::mutate(checkdate3 = CheckDate + lubridate::days(3)) %>%
  dplyr::mutate(begin = checkdate3 - lubridate::days(30)) %>%
  dplyr::rename("bandid" = BandID)
glimpse(nests.1)

#' Creates an empty dataframe with the following information 
#' Transmitter: Character (TransmitterID)
#' nestid: Character (Unique Nest Identifier)
#' startI: Date Object (Start of Incubation)
#' endI: Date Object (End of Incubation)
nest.attemps.id<-list()

#' Initialize an empty list to store `df.prop.15.complete` for each bird
prop_15_complete_list <- list()

#' List all CSV files in the directory using list.files
files <- list.files("E:/ACC/Draft2/", pattern = "\\.csv$", full.names = TRUE)

#' Loop through each file and extract incubation start and end dates
#' For the length of files
#' i is indexed as the length of individuals (Csvs)
#' j is indexed as the length of unique days
for (i in 1:length(files)) {
  
  nest.attemps.id[[i]]<-data.frame(band=character(),
                                   nestid=character(),
                                   startI=as.Date(character()),
                                   endI=as.Date(character()),
                                   stringsAsFactors=FALSE)
  
  bandid_data <- fread(files[i])
  nest.ID.match<-nests.1[which(nests.1$bandid ==bandid_data$bandid[1]),]
  
  for (n in 1:nrow(nest.ID.match)){
    
    if (nrow(nest.ID.match)==1){
      sub<- bandid_data[bandid_data$date >= nest.ID.match[n] & 
                          bandid_data$date <= nest.ID.match$checkdate3[n], ]
    }else{                    
      sub<- bandid_data[bandid_data$date >= max(nest.ID.match$checkdate3[n-1],nest.ID.match$begin[n]) & 
                          bandid_data$date <= nest.ID.match$checkdate3[n], ]                    
    }                   
    
    dates<-unique(sub$date)
    
    prop.15<-rep(0,length(dates))
    
    for (j in 1:length(unique(dates))){
      tmp<-sub[sub$date==dates[[j]],]
      prop.15[j]<-sum(tmp$z.sd<15)/nrow(tmp)
    }
    
    df.prop.15<-data.frame(band=sub$bandid[1],
                           nestid=nest.ID.match$NestID[n],
                           prop.15=prop.15,
                           date=as.Date(dates))
    
    df.prop.15.complete<- complete(df.prop.15, date=seq.Date(max(nest.ID.match$checkdate3[n-1],
                                                                 nest.ID.match$begin[n]), 
                                                             nest.ID.match$checkdate3[n], 
                                                             by="day"))
    
    df.prop.15.complete$startI<-rep(0,nrow(df.prop.15.complete))
    df.prop.15.complete$endI<-rep(0,nrow(df.prop.15.complete))
    
    for (k in 1:(nrow(df.prop.15.complete)-2)){
      df.prop.15.complete$startI[k][df.prop.15.complete$prop.15[k]>=0.85 & 
                                      df.prop.15.complete$prop.15[k+1]>=0.85 & 
                                      df.prop.15.complete$prop.15[k+2]>=0.85]<-1
    }
    startI<-min(which(df.prop.15.complete$startI==1))  
    if(startI=="Inf"){
      for (m in 1:nrow(df.prop.15.complete)-1){
        df.prop.15.complete$endI[m][df.prop.15.complete$prop.15[m]<=0.85 & 
                                      df.prop.15.complete$prop.15[m+1]<=0.85] <-1
      }} else{
        for (m in startI:nrow(df.prop.15.complete)-1){
          df.prop.15.complete$endI[m][(df.prop.15.complete$prop.15[m] <= 0.85 & 
                                         df.prop.15.complete$prop.15[m+1] <= 0.85)| 
                                        df.prop.15.complete$prop.15[m] == 1] <- 1
        } }
    endI<-min(which(df.prop.15.complete$endI==1))
    
    nest.attemps.id[[i]][n,]$band<-as.character(df.prop.15$band[n])
    nest.attemps.id[[i]][n,]$nestid<-as.character(df.prop.15$nestid[n])
    nest.attemps.id[[i]][n,]$startI<-as.Date(df.prop.15.complete$date[startI])
    nest.attemps.id[[i]][n,]$endI<-as.Date(df.prop.15.complete$date[endI])
    
  }
    
    prop_15_complete_list[[paste("bird_",i, "_nest_", n, sep = "")]] <- df.prop.15.complete
    
  }

#' Create object with all nest attempts
#' Convert from a list to a dataframe 
nest.attemps.df <- do.call(rbind, nest.attemps.id)
nest.attemps.df <- nest.attemps.df %>%
  dplyr::rename(bandid = "band") %>%
  dplyr::rename("NestID" = nestid)
nest.attemps.df <- merge(nest.attemps.df, 
                         nests[, c("NestID", "CheckDate", "NestFate")], 
                         by = "NestID", all.x = TRUE)

#' View all nests with NAs
na_nests <- nest.attemps.df[is.na(nest.attemps.df$startI) | is.na(nest.attemps.df$endI), ]
prop_15_complete.df<- do.call(rbind, prop_15_complete_list)
na.nests.prop15 <- prop_15_complete.df %>%
  dplyr::filter(nestid %in% na_nests$NestID)

#' Filter NA nests out of data
filtered_nest.attemps.all.df <- tidyr::drop_na(nest.attemps.df)

#' Constrain data to only include nests in which the start of incubation is not equal to the end
#' The CheckDate is greater than or equal to the end of incubation 
#' The end of incubation is greater than the start of incubation (Didn't have observations of this)
filtered_nest.attemps.check.df <- filtered_nest.attemps.all.df %>%
  dplyr::filter(startI != endI) %>%
  dplyr::filter(CheckDate >= endI) %>%
  dplyr::filter(endI >= startI)


################################################################################
## Output Data 

write.table(filtered_nest.attemps.check.df,
            "Data Management/Csvs/Processed/IncubationDates/20250131_NestAttempts_allbirds.csv",
            sep = ",",
            row.names = FALSE, 
            col.names = TRUE)


