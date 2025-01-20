#'---
#' title: Habitat selection of female wild turkeys during pre-nesting (an SSF analysis)
#' author: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output: 
#'   html_document: 
#'     toc: true
#'---
#'

######################
## Load Packages 

library(data.table)
library(matrixStats)
library(lubridate)
library(tidyr)

###############################################
## Data Prep- Process Csv Files for Each Hen

#' This loop writes a csv file for each bird filtered by bandid that has daily z-axis sd calculated during daylight hours
#' The issue is it doesn't filter properly by bandid and when it iterates, it'll copy over the same data and name the file differently (Checkpoint is files are all the same size)
#' I took the timezone portion from the second part of the script and added it to the loop, this saves time with having to read in files again
#' Other than the issue with subsetting properly this loop works
#' I will run this over the weekend and get separate csvs for each unique band id, I will probably need help with the identifying nest attempts (bottom of script) next week

#' Get working directory
working.dir <- getwd()

#' List all CSV files in the directory using list.files
file_paths <- list.files(path = "E:/ACC/Raw", pattern = "\\.csv$", full.names = TRUE)

#' Loop through each file
for (file in file_paths) {
  
  #' Read in the large CSV file
  temp_data <- fread(file)
  
  #' Get unique bandid values within the file
  unique_bandids <- unique(temp_data$bandid)
  
  #' Loop through each unique bandid
  for (i in 1:length(unique_bandids)) {
    
    #' Get the current bandid
    bandid.tmp <- unique_bandids[i]
    
    #' Filter the data for the current bandid
    bandid_data <- subset(temp_data, bandid == bandid.tmp)
    
    #' Extract a substring from the transmitter ID column, characters in position 4-8
    bandid_data$transmitter <- as.factor(substr(bandid_data$transmitter, 4, 8))
    
    #' Convert from UTC to America/New York time zone
    bandid_data$timeofday <- ymd_hms(paste(bandid_data$date, bandid_data$time), tz = "UTC")
    bandid_data$timeofday <- with_tz(bandid_data$timeofday, "America/New_York")
    
    #' Constrain dataset to daylight hours (6 AM to 7 PM)
    bandid_data <- bandid_data[hour(bandid_data$timeofday) >= 6 & hour(bandid_data$timeofday) <= 19, ]
    
    #' Calculate rowwise standard deviation for the z-axis columns
    bandid_data[, z.sd := rowSds(as.matrix(.SD)), .SDcols = seq(11, 128, 3)]
    
    #' Create a new data table with the required columns
    bandid_data <- data.table(transmitter = bandid_data$transmitter,
                              z.sd = bandid_data$z.sd,
                              date = bandid_data$date,
                              bandid = bandid_data$bandid)
    
    #' Construct the output file name for this bandid
    output_file <- paste0("E:/ACC/processed/bandid_", bandid.tmp, "_data.csv")
    
    #' Check if the file exists already and modify the name if needed
    counter <- 1
    original_output_file <- output_file
    
    #' Loop to increment counter until a unique file name is found
    while (file.exists(output_file)) {
      output_file <- paste0(tools::file_path_sans_ext(original_output_file), "_", 
                            counter,
                            ".",
                            tools::file_ext(original_output_file))
      counter <- counter + 1
    }
    
    #' Write the filtered data to a separate CSV file
    fwrite(bandid_data, output_file, path = "E:/ACC/")
    
    
    #' Checkpoint
    message("Finished processing bandid: ", bandid.tmp, " in file: ", file)
    
    #' Clean up memory
    rm(bandid.tmp, bandid_data)
    gc()  
    
  }
  
  #' Clean up memory after processing the entire file
  rm(temp_data)
  gc()  
}


#######################################
## Data Prep-Consolidate Nesting Data 

#' Read in nest csv 
nests <- read.csv("Data Management/Csvs/Raw/nests_raw_FEB.csv") %>%
  dplyr::filter(nestfound == "Y")

#' Format days so that there are two digits per each day
nests$checkday <- stringr::str_pad(nests$checkday, width = 2, pad = "0")

#' Format months so that there are two digits per each month
nests$checkmo <- stringr::str_pad(nests$checkmo, width = 2, pad = "0")

#' Paste day, month, and year columns together to create checkdate
nests$checkdate <- paste0(nests$checkyr, 
                          nests$checkmo, 
                          nests$checkday)

#' Remove any rows that contain NAs from the nests df
nests <- nests %>%tidyr::drop_na()

#' Format checkdate column as a date object
#' Checkdate3 is 3 days from the checkdate
nests$checkdate <- as.Date(nests$checkdate, format = "%Y%m%d")
nests$checkdate3 <- nests$checkdate + days(3)

#' 50 days prior to the nest being checked 
#' Truncate ACC data for each hen 50 days before the check date
nests$begin<-nests$checkdate - days(50)

unique(nests$nestid)

###################################################
## Process Data-Get Incubation Start and End Dates

#' Get incubation start and end dates using daily z-axis standard deviation calculations
#' Will just put these into the nests raw file 


#' Creates an empty dataframe with the following information 
#' Transmitter: Character (TransmitterID)
#' nestid: Character (Unique Nest Identifier)
#' startI: Date Object (Start of Incubation)
#' endI: Date Object (End of Incubation)
nest.attemps.id<-list()

#' Initialize an empty list to store `df.prop.15.complete` for each bird
prop_15_complete_list <- list()

#' List all CSV files in the directory using list.files
files <- list.files("E:/ACC/Processed/zaxis_calcs/", pattern = "\\.csv$", full.names = TRUE)

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
  nest.ID.match<-nests[which(nests$bandid ==bandid_data$bandid[1]),]
  
  for (n in 1:nrow(nest.ID.match)){
    
    if (nrow(nest.ID.match)==1){
      sub<- bandid_data[bandid_data$date >= nest.ID.match$begin[n] & 
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
                           nestid=nest.ID.match$nestid[n],
                           prop.15=prop.15,
                           date=as.Date(dates))
    
    df.prop.15.complete<- complete(df.prop.15, date=seq.Date(max(nest.ID.match$checkdate3[n-1],
                                                                 nest.ID.match$begin[n]), 
                                                             nest.ID.match$checkdate3[n], 
                                                             by="day"))
    
    df.prop.15.complete$startI<-rep(0,nrow(df.prop.15.complete))
    df.prop.15.complete$endI<-rep(0,nrow(df.prop.15.complete))
    
    for (k in 1:(nrow(df.prop.15.complete)-2)){
      df.prop.15.complete$startI[k][df.prop.15.complete$prop.15[k]>=0.8 & 
                                      df.prop.15.complete$prop.15[k+1]>=0.8 & 
                                      df.prop.15.complete$prop.15[k+2]>=0.8]<-1
    }
    startI<-min(which(df.prop.15.complete$startI==1))  
    if(startI=="Inf"){
      for (m in 1:nrow(df.prop.15.complete)-1){
        df.prop.15.complete$endI[m][df.prop.15.complete$prop.15[m]<=0.8 & 
                                      df.prop.15.complete$prop.15[m+1]<=0.8] <-1
      }} else{
        for (m in startI:nrow(df.prop.15.complete)-1){
          df.prop.15.complete$endI[m][(df.prop.15.complete$prop.15[m] <= 0.8 & 
                                         df.prop.15.complete$prop.15[m+1] <= 0.8)| 
                                        df.prop.15.complete$prop.15[m] == 1] <- 1
        } }
    endI<-min(which(df.prop.15.complete$endI==1))  
    
    
    if (is.infinite(endI)) {
      checkdate <- nest.ID.match$checkdate[n]  
      endI <- which(df.prop.15.complete$date == checkdate)  
    }
    
    nest.attemps.id[[i]][n,]$band<-as.character(df.prop.15$band[n])
    nest.attemps.id[[i]][n,]$nestid<-as.character(df.prop.15$nestid[n])
    nest.attemps.id[[i]][n,]$startI<-as.Date(df.prop.15.complete$date[startI])
    nest.attemps.id[[i]][n,]$endI<-as.Date(df.prop.15.complete$date[endI])
    
  }
    
    #' Store df.prop.15.complete for each bird in the list
    prop_15_complete_list[[paste("bird_",i, "_nest_", n, sep = "")]] <- df.prop.15.complete
    
  }

#' Create dataframe
nest.attemps.df <- do.call(rbind, nest.attemps.id)

#' Rename column
nest.attemps.df <- nest.attemps.df %>%
  dplyr::rename(bandid = "band") 

#' Merge 'nests' with 'nest.attemps.df' using 'nestid'
nest.attemps.df <- merge(nest.attemps.df, 
                         nests[, c("nestid", "checkdate", "nestfate")], 
                         by = "nestid", all.x = TRUE)

#' Filter nest.attempts.df for rows with NA in startI or endI
na_nests <- nest.attemps.df[is.na(nest.attemps.df$startI) | is.na(nest.attemps.df$endI), ]

#' Create dataframe
prop_15_complete.df<- do.call(rbind, prop_15_complete_list)


#' Filter 'nests' to only include rows where 'checkdate' matches any 'bandid' in 'nest.attemps.df'
na.nests.prop15 <- prop_15_complete.df %>%
  dplyr::filter(nestid %in% na_nests$nestid)

#' Drop all NA values
filtered_nest.attemps.all.df <- tidyr::drop_na(nest.attemps.df)

#' Filter the dataframe to only include rows where the 'checkdate' is within 14 days of 'endI'
filtered_nest.attemps.14.df <- nest.attemps.df %>%
  dplyr::filter(checkdate <= (endI + 14) & checkdate >= endI) 


############################
## Output Data 

write.table(filtered_nest.attemps.df, "NestAttempts_allbirds.csv", sep = ",", row.names = FALSE, col.names = TRUE)
write.table(filtered_nest.attemps.df, "NestAttempts_allbirdsfiltered.csv", sep = ",", row.names = FALSE, col.names = TRUE)
