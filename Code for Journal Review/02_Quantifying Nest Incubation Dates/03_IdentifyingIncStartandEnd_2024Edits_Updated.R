#---
#' title: Incubation Start and End Derivation 
#' author: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output: 
#'   html_document: 
#'     toc: true
#---
#' **Purpose**: This script identifies incubation behavior using daily average standard deviation calculations 
#' **Last Updated**: 1/31/2025

################################################################################
## Load Packages 

library(data.table)
library(matrixStats)
library(lubridate)
library(tidyr)
library(tidyverse)

################################################################################
## Data Prep- Process Csv Files for Each Hen

# This loop writes a csv file for each bird filtered by bandid that has daily z-axis sd calculated during daylight hours
# The issue is it doesn't filter properly by bandid and when it iterates, it'll copy over the same data and name the file differently (Checkpoint is files are all the same size)
# I took the timezone portion from the second part of the script and added it to the loop, this saves time with having to read in files again
# Other than the issue with subsetting properly this loop works
# I will run this over the weekend and get separate csvs for each unique band id, I will probably need help with the identifying nest attempts (bottom of script) next week

working.dir <- getwd()
file_paths <- list.files(path = "E:/ACC/Raw/PA/", pattern = "\\.csv$", full.names = TRUE)


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
    output_file <- paste0("E:/ACC/Processed/PA/Updated/", bandid.tmp, "_data.csv")
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

# Read in cleaned nest csv 
# Just use the used nests
nests <- read.csv("Data Management/Csvs/Processed/Nests/Nests/Pennsylvania/20250629_CleanedNests_2022_2023_2024.csv") %>%
  dplyr::filter(Case == "1")


################################################################################
## Sanity Check: Assess if there are missing bandids that didn't download

# ids <- (unique(nests$bandid)) %>% as.data.frame()
# files <- (list.files("E:/ACC/Draft2/", pattern = "\\.csv$", full.names = TRUE))
# bandids <- str_extract(files, "(?<=bandid_)\\d+") %>% as.data.frame()
# missing_bandids <- setdiff(ids$., bandids$.)
# nests.missing <- nests %>%
#   dplyr::filter(bandid %in% missing_bandids) 
# 
# # Use this file to download the rest of the bandids
# write.csv(nests.missing, "Calculating Nest Incubation Dates/20250124_Nests.Missing.csv")

################################################################################
## Process Data-Get Incubation Start and End Dates

# Get incubation start and end dates using daily z-axis standard deviation calculations
# Make sure that only nest attempts where the nest was found are included
# Create a begin column which truncates the ACC data to 30 days prior to the checkdate + 3 days
nests.1 <- nests %>%
  mutate(CheckDate = as.Date(CheckDate, format = "%m/%d/%Y")) %>%
  mutate(checkdate3 = CheckDate + 3,
         begin = checkdate3 - 30) %>%
  rename("bandid" = BandID)
glimpse(nests.1)

# Create an empty list to store incubation attempt data for each individual bird
nest.attemps.id <- list()

# Initialize an empty list to store dataframes with complete prop.15 time series for each bird-nest combo
prop_15_complete_list <- list()

# Get list of CSV files from directory (excluding ones ending in _1.csv)
files <- list.files("E:/ACC/Processed/PA/Updated/",
                    pattern = "\\.csv$",          # Only .csv files
                    full.names = TRUE)            # Include full file paths

# Exclude files with "_1" in their names (likely duplicates or not of interest)
# files <- files[!str_detect(basename(files), "_1")]

# Loop over each file, where each file corresponds to an individual bird
for (i in 1:length(files)) {
  
  # Initialize a dataframe to hold incubation attempts for the current individual
  nest.attemps.id[[i]] <- data.frame(
    band = character(),                   # Bird ID
    nestid = character(),                 # Nest ID
    startI = as.Date(character()),        # Start of incubation
    endI = as.Date(character()),          # End of incubation
    stringsAsFactors = FALSE
  )
  
  # Read the current individual's data
  bandid_data <- fread(files[i])
  
  # Match the individual's nest data from an external reference dataframe `nests.1`
  nest.ID.match <- nests.1[which(nests.1$bandid == bandid_data$bandid[1]), ]
  
  # Loop over each nest attempt by this individual
  for (n in 1:nrow(nest.ID.match)) {
    
    # If only one nesting attempt exists
    if (nrow(nest.ID.match) == 1) {
      sub <- bandid_data[bandid_data$date >= nest.ID.match$begin[n] & 
                           bandid_data$date <= nest.ID.match$checkdate3[n], ]
    } else {
      # If multiple nesting attempts, avoid overlap using previous checkdate3
      sub <- bandid_data[bandid_data$date >= max(nest.ID.match$checkdate3[n - 1], nest.ID.match$begin[n]) & 
                           bandid_data$date <= nest.ID.match$checkdate3[n], ]
    }
    
    # Extract unique days of observation
    dates <- unique(sub$date)
    prop.15 <- rep(0, length(dates))  # Initialize proportion <15 vector
    
    # Loop through each date and calculate proportion of z.sd < 15
    for (j in 1:length(dates)) {
      tmp <- sub[sub$date == dates[[j]], ]
      prop.15[j] <- sum(tmp$z.sd < 15) / nrow(tmp)
    }
    
    # Create a dataframe with proportion and date info
    df.prop.15 <- data.frame(
      band = sub$bandid[1],
      nestid = nest.ID.match$NestID[n],
      prop.15 = prop.15,
      date = as.Date(dates)
    )
    
    # Generate a complete date sequence for the nesting period
    complete_dates <- seq.Date(
      from = max(nest.ID.match$checkdate3[n - 1], nest.ID.match$begin[n]),
      to = nest.ID.match$checkdate3[n],
      by = "day"
    )
    
    # Fill in missing dates (if any) in the dataframe
    df.prop.15.complete <- complete(df.prop.15, date = complete_dates)
    
    # Initialize columns to mark start and end of incubation
    df.prop.15.complete$startI <- 0
    df.prop.15.complete$endI <- 0
    
    # Incubation Start Detection: 3 consecutive days with prop.15 ≥ 0.85 
    for (k in 1:(nrow(df.prop.15.complete) - 2)) {
      if (!is.na(df.prop.15.complete$prop.15[k]) &&
          !is.na(df.prop.15.complete$prop.15[k + 1]) &&
          !is.na(df.prop.15.complete$prop.15[k + 2]) &&
          df.prop.15.complete$prop.15[k] >= 0.85 &&
          df.prop.15.complete$prop.15[k + 1] >= 0.85 &&
          df.prop.15.complete$prop.15[k + 2] >= 0.85) {
        df.prop.15.complete$startI[k] <- 1  # Mark potential start
      }
    }
    
    # Get the first occurrence of a marked start
    startI <- suppressWarnings(min(which(df.prop.15.complete$startI == 1)))
    
    ##  Incubation End Detection: 3 consecutive days with prop.15 ≤ 0.80 
    if (is.finite(startI)) {
      # Clip dataframe to the nesting window
      df_clip <- df.prop.15.complete[df.prop.15.complete$date <= nest.ID.match$checkdate3[n], ]
      df_clip$endI <- 0  # Initialize endI column
      endI <- NA  # Initialize end index
      
      for (m in startI:(nrow(df_clip) - 2)) {
        if (!is.na(df_clip$prop.15[m]) &&
            !is.na(df_clip$prop.15[m + 1]) &&
            !is.na(df_clip$prop.15[m + 2]) &&
            df_clip$prop.15[m] <= 0.80 &&
            df_clip$prop.15[m + 1] <= 0.80 &&
            df_clip$prop.15[m + 2] <= 0.80) {
          df_clip$endI[m] <- 1  # Mark potential end
          endI <- m             # Save the index
          break                 # Stop at first match
        }
      }
      
      # If no end detected, default to last date in window
      if (is.na(endI)) {
        endI <- which(df_clip$date == nest.ID.match$checkdate3[n])
        df_clip$endI[endI] <- 1
      }
      
      # Save results in master list for the current bird
      nest.attemps.id[[i]][n, ]$band <- as.character(df.prop.15$band[1])
      nest.attemps.id[[i]][n, ]$nestid <- as.character(df.prop.15$nestid[1])
      nest.attemps.id[[i]][n, ]$startI <- as.Date(df_clip$date[startI])
      nest.attemps.id[[i]][n, ]$endI <- as.Date(df_clip$date[endI])
      
      # Save complete dataframe with marked start/end to output list
      prop_15_complete_list[[paste("bird_", i, "_nest_", n, sep = "")]] <- df_clip
      
    } else {
      # If no incubation start was detected, store NAs
      nest.attemps.id[[i]][n, ]$band <- as.character(df.prop.15$band[1])
      nest.attemps.id[[i]][n, ]$nestid <- as.character(df.prop.15$nestid[1])
      nest.attemps.id[[i]][n, ]$startI <- NA
      nest.attemps.id[[i]][n, ]$endI <- NA
      
      # Still save the completed time series without start/end indicators
      prop_15_complete_list[[paste("bird_", i, "_nest_", n, sep = "")]] <- df.prop.15.complete
    }
  }
}

# Create object with all nest attempts
# Convert from a list to a dataframe 
nest.attemps.df <- do.call(rbind, nest.attemps.id)
nest.attemps.df <- nest.attemps.df %>%
  dplyr::rename(bandid = "band") %>%
  dplyr::rename("NestID" = nestid) %>%
  dplyr::distinct(NestID, .keep_all = TRUE)

# Merge incubation ranges with nest data  
nest.attemps.df <- merge(nest.attemps.df, 
                         nests[, c("NestID", "CheckDate", "NestFate")], 
                         by = "NestID", all.x = TRUE)

# View all nests with NAs
na_nests <- nest.attemps.df[is.na(nest.attemps.df$startI) | is.na(nest.attemps.df$endI), ]
prop_15_complete.df<- do.call(rbind, prop_15_complete_list)
na.nests.prop15 <- prop_15_complete.df %>%
  dplyr::filter(nestid %in% na_nests$NestID)

# Filter NA nests out of data
filtered_nest.attemps.all.df <- tidyr::drop_na(nest.attemps.df)

# Make date formats match
filtered_nest.attemps.all.df$CheckDate <- as.Date(filtered_nest.attemps.all.df$CheckDate, format = "%m/%d/%Y")

# Constrain data to only include nests in which the start of incubation is not equal to the end
# The CheckDate is greater than or equal to the end of incubation 
# The end of incubation is greater than the start of incubation (Didn't have observations of this)
# The end of incubation must be within 14 days of the day the nest was checked
filtered_nest.attemps.check.df <- filtered_nest.attemps.all.df %>%
  dplyr::mutate(endI = dplyr::if_else(endI > CheckDate, CheckDate, endI)) %>%
  dplyr::filter(startI != endI) %>%
  dplyr::filter(CheckDate >= endI) %>%
  dplyr::filter(endI >= startI) 

################################################################################
## Output Data 

# Output nesting data for each state 
# write.table(filtered_nest.attemps.check.df,
#             "Data Management/Csvs/Processed/Incubation Dates/Pennsylvania/20250709_NestAttempts_allbirds_PA_Ready.csv",
#             sep = ",",
#             row.names = FALSE, 
#             col.names = TRUE)

################################################################################
###############################################################################X

