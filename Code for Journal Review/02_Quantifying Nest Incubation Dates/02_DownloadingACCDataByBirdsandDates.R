#---
# title: Downloading ACC Data For Specified Birds and Dates
# author: "V. Winter, K. Smelter
# date: "`r format(Sys.time(), '%d %B, %Y')`"
# output:
#   html_document: 
#     toc: true
#---
#  
# **Purpose**: This script downloads ACC data needed to obtain infomation about start and termination of incubation

################################################################################
## Load Packages and Data

# Initialize workspace
rm(list = ls())
gc()

# load in packages
library(dplyr)
library(lubridate)
library(dplyr)
library(DBI)
library(RPostgres)

################################################################################
## Connect to Turkey DB

# Each state has its own DB
# PA: PAturkey
# MD: MDturkey
# NJ: NJturkey
# Password is the same for each 

# Load in SQL ----
DBI::dbDriver('Postgres')
drv <- dbDriver("Postgres")
db <- "NJturkey"  # Switch database for each state 
host_db <- "localhost"
db_port <- 5432
db_user <- 'postgres'
db_password <- 'PAturk3y'

# Connect to the drive
conn <- dbConnect(drv, dbname = db, host = host_db, 
                  port = db_port, user = db_user, password = db_password)

# Check connection and tables
dbListTables(conn)

################################################################################
## Data Prep

# Load the RDS file w. the birds IDs and date ranges
nests <- readRDS("Nest_data/Nesting_KJS/NJ/20260109_NestingRanges_2025_NJ.rds")
head(nests$BandID)

# Set up unique bird IDs and date ranges
# Download 25 individuals at a time due to large sizes
bird_ids <- unique(nests$BandID) [1:25]

#... filter dates that match the bird IDs specified above
daterange <- nests %>% 
  dplyr::filter(BandID %in% bird_ids) 

# Define start and end dates for the loop
start_dates <- daterange$startsearchdate
end_dates <- daterange$endsearchdate

################################################################################
## Query to get the total number of observations for an individual

# dbGetQuery(conn, "
#       SELECT COUNT(*)
#       FROM accelerometer
#      WHERE bandid = '1013';   
#    ")

################################################################################
## Query to get minimum and maximum dates from each bird in database

# This was used to check whether ACC data existed in Maryland and New Jersey in 2025
# 1/9/25: No ACC data yet for 2025 past March 1st

# # Get unique list
# nj_birds <- unique(as.character(daterange$BandID))
# bird_list_sql <- paste0("'", nj_birds, "'", collapse = ",")
# 
# # Query to see if data exists from 2025 birds
# nj_accel_range <- dbGetQuery(
#   conn,
#   paste0("
#     SELECT bandid,
#            MIN(date) AS min_dt,
#            MAX(date) AS max_dt
#     FROM accelerometer
#     WHERE bandid IN (", bird_list_sql, ")
#     GROUP BY bandid
#     ORDER BY bandid;
#   ")
# )
# 
# # Get range of values for each bird
# nj_accel_range <- nj_accel_range %>%
#   mutate(
#     min_dt = as.Date(min_dt),
#     max_dt = as.Date(max_dt),
#     accel_year_min = format(min_dt, "%Y"),
#     accel_year_max = format(max_dt, "%Y"),
#     overlaps_2025 = max_dt >= as.Date("2025-03-01") # No ACC data prior to March 1st in 2025
#   )
# 
# View(nj_accel_range)

# Write csv 
# write.csv(
#   nj_accel_range,
#   "D:/PA_TurkeyDataProcessing/Nest_data/Nesting_KJS/NJ/NJ_accelerometer_date_ranges.csv",
#   row.names = FALSE
# )

################################################################################
## Run loop to pull ACC data

# Define output file path
output_file_path <- "D:/PA_TurkeyDataProcessing/Nest_data/Nesting_KJS/NJ/20260109_nestingbirds_1_30_KJS.csv"   # Change csv name after each run to match bird range

# Initialize a variable to track if headers have been written in CSV
headers_written <- FALSE

# Loop over each bird ID and date range, and write each result directly to CSV
start <- Sys.time()
for (bird_id in bird_ids) {
  for (i in seq_along(start_dates)) {
    # Build the SQL query
    query <- paste0("
      SELECT *
      FROM accelerometer
      WHERE date >= '", start_dates[i], "' AND date <= '", end_dates[i], "'
      AND bandid = '", bird_id, "';
    ")
    
    # Execute the query and store the result
    result <- dbGetQuery(conn, query)
    
    # Check if the result is not empty
    if (nrow(result) > 0) {
      # Clean the result: remove rows with NA in "studyarea" and create "datetime"
      result <- result[!is.na(result$studyarea),]
      result$datetime <- as.POSIXct(strptime(paste(result$date, result$time), 
                                             "%Y-%m-%d %H:%M:%S", tz="UTC"))
      
      # Order rows by "datetime"
      result <- result[order(result$datetime),]
      
      # Remove "datetime" column before saving
      result <- subset(result, select = -c(datetime))
      
      # Write data to CSV
      if (!headers_written) {
        # Write the first chunk with headers
        write.csv(result, file = output_file_path, row.names = FALSE)
        headers_written <- TRUE  # Set to TRUE after writing headers
      } else {
        # Append the data without headers for subsequent chunks
        write.table(result, file = output_file_path, append = TRUE, sep = ",", 
                    col.names = FALSE, row.names = FALSE)
      }
    }
  }
}
end <- Sys.time()

# Disconnect from the database after processing
dbDisconnect(conn)

# Processing time 
print(paste("Total time:" (end - start)))

# Confirm the file has been saved
cat("Data has been saved to:", output_file_path, "\n")

# Test read in
# library(data.table)
# df <- fread(output_file_path) # use this to read in large csv's

################################################################################
###############################################################################X