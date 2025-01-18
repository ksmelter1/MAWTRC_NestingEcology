rm(list = ls())
gc()

# load in packages
library(dplyr)
library(lubridate)
library(dplyr)
library(DBI)
library(RPostgres)

# Load in SQL ----
DBI::dbDriver('Postgres')
drv <- dbDriver("Postgres")
db <- "PAturkey"  
host_db <- "localhost"
db_port <- 5432
db_user <- 'postgres'
db_password <- 'PAturk3y'

# Connect to the drive
conn <- dbConnect(drv, dbname = db, host = host_db, 
                  port = db_port, user = db_user, password = db_password)

# Check connection and tables
dbListTables(conn)

# Load the RDS file w. the birds IDs and date ranges
# All dates are date objects
nests <- readRDS("Nest_data/20250103_NestingDates_KS.rds")
str(nests$endsearchdate)
str(nests$startsearchdate)
str(nests$checkdate)

# Set up unique bird IDs and date range
bird_ids <- unique(nests$bandid)[1]

# Range of dates for the birds above
# Change dates to characters
daterange <- nests %>%
  dplyr::filter(bandid %in% bird_ids) 

# Define start and end dates for the loop
start_dates <- daterange$startsearchdate
end_dates <- daterange$endsearchdate

# Define output file path
output_file_path <- "D:/PA_TurkeyDataProcessing/Nest_data/20241109_PA_nesting_4255_birds.csv"

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
library(data.table)
df <- fread(output_file_path) # use this to read in large csv's