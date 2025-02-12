
#'---
#' title: Calculating Start and End Incubation Dates
#' author: K. Smelter, F. Buderman, V. Winter, K. Lamp
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output: 
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script downloads ACC data off the VM by individuals and dates
#' **Last Updated**: 1/31/2025


rm(list = ls())
gc()

################################################################################
## load in packages

library(dplyr)
library(lubridate)
library(dplyr)
library(DBI)
library(RPostgres)

################################################################################
## Connect to SQL Turkey Database

DBI::dbDriver('Postgres')
drv <- dbDriver("Postgres")
db <- "PAturkey"  
host_db <- "localhost"
db_port <- 5432
db_user <- 'postgres'
db_password <- 'PAturk3y'


conn <- dbConnect(drv, dbname = db, host = host_db, 
                  port = db_port, user = db_user, password = db_password)


dbListTables(conn)

################################################################################
## Read in Nest Files

nests <- readRDS("Nest_data/Nesting_KJS/20250103_NestingDates_KS.rds") 

bird_ids <- unique(nests$bandid)  

daterange <- nests %>% 
  dplyr::filter(bandid %in% bird_ids) 

start_dates <- daterange$startsearchdate
end_dates <- daterange$endsearchdate

output_file_path <- "D:/PA_TurkeyDataProcessing/Nest_data/Nesting_KJS/20250124_RetrievedMissingBirds.csv"

headers_written <- FALSE

################################################################################
## Loop to Download ACC Data by Individual and Date 


start <- Sys.time()
for (bird_id in bird_ids) {
  for (i in seq_along(start_dates)) {
    #' Build the SQL query
    query <- paste0("
      SELECT *
      FROM accelerometer
      WHERE date >= '", start_dates[i], "' AND date <= '", end_dates[i], "'
      AND bandid = '", bird_id, "';
    ")
    
    #' Execute the query and store the result
    result <- dbGetQuery(conn, query)
    
    #' Check if the result is not empty
    if (nrow(result) > 0) {
      #' Clean the result: remove rows with NA in "studyarea" and create "datetime"
      result <- result[!is.na(result$studyarea),]
      result$datetime <- as.POSIXct(strptime(paste(result$date, result$time), 
                                             "%Y-%m-%d %H:%M:%S", tz="UTC"))
      
      #' Order rows by "datetime"
      result <- result[order(result$datetime),]
      
      #' Remove "datetime" column before saving
      result <- subset(result, select = -c(datetime))
      
      #' Write data to CSV
      if (!headers_written) {
        #' Write the first chunk with headers
        write.csv(result, file = output_file_path, row.names = FALSE)
        headers_written <- TRUE  
      } else {
        #' Append the data without headers for subsequent chunks
        write.table(result, file = output_file_path, append = TRUE, sep = ",", 
                    col.names = FALSE, row.names = FALSE)
      }
    }
  }
}
end <- Sys.time()

#' Disconnect from the database after processing
dbDisconnect(conn)

#' Processing time 
print(paste("Total time:" (end - start)))

#' Confirm the file has been saved
cat("Data has been saved to:", output_file_path, "\n")

#' Test read in
#' use this to read in large csv's
library(data.table)
df <- fread(output_file_path) 
