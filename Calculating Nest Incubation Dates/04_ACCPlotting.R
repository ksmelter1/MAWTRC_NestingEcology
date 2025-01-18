
# Part for ACC Plotting R script -------------

#####################################################################X
## Create a new dataframe to reformat for acc plotting R script
data.accplot <- combined_results

## Subset by study area (IF WANTED, OTHERWISE COMMENT OUT)
# data.accplot <- subset(data.accplot, studyarea == "2D")

## Subset the data to eastern time dates and times
data.accplot$datetime <- strptime(paste(data.accplot$date, data.accplot$time), "%Y-%m-%d %H:%M:%S", tz="UTC")
data.accplot$datetime <- as.POSIXct(data.accplot$datetime)
data.accplot$datetime <- with_tz(data.accplot$datetime, "America/New_York")

# Optional: Further filter the data by date range
# data.accplot <- subset(data.accplot, datetime >= "2024-04-06 00:00:00" & datetime <= "2024-04-10 23:59:59") 

## Remove unneeded columns
data.accplot <- subset(data.accplot, select = -c(sex, age, studyarea, state, datetime))

## Add column for day of the week
data.accplot$day <- weekdays(data.accplot$date)

## Reorder the columns to match the order in the raw acc text file
data.accplot <- data.accplot %>% relocate(day, .before = time)

## Reformat the date field to match the format in raw acc text file
data.accplot$date <- format(data.accplot$date, format = "%d.%m.%Y")

## Split dataframe by transmitter and save each as a separate txt file
ids <- unique(data.accplot$transmitter)

for(i in seq_along(ids)){
  df <- subset(data.accplot, transmitter == ids[i])
  tagname <- ids[i]
  write.table(df, paste0("D:/PA_TurkeyDataProcessing/SQL_ACC_Dwnld/", tagname, "_acc.txt"), sep = ',', row.names = FALSE, col.names = FALSE)
}

# Confirm that the files have been saved
cat("ACC plot data files have been saved.\n")


# For when you want one bird ----
# query database for records between a set of dates ## REMEMBER THIS IS IN UTC NOT EASTERN TIME
q1 <- paste0("
            SELECT *
            FROM accelerometer
            WHERE date >= '2022-04-06' AND date <= '2022-05-26;
            ")

data.sql <- dbGetQuery(conn,q1)

dbDisconnect(conn)  #close database  
