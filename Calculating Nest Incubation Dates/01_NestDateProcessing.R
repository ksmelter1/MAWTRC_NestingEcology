###################################################
## Download ACC Data from Specified Date Ranges ##
###################################################

rm(list = ls())
gc()

# load in packages
library(dplyr)

# Load in data ---
nsts <- read.csv("Nest_data/nests_raw_FEB.csv")
colnames(nsts)

# Clean data -----
# We want birds with (nestfound=Y)
nsts_clean <- nsts %>% 
  # Filter out where nests were found
  filter(nestfound %in% "Y") %>%
  dplyr::rename("rawstartsearch" = startsearchdate)


# Great, now we want to find a lists of bird IDs and start date
ids <- unique(nsts_clean$nestid)
ids
dates <- unique(nsts_clean$checkdate)
dates

# Ensuring that checkdate is parsed correctly before adding/subtracting days
# Startsearchdate is the checkdate + 3 days
# Endsearchdate is the checkdate - 50 days 
date_id <- nsts_clean %>% 
  dplyr::select(bandid, nestid, checkdate) %>% 
  dplyr::mutate(checkdate = lubridate::mdy(trimws(checkdate)),  # Removing any extra spaces
                endsearchdate = checkdate - lubridate::days(50),  # Subtract 50 days for endsearchdate
                startsearchdate = checkdate + lubridate::days(3)) # Add 3 days for startsearchdate
head(date_id)

write.csv(date_id, "Nest_data/20250103_dates_KS.csv")
saveRDS(date_id, "Nest_data/20250103_NestingDates_KS.rds")