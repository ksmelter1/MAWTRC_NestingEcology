#---
# title: Capture Data Management
# author: "K. Smelter"
# date: "`r format(Sys.time(), '%d %B, %Y')`"
# output:
#   html_document: 
#     toc: true
#---
#  
# **Purpose**: This script creates an output csv containing hen capture data and a shapefile for hen capture locations
#
################################################################################

# Load in capture csv 
caps <- read_csv("Data Management/Csvs/Raw/Captures/captures_pa.csv")
caps

# Filter data to contain only hen data
caps.f <- caps %>%
  dplyr::filter(sex=="F") %>%
  dplyr::filter(captyr != "2025") 

# Create a capture date column
# This data is now clean
  caps.f <- caps.f %>%
  dplyr::mutate(CaptureDate = lubridate::make_date(year = captyr, month = captmo, day = captday)) %>%
  dplyr::rename(BandID = bandid,
                Lat = lat,
                Long = long,
                WMU = studyarea,
                Age = age) %>%
    dplyr::select(BandID, CaptureDate, WMU, Lat, Long, Age) %>%
    dplyr::select(BandID, Age, WMU, CaptureDate, everything())

# Output capture data  
write_csv(caps.f, "Data Management/Csvs/Processed/20250629_PAHenCaptures_2022_2023_2024.csv")

################################################################################
## Create Map of Capture Sites

# Selecting columns of interest for site map
# Filtering distinct values for study area map
caps.site.map <- dplyr::select(caps.f, sitename, lat, long, studyarea) %>%
  distinct()

# Convert capture sites to sf object
caps.sf <- st_as_sf(caps.site.map, coords= c("long", "lat"), crs=4326, remove=F)
mapview(caps.sf)

# Write shapefile for distinct capture sites
st_write(caps.sf, "capturelocations.shp")
