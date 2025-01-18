#'---
#' title: Nest and Vegetation Survey Data Management
#' author: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'
#+ include = FALSE
knitr::opts_chunk$set(warning = FALSE, message = FALSE, cache = TRUE)
#'  
#' **Purpose**: This script shows the amount of hens that nested across years and organizes nest data and shapefiles for analysis

#####################
## Load Packages ##
####################

require(dplyr)
require(ggplot2)

###################
## Load in Data ##
###################

#' Load in raw nest csv file with data from 2022-2023
#' This is just nests, no veg sampling or no nest values
nests.raw <- read.csv("Data Management/Csvs/raw data/nests_raw.csv")

#' Load in raw nest vegetation csv file with data from 2022-2023 surveys
#' This file also contains original nest data 
#' For nest fates we will have to merge the two 
veg.raw <- read.csv("Data Management/Csvs/raw data/nest_vegetation_raw.csv")

#######################################################
## Check Nesting Attempts per Each Hen Across Years ##
#######################################################

#' Create a dataframe and remove "No Nest" values as nest success models are reliant upon the nest being found
nest_df<- nests.raw%>% dplyr::mutate(col1=stringr::str_detect(nestid,"NoNest"))%>%
  dplyr::filter(!col1)%>%
  dplyr::select(-col1)

#' Subset nest dataframe by grouping it by nest check year, then keeping columns BandID and NestID
nest.year <- nest_df %>%
  dplyr::group_by(checkyr)%>%
  dplyr::select(checkyr, bandid, nestid, lat_n, long_n,wmu_n, checkday, nestfate)

#' Using dplyr to filter birds that nested in both 2022 and 2023
nests.both <- nest.year %>% dplyr::select(bandid, checkyr)%>%
  dplyr::mutate(x=1) %>%
  distinct() %>%
  tidyr::pivot_wider( names_from = checkyr, values_from = x, 
                      values_fill = 0 ) %>%
  dplyr::filter(`2022`==1 & `2023`==1) %>%
  pull(bandid) 

#' Create dataframe 
nests.year.both <- dplyr::filter(nest.year,bandid %in% nests.both)
unique(nests.year.both$bandid)

## 16 hens nested across 2022 and 2023
## We determined this was not as large of a sample as we needed to look into the dominant hen hypothesis

#############################################
## Clean Nesting Data from 2022 and 2023 ##
#############################################

#' Fill empty cells with 99 as a proxy for NA
#' These are in columns that are structured as characters to not throw off math
nest_df$nestsubfate1[nest_df$nestsubfate1== ""] <- as.character("99")
nest_df$nestsubfate2[nest_df$nestsubfate2 == ""] <- as.character("99")
nest_df$unhatchids[nest_df$unhatchids == ""] <- as.character("99")
nest_df$nestcomments[nest_df$nestcomments == ""] <- as.character("99")

#' Create cleaned nest df and remove unhatched ID and nest found columns as nest success models rely upon
#' the assumption that the nest was found 
#' Remove nests that have longitude values greater than zero
threshold_long <- 0
nests_clean <- nest_df %>%
  dplyr::select(-unhatchids, -nestfound) %>%
  dplyr::filter(nest_df$long_n <= threshold_long)
  
#####################################################
## Clean Nest Vegetation Data from 2022 and 2023 ##
#####################################################

#' Subset veg dataframe 
veg_df <- veg.raw %>%
  dplyr::select( -county_v, -township_v, -landownership_v, 
                 -transmitterid_v)

#' Fill empty cells with 99 as a proxy for NA
veg_df$vegcomments[veg_df$vegcomments== ""] <- as.character("99")
veg_df$woodytype1[veg_df$woodytype1== ""] <- as.character("99")
veg_df$woodytype2[veg_df$woodytype2== ""] <- as.character("99")
veg_df$guardvo[veg_df$guardvo== ""] <- as.character("99")

#' Remove no nest values from nestid column 
veg_clean <- veg_df %>% dplyr::mutate(col1=stringr::str_detect(nestid_v,"NoNest"))%>%
  dplyr::filter(!col1)%>%
  dplyr::select(-col1)%>%
  dplyr::filter(long_v <= threshold_long)%>%
  dplyr::filter(lat_v !=70.72393)%>%
  dplyr::group_by(plottype)
  
  veg_clean$case <- ifelse(veg_clean$plottype =="North"|
                             veg_clean$plottype == "South"|
                             veg_clean$plottype == "East"|
                             veg_clean$plottype == "West", "0","1")

#' Create unique identifier column 
veg_clean$uniqueID <- with(veg_clean, paste0(nestid_v, plottype))

#' Write cleaned vegetation csv
write.csv(veg_clean, "nests.veg_22_23.csv")
  
  
##########################
## Create Sf Dataframes ##
##########################

#' Create spatial dataframe for merged nest and vegetation datasets  
nests.veg.sf <- veg_clean %>%
  st_as_sf(coords = c("long_v", "lat_v"), crs = 4326, remove = FALSE) %>%
  st_transform(9822) #Albers equal area conic epsg code to highlight region
#' check
#mapview(nests.veg.sf)

#' Create spatial dataframe for nests
nests.only.sf <- nests_clean %>%
  st_as_sf(coords = c("long_n", "lat_n"), crs = 4326, remove = FALSE) %>%
  st_transform(9822) #Albers equal area conic epsg code to highlight region
#check
mapview(nests.only.sf)


################################
## Write csvs and shapefiles ##
################################

#' Shapefiles
#' Cleaned nest shapefile
sf::st_write(nests.sf, ".", layer = "nests_cleaned2022-2023.shp", 
             driver = "ESRI Shapefile", append = T )
#' Cleaned vegetation survey shapefile
sf::st_write(veg.sf, ".", layer = "veg_cleaned2022-2023.shp",  
             driver = "ESRI Shapefile", append = T )
#' Cleaned merged vegetation and nest data
sf::st_write(nests.veg.sf, ".", layer = "nest.veg.locations_cleaned.shp", 
             driver = "ESRI Shapefile", append = T)


#' Csvs
#' Cleaned nest csv
#write.csv(nests_clean, "nests_22_23_clean.csv")
#' Cleaned vegetation survey csv
#write.csv(veg_clean, "veg_22_23_clean.csv")
#' Cleaned merged vegetation and nest data
#write.csv(nest.veg.df, "nests.veg_22_23_clean.csv")

################################################################################
################################################################################

###################
## Merge Data ##
###################


#' If else statement to consolidate North, South, East and West to Available locations
#' the straight line means or 
# veg_clean$case <- ifelse(veg_clean$plottype =="North"|
#          veg_clean$plottype == "South"|
#          veg_clean$plottype == "East"|
#          veg_clean$plottype == "West", "Available","Used")

veg_clean$case <- as.character(veg_clean$case)

#' Rename column for merge                           
veg_clean <- veg_clean %>%
  dplyr::rename("nestid"=nestid_v)
  
#' Check values to make sure they match
table(veg_clean$plottype)
table(veg_clean$case)


#' Merge cleaned nest data with cleaned veg data
#' Create analysis fate column where nests that reach hatch get a 1 and everything else gets a 0
#' We define a nest as successful if the nest reaches hatch 
nests.veg <- merge(nests_clean,
                   veg_clean,
                   by= "nestid") %>%
                   dplyr::select(case, 
                                 uniqueID, checkyr, 
                                 checkmo, checkday,wmu_v,
                                 lat_v, long_v,
                                 averagevo, averagemaxvo,
                                 bandid.x, percgrassforb, percwoody, nestfate,
                                 litterht, nestsubfate1, nestsubfate2, eggshatched,
                                 eggsunhatched) %>%
                  dplyr::rename("birdid"= bandid.x) %>%
                  dplyr::rename("nestid"=uniqueID)%>%
                  dplyr::rename("lat"=lat_v)%>%
                dplyr::rename("long"=long_v)%>%
                dplyr::rename("wmu"=wmu_v)%>%
                dplyr::rename("year"=checkyr)%>%
                dplyr::rename("day_checked"=checkday)%>%
                dplyr::rename("month_checked"= checkmo)%>%
                dplyr::mutate(analysisfate = ifelse(nestfate == "Hatched", 1, 0)) %>%
                dplyr::arrange(nestid, case)

##########################
## Create Date Columns ##
##########################

#' Clean up issue with zeros in dates column
nests.veg$day_checked<-ifelse(nchar(nests.veg$day_checked)==1,paste(0,nests.veg$day_checked,sep=""),
                           nests.veg$exposure)
nests.veg$month_checked<-ifelse(nchar(nests.veg$month_checked)==1,paste(0,nests.veg$month_checked,sep=""),nests.veg$month_checked)

#' Build dataframe with information we need
nests.veg$checkdate <- paste0(nests.veg$year, 
                              nests.veg$month_checked, 
                              nests.veg$day_checked) 

#' Format timestamps correctly in year, month, day format and create start date column 
nests.veg$checkdate<-as.Date(as.character(nests.veg$checkdate),
                               format="%Y%m%d")
#' Subset dataframe
nests.veg <- dplyr::select(nests.veg, -checkmo, -checkday)

################################################################################