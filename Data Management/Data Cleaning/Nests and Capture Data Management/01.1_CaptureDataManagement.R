#'---
#' title: Capture Data Management
#' author: "K. Smelter"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'
#+ include = FALSE
knitr::opts_chunk$set(warning = FALSE, message = FALSE, cache = TRUE)
#'  
#' **Purpose**: This script creates an output csv containing hen capture data and a shapefile for hen capture locations
#'

#' Load in capture csv 
caps <- read.csv("Csvs/raw data/captures.csv")

#' Filter data to contain only hen data
caps.f <- dplyr::filter(caps, sex=="F")

#' Write csv
write.csv(caps.f, "HenCaptures.csv")

#' Selecting columns of interest for site map
#' Filtering distinct values for study area map
caps.site.map <- dplyr::select(caps.f, sitename, lat, long, studyarea) %>%
  distinct()

#' Convert capture sites to sf object
caps.sf <- st_as_sf(caps.site.map, coords= c("long", "lat"), crs=4326, remove=F)
mapview(caps.sf)

#' Write shapefile for distinct capture sites
st_write(caps.sf, "capturelocations.shp")
