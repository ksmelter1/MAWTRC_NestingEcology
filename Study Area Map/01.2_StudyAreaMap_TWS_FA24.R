######################
## Study Area Map ##
######################
#'---
#' title: Pennsylvania Study Area Map
#' author: "K. Smelter
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'
#+ include = FALSE
knitr::opts_chunk$set(warning = FALSE, message = FALSE, cache = TRUE)
#'  
#' **Purpose**: This script produces a study area map for the Pennsylvania Wild Turkey Project


#' Vector of package names
packages <- c("dplyr",
              "ggspatial",
              "ggplot2",
              "mapview",
              "sf",
              "rmapshaper")

#' Function to load a package or install it if not already installed
load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

#' Apply the function to each package 
lapply(packages, load_packages)

#' Set font to Times New Roman
windowsFonts(Times=windowsFont("Times New Roman"))

#' Use the tigris package to filter through states and remove GEOIDs above 60
#' GEOIDs above 60 are territories and islands so they are being removed for scaling
st <- tigris::states() %>%
  dplyr::filter(GEOID < "60") %>% 
  tigris::shift_geometry()

#' Transform to albers equal area conic projection, epsg code is 5070
st <- st_transform(st, 4326)

#' Grab outline of PA
pa.outline <- subset(st, st$NAME=="Pennsylvania")
plot(st_geometry(pa.outline))

#' Pennsylvania
pa.wmus <- st_read("Data Management/Shapefiles/Pennsylvania WMUs/PGC_BNDWildlifeManagementUnits2021.shp") %>%
  st_transform(4326)%>%
  rmapshaper::ms_simplify(keep = 0.002, weighting = 1)

#' The hen study takes place here in PA in 2D, 3D, 4D, and 5C
study.area <- subset(pa.wmus, WMU_ID=="2D"| WMU_ID=="3D"| WMU_ID=="4D"| WMU_ID=="5C")

#' Read in capture sites shapefile
capture.sites <- st_read("Data Management/Shapefiles/capture sites/capturelocations.shp")

#' Map data
map_data <- st_zm(pa.wmus, drop=TRUE, what="ZM")
ggplot()+
  geom_sf(data = map_data)

map_data <- map_data %>% st_transform(4326)

#' Visualize map
map <- ggplot() +
  geom_sf(data = pa.wmus, color = "gray30", fill="grey90") +
  geom_sf(data = map_data, color = "black") +
  geom_sf(data = capture.sites, shape=17, size = 3, colour="black", show.legend = F) +
  theme_void() 
#' View map
map
