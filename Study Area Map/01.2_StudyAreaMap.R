######################
## Study Area Map ##
######################
##'---
#' title: Study Area Map for "A Conditional Behavioral Process Model for Wild Turkeys in Pennsylvania"
#' authors: "K. Smelter, F. Buderman
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'
#+ include = FALSE
knitr::opts_chunk$set(warning = FALSE, message = FALSE, cache = TRUE)
#'  
#' **Purpose**: This script creates a study area map for "A Framework for Analyzing Wild Turkey Summer Sighting Data" Manuscript


#####################
## Load Packages ##
####################

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


############################
## Shapefile Management ##
############################

#' Read in wildlife management unit shapefiles
#' Transform all shapefiles to Albers (5070)
#' st_zm removes the z dimensions
#' ms_simplify: Uses mapshaper to simplify polygons

#' Pennsylvania
pa.wmus <- st_read("Data Management/Shapefiles/Pennsylvania WMUs/PGC_BNDWildlifeManagementUnits2021.shp") %>%
  st_transform(5070)%>%
  rmapshaper::ms_simplify(keep = 0.002, weighting = 1)

#' The hen study takes place here in PA in 2D, 3D, 4D, and 5C in PA
pa.study.area <- subset(pa.wmus, WMU_ID=="2D"| WMU_ID=="3D"| WMU_ID=="4D"| WMU_ID=="5C")

#' Read in Pennsylvania capture sites shapefile
capture.sites <- st_read("Data Management/Shapefiles/capture sites/capturelocations.shp")

#' New Jersey
nj.thas <-  st_read("Data Management/Shapefiles/New Jersey Turkey Hunting Areas/NJ_THA.shp") %>%
  st_transform(5070) %>%
  st_zm()%>%
  rmapshaper::ms_simplify(keep = 0.002, weighting = 1)

#' Check
mapview(nj.thas)

#' Maryland
md.wmus <- st_read("Data Management/shapefiles/Maryland Turkey Regions/turkey_regions.shp") %>%
  st_set_crs(26985) %>%
  st_transform(5070) %>%
  dplyr::rename("Region1"=Region) %>%
  dplyr::group_by(Region1) %>%
  summarise() %>%
  rmapshaper::ms_simplify(keep=0.002, weighting = 1)

#' Check
mapview(md.wmus)

################
## Outlines ##
################

#' Grab outline of USA states
#' Shapefile from Alberto that crops out the great lakes
sa.outline <- st_read("Data Management/shapefiles/USA_no_greatlakes/USA_adm1_without_great_lakes.shp")%>%
  st_transform(5070)%>%
  rmapshaper::ms_simplify(keep = 0.002, weighting = 1)

#' Check structure
str(sa.outline)

#' Crop map to study area by using st_bbox
box_zoom = (st_bbox(c(xmin= 822757.4, ymin= 1776761, xmax= 1891072, ymax= 2558539))) %>%
  st_as_sfc() %>% st_set_crs(5070)

#' Check
ggplot()+
  geom_sf(data= sa.outline)+
  geom_sf(data= box_zoom)


#' Crop map to our study area
zoomed <- st_intersection(sa.outline, box_zoom)
plot(st_geometry(zoomed))


#' Combine shapefiles into one sf object
combined_sf <- bind_rows(pa.wmus, md.wmus, nj.thas)

#' Map data
map_data <- st_zm(pa.wmus, drop=TRUE, what="ZM")
ggplot()+
  geom_sf(data = map_data)


map <- ggplot() +
  geom_sf(data = zoomed, color = "black", fill="white", linewidth= 0.4)+
  geom_sf(data = combined_sf, color = "grey48", fill="grey90")  +
  geom_sf(data = zoomed, color = "black", fill="transparent", linewidth= 0.4)+
  geom_sf(data = capture.sites, color = "black", fill="transparent", linewidth= 0.4)+
  theme_classic() +
  theme(axis.text = element_text(size = 11, family="sans"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  annotation_north_arrow(location = 'br',
                         style = north_arrow_nautical(text_family = "sans"),
                         which_north = "false") +
  annotation_scale(style = "bar",
                   pad_x = unit(0.05, "in"),
                   pad_y = unit(0.05, "in"),
                   text_family = "sans")
#' Plot map
map

#' View map
map
