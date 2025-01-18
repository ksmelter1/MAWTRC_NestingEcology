
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
#' **Purpose**: This script creates a study area map for "A Conditional Behavioral Process Model for Wild Turkeys in Pennsylvania" Manuscript


#####################
## Load Packages ##
####################

#' Vector of package names
packages <- c("dplyr",
              "ggspatial",
              "ggplot2",
              "mapview",
              "sf",
              "rmapshaper",
              "DBI",
              "RODBC")

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
################################################################################

########################
## Capture site data ##
########################

#################
## Pennsylvania

# Read in tables from access database
db.pa <- file.path("Databases/TurkeyDB.accdb")
ch.pa <- odbcConnectAccess2007(db.pa) # open channel to database
dcapt.pa <- sqlFetch(ch.pa, "captures", as.is = TRUE)
close(ch.pa) # close channel to database

###############
## Maryland

# Read in tables from access database
db.md <- file.path("Databases/MD_TurkeyDB_080124_CRL.accdb")
ch.md <- odbcConnectAccess2007(db.md) # open channel to database
dcapt.md <- sqlFetch(ch.md, "captures", as.is = TRUE)
close(ch.md) # close channel to database

#' Create an object with unique captures
#' Select columns we need
#' Filter by distinct values
md.unique.caps <- dcapt.md %>%
  dplyr::select(lat, long, sitename, captyr) %>%
  dplyr::filter(sitename != "AKRIDGE WIDGEON PT FARM")

#' Create sf object
#' Project to Albers
md.caps.sf <- st_as_sf(md.unique.caps, coords = c("long", "lat"), crs= 4326) %>%
  st_transform(5070)

#' Check
mapview(md.caps.sf)

################
## New Jersey

# Read in tables from access database
db.nj <- file.path("Databases/NJ_TurkeyDB080724_KJL.accdb")
ch.nj <- odbcConnectAccess2007(db.nj) # open channel to database
dcapt.nj <- sqlFetch(ch.nj, "captures", as.is = TRUE)
close(ch.nj) # close channel to database

#' Create an object with unique captures
#' Select columns we need
#' Filter by distinct values
nj.unique.caps <- dcapt.nj %>%
  dplyr::select(lat, long, sitename, captyr) %>%
  dplyr::filter(lat > 14.00252)

#' Create sf object
#' Project to Albers
nj.caps.sf <- st_as_sf(nj.unique.caps, coords = c("long", "lat"), crs= 4326) %>%
  st_transform(5070)

#' Check
mapview(nj.caps.sf)

########################
## Pennsylvania

#' Read in Pennsylvania capture sites shapefile
pa.caps.sf <- st_read("Data Management/Shapefiles/capture sites/capturelocations.shp") %>%
  st_transform(5070)

#' Merge together
capture.sites <- bind_rows(nj.caps.sf, md.caps.sf, pa.caps.sf)
################################################################################

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


#################
## New Jersey ##
#################

#' Read in shapefile data 
nj.thas <-  st_read("Data Management/Shapefiles/New Jersey Turkey Hunting Areas/NJ_THA.shp") %>%
  st_transform(5070) %>%
  st_zm()%>%
  rmapshaper::ms_simplify(keep = 0.002, weighting = 1)

#' Check
mapview(nj.thas)

#' Consolidate turkey hunting areas
nj.thas$THA <- ifelse(nj.thas$THA<11,"North",
                      ifelse(nj.thas$THA>10 & nj.thas$THA<15,"Central","South")) 



#' Summarise allows you to filter by a group's geometries
nj.thas.ready <-nj.thas %>% 
  dplyr::dplyr::group_by(THA) %>%
  summarise()

# Remove the first row
nj.thas.ready <- nj.thas.ready[-1, ]

#' Check
mapview(nj.thas.ready, zcol="THA")
mapview(nj.thas)+mapview(nj.caps.sf)

################
## Maryland ##
################

#' Maryland wmus (Thanks Jake!)
#' Project to Albers
#' Read in shapefile 
md.wmus <- st_read("Data Management/Shapefiles/Maryland Turkey Regions/Sa.map/md_wmu_boundaries.gpkg") %>%
  st_transform(5070) %>%
  rmapshaper::ms_simplify(keep = 0.007, weighting = 1) %>%
  dplyr::rename("geometry"=geom)

#' Check
mapview(md.wmus)

#' Filter the rows to combine
rows_2combine <- md.wmus %>% 
  filter(wmu %in% c("Upper Eastern Shore", "Lower Eastern Shore"))

#' Combine the filtered rows into a new geometry
combined_row <- rows_2combine %>%
  summarise(geometry = st_union(geometry), wmu = "Eastern Shore Combined")

#' Remove the original rows and bind the new combined row
md.wmus.ready <- md.wmus %>%
  filter(!wmu %in% c("Upper Eastern Shore", "Lower Eastern Shore")) %>%
  bind_rows(combined_row) 

#' Check
mapview(md.wmus.ready)
mapview(md.wmus.ready)+mapview(md.caps.sf)

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
#' Will need to adjust a little bit to cut it east
box_zoom = (st_bbox(c(xmin= 903751, ymin= 1776600, xmax= 1871955, ymax= 2508599))) %>%
  st_as_sfc() %>% st_set_crs(5070)

#' Check
ggplot()+
  geom_sf(data= sa.outline)+
  geom_sf(data= box_zoom)


#' Crop map to our study area
zoomed <- st_intersection(sa.outline, box_zoom)
plot(st_geometry(zoomed))


#' Combine shapefiles into one sf object
combined_sf <- bind_rows(pa.wmus, md.wmus.ready, nj.thas.ready)

#' Map data
map_data <- st_zm(combined_sf, drop=TRUE, what="ZM")
ggplot()+
  geom_sf(data = map_data)

#' Visualize map
map <- ggplot() +
  geom_sf(data = zoomed, color = "black", fill="white", linewidth= 0.4)+
  geom_sf(data = combined_sf, color = "grey48", fill="grey90")  +
  geom_sf(data = zoomed, color = "black", fill="transparent", linewidth= 0.4)+
  geom_sf(data = capture.sites, color = "black", shape = 17, fill="transparent", linewidth= 0.4)+
  theme_classic() +
  theme(axis.text = element_text(size = 11, family="sans"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  annotation_north_arrow(location = 'br',
                         style = north_arrow_nautical(text_family = "sans"),
                         which_north = "false") +
  annotation_scale(style = "bar",
                   pad_x = unit(0.03, "in"),
                   pad_y = unit(0.03, "in"),
                   text_family = "sans")
#' Plot map
map
