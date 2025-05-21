##'---
#' title: Study Area Map for Nesting Ecology Manuscript
#' authors: "K. Smelter, F. Buderman
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'
#+ include = FALSE
#'  
#' **Purpose**: This script creates a study area map for the nesting ecology manuscript
#' **Last Updated**: 5/5/2025

################################################################################
## Load Packages

library("sf")
library("ggplot2")
library("terra")
library("mapview")
library("dplyr")
library("ggspatial")

################################################################################
## US Map

# Us map shapefile without great lakes
us <- st_read("Data Management/Shapefiles/USA_no_greatlakes/USA_adm1_without_great_lakes.shp")
us <- st_transform(us, crs = 5070)

################################################################################
## Pennsylvania

# Pennsylvania WMU shapefile
# Project to Albers
pennsylvania <- st_read("Data Management/Shapefiles/Pennsylvania/Pennsylvania WMUs/PGC_BNDWildlifeManagementUnits2021.shp")
pennsylvania <- st_transform(pennsylvania, crs = 5070) 

# The hen study takes place here in PA in 2D, 3D, 4D, and 5C in PA
pa.study.area <- subset(pennsylvania, WMU_ID=="2D"| WMU_ID=="3D"| WMU_ID=="4D"| WMU_ID=="5C")


# Create WMU_ID column for labelling
pa.study.area$WMU_ID_Map <- ifelse(pa.study.area$WMU_ID == "2D", "PA 2D",
                               ifelse(pa.study.area$WMU_ID == "3D", "PA 3D",
                                      ifelse(pa.study.area$WMU_ID == "4D", "PA 4D",
                                             ifelse(pa.study.area$WMU_ID == "5C", "PA 5C",NA))))


################################################################################
## New Jersey

# New Jersey Turkey hunting area shapefile
# Transform to Albers and reclassify THAs
newjersey <- st_read("Data Management/Shapefiles/New Jersey/New Jersey Turkey Hunting Areas/NJ_THA.shp") %>%
  st_zm() 
newjersey <- st_transform(newjersey, crs = 5070)
newjersey <- st_intersection(newjersey, us[us$NAME_1=="New Jersey",])
newjersey$THA <- ifelse(newjersey$THA<11,"North",
                      ifelse(newjersey$THA>10 & newjersey$THA<15,"Central","South")) 

# Use summarise to display new THAs on map
# Check work with ggplot
nj.thas.ready <-newjersey %>% 
  dplyr::group_by(THA) %>%
  summarise()
ggplot() + geom_sf(data = pennsylvania) + geom_sf(data = nj.thas.ready) 

# The New Jersey study area consists of NJ North and NJ South
nj.study.area <- subset(nj.thas.ready,THA == "South")

# Create WMU_ID column for labelling
nj.study.area$WMU_ID_Map <- ifelse(nj.study.area$THA == "South", "NJ South", NA)
  
################################################################################
## Maryland

# Read in Maryland geopackage with Turkey hunting zones
# Check work with MD, PA, and NJ
maryland <- st_read("Data Management/Shapefiles/Maryland/Maryland Turkey Regions/Sa.map/md_wmu_boundaries.gpkg") 
ggplot() + geom_sf(data = pennsylvania) + 
  geom_sf(data = nj.thas.ready) + geom_sf(data = maryland)

# The Maryland study takes place in two area MD West and MD East
md.study.area <- subset(maryland, wmu == "Western" | wmu == "Lower Eastern Shore") 

# Create WMU_ID column for labelling
md.study.area$WMU_ID_Map <- ifelse(md.study.area$wmu == "Western", "MD West",
                               ifelse(md.study.area$wmu == "Lower Eastern Shore", "MD East", NA))

  
################################################################################
## Create Map

# Display PA, MD, and NJ
ggplot() + 
  geom_sf(data = pennsylvania) + 
  geom_sf(data = nj.thas.ready) + 
  geom_sf(data = maryland) 

# Create object with all three states and create bounding box
all <- st_union(pennsylvania, nj.thas.ready, maryland)
st_bbox(all)

# Create bounding box
box_zoom = (st_bbox(c(xmin= 1113751, ymin= 1776600, xmax= 1871955, ymax= 2508599))) %>%
  st_as_sfc() %>% st_set_crs(5070)
us_crop <- st_crop(us, box_zoom)

# Create state-only boundaries by dissolving internal WMUs/THAs
pa_boundary <- st_union(pennsylvania)
nj_boundary <- st_union(nj.thas.ready)
md_boundary <- st_union(maryland)


# Create object with all other US states
us_crop2 <- us_crop[!us_crop$NAME_1 %in% c("Maryland",
                                           "Pennsylvania","New Jersey"),]
sa_map <- ggplot() +
  
  geom_sf(data = pa_boundary, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = nj_boundary, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = md_boundary, fill = NA, color = "black", linewidth = 0.6) +
  
  # State boundaries for PA, NJ, and MD
  geom_sf(data = pennsylvania, fill = "grey90", color = "black", linewidth = 0.0001) +
  geom_sf(data = nj.thas.ready, fill = "grey90", color = "black", linewidth = 0.001) +
  geom_sf(data = maryland, fill = "grey90", color = "black", linewidth = 0.001) +
  
  # Study area boundaries (thinner)
  geom_sf(data = pa.study.area, fill = "darkgray", color = "black", linewidth = 0.0001) +
  geom_sf(data = md.study.area, fill = "darkgray", color = "black", linewidth = 0.0001) +
  geom_sf(data = nj.study.area, fill = "darkgray", color = "black", linewidth = 0.0001) +
  
  # Labels
  geom_sf_label(data = pa.study.area, aes(label = WMU_ID_Map), size = 4, fontface = "bold", label.size = 0.3) +
  geom_sf_label(data = md.study.area, aes(label = WMU_ID_Map), size = 4, fontface = "bold", label.size = 0.3) +
  geom_sf_label(data = nj.study.area, aes(label = WMU_ID_Map), size = 4, fontface = "bold", label.size = 0.3) +
  
  # Map styling
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  annotation_north_arrow(location = 'br',
                         style = north_arrow_orienteering,
                         height = unit(1.25, "cm"), 
                         width = unit(1.25, "cm")) +
  annotation_scale(style = "bar",
                   pad_x = unit(0.05, "in"),
                   pad_y = unit(0.05, "in"))

plot(sa_map)



