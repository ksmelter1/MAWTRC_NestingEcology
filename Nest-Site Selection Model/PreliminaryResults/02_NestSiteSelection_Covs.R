
#'---
#' title: Nest-site selection of wild turkeys in Pennsylvania (an SSF analysis)
#' author: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output: *InsertDate*_Covs.RData
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script calculates the distance from each 30m x 30m raster cell to the nearest road and extracts elevation data from DEM
#' **Last Updated**: 1/24/25

#####################
## Load Packages 

packages <- c("sf",
              "amt",
              "tigris",
              "terra",
              "mapview",
              "tidyverse",
              "stars")

load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

lapply(packages, load_packages)

########################################
## Create Distance from Road Rasters 

#' Read in Pennsylvania state roads shapefile and project to Albers
#' Create a column for road ID that labels the road as primary
#' Select the Road ID column
pa.roads.prim <- st_read("Data Management/Shapefiles/roads/Pennsylvania/Primary Roads/PaStateRoads2023_10.shp") %>%
  st_transform(5070) %>%
  dplyr::mutate("Road_ID"="Primary") %>%
  dplyr::select(Road_ID) 

#' Read in NLCD
pa.nlcd <- terra::rast("Data Management/Rasters/nlcd/pa.nlcd.tif")

#' Create skeleton raster
#' To call the raster just type in the name in the console 
pa.nlcd
#' Format skeleton raster dimensions based on Maine NLCD
r <- terra::rast(nrow= 11063, ncol= 17330, xmin= 1261935, ymin= 1962645, nlyr=1,
                 xmax= 1781835, ymax= 2294535, crs= "epsg:5070")

#####################
## Primary Roads 

#' Rasterize road layer
pa.roads.prim.rast <- terra::rasterize(pa.roads.prim, r, fun=min)

#' Calculate distance to each 30m x 30m cell
dist.prim <- distance(pa.roads.prim.rast)

dist.prim <- crop(dist.prim, pa.nlcd)
dist.prim <- mask(dist.prim, pa.nlcd)

#' Check raster crs
raster_crs <- crs(dist.prim)

#' Save raster
 #writeRaster(dist.prim, "paroadrast.prim.tiff")

#' Load Raster
dist.prim<- terra::rast("Data Management/Rasters/PA Roads/paroadrast.prim.tiff")

########################
## Secondary Roads 

#' source: https://gis.stackexchange.com/questions/310489/calculating-euclidian-distance-in-r-between-lines-and-points
#' Read in secondary roads shapefile
pa.roads.sec <- st_read("Data Management/Shapefiles/roads/Pennsylvania/Secondary Roads/PaLocalRoads2023_10.shp") %>%
  sf::st_transform(5070) %>%
  dplyr::mutate("Road_ID"="Secondary") %>%
  dplyr::select(Road_ID, ID) 

#' Rasterize road layer
pa.roads.sec.rast <- terra::rasterize(pa.roads.sec, r, fun=min)

#' Calculate distance to each 30m x 30m cell
dist.sec <- distance(pa.roads.sec.rast)

dist.sec<- crop(dist.sec, pa.nlcd)
dist.sec <- mask(dist.sec, pa.nlcd)

#' Check raster crs
raster_crs <- crs(dist.sec)

#' Save raster
#writeRaster(dist.sec, "paroadrast.sec.tiff")

#' Load Raster
dist.sec<- terra::rast("Data Management/Rasters/PA Roads/paroadrast.sec.tiff")

###############################
## Extract Distances to Roads

#' Read in nest data with landcover
pa.nests.landcov.sf <- pa.nests.landcov

#' Extract point value at each nest
dist.prim.out <-terra::extract(dist.prim, pa.nests.landcov.sf) %>%
  dplyr::select(-ID)%>%
  dplyr::rename("primary" = layer)
dist.sec.out <-terra::extract(dist.sec, pa.nests.landcov.sf) %>%
  dplyr::select(-ID) %>%
  dplyr::rename("secondary" = layer)

pa.nests.landcov.sf.roads <- cbind(pa.nests.landcov.sf, dist.prim.out)
pa.nests.landcov.sf.roads <- cbind(pa.nests.landcov.sf.roads, dist.sec.out) 

pa.nests.covs <- cbind(pa.nests.landcov.sf, pa.nests.landcov.sf.roads) %>%
  dplyr::select(landuse, 
                Case,
                NestID, 
                BandID, 
                PercGrassForb,
                PercWoody, 
                PercLitter,
                PercBare,
                PercBoulder,
                AvgMaxVO, 
                primary,
                secondary,
                StemCount,
                WoodyType1,
                AvgVO,
                PercFern, 
                GuardHt,
                HtWoody,
                HtFern,
                HtGrassForb,
                HtLitter, 
                HtBoulder
                ) 


#' Convert land cover classifications to a categorical variable and create separate columns
pa.nests.covs$Agriculture <- ifelse(pa.nests.covs$landuse == "Agriculture", 1, 0)
pa.nests.covs$Developed <- ifelse(pa.nests.covs$landuse == "Developed", 1, 0)
pa.nests.covs$Deciduous <- ifelse(pa.nests.covs$landuse == "Deciduous Forest", 1, 0)
pa.nests.covs$Evergreen <- ifelse(pa.nests.covs$landuse == "Evergreen Forest", 1, 0)
pa.nests.covs$Mixed <- ifelse(pa.nests.covs$landuse == "Mixed Forest", 1, 0)
pa.nests.covs$Grassland <- ifelse(pa.nests.covs$landuse == "Grassland/Shrub", 1, 0)
pa.nests.covs$Water <- ifelse(pa.nests.covs$landuse == "Water", 1, 0)


########################################
## Save RData file *InsertDate_Covs

################################################################################
################################################################################

