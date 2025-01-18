
#'---
#' title: Nest-site selection of wild turkeys in Pennsylvania (an SSF analysis)
#' author: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'  
#' **Purpose**: This script calculates the distance from each 30m x 30m raster cell to the nearest road and extracts elevation data from DEM

#####################
## Load Packages ##
#####################

#' Vector of package names
packages <- c("sf",
              "amt",
              "tigris",
              "terra",
              "mapview",
              "tidyverse",
              "stars")

#' Function to load a package or install it if not already installed
load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

#' Apply the function to each package name
lapply(packages, load_packages)

########################################
## Create Distance from Road Rasters ##
########################################

#' Read in Pennsylvania state roads shapefile and project to Albers
#' Create a column for road ID that labels the road as primary
#' Select the Road ID column
pa.roads.prim <- st_read("Data Management/Shapefiles/roads/Pennsylvania/Primary Roads/PaStateRoads2023_10.shp") %>%
  st_transform(5070) %>%
  dplyr::mutate("Road_ID"="Primary") %>%
  dplyr::select(Road_ID) 

#' Read in NLCD
pa.nlcd <- terra::rast("Data Management/Rasters/nlcd/paNLCD.tiff")

#' Create skeleton raster
#' To call the raster just type in the name in the console 
pa.nlcd
#' Format skeleton raster dimensions based on Maine NLCD
r <- terra::rast(nrow= 11063, ncol= 17330, xmin= 1261935, ymin= 1962645, nlyr=1,
                 xmax= 1781835, ymax= 2294535, crs= "epsg:5070")

#####################
## Primary Roads ##
####################

#' Rasterize road layer
pa.roads.prim.rast <- terra::rasterize(pa.roads.prim, r, fun=min)

#' Calculate distance to each 30m x 30m cell
dist.prim <- distance(pa.roads.prim.rast)

dist.prim <- crop(dist.prim, pa.nlcd)
dist.prim <- mask(dist.prim, pa.nlcd)

#' Check raster crs
raster_crs <- crs(dist.prim)

#' Save raster
 writeRaster(dist.prim, "paroadrast.prim.tiff")

#' Load Raster
dist.prim<- terra::rast("Data Management/Rasters/PA Roads/paroadrast.prim.tiff")

########################
## Secondary Roads ##
#######################

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
plot(dist.sec)

###############################
## Extract Distances to Roads

#' Read in nest data with landcover
pa.nests.landcov.sf <- st_read("Data Management/Shapefiles/Nest-Site Selection/pa.nests.landcov.sf.shp")
class(pa.nests.landcov.sf)

#' Extract point value at each nest
dist.prim.out <-terra::extract(dist.prim, pa.nests.landcov.sf) %>%
  dplyr::select(-ID)%>%
  dplyr::rename("primary" = layer)
dist.sec.out <-terra::extract(dist.sec, pa.nests.landcov.sf) %>%
  dplyr::select(-ID) %>%
  dplyr::rename("secondary" = layer)

pa.nests.landcov.sf.roads <- cbind(pa.nests.landcov.sf, dist.prim.out)
pa.nests.landcov.sf.roads <- cbind(pa.nests.landcov.sf.roads, dist.sec.out) 
  
#################
## Elevation ##
################

#' Read in elevation raster
DEM.proj <- terra::rast("Data Management/Rasters/DEM/paDEM.tiff") 
DEM.crs <- terra::project(DEM.proj, pa.nlcd, method= "bilinear")

#' Extract elevation data
pa.nests.elev <- terra::extract( DEM.crs, pa.nests.landcov.sf) %>%
  dplyr::rename("elev" = Layer_1) %>%
  dplyr::select(-ID)

pa.nests.covs <- cbind(pa.nests.elev, pa.nests.landcov.sf.roads) %>%
  dplyr::rename("landuse" = Class) %>%
  dplyr::rename("nest_id" = nestd_v) %>%
  dplyr::select(landuse, 
                case,
                elev, 
                nest_id, 
                bandid, 
                prcgrss,
                percwdy, 
                avrgmxv, 
                uniqueid,
                primary,
                secondary)

#' Convert land cover classifications to a categorical variable and create separate columns
pa.nests.covs$Agriculture <- ifelse(pa.nests.covs$landuse == "Agriculture", 1, 0)
pa.nests.covs$Developed <- ifelse(pa.nests.covs$landuse == "Developed", 1, 0)
pa.nests.covs$Deciduous <- ifelse(pa.nests.covs$landuse == "Deciduous Forest", 1, 0)
pa.nests.covs$Evergreen <- ifelse(pa.nests.covs$landuse == "Evergreen Forest", 1, 0)
pa.nests.covs$Mixed <- ifelse(pa.nests.covs$landuse == "Mixed Forest", 1, 0)
pa.nests.covs$Grassland<- ifelse(pa.nests.covs$landuse == "Grassland/Shrub", 1, 0)

#' Save RData file with everything we need in top right

################################################################################
################################################################################

############################################
## Distance to Primary Road Calculations ##
############################################

# pa.nests.landcov <- pa.nests.landcov %>%
#   dplyr::select(-geometry) %>%
#   st_as_sf(., coords = c("long_v", "lat_v"), crs = 4326) %>%
#   st_transform(5070)
# 
# mapview(pa.nests.landcov.sf)


#' #### Franny's code ####
#' 
#' #' Create random steps vector 
#' #' The length of the list is the number of unique IDs in the random_steps object
#' #' Should be 180
#' pa.nests.landcov.list<-vector(mode = "list", length = length(unique(pa.nests.landcov$uniqueid)))
#' for (i in 1:length(unique(pa.nests.landcov$uniqueid))){
#'   pa.nests.landcov.list[[i]]<-pa.nests.landcov[which(pa.nests.landcov$uniqueid==unique(pa.nests.landcov$uniqueid)[i]),]
#' }
#' 
#' #' Read in Pennsylvania state roads shapefile and project to Albers
#' #' Create a column for road ID that labels the road as primary
#' #' Select the Road ID column
#' pa.roads.prim <- st_read("Data Management/Shapefiles/roads/Pennsylvania/Primary Roads/PaStateRoads2023_10.shp") %>%
#'   st_transform(5070) %>%
#'   dplyr::mutate("Road_ID"="Primary") %>%
#'   dplyr::select(Road_ID) 
#' 
#' 
#' #' Create an empty vector to store minimum distances from the nearest road
#' min.road.dist.prim<-vector(mode = "list", length = length(pa.nests.landcov.list))
#' 
#' #' Create processing time object
#' ptm <- proc.time()
#' 
#' #' For loop to calculate distance to the nearest primary road 
#' for (i in 1: length( pa.nests.landcov.list)) {
#'   
#'   #' Create step.coords object which contains all x and y coordinates from the random steps object
#'   nest.coords<-data.frame(x= pa.nests.landcov.list[[i]]$lat_v,y= pa.nests.landcov.list[[i]]$long_v)
#'   
#'   #' Buffer the mean.coord by 10,000 meters
#'   mean_buffer.prim = st_buffer(pa.nests.landcov, 20000)
#'   
#'   #' Stratify roads object to the roads that are within the buffer
#'   roads.in.buffer.prim<-st_intersection(pa.roads.prim,  pa.nests.landcov)
#'   
#'   #' Rasterize the roads in the buffer
#'   roads.raster.prim<-as((stars::st_rasterize(roads.in.buffer.prim)),"Raster")
#'   
#'   #' Remove all NA's from the road raster
#'   road.coords.prim<-data.frame(xyFromCell(roads.raster.prim,which(!is.na(values(roads.raster.prim)))))
#'   
#'   #' Nested loop
#'   #' For each row in step.coords
#'   for(j in 1:nrow(nest.coords)){
#'     
#'     #' Extract the minimum euclidean distance in step.coords, road.coords.prim and store in the min.road.dist.prim
#'     min.road.dist.prim[[i]][j]<-min(dist(nest.coords[j,],road.coords.prim,method="euclidean"))
#'   }
#'   
#'   print(i)
#' }
#' proc.time() - ptm 
#' 
#' #' Create distance to primary road object
#' dist.to.prim.road.out <- do.call(c,min.road.dist.prim)
#' dist.to.prim.road.out.df <- as.data.frame(dist.to.prim.road.out)
#' 
#' #' ###############################################
#' #' ## Distance to Secondary Road Calculations ##
#' #' ##############################################
#' 
#' #' source: https://gis.stackexchange.com/questions/310489/calculating-euclidian-distance-in-r-between-lines-and-points
#' #' Read in secondary roads shapefile
#' pa.roads.sec <- st_read("Data Management/Shapefiles/roads/Pennsylvania/Secondary Roads/PaLocalRoads2023_10.shp") %>%
#'   sf::st_transform(5070) %>%
#'   dplyr::mutate("Road_ID"="Secondary") %>%
#'   dplyr::select(Road_ID, ID) 
#' 
#' random.steps.list<-vector(mode = "list", length = length(unique(random_steps$id)))
#' for (i in 1:length(unique(random_steps$id))){
#'   random.steps.list[[i]]<-random_steps[which(random_steps$id==unique(random_steps$id)[i]),]
#' }
#' 
#' #' Create an empty vector to store minimum distances from the nearest road
#' min.road.dist.sec<-vector(mode = "list", length = length( random.steps.list))
#' 
#' #' Create processing time object
#' ptm <- proc.time()
#' 
#' #' For loop to calculate distance to the nearest primary road 
#' for (i in 1: length(random.steps.list)) {
#'   
#'   #' Create step.coords object which contains all x and y coordinates from the random steps object
#'   step.coords<-data.frame(x=random.steps.list[[i]]$x2_,y=random.steps.list[[i]]$y2_)
#'   
#'   #' Create mean.coord.prim object which contains the mean x and y coordinates for each bird
#'   mean.coord.sec <-data.frame(x=mean(step.coords$x,na.rm=T),y=mean(step.coords$y,na.rm=T))
#'   
#'   #' Create mean coord as a spatial object
#'   mean.coord.sf.sec<-st_as_sf(mean.coord.sec, coords = c("x","y"), crs=5070)
#'   
#'   #' Buffer the mean.coord by 10,000 meters
#'   mean_buffer.sec = st_buffer(mean.coord.sf.sec, 20000)
#'   
#'   
#'   #' Stratify roads object to the roads that are within the buffer
#'   roads.in.buffer.sec<-st_intersection(pa.roads.sec,  mean_buffer.sec)
#'   
#'   #' Rasterize the roads in the buffer
#'   roads.raster.sec<-as((stars::st_rasterize(roads.in.buffer.sec)),"Raster")
#'   
#'   #' Remove all NA's from the road raster
#'   road.coords.sec<-data.frame(xyFromCell(roads.raster.sec,which(!is.na(values(roads.raster.sec)))))
#'   
#'   #' Nested loop
#'   #' For each row in step.coords
#'   for(j in 1:nrow(step.coords)){
#'     
#'     #' Extract the minimum euclidean distance in step.coords, road.coords.prim and store in the min.road.dist.prim
#'     min.road.dist.sec[[i]][j]<-min(dist(step.coords[j,],road.coords.sec,method="euclidean"))
#'   }
#'   
#'   print(i)
#' }
#' proc.time() - ptm 
#' 
#' #' Create distance to primary road object
#' dist.to.sec.road.out <- do.call(c,min.road.dist.sec)
#' dist.to.sec.road.out.df <- as.data.frame(dist.to.sec.road.out)
#' 
