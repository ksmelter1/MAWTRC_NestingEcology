
#'---
#' title: Habitat selection of wild turkeys during incubation (an SSF analysis)
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
#' **Purpose**: This script creates a 30m x 30m cell distance from road raster and we extract distance from road and elevation from the NED

#####################
## Load Packages ##
#####################

#' Vector of package names
packages <- c("sf",
              "amt",
              "tigris",
              "dplyr",
              "terra",
              "mapview",
              "tidyr")

#' Function to load a package or install it if not already installed
load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

#' Apply the function to each package name
lapply(packages, load_packages)


#' Visually checked and the projections look good
load("Data Management/RData/Individual-Specific Movement Process/RData Files/RandomSteps_Prep.RData") 

#' Convert the random steps object to an sf object
#' 5070 is Albers equal area conic projection
random_steps.sf <- random_steps %>%
  st_as_sf(coords = c("x2_", "y2_") , crs=5070)

############################################
## Distance to Primary Road Calculations ##
############################################

#source: https://gis.stackexchange.com/questions/310489/calculating-euclidian-distance-in-r-between-lines-and-points

#' Read in primary roads shapefile
pa.roads.prim <- st_read("Data Management/Shapefiles/roads/Pennsylvania/Primary Roads/PaStateRoads2023_10.shp") %>%
  sf::st_transform(5070) %>%
  dplyr::mutate("Road_ID"="Primary") %>%
  dplyr::select(Road_ID, NLF_ID) 
  
#Empty lists to store final outputs
random_steps.prim.out <- list()
prim.road.out <- list()

#Bounding box of nests in PA as a spatial object
#st_as_sfc is class 3 method for bbox
random_steps.bbox <- sf::st_bbox(random_steps.sf) %>% st_as_sfc()

#Bounding box of roads in PA as a spatial object
#st_as_sfc is class 3 method for bbox
pa.roads.prim.bbox <- sf::st_bbox(pa.roads.prim) %>% st_as_sfc() 

#For loop to extract the distance from each nest to the nearest road
#Takes me about 40 minutes to run so I've included RData results below
for (i in 1: nrow(random_steps.sf)) {
  
  #join the nests and roads
  #st_nearest feature: get index of nearest feature
  #Spatial join of pa.nests and pa.roads
  #[[]] first entire object within list, indexing of lists
  neardist <- sf::st_join(
    random_steps.sf[i,],
    pa.roads.prim, 
    join= sf::st_nearest_feature
  ) 
  
  #find the distance of each nest to the nearest road
  #we know this is in meters by units(dist.to.road)
  dist.to.prim.road <- sf::st_distance(random_steps.sf[i,],
                                  pa.roads.prim %>%
                                    dplyr::filter(NLF_ID==neardist$NLF_ID))
  
  #create a df with things we need
  dist.to.road.prim.df <- data.frame(id = (random_steps.sf$id)[i],
                                case_ = (random_steps.sf$case_)[i],
                                step_id_= (random_steps.sf$step_id_)[i],
                                dist.to.prim.road = as.numeric(dist.to.prim.road))
  
  #put the df in output
  prim.road.out[[i]] <- dist.to.road.prim.df 
  
  #Printing i here means that for every time the loop runs the row is recorded in the console
  print(i)
  
}

#Create dist2road object
dist.to.prim.road.out <- do.call(rbind,(prim.road.out))

#Write csv
#write.csv(dist2road, "Csvs/dist2road.out.csv")

###############################################
## Distance to Secondary Road Calculations ##
##############################################

#source: https://gis.stackexchange.com/questions/310489/calculating-euclidian-distance-in-r-between-lines-and-points

#' Read in secondary roads shapefile
pa.roads.sec <- st_read("Data Management/Shapefiles/roads/Pennsylvania/Secondary Roads/PaLocalRoads2023_10.shp") %>%
  sf::st_transform(5070) %>%
  dplyr::mutate("Road_ID"="Secondary") %>%
  dplyr::select(Road_ID, ID) 

#Empty lists to store final outputs
random_steps.sec.out <- list()
sec.road.out <- list()

#Bounding box of nests in PA as a spatial object
#st_as_sfc is class 3 method for bbox
random_steps.bbox <- sf::st_bbox(random_steps.sf) %>% st_as_sfc()

#Bounding box of roads in PA as a spatial object
#st_as_sfc is class 3 method for bbox
pa.roads.sec.bbox <- sf::st_bbox(pa.roads.sec) %>% st_as_sfc() 

#For loop to extract the distance from each nest to the nearest road
#Takes me about 40 minutes to run so I've included RData results below
for (i in 1: nrow(random_steps.sf)) {
  
  #join the nests and roads
  #st_nearest feature: get index of nearest feature
  #Spatial join of pa.nests and pa.roads
  #[[]] first entire object within list, indexing of lists
  neardist <- sf::st_join(
    random_steps.sf[i,],
    pa.roads.sec, 
    join= sf::st_nearest_feature
  ) 
  
  #find the distance of each nest to the nearest road
  #we know this is in meters by units(dist.to.road)
  dist.to.sec.road <- sf::st_distance(random_steps.sf[i,],
                                       pa.roads.sec %>%
                                         dplyr::filter(ID==neardist$ID))
  
  #create a df with things we need
  dist.to.road.sec.df <- data.frame(id = (random_steps.sf$id)[i],
                                     case_ = (random_steps.sf$case_)[i],
                                     step_id_= (random_steps.sf$step_id_)[i],
                                     dist.to.sec.road = as.numeric(dist.to.sec.road))
  
  #put the df in output
  sec.road.out[[i]] <- dist.to.road.sec.df 
  
  #Printing i here means that for every time the loop runs the row is recorded in the console
  print(i)
  
}

#Create dist2road object
dist.to.sec.road.out <- do.call(rbind,(sec.road.out))


#' #################################
#' ## Distance from Road Raster ##
#' #################################
#' 
#' #' Not quite as accurate as exact distances calculated in terra
#' #' Great for a predictive surface if results are significant
#' #' Code runs 
#' 
#' #' Read in roads shapefile
#' pa.roads <- st_read("Data Management/Shapefiles/roads/pa.roads.shp") %>%
#'   sf::st_transform(5070) %>%
#'   dplyr::select(ID) %>%
#'   tidyr::drop_na(ID)
#' 
#' #' NLCD without Lake Erie, thanks Jake 
#' pa<- sf::st_read("Data Management/Shapefiles/PaMunicipalities2024_03/PaMunicipalities2024_03.shp") %>%
#'   st_transform(5070)
#' 
#' #' Create skeleton raster
#' #' To call the raster just type in the name in the console 
#' r <- terra::rast(nrow= 11063, ncol= 17117, xmin= 1268325 , ymin= 1962645 , nlyr=1,
#'           xmax= 1781835 , ymax= 2294535, crs= "epsg:5070")
#' 
#' #' Rasterize road layer
#' pa.roads.rast <- terra::rasterize(pa.roads, r, fun=min)
#' 
#' #' Calculate distance to each 30m x 30m cell
#' dist <- distance(pa.roads.rast)
#' 
#' dist <- crop(dist, pa)
#' dist <- mask(dist, pa)
#' 
#' #' Check raster crs
#' raster_crs <- crs(dist)
#' 
#' #' Save raster
#' #writeRaster(dist, "paroadrast.tiff")
#' 
#' ##########################
#' ## Use Created Raster ##
#' ##########################
#' 
#' #' Read in raster
#' dist <- terra::rast("Data Management/Rasters/PA Roads/paroadrast.tiff")
#' 
#' #' Check
#' glimpse(random_steps)
#' 
#' #' Transform to Albers
#' dat_2 <- random_steps %>% st_as_sf(coords = c("x2_", "y2_"), crs = 5070) 
#' 
#' #' Check projection using a subset of data
#' dat_2.subset <- dat_2 %>% dplyr::filter(id=="8579_2023_1_1")
#' mapview(dat_2.subset)
#' 
#' #' Transform the random_steps to the raster crs
#' dat_2 <- st_transform(dat_2, crs = raster_crs)
#' mapview(dat_2.subset)
#' 
#' #' Extract distance from road outputs from raster
#' #' Min is the minimum distance to the nearest road
#' dist.outputs <- terra::extract(dist, dat_2, fun=min)
#' 
#' #' Bind columns back together 
#' dat_2 <- cbind(dist.outputs, dat_2) %>%
#'   dplyr::rename("distfromroad"= layer) %>%
#'   dplyr::select(-ID)

###################################
## Vectorize Distance from Road ##
###################################

#' Jake's suggestion
#' Code runs 

# test <- data.frame()
# for (i in 1:nrow(dat_2)){
#   
#   tmpbuffer <- buffer(vect(dat_2[i,]),width=2000)
#   tmpcrop <-crop(vect(pa.roads), tmpbuffer)
# 
# #Use terra nearest function
# #Vectorize data for faster processing
# distances <- terra::nearest(vect(dat_2[i, "geometry"]), 
#                             tmpcrop)
# distances= st_as_sf(distances, crs=5070)
# test= rbind(test, distances$distance)
# }
# 
# 
#' Check, thanks Jake
# mapview(dat_2[1:10,]) + mapview(st_as_sf(buffer(vect(dat_2[1:10,]), width=test$X208.842885201443)))
