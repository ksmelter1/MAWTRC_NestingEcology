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
