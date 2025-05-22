#'---
#' title: Nest-site selection of wild turkeys in Pennsylvania (an SSF analysis)
#' author: "K. Smelter
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#'---
#'  
#' **Purpose**: This script obtains NLCD values for used and available nests in each state
#' **Last Updated**: 5/12/2025


################################################################################
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

#' Load in landcover
load("Data Management/RData/Nest-Site Selection/Covs/Multi-State/Manuscript/01_Covs_AllStates_Landcover.RData")

################################################################################
## Create Distance from Road Rasters 

#' Download local and state road data from all counties in each of the three states
# pa.roads <- tigris::roads(state = "Pennsylvania", 
#                           county = c("Adams", "Allegheny", "Armstrong", "Beaver", "Bedford", 
#                                      "Berks", "Blair", "Bradford", "Bucks", "Butler", "Cambria", 
#                                      "Cameron", "Carbon", "Centre", "Chester", "Clarion", 
#                                      "Clearfield", "Clinton", "Columbia", "Crawford", "Cumberland", 
#                                      "Dauphin", "Delaware", "Elk", "Erie", "Fayette", "Forest", 
#                                      "Franklin", "Fulton", "Greene", "Huntingdon", "Indiana", 
#                                      "Jefferson", "Juniata", "Lackawanna", "Lancaster", 
#                                      "Lawrence", "Lebanon", "Lehigh", "Luzerne", "Monroe", 
#                                      "Montgomery", "Montour", "Northampton", "Northumberland", 
#                                      "Perry", "Philadelphia", "Pike", "Potter", "Schuylkill", 
#                                      "Snyder", "Somerset", "Sullivan", "Susquehanna", "Tioga", 
#                                      "Union", "Venango", "Warren", "Washington", "Wayne", 
#                                      "Westmoreland", "Wyoming", "York"))

nj.roads <- tigris::roads(state = "New Jersey", 
                          county = c("Atlantic", "Bergen", "Burlington", "Camden", "Cape May", 
                                     "Cumberland", "Essex", "Gloucester", "Hudson", "Hunterdon", 
                                     "Mercer", "Middlesex", "Monmouth", "Morris", "Ocean", 
                                     "Passaic", "Salem", "Somerset", "Sussex", "Union", "Warren"))

md.roads <- tigris::roads(state = "Maryland", 
                          county = c("Allegany", "Anne Arundel", "Baltimore County", "Baltimore City",
                                     "Calvert", "Caroline", "Carroll", 
                                     "Cecil", "Charles", "Dorchester", "Frederick", 
                                     "Garrett", "Harford", "Howard", "Kent", "Montgomery", 
                                     "Prince George's", "Queen Anne's", "St. Mary's", 
                                     "Somerset", "Talbot", "Washington", "Wicomico", 
                                     "Worcester"))

# table(pa.roads$RTTYP)
table(nj.roads$RTTYP)
table(md.roads$RTTYP)


#' Define Albers Equal Area CRS (EPSG:5070)
albers_crs <- 5070

#' Transform road data to Albers CRS for PA, NJ, and MD roads
# pa.roads <- st_transform(pa.roads, crs = albers_crs)
nj.roads <- st_transform(nj.roads, crs = albers_crs)
md.roads <- st_transform(md.roads, crs = albers_crs)


#' Remove NA values for road type
#' If else statement to group by road type
# pa.roads <- pa.roads %>% drop_na(RTTYP)
# pa.roads$Type <- ifelse(pa.roads$RTTYP == "M" | pa.roads$RTTYP == "C", "Secondary","Primary")
# pa.roads.prim <- pa.roads %>%
#   dplyr::filter(Type == "Primary") 
# pa.roads.sec <- pa.roads %>%
#   dplyr::filter(Type == "Secondary")

#' Remove NA values for road type
#' If else statement to group by road type
nj.roads <- nj.roads %>% drop_na(RTTYP)
nj.roads$Type <- ifelse(nj.roads$RTTYP == "M" | nj.roads$RTTYP == "C","Secondary","Primary")
nj.roads.prim <- nj.roads %>%
  dplyr::filter(Type == "Primary")
nj.roads.sec <- nj.roads %>%
  dplyr::filter(Type == "Secondary")

#' Remove NA values for road type
#' If else statement to group by road type
md.roads <- md.roads %>% drop_na(RTTYP)
md.roads$Type <- ifelse(md.roads$RTTYP == "M" | md.roads$RTTYP =="C", "Secondary","Primary")
md.roads.prim <- md.roads %>%
  dplyr::filter(Type == "Primary")
md.roads.sec <- md.roads %>%
  dplyr::filter(Type == "Secondary")


################################################################################
## Load in NLCD Rasters for use as Templates

#' Read in NLCD
pa.nlcd <- terra::rast("Data Management/Rasters/NLCD/Pennsylvania/pa.nlcd.tif")

#' Create skeleton raster
pa.nlcd

#' Format skeleton raster dimensions based on NLCD
pa.r <- terra::rast(nrow= 11063, ncol= 17330, xmin= 1261935, ymin= 1962645, nlyr=1,
                 xmax= 1781835, ymax= 2294535, crs= "epsg:5070")

#' Read in NLCD
nj.nlcd <- terra::rast("Data Management/Rasters/NLCD/New Jersey/nj.nlcd.tif")

#' Create skeleton raster
nj.nlcd

#' Format skeleton raster dimensions based on NLCD
nj.r <- terra::rast(nrow= 9556, ncol= 4101, xmin= 1725345  , ymin= 1949775, nlyr=1,
                     xmax= 1848375, ymax= 2236455, crs= "epsg:5070")
#' Read in NLCD
md.nlcd <- terra::rast("Data Management/Rasters/NLCD/Maryland/md.nlcd.tif")

#' Create skeleton raster
md.nlcd

#' Format skeleton raster dimensions based on NLCD
md.r <- terra::rast(nrow= 6975, ncol= 13535, xmin= 1396635, ymin= 1828485, nlyr=1,
                  xmax= 1802685, ymax= 2037735, crs= "epsg:5070")
  
  
################################################################################
## Pennsylvania

#' This code uses the tigris package
#' I am going to move forward with PASDA

#' Rasterize road layer
#' pa.roads.rast <- terra::rasterize(pa.roads, pa.r, fun = min)
#' pa.dist <- distance(pa.roads.rast)
#' pa.dist <- crop(pa.dist, pa.nlcd)
#' pa.dist <- mask(pa.dist, pa.nlcd)
#' raster_crs <- crs(pa.dist)
#' plot(pa.dist)
#' writeRaster(pa.dist, "Data Management/Rasters/Roads/Pennsylvania/PaRoadRast.tiff",
#'             overwrite = T)
#' 
#' #' Rasterize road layer
#' pa.roads.prim.rast <- terra::rasterize(pa.roads.prim, pa.r, fun=min)
#' pa.dist.prim <- distance(pa.roads.prim.rast) #' Calculate distance to each 30m x 30m cell
#' pa.dist.prim <- crop(pa.dist.prim, pa.nlcd)
#' pa.dist.prim <- mask(pa.dist.prim, pa.nlcd)
#' raster_crs <- crs(pa.dist.prim)
#' plot(pa.dist.prim)
#' writeRaster(pa.dist.prim, "Data Management/Rasters/Roads/Pennsylvania/PaRoadRast.Prim.tiff",
#'             overwrite = T)
#' pa.dist.prim<- terra::rast("Data Management/Rasters/Roads/Pennsylvania/PaRoadRast.Prim.tiff")
#' 
#' pasda.dist.prim <- terra::rast("Data Management/Rasters/PA Roads/paroadrast.prim.tiff")
#' 
#' #' Rasterize road layer
#' pa.roads.sec.rast <- terra::rasterize(pa.roads.sec, pa.r, fun=min)
#' pa.dist.sec <- distance(pa.roads.sec.rast) #' Calculate distance to each 30m x 30m cell
#' pa.dist.sec <- crop(pa.dist.sec, pa.nlcd)
#' pa.dist.sec <- mask(pa.dist.sec, pa.nlcd)
#' raster_crs <- crs(pa.dist.sec)
#' writeRaster(pa.dist.sec, "Data Management/Rasters/Roads/Pennsylvania/PaRoadRast.Sec.tiff", 
#'             overwrite = T)
#' pa.dist.sec<- terra::rast("Data Management/Rasters/Roads/Pennsylvania/PaRoadRast.Sec.tiff")
#' 
#' pasda.dist.sec <- terra::rast("Data Management/Rasters/PA Roads/paroadrast.sec.tiff")


################################################################################
## New Jersey

#' Rasterize road layer
nj.roads.rast <- terra::rasterize(nj.roads, nj.r, fun = min)
nj.dist <- distance(nj.roads.rast) #' Calculate distance to each 30m x 30m cell
nj.dist <- crop(nj.dist, nj.nlcd)
nj.dist <- mask(nj.dist, nj.nlcd)
raster_crs <- crs(nj.dist)
plot(nj.dist)
writeRaster(nj.dist, "Data Management/Rasters/Roads/New Jersey/NjRoadRast.tiff",
            overwrite = T)


#' Rasterize road layer
nj.roads.prim.rast <- terra::rasterize(nj.roads.prim, nj.r, fun=min)
nj.dist.prim <- distance(nj.roads.prim.rast) #' Calculate distance to each 30m x 30m cell
nj.dist.prim <- crop(nj.dist.prim, nj.nlcd)
nj.dist.prim <- mask(nj.dist.prim, nj.nlcd)
raster_crs <- crs(nj.dist.prim)
plot(nj.dist.prim)
writeRaster(nj.dist.prim, "Data Management/Rasters/Roads/New Jersey/NjRoadRast.Prim.tiff",
            overwrite = T)
nj.dist.prim<- terra::rast("Data Management/Rasters/Roads/New Jersey/NjRoadRast.Prim.tiff")

#' Rasterize road layer
nj.roads.sec.rast <- terra::rasterize(nj.roads.sec, nj.r, fun=min)
nj.dist.sec <- distance(nj.roads.sec.rast) #' Calculate distance to each 30m x 30m cell
nj.dist.sec <- crop(nj.dist.sec, nj.nlcd)
nj.dist.sec <- mask(nj.dist.sec, nj.nlcd)
raster_crs <- crs(nj.dist.sec)
plot(nj.dist.sec)
writeRaster(nj.dist.sec, "Data Management/Rasters/Roads/New Jersey/NjRoadRast.Sec.tiff",
            overwrite = T)
nj.dist.sec<- terra::rast("Data Management/Rasters/Roads/New Jersey/NjRoadRast.Sec.tiff")


################################################################################
## Maryland

#' Rasterize road layer
md.roads.rast <- terra::rasterize(md.roads, md.r, fun = min)
md.dist <- distance(md.roads.rast) #' Calculate distance to each 30m x 30m cell
md.dist <- crop(md.dist, md.nlcd)
md.dist <- mask(md.dist, md.nlcd)
raster_crs <- crs(md.dist)
plot(md.dist)
writeRaster(md.dist, "Data Management/Rasters/Roads/Maryland/MdRoadRast.tiff",
            overwrite = T)
md.dist <- terra::rast("Data Management/Rasters/Roads/Maryland/MdRoadRast.tiff")

#' Rasterize road layer
md.roads.prim.rast <- terra::rasterize(md.roads.prim, md.r, fun=min)
md.dist.prim <- distance(md.roads.prim.rast) #' Calculate distance to each 30m x 30m cell
md.dist.prim <- crop(md.dist.prim, md.nlcd)
md.dist.prim <- mask(md.dist.prim, md.nlcd)
raster_crs <- crs(md.dist.prim)
plot(md.dist.prim)
writeRaster(md.dist.prim, "Data Management/Rasters/Roads/Maryland/MdRoadRast.Prim.tiff",
            overwrite = T)
md.dist.prim<- terra::rast("Data Management/Rasters/Roads/Maryland/MdRoadRast.Prim.tiff")

#' Rasterize road layer
md.roads.sec.rast <- terra::rasterize(md.roads.sec, md.r, fun=min)
md.dist.sec <- distance(md.roads.sec.rast) #' Calculate distance to each 30m x 30m cell
md.dist.sec <- crop(md.dist.sec, md.nlcd)
md.dist.sec <- mask(md.dist.sec, md.nlcd)
raster_crs <- crs(md.dist.sec)
plot(md.dist.sec)
writeRaster(md.dist.sec, "Data Management/Rasters/Roads/Maryland/MdRoadRast.Sec.tiff",
            overwrite = T)
md.dist.sec<- terra::rast("Data Management/Rasters/Roads/Maryland/MdRoadRast.Sec.tiff")


################################################################################
## Extract Distances to Roads- Pennsylvania

#' PASDA Raster
pasda.dist.prim <- terra::rast("Data Management/Rasters/Roads/Pennsylvania/PaRoadRast.Prim.tiff")
pasda.dist.sec <- terra::rast("Data Management/Rasters/Roads/Pennsylvania/PaRoadRast.Sec.tiff")

#' Read in nest data with landcover
pa.nests.landcov.sf <- pa.nests.landcov

#' Extract point value at each nest
dist.prim.out <-terra::extract(pasda.dist.prim, pa.nests.landcov.sf) %>%
  dplyr::select(-ID)%>%
  dplyr::rename("primary" = layer)
dist.sec.out <-terra::extract(pasda.dist.sec, pa.nests.landcov.sf) %>%
  dplyr::select(-ID) %>%
  dplyr::rename("secondary" = layer)
# dist.out <- terra::extract(pa.dist, pa.nests.landcov) %>%
#   dplyr::rename("NearestRoad" = layer)

pa.nests.landcov.sf.roads <- cbind(pa.nests.landcov.sf, dist.prim.out)
pa.nests.landcov.sf.roads <- cbind(pa.nests.landcov.sf.roads, dist.sec.out) 
# pa.nests.landcov.sf.roads <- cbind(pa.nests.landcov.sf.roads, dist.out)

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
                HtBoulder,
                PlotType
                ) 


#' Convert land cover classifications to a categorical variable and create separate columns
pa.nests.covs$Pasture <- ifelse(pa.nests.covs$landuse == "Pasture", 1, 0)
pa.nests.covs$Crop <- ifelse(pa.nests.covs$landuse == "Crop", 1, 0)
pa.nests.covs$Developed <- ifelse(pa.nests.covs$landuse == "Developed", 1, 0)
pa.nests.covs$Deciduous <- ifelse(pa.nests.covs$landuse == "Deciduous Forest", 1, 0)
pa.nests.covs$Evergreen <- ifelse(pa.nests.covs$landuse == "Evergreen Forest", 1, 0)
pa.nests.covs$Mixed <- ifelse(pa.nests.covs$landuse == "Mixed Forest", 1, 0)
pa.nests.covs$Grassland <- ifelse(pa.nests.covs$landuse == "Grassland/Shrub", 1, 0)
pa.nests.covs$Wetland <- ifelse(pa.nests.covs$landuse == "Wetland", 1, 0)
pa.nests.covs$Water <- ifelse(pa.nests.covs$landuse == "Open Water", 1, 0)


################################################################################
## Extract Distances to Roads- New Jersey

#' Read in nest data with landcover
nj.nests.landcov <- nj.nests.landcov

#' Extract point value at each nest
dist.prim.out <-terra::extract(nj.dist.prim, nj.nests.landcov) %>%
  dplyr::select(-ID)%>%
  dplyr::rename("primary" = layer)
dist.sec.out <-terra::extract(nj.dist.sec, nj.nests.landcov) %>%
  dplyr::select(-ID) %>%
  dplyr::rename("secondary" = layer)
dist.out <- terra::extract(nj.dist, nj.nests.landcov) %>%
  dplyr::rename("NearestRoad" = layer)

nj.nests.landcov.sf.roads <- cbind(nj.nests.landcov, dist.prim.out)
nj.nests.landcov.sf.roads <- cbind(nj.nests.landcov.sf.roads, dist.sec.out) 
nj.nests.landcov.sf.roads <- cbind(nj.nests.landcov.sf.roads, dist.out)

nj.nests.covs <- cbind(nj.nests.landcov, nj.nests.landcov.sf.roads) %>%
  dplyr::select(landuse, 
                Case,
                NestID,
                primary,
                secondary)

nj.nests.covs$BandID <- str_sub(nj.nests.covs$NestID, 1, 4)

#' Convert land cover classifications to a categorical variable and create separate columns
nj.nests.covs$Pasture <- ifelse(nj.nests.covs$landuse == "Pasture", 1, 0)
nj.nests.covs$Crop <- ifelse(nj.nests.covs$landuse == "Crop", 1, 0)
nj.nests.covs$Developed <- ifelse(nj.nests.covs$landuse == "Developed", 1, 0)
nj.nests.covs$Deciduous <- ifelse(nj.nests.covs$landuse == "Deciduous Forest", 1, 0)
nj.nests.covs$Evergreen <- ifelse(nj.nests.covs$landuse == "Evergreen Forest", 1, 0)
nj.nests.covs$Mixed <- ifelse(nj.nests.covs$landuse == "Mixed Forest", 1, 0)
nj.nests.covs$Grassland <- ifelse(nj.nests.covs$landuse == "Grassland/Shrub", 1, 0)
nj.nests.covs$Wetland <- ifelse(nj.nests.covs$landuse == "Wetland", 1, 0)
nj.nests.covs$Water <- ifelse(nj.nests.covs$landuse == "Open Water", 1, 0)


################################################################################
## Extract Distances to Roads- Maryland

#' Extract point value at each nest
dist.prim.out <-terra::extract(md.dist.prim, md.nests.landcov) %>%
  dplyr::select(-ID)%>%
  dplyr::rename("primary" = layer)
dist.sec.out <-terra::extract(md.dist.sec, md.nests.landcov) %>%
  dplyr::select(-ID) %>%
  dplyr::rename("secondary" = layer)
dist.out <- terra::extract(md.dist, md.nests.landcov) %>%
  dplyr::rename("NearestRoad" = layer)

md.nests.landcov.sf.roads <- cbind(md.nests.landcov, dist.prim.out)
md.nests.landcov.sf.roads <- cbind(md.nests.landcov.sf.roads, dist.sec.out) 
md.nests.landcov.sf.roads <- cbind(md.nests.landcov.sf.roads, dist.out)

md.nests.covs <- cbind(md.nests.landcov, md.nests.landcov.sf.roads) %>%
  dplyr::select(landuse, 
                Case,
                NestID,
                primary,
                secondary)

md.nests.covs$BandID <- str_sub(md.nests.covs$NestID, 1, 4)

#' Create dummy variables for land cover classifications
md.nests.covs$Pasture <- ifelse(md.nests.covs$landuse == "Pasture", 1, 0)
md.nests.covs$Crop <- ifelse(md.nests.covs$landuse == "Crop", 1, 0)
md.nests.covs$Developed <- ifelse(md.nests.covs$landuse == "Developed", 1, 0)
md.nests.covs$Deciduous <- ifelse(md.nests.covs$landuse == "Deciduous Forest", 1, 0)
md.nests.covs$Evergreen <- ifelse(md.nests.covs$landuse == "Evergreen Forest", 1, 0)
md.nests.covs$Mixed <- ifelse(md.nests.covs$landuse == "Mixed Forest", 1, 0)
md.nests.covs$Grassland <- ifelse(md.nests.covs$landuse == "Grassland/Shrub", 1, 0)
md.nests.covs$Wetland <- ifelse(md.nests.covs$landuse == "Wetland", 1, 0)
md.nests.covs$Water <- ifelse(md.nests.covs$landuse == "Open Water", 1, 0)


################################################################################
## Save RData file *InsertDate_Covs

################################################################################
###############################################################################X