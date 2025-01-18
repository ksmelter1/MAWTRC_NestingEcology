################################################
## Nest-Site Selection of Wild Turkeys in PA ##
################################################
## Kyle J. Smelter
## 3/26/24
## Applied Spatial Ecology Class Project 

#######################
## Load Packages ##
#######################

require(tigris)
require(AICcmodavg)
require(FedData)
require(dplyr)
require(survival)
require(terra)
require(sf)
require(exactextractr)
require(mapview)
require(tidyr)
require(zoo)
require(broom)
require(ggplot2)
require(amt)
require(landscapemetrics)

##########################
## Download NLCD Data ##
##########################

#Use the tigris package to filter through states and remove GEOIDs above 60
#GEOIDs above 60 are territories and islands so they are being removed for scaling
st <- tigris::states() %>%
  dplyr::filter(GEOID < "60") %>% 
  tigris::shift_geometry()

#Transform to albers equal area conic projection, epsg code is 5070
st <- st_transform(st, 5070)

#Grab outline of PA
pa.outline <- subset(st, st$NAME=="Pennsylvania")

#Grab outline of the Mid-Atlantic Region 
#mar.outline <- subset(st, st$NAME== "Pennsylvania"| st$NAME== "Maryland"| st$NAME == "New Jersey") 

#Obtain Pennsylvania NLCD raster
pa.nlcd <- FedData::get_nlcd(template= pa.outline, year = 2019, 
                             label = 'pa', 
                             force.redo = T)

#Obtain study area NLCD information 
#mar.nlcd <- FedData::get_nlcd(template= mar.outline, year = 2019, label = 'mar', force.redo = T)

#Write .tiff files 
#Stored in onedrive folder with rasters
#raster::writeRaster(mar.nlcd, "MAR_NLCD.tiff")
#raster::writeRaster(pa.nlcd, "PA.NLCD.tiff")

################################################################################
################################################################################

##########################
## Study Area Analysis ##
##########################

#Incomplete
#Landscapemetrics analysis doesn't run
#Goal was to get study area statistics for landcover

#Load PA Wildlife Management Units Shapefile from PASDA
pa.wmus <- st_read("shapefiles/WMUs/PGC_BNDWildlifeManagementUnits2021.shp") %>%
  sf::st_transform(5070)

#The hen study takes place here in PA in 2D, 3D, 4D, and 5C
study.area <- subset(pa.wmus, WMU_ID=="2D"| WMU_ID=="3D"| WMU_ID=="4D"| WMU_ID=="5C")

pa.wmus.nlcd <- FedData::get_nlcd(template= study.area, year = 2019, 
                                  label = 'pa', 
                                  force.redo = T)

# m <- c(0, 24, 1, 30,31, 2, 39, 52, 3, 70, 95, 4)
# rclmat <- matrix(m, ncol=3, byrow=TRUE)
# pa.wmus.rc<- classify(pa.wmus.nlcd, rclmat)
# studyarea.rast <- crop(pa.wmus.rc,study.area, mask=TRUE)
# plot(studyarea.rast)

#Calculate the landscape statistics for available metrics
# wmu.l.level= calculate_lsm(studyarea.rast, level = "landscape")
# wmu.l.level

# #Calculate the class level statistics for available metrics
# wmu.c.level= calculate_lsm(studyarea.rast, level = "Class")
# wmu.c.level

#Calculate the patch level statistics for available metrics
# wmu.p.level= calculate_lsm(studyarea.rast, level = "Patch")
# wmu.p.level

###########################
## Landcover Analysis ##
###########################

#Will be able to extract nest locations from entire region so no separate code like distance to road below

#Read in raster with NLCD data from the Mid-Atlantic Region 
#May save space later and not write these out 
#mar.nlcd <- terra::rast("rasters/nlcd/MAR_NLCD.tiff")
#pa.nlcd <- terra::rast("rasters/nlcd/PA.NLCD.tiff")

#1=Developed 
#2=Barren/Rock 
#3= Deciduous Forest
#4= Evergreen Forest
#5= Mixed Forest
#6= Shrub/Scrub
#7= Agriculture
# reclassify the values into 7 groups all values between 0 and 20 equal 1, etc. while removing water (11)
m <- c(0, 24, 1, 30,31, 2, 39,41, 3, 41.5, 42, 4, 42.5, 43, 5, 44, 52, 6, 70, 95, 7)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- terra::classify(pa.nlcd, rclmat)

landuse <- rc


#Read in nests shapefile
#Transform to albers (Hoping to later combine these nests with others from MD and NJ)
#Rename unique identifier column 
pa.nests <- st_read("shapefiles/nests/nests and veg/nest.veg.locations_cleaned1.shp") %>%
  sf::st_transform(5070) %>%
  dplyr::rename("uniqueid"=uniquID) 

#Buffer nests by 80m to simulate nesting home range (Badyaev 1995)
nests_buffered <- sf::st_buffer(pa.nests, 90)
mapview(nests_buffered) #Check

#Create dataframe that has buffered nest values 
#Extract land cover type values using the exact extractr in R
nests.buffered.values <- exactextractr::exact_extract(rc, st_as_sf(nests_buffered))
head(nests.buffered.values[[1]])#Percent in each category (buffers 1-6)
#Create a table of land cover types within each buffer
et=lapply(nests.buffered.values,table)

#Create proportions list object using for loop
prop <- list()
for(i in 1:length(nests.buffered.values)[1] ){
  prop[[i]] <- round((margin.table(et[[i]],1)/margin.table(et[[i]])),digits = 6)
}

M <- coredata(do.call(cbind, lapply(prop, zoo)))
colnames(M) <- NULL
#Transpose matrix so land cover become separate columns of data
matrix <- t(M)
#Now convert the matrix to a data frame so it is easier to manipulate
dfland <- as.data.frame(matrix)
head(dfland) #Check
#Assing column names to land cover
colnames(dfland) <- c("Developed","Barren","Deciduous Forest", "Evergreen Forest", "Agriculture")
head(dfland)

#Function to fill dataframe NA values with 0.000001
fill_NA_with_value <- function(df, value = 0.000001) {
  df[is.na(dfland)] <- 0.000001
  return(df)
}

# Apply the function to fill NA values with 0.0001
dfland.final <- fill_NA_with_value(dfland)

#Write csv
#write.csv(dfland.final, "dfland.final.csv")

#The barren is really surprising to me. I haven't seen that in the turkey literature
#The field crew in 4D hinted at this when I went out trapping though

#Create nests dataframe with proportions
pa.nests.prop <- cbind(pa.nests, dfland.final)

########################
##  Road Shapefiles ##
#######################

#PA State Roads
#WGS to Albers
#pa.state.roads <- sf::st_read("shapefiles/roads/PaStateRoads2023_10.shp") %>%
  #st_transform(5070)

#PA local roads
#WGS to Albers
#pa.local.roads <- sf::st_read("shapefiles/roads/PaLocalRoads2023_10.shp") %>%
  #st_transform(5070)

#Merge shapefiles
#Create PA roads
#pa.roads <-dplyr::bind_rows(pa.state.roads, pa.local.roads)

#Write new shapefile
#st_write(pa.roads, "pa.roads.shp")

###########################################
## Nest Distance from Road Calculations ##
###########################################

#source: https://gis.stackexchange.com/questions/310489/calculating-euclidian-distance-in-r-between-lines-and-points

#Read in roads shapefile
#Drop all roads that have an NA value for ID
#Will check outputs for this later
#Big file
pa.roads <- st_read("shapefiles/roads/pa.roads.shp") %>%
  sf::st_transform(5070) %>%
  dplyr::select(ID) %>%
  tidyr::drop_na(ID)

#Empty lists to store final outputs
pa.nests.out <- list()
road.out <- list()

#Bounding box of nests in PA as a spatial object
#st_as_sfc is class 3 method for bbox
pa.nests_bbox <- sf::st_bbox(pa.nests) %>% st_as_sfc()

#Bounding box of roads in PA as a spatial object
#st_as_sfc is class 3 method for bbox
pa.roads_bbox <- sf::st_bbox(pa.roads) %>% st_as_sfc() 

#For loop to extract the distance from each nest to the nearest road
#Takes me about 40 minutes to run so I've included RData results below
for (i in 1: nrow(pa.nests)) {
  
  #join the nests and roads
  #st_nearest feature: get index of nearest feature
  #Spatial join of pa.nests and pa.roads
  #[[]] first entire object within list, indexing of lists
  neardist <- sf::st_join(
    pa.nests[i,],
    pa.roads, 
    join= sf::st_nearest_feature
  ) 
  
  #find the distance of each nest to the nearest road
  #we know this is in meters by units(dist.to.road)
  dist.to.road <- sf::st_distance(pa.nests[i,],
                                  pa.roads %>%
                                    dplyr::filter(ID==neardist$ID))
  
  #create a df with things we need
  dist.to.road.df <- data.frame(nestid = (pa.nests$uniqueid)[i],
                                plottype = (pa.nests$plottyp)[i],
                                bandid= (pa.nests$bandid)[i],
                                case= (pa.nests$case)[i],
                                NEARDIST = as.numeric(dist.to.road))
  
  #put the df in output
  road.out[[i]] <- dist.to.road.df 
  
  #Printing i here means that for every time the loop runs the row is recorded in the console
  print(i)
  
}

#Create dist2road object
dist2road <- do.call(rbind,(road.out))

#Write csv
#write.csv(dist2road, "Csvs/dist2road.out.csv")
################################################################################

####################
## Organize Data ##
####################

#Load this to avoid running the for loop
#Contains data for everything above 
data <-load("RData/ClassProject.RData")

#Create final.data object for models
#Scale covariates
final.data <- pa.nests.prop %>%
  dplyr::select(nestd_v, bandid, wmu_v, avrgmxv, grssfrb, percwdy,
                case, uniqueid, Developed, Barren, Deciduous.Forest, Evergreen.Forest,
                Agriculture,geometry)  %>%
  mutate(Developed.Z= scale(Developed)[,1])%>%
  mutate(Deciduous.Forest.Z= scale(Deciduous.Forest)[,1])%>%
  mutate(Evergreen.Forest.Z= scale(Evergreen.Forest)[,1]) %>%
  mutate(Agriculture.Z= scale(Agriculture)[,1]) %>%
  mutate(Barren.Z= scale(Barren)[,1]) %>%
  mutate(Visual.Obstruction=scale(avrgmxv)[,1]) %>%
  mutate(GrassForb=scale(grssfrb)[,1]) %>%
  mutate(PercentWoody=scale(percwdy)[,1])

#Change case to numeric
final.data$case <- as.numeric(final.data$case)

#Bring in Distance from Road
final.data$NEARDIST <- dist2road$NEARDIST 
final.data$NEARDIST <- scale(final.data$NEARDIST)

#Drop all rows with NA values
final.data <- tidyr::drop_na(final.data)

###############################
## Step Selection Analysis ##
###############################

#Fit SSFs using conditional logistic regression in the survival package 
#Used nests=1, Stratified randomly sampled available nests=0
#Nest-site metrics, landscape covariates, and distance from road
#May add a distance to early succession covariate
#I only included one nest-site metric per model (Univariate) to avoid collinearity

#Distance from Road
fit1=clogit(case~ NEARDIST + 
                 strata(nestd_v),
               data = final.data)

#Visual Obstruction
fit2=clogit(case~ Visual.Obstruction + 
         strata(nestd_v), 
       data = final.data)

#Developed
fit3=clogit(case~ Developed.Z + 
              strata(nestd_v), 
            data = final.data)

#Barren 
fit4=clogit(case~ Barren.Z +
              strata(nestd_v), 
            data = final.data)

#Forest 
fit5=clogit(case~ Deciduous.Forest.Z +
            strata(nestd_v), 
            data = final.data)

#Forest + Distance from Road
fit6=clogit(case~ Evergreen.Forest.Z + 
            strata(nestd_v), 
            data = final.data)
#GrassForb
fit7=clogit(case~ GrassForb +
               strata(nestd_v), 
             data = final.data)

#Percent Woody
fit8=clogit(case~ PercentWoody +
               strata(nestd_v), 
             data = final.data)

#Check AIC values
AIC(fit2,fit3,fit4,fit5,fit6, fit7, fit8)

#Create AIC Table 
mynames <- paste("fit", as.character(1:16), sep = "")
AIC <- aictab(list(fit1,fit2,fit3,fit4,fit5,fit6, fit7, fit8, fit9, fit10, fit11,
                      fit12, fit13, fit14, fit15, fit16), 
                 modnames = mynames)
print(AIC, LL = FALSE)

#write.csv(AIC, "ModelSelectionClass.csv")

############################
## Extract Coefficients ##
############################

#I only extracted betas from the case~covariate models
Dist2Road <- broom::tidy(fit1$coefficients) %>% 
  filter( names== "NEARDIST") %>%
  bind_cols(confint(fit1)) %>%
  rename(est=x, l95=`2.5 %`, u95=`97.5 %`) 

VisOb <- broom::tidy(fit3$coefficients) %>% 
  filter(names == "Visual.Obstruction") %>%
  bind_cols(confint(fit3)) %>%
  rename(est=x, l95=`2.5 %`, u95=`97.5 %`)

Barren <- broom::tidy(fit8$coefficients) %>% 
  filter(names == "Barren.Z") %>%
  bind_cols(confint(fit8)) %>%
  rename(est=x, l95=`2.5 %`, u95=`97.5 %`)

Developed <- broom::tidy(fit5$coefficients) %>% 
  filter(names== "Developed.Z") %>%
  bind_cols(confint(fit5)) %>%
  rename(est=x, l95=`2.5 %`, u95=`97.5 %`)

Forest <- broom::tidy(fit9$coefficients) %>% 
  filter(names == "Forested.Z") %>%
  bind_cols(confint(fit9)) %>%
  rename(est=x, l95=`2.5 %`, u95=`97.5 %`)

Agriculture <- broom::tidy(fit12$coefficients) %>% 
  filter(names == "Agriculture.Z") %>%
  bind_cols(confint(fit12)) %>%
  rename(est=x, l95=`2.5 %`, u95=`97.5 %`)

PercentWoody <- broom::tidy(fit16$coefficients) %>% 
  filter(names == "PercentWoody") %>%
  bind_cols(confint(fit16)) %>%
  rename(est=x, l95=`2.5 %`, u95=`97.5 %`)

GrassForb <- broom::tidy(fit15$coefficients) %>% 
  filter(names == "GrassForb") %>%
  bind_cols(confint(fit15)) %>%
  rename(est=x, l95=`2.5 %`, u95=`97.5 %`)
  
  ggdf <- bind_rows(Dist2Road,VisOb,Barren, Developed, Forest, Agriculture, 
                    PercentWoody, GrassForb) 
  
  #changing data to factor in order to rename names from rows using levesls
  ggdf$names<-as.factor(ggdf$names)
  levels(ggdf$names)[levels(ggdf$names) == 'NEARDIST'] <- 'Distance from Road'
  levels(ggdf$names)[levels(ggdf$names) == 'Agriculture.Z'] <- 'Percent Agriculture'
  levels(ggdf$names)[levels(ggdf$names) == 'Barren.Z'] <- 'Percent Barren'
  levels(ggdf$names)[levels(ggdf$names) == 'Developed.Z'] <- 'Percent Developed'
  levels(ggdf$names)[levels(ggdf$names) == 'Forested.Z'] <- 'Percent Forest'
  levels(ggdf$names)[levels(ggdf$names) == 'Visual.Obstruction'] <- 'Visual Obstruction'
  levels(ggdf$names)[levels(ggdf$names) == 'GrassForb'] <- 'Grass Forb'
  levels(ggdf$names)[levels(ggdf$names) == 'PercentWoody'] <- 'Percent Woody'
  
  ##########################
  ## Make Plot of Betas ##
  ##########################

  covs <-ggdf %>% 
    ggplot(mapping = aes(x = names, y = est)) +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_point() +
    geom_errorbar(aes(ymin = l95, ymax = u95)
                  ,position=position_dodge2(0.9,reverse = TRUE))+
    theme_classic() +
    xlab("Covariate") +
    ylab("Beta Estimate") +
    coord_flip() +
    scale_x_discrete(limits=rev)
covs  

#########################
## Prediction Plots ##
#########################

#https://stackoverflow.com/questions/47080625/plotting-predictions-from-a-logistic-regression

#Visual Obstruction
Visob.Predict <- ggplot(final.data, aes(x=Visual.Obstruction, y=case)) +
  stat_smooth(method="glm",
  method.args = list(family="binomial"), se=TRUE,
                 fullrange=TRUE) +
  theme_classic()+
    labs(x="Visual Obstruction", y="Probability of Use")+
  ylim(0,1)+
    expand_limits(x=10)
Visob.Predict

#Grass Forb
GrassForb.Predict <- ggplot(final.data, aes(x= GrassForb, y=case)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  theme_classic()+
  labs(x="Grass Forb", y="Probability of Use")+
  ylim(0,1)+
  expand_limits(x=20)
GrassForb.Predict

#Percent Woody
PercWoody.Predict <- ggplot(final.data, aes(x= PercentWoody, y=case)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  theme_classic()+
  labs(x="Percent Woody", y="Probability of Use")+
  ylim(0,1)+
  expand_limits(x=50)
PercWoody.Predict

#Distance from Road
Dist2Road.Predict <- ggplot(final.data, aes(x= NEARDIST, y=case)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  theme_classic()+
  labs(x="Distance from Road", y="Probability of Use")+
  ylim(0,1)+
  expand_limits(x=80)
Dist2Road.Predict

#Agriculture
Agriculture.Z.Predict <- ggplot(final.data, aes(x= Agriculture.Z, y=case)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  theme_classic()+
  labs(x="Percent Agriculture", y="Probability of Use")+
  ylim(0,1)+
  expand_limits(x=100)
Agriculture.Z.Predict

#Forest
Forest.Z.Predict <- ggplot(final.data, aes(x= Forested.Z, y=case)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  theme_classic()+
  labs(x="Percent Forest", y="Probability of Use")+
  ylim(0,1)+
  expand_limits(x=100)
Forest.Z.Predict 

#Developed
Developed.Z.Predict <- ggplot(final.data, aes(x= Developed.Z, y=case)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  theme_classic()+
  labs(x="Percent Developed", y="Probability of Use")+
  ylim(0,1)+
  expand_limits(x=100)
Developed.Z.Predict 

#Barren
Barren.Z.Predict <- ggplot(final.data, aes(x= Barren.Z, y=case)) +
  stat_smooth(method="glm",
              method.args = list(family="binomial"), se=TRUE,
              fullrange=TRUE) +
  theme_classic()+
  labs(x="Percent Barren", y="Probability of Use")+
  ylim(0,1)+
  expand_limits(x=100)
Barren.Z.Predict 

#Odds Ratios
odds.data <- pa.nests.prop %>%
  dplyr::select(avrgmxv, case, nestd_v, percwdy, grssfrb)

odds.data$case <- as.numeric(odds.data$case)

fit17=clogit(case~ avrgmxv +
                 strata(nestd_v), 
               data = odds.data)

summary(fit17)

 fit18=clogit(case~ percwdy +
                      strata(nestd_v), 
                    data = odds.data)
 
 summary(fit18)

 
 fit19= clogit(case~ grssfrb +
                 strata(nestd_v), 
               data = odds.data)
 
 summary(fit19)
 