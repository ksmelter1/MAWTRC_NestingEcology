
#'---
#' title: Habitat selection of female wild turkeys during pre-incubation (an SSF analysis)
#' author: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output: RandomSteps_Prep(R Workspace)
#'   html_document: 
#'     toc: true
#'---
knitr::opts_chunk$set(warning = FALSE, message = FALSE, cache = TRUE)
#'  
#' **Purpose**: This script checks on the variability of NLCD covs across 4 scales 
#'

################################
## Step Length Distributions ##
################################

#' Create random subset of data from hen 8262 during her first nesting attempt in 2022
#' Remove NonForest cover type
random_steps_subset <- random_steps[
  (random_steps$id %in% c("8262_2022_1_1")) &
    (random_steps$landuse != "NonForest"), 
]

#' Maximum distance moved by 8262 during 30 minute interval
#' 599.27 m
max(random_steps_subset$sl_)

##############################################
## Checking Proportions at the point scale ##
##############################################

#' Convert the landuse and case columns to factors
random_steps_subset$landuse <- as.factor(random_steps_subset$landuse)
random_steps_subset$case_ <- as.factor(random_steps_subset$case_)

#'Creating a frequency table that combines landuse and case
table(random_steps_subset$landuse, random_steps_subset$case_)

###########################################
## Checking Proportions at Patch Scale ##
##########################################

random_steps.sf.patch <- st_as_sf(random_steps_subset, coords = c("x2_", "y2_"), crs= 5070) %>%
  st_buffer(200)


#' Use terra to extract the proportion of land cover categories within each buffer
#' Use mutate to create new columns for each land cover category
#' Developed, Forest, Agriculture, Other
#' Check variability between used and available steps
landcov_table <- terra::extract(x = pa.nlcd, y = vect(random_steps.sf.patch), fun = "table")%>%
  dplyr::mutate(pct_forest = Forest/(Forest+Shrub+Developed+Agriculture)) %>%
  dplyr::mutate(pct_developed = Developed/(Developed+Shrub+Forest+Agriculture)) %>%
  dplyr::mutate(pct_agriculture= Agriculture/(Agriculture+Shrub+Forest+Developed)) %>%
  dplyr::mutate(pct_shrub = Shrub/(Shrub+Forest+Developed+Agriculture))

#' Cbind hen subset and landcover table
random_steps.sf.landcov.patch <- cbind(random_steps_subset, landcov_table)

#' Change percentages to numeric
random_steps.sf.landcov.patch$pct_forest <- as.numeric(random_steps.sf.landcov.patch$pct_forest)
random_steps.sf.landcov.patch$pct_agriculture <- as.numeric(random_steps.sf.landcov.patch$pct_agriculture)
random_steps.sf.landcov.patch$pct_developed <- as.numeric(random_steps.sf.landcov.patch$pct_developed)

#' Create two classes for proportion forest
prop.cat = cut(
  random_steps.sf.landcov.patch$pct_forest, 
  breaks = seq(0, 1, by = 0.5),
  labels = c("0-0.5","0.5-1"),
  include.lowest = TRUE,
  right = FALSE
)

table(prop.cat, random_steps.sf.landcov.patch$case_)

################################################
## Checking Proportions at Home Range Scale ##
###############################################

random_steps.sf.hr <- st_as_sf(random_steps_subset, coords = c("x2_", "y2_"), crs= 5070) %>%
  st_buffer(500)


#' Use terra to extract the proportion of land cover categories within each buffer
#' Use mutate to create new columns for each land cover category
#' Developed, Forest, Agriculture, Other
#' Check variability between used and available steps
landcov_table.hr <- terra::extract(x = pa.nlcd, y = vect(random_steps.sf.hr), fun = "table")%>%
  dplyr::mutate(pct_forest = Forest/(Forest+Shrub+Developed+Agriculture)) %>%
  dplyr::mutate(pct_developed = Developed/(Developed+Shrub+Forest+Agriculture)) %>%
  dplyr::mutate(pct_agriculture= Agriculture/(Agriculture+Shrub+Forest+Developed)) %>%
  dplyr::mutate(pct_shrub = Shrub/(Shrub+Forest+Developed+Agriculture))

random_steps.sf.landcov.hr<- cbind(random_steps_subset, landcov_table.hr)

#' Change percentages to numeric
random_steps.sf.landcov.hr$pct_forest <- as.numeric(random_steps.sf.landcov.hr$pct_forest)
random_steps.sf.landcov.hr$pct_agriculture <- as.numeric(random_steps.sf.landcov.hr$pct_agriculture)
random_steps.sf.landcov.hr$pct_developed <- as.numeric(random_steps.sf.landcov.hr$pct_developed)

#' Create two classes for proportion forest
#' >0.5 and <0.5
prop.cat = cut(
  random_steps.sf.landcov.hr$pct_forest, 
  breaks = seq(0, 1, by = 0.5),
  labels = c("0-0.5","0.5-1"),
  include.lowest = TRUE,
  right = FALSE
)

#' Visualize by case
table(prop.cat, random_steps.sf.landcov.hr$case)


################################################
## Checking Proportions at Landscape Scale ##
###############################################

random_steps.sf.land <- st_as_sf(random_steps_subset, coords = c("x2_", "y2_"), crs= 5070) %>%
  st_buffer(1000)

#' Use terra to extract the proportion of land cover categories within each buffer
#' Use mutate to create new columns for each land cover category
#' Developed, Forest, Agriculture, Other
#' Check variability between used and available steps
landcov_table.land <- terra::extract(x = pa.nlcd, y = vect(random_steps.sf.land), fun = "table")%>%
  dplyr::mutate(pct_forest = Forest/(Forest+Shrub+Developed+Agriculture)) %>%
  dplyr::mutate(pct_developed = Developed/(Developed+Shrub+Forest+Agriculture)) %>%
  dplyr::mutate(pct_agriculture= Agriculture/(Agriculture+Shrub+Forest+Developed)) %>%
  dplyr::mutate(pct_shrub = Shrub/(Shrub+Forest+Developed+Agriculture))

random_steps.sf.landcov.land<- cbind(random_steps_subset, landcov_table.land)

#' Change percentages to numeric
random_steps.sf.landcov.land$pct_forest <- as.numeric(random_steps.sf.landcov.land$pct_forest)
random_steps.sf.landcov.land$pct_agriculture <- as.numeric(random_steps.sf.landcov.land$pct_agriculture)
random_steps.sf.landcov.land$pct_developed <- as.numeric(random_steps.sf.landcov.land$pct_developed)

#' Create two classes for proportion forest
prop.cat = cut(
  random_steps.sf.landcov.land$pct_forest, 
  breaks = seq(0, 1, by = 0.5),
  labels = c("0-0.5","0.5-1"),
  include.lowest = TRUE,
  right = FALSE
)

table(prop.cat, random_steps.sf.landcov.land$case_)

################################################################################
################################################################################