#'---
#' title: Nest-site selection of female wild turkeys in Pennsylvania (an SSF analysis)
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
#' **Purpose**: This script creates a Bayesian conditional logistic regression model for nest-site selection in JAGs using the gathered covariates 

#####################
## Load Packages ##
####################

#' Vector of package names
packages <- c("R2jags",
              "MCMCvis",
              "dplyr",
              "tidyverse")

#' Function to load a package or install it if not already installed
load_packages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    require(package_name, character.only = TRUE)
  }
}

#' Apply the function to each package name
lapply(packages, load_packages)

########################
## Data Preparation ##
########################

#' Load in RData
load("Data Management/RData/Nest-Site Selection/Covs/02_Covs.RData")

#' Read in captures csv
captures <- read_csv("Data Management/Csvs/raw data/captures.csv") 
captures

#' Filter data to include only hens 
captures <- captures %>%
  dplyr::filter(sex == "F") %>%
  dplyr::select(bandid, sex, age, weight)
captures

#' Create nest.data object for models
#' Scale continuous predictors
nest.data <- pa.nests.covs 

#' Merge data
nest.data <- merge(captures, nest.data, by = "bandid")

#' Add in columns
nest.data$basal <- pa.nests$stemcnt
nest.data$grssfrb <- pa.nests$grssfrb

#' Select columns of interest
nest.data <- nest.data %>%
  dplyr::select(nest_id, bandid, avrgmxv, grssfrb, percwdy,
                case, Developed, Deciduous, Mixed, Evergreen, Agriculture, elev,
                primary, secondary, basal, weight, age)

#' Scale predictors
nest.data <- nest.data %>%
  dplyr::mutate(basal = scale(basal)) %>%
  dplyr::mutate(percwdy = scale(percwdy)) %>%
  dplyr::mutate(grssfrb = scale(grssfrb)) %>%
  dplyr::mutate(avgvisob = scale(avrgmxv)) %>%
  dplyr::mutate(elev = scale(elev)) %>%
  dplyr::mutate(primary = scale(primary)) %>%
  dplyr::mutate(secondary = scale(secondary)) %>%
  dplyr::mutate(weight = scale(weight)) 


str(nest.data)
glimpse(nest.data)

#' Change structure from a matrix to numeric (Scaled covariates)
nest.data <- nest.data %>%
  dplyr::mutate(percwdy = as.numeric(percwdy)) %>%
  dplyr::mutate(grssfrb = as.numeric(grssfrb)) %>%
  dplyr::mutate(avgvisob = as.numeric(avgvisob)) %>%
  dplyr::mutate(elev = as.numeric(elev)) %>%
  dplyr::mutate(primary = as.numeric(primary)) %>%
  dplyr::mutate(secondary = as.numeric(secondary)) %>%
  dplyr::mutate(basal = as.numeric(basal)) %>%
  dplyr::mutate(weight = as.numeric(weight)) %>%
  dplyr::mutate(Deciduous = as.numeric(Deciduous)) %>%
  dplyr::mutate(Mixed = as.numeric(Mixed)) %>%
  dplyr::mutate(Evergreen= as.numeric(Evergreen)) %>%
  dplyr::mutate(Agriculture = as.numeric(Agriculture)) %>%
  dplyr::mutate(Developed = as.numeric(Developed))

glimpse(nest.data)
str(nest.data)

#' Change case to numeric
nest.data$case <- as.numeric(nest.data$case)

#' Order df by nestid
nest.data <- nest.data[order(nest.data$nest_id),]

#' Change Nest_ID_V to numeric 
#' First must remove underlines
nest.data$nest_id <- gsub("_", "", nest.data$nest_id)
nest.data$nest_id<- as.numeric(nest.data$nest_id)

#' Create str_id column 
#' This allows the loop to iterate through the steps associated with each bird 
#' cur_group_id() gives a unique numeric identifier for the current group.
nest.data <- nest.data %>%
  group_by(nest_id) %>%
  mutate(str_ID=cur_group_id())

#' Check
str(nest.data)
glimpse(nest.data)
summary(nest.data)

#' Function to replace all NAs with 0 in a dataframe
#' Since we scaled the data, the mean is zero so we are filling all values as the mean
fill_na_with_zero <- function(df) {
  df[is.na(df)] <- 0
  return(df)
}

#' Apply function
nest.data.ready <- fill_na_with_zero(nest.data)

summary(nest.data.ready)

#' Bring in age 
nest.data.ready$age <- ifelse(nest.data.ready$age == "A", 1, 
                              ifelse(nest.data.ready$age == "J", 0, NA))


# Identify if the individual was initially a juvenile
# nest.data.ready$initial_age <- ifelse(nest.data.ready$age == "J", 1, 0)
# 
# Create a variable for the time since banding (if `year` is available)
# nest.data.ready$years_since_banding <- nest.data.ready$year - nest.data.ready$year


# nest.data.ready$age <- ifelse(nest.data.ready$initial_age == 1 & nest.data.ready$years_since_banding > 1, "A", nest.data.ready$age)

# Recreate dummy variables based on the updated age
# nest.data.ready$age.A <- ifelse(nest.data.ready$age == "A", 1, 0)
# nest.data.ready$age.J <- ifelse(nest.data.ready$age == "J", 1, 0)

#' Check
summary(nest.data.ready)
str(nest.data.ready$age)

######################################################
## JAGS Models

####################
## No Developed ##
####################
#' We keep developed out of the model because it is correlated with agriculture

nest.selection.subset.mod <- "  
 model{
  for (i in 1:I){
   case[i]~dpois(lambda[i])
   log(lambda[i])<-primary[i]*beta.primary +
   secondary[i]*beta.secondary +
   deciduous[i]*beta.deciduous_forest + 
   mixed[i]*beta.mixed_forest + 
   evergreen[i]*beta.evergreen_forest + 
   agriculture[i]*beta.agriculture + 
   elevation[i] *beta.elevation + 
   visob[i]* beta.visob +
   grassforb[i]*beta.grassforb +
   percwoody[i]*beta.percwoody + 
   basalarea[i]*beta.basalarea +
   age[i] * beta.age
   alpha[str_ID[i]] 
  }
  
  #Priors
  beta.deciduous_forest~dnorm(0,0.0001)  #Deciduous Forest
  beta.mixed_forest~dnorm(0,0.0001)      #Mixed Forest
  beta.evergreen_forest~dnorm(0,0.0001)  #Mixed Forest
  beta.agriculture~dnorm(0,0.0001)       #Agriculture
  beta.primary~dnorm(0,0.0001)           #Primary
  beta.secondary~dnorm(0,0.0001)         #Secondary
  beta.elevation~dnorm(0,0.0001)         #Elevation
  beta.visob~dnorm(0,0.0001)             #Average Visual Obstruction
  beta.grassforb~dnorm(0,0.0001)         #Percent Grass Forb
  beta.percwoody~dnorm(0,0.0001)         #Percent Woody
  beta.age~dnorm(0,0.0001)             #Age Adult (Juvenile is the reference level)
  beta.basalarea~dnorm(0,0.0001)         #Basal Area
  
  for (k in 1:K){
   alpha[k]~dnorm(0,0.000001)
  }
 }
"
#' Data list for JAGS
jags_data <- list(
  case = nest.data.ready$case,
  primary = nest.data.ready$primary,
  secondary = nest.data.ready$secondary,
  mixed = nest.data.ready$Mixed,
  evergreen = nest.data.ready$Evergreen,
  deciduous = nest.data$Deciduous,
  agriculture = nest.data.ready$Agriculture,
  elevation = nest.data.ready$elev,
  visob = nest.data.ready$avgvisob,
  grassforb = nest.data.ready$grssfrb,
  percwoody = nest.data.ready$percwdy,
  age = nest.data.ready$age,
  basalarea = nest.data.ready$basal,
  str_ID = c(nest.data.ready$str_ID),
  I = nrow(nest.data.ready),
  K =length(unique(nest.data.ready$str_ID))
)

#' Initialize JAGS
jags_model <- jags.model(textConnection(nest.selection.subset.mod), 
                         data = jags_data,
                         n.chains = 3)

#' Burn-in and sampling from posterior distribution
jags_samples.nss.com1 <- coda.samples(jags_model, c("beta.deciduous_forest",
                                                    "beta.mixed_forest",
                                                    "beta.evergreen_forest",
                                                    "beta.agriculture",
                                                    "beta.primary",
                                                    "beta.secondary",
                                                    "beta.elevation",
                                                    "beta.visob",
                                                    "beta.grassforb",
                                                    "beta.percwoody",
                                                    "beta.age.A",
                                                    "beta.basalarea"), 
                                      n.iter = 200000, 
                                      n.burnin = 30000,
                                      n.thin = 5)
#' View traceplots
MCMCtrace(jags_samples.nss.com1, pdf=F)
