#'---
#' title: Nest-site selection of female wild turkeys in Pennsylvania (an SSF analysis)
#' author: "K. Smelter, F. Buderman"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document: 
#'     toc: true
#'---
#'
#'  
#' **Purpose**: This script creates a Bayesian conditional logistic regression model for nest-site selection in JAGs using the gathered covariates 

#####################
## Load Packages ##
####################

#' Vector of package names
packages <- c("matrixStats",
              "R2jags",
              "jagsUI",
              "MCMCvis",
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

#' Read in nests
pa.nests <- read_csv("Data Management/Csvs/processed data/Nest-Site Selection/nest.landcov.csv")
pa.nests

#' Read in captures csv
captures <- read_csv("Data Management/Csvs/raw data/captures.csv") 
captures

#' Filter data to include only hens 
captures <- captures %>%
  dplyr::filter(sex == "F") %>%
  dplyr::select(bandid, sex, age, weight, captyr)
captures

#' Create nest.data object for models
#' Scale continuous predictors
nest.data <- pa.nests.covs 

#' Merge data
nest.data <- merge(captures, nest.data, by = "bandid")

#' Add in columns
nest.data$basal<- pa.nests$stemcnt
nest.data$grssfrb <- pa.nests$grssfrb
nest.data$woodytyp <- pa.nests$wdytyp1

###############################
## Prep Woody Type Predictors

nest.data$woodytyp <- ifelse(nest.data$woodytyp %in% c("Decid. Eric", "decid. ericaceous", "Deciduous ericaceous", "evergreen shrub", "Evrgrn Shrub", "grpvns dwn/dead wood", "Grpvines dwn/dead wood"), yes = "Native",
                                                     no = ifelse(nest.data$woodytyp %in% c("Non-native invasive", "Non-native Invasive"), yes = "Invasive",
                                                                 no = "Other"))

#' Create dummy variables for microscale classification
nest.data$Invasive <- ifelse(nest.data$woodytyp == "Invasive", 1, 0)
nest.data$Native <- ifelse(nest.data$woodytyp == "Native", 1, 0)
nest.data$Other <- ifelse(nest.data$woodytyp == "Other", 1, 0)

View(nest.data)

#' (?<=_): Only match is there is an underscore immediately before the number we are trying to extract
#' (\\d{4}): Match exactly 4 digits
#' (?=_) : Only match if the four digits are followed 
nest.data <- nest.data %>%
  dplyr::mutate(nestyr = str_extract(nest_id, "(?<=_)(\\d{4})(?=_)"))
glimpse(nest.data)

nest.data$nestyr <- as.numeric(nest.data$nestyr)
nest.data$captyr <- as.numeric(nest.data$captyr)

#' Create a years since capture column
nest.data <- nest.data %>%
  dplyr::mutate(yrsincecap = nestyr-captyr)
glimpse(nest.data)

#' Select columns of interest
nest.data <- nest.data %>%
  dplyr::select(nest_id, bandid, avrgmxv, grssfrb, percwdy,
                case, Developed, Deciduous, Mixed, Evergreen, Agriculture, elev,
                primary, secondary, basal, weight, age, captyr, nestyr, yrsincecap,
                Native, Invasive, Other, Grassland)

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
#' Since we scaled the data, the mean is zero
fill_na_with_zero <- function(df) {
  df[is.na(df)] <- 0
  return(df)
}

#' Apply the function
nest.data.ready <- fill_na_with_zero(nest.data)
summary(nest.data.ready)

#' Assign juvenile as the reference level
nest.data.ready$age <- ifelse(nest.data.ready$age == "A", 1, 
                              ifelse(nest.data.ready$age == "J", 0, NA))


#' Dealing with scaling age ad hoc
#' If the bird is an adult and the years since capture is >1 assign it as an adult
#' If not keep the age as juvenile because turkeys will nest the first year as a juvenile
nest.data.ready$age <- ifelse(nest.data.ready$age == "A" & nest.data.ready$yrsincecap >= 1, "A", nest.data.ready$age)
glimpse(nest.data.ready)

summary(nest.data.ready)
str(nest.data.ready$age)

###############################################
## JAGS Nest-site Selection Model

JAGS.mod <- "  
 model{
  for (i in 1:I){
   use[i]~dpois(lambda[i])

   log(lambda[i])<-inprod(beta,X[i,])+alpha[str_ID[i]]
  }
  
  #Priors
  for(j in 1:J){
   beta[j]~dnorm(0,0.0001)
  }
  
  for(k in 1:K){
   alpha[k]~dnorm(0,0.000001)
  }
  
 }
"
#' Model parameters
#' Leave Grasland/Shrub and Developed out of model
X <- cbind(
  rep(1, nrow(nest.data.ready)),  # Intercept (1)
  nest.data.ready$primary,         # Distance to Primary Road
  nest.data.ready$secondary,       # Distance to Secondary Road
  nest.data.ready$Mixed,           # Mixed Forest
  nest.data.ready$Evergreen,       # Evergreen Forest
  nest.data.ready$Deciduous,       # Deciduous Forest
  nest.data.ready$Agriculture,     # Agriculture
  nest.data.ready$elev,            # Elevation
  nest.data.ready$avgvisob,        # Visual Obstruction
  nest.data.ready$grssfrb,         # Percent Grass/Forb
  nest.data.ready$percwdy,         # Percent Woody
  nest.data.ready$basal,           # Basal Area
  nest.data.ready$Invasive,        # Invasive
  nest.data.ready$Native,          # Native
  nest.data.ready$Grassland        # Grassland/Shrub
)

#' Data prep for JAGS
jags_data <- list(use = nest.data.ready$case,
                  X=X,
                  str_ID=as.numeric(nest.data.ready$str_ID),
                  I=nrow(nest.data.ready),
                  J=ncol(X),
                  K=length(unique(nest.data.ready$str_ID)))

#' Initialize JAGS
jags_model <- jags.model(textConnection(JAGS.mod), 
                         data = jags_data,
                         n.chains = 1,
                         n.adapt=2000)

#' Sample from posterior distribution using parallel
jags_samples <- coda.samples(jags_model, c("beta"),
                             n.iter = 40000, 
                             n.burnin = 10000,
                             n.thin = 1,
                             do.parallel = T)
colMeans(jags_samples[[1]])
colSds(jags_samples[[1]])

#' View traceplots
MCMCtrace(jags_samples, pdf = F)

#' Convert mcmc.list to matrix
samples_matrix <- as.matrix(jags_samples)

#' Convert the matrix to a data frame
samples_df <- as.data.frame(samples_matrix) 

#' View the first few rows of the samples data frame
head(samples_df)

#' Check the actual column names to confirm their structure
colnames(samples_df)

#' Create a vector of new names
new_names <- c("Intercept", "Distance to Primary Road", "Distance to Secondary Road", 
               "Mixed Forest", "Evergreen Forest", "Deciduous Forest", "Agriculture", 
               "Elevation", "Visual Obstruction", "Percent Grass/Forb", 
               "Percent Woody Vegetation", "Basal Area", "Invasive Woody Vegetation", 
               "Native Woody Vegetation", "Grassland/Shrub")

#' Rename the columns based on position
colnames(samples_df) <- new_names

#' View the renamed data frame
head(samples_df)

#' Create RDS file of model outputs
saveRDS(jags_samples, "01_NestModelSamples.Raw.RDS_1_10")
saveRDS(samples_df, "02_NestModelSamples.RDS_1_10")

################################################################################
################################################################################
